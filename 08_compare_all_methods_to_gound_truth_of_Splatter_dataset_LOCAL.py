#!/usr/bin/env python
# coding: utf-8
import os
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
import pickle5 as pickle
import multiprocessing as mp
#from funcs import *
import time
import json


def method_in_2_lines(method):
    return method.replace(' (', '\n(')


def compute_sum_of_diff(clustering1, clustering2, nclus1, nclus2, gene_list, ind1, ind2):
    d = 0
    for ig1, g1 in enumerate(gene_list[ind1:ind2]):
        if ig1 % 50 == 0:
            print(ig1, '/', ind2-ind1)
        for ig2, g2 in enumerate(gene_list[ind1+ig1+1:]):
            #print(ind1+ig1, ind1+ig1+1+ig2)
            E1 = 1 - np.linalg.norm(clustering1[ind1+ig1]-clustering1[ind1+ig1+1+ig2], ord=1)/nclus1
            E2 = 1 - np.linalg.norm(clustering2[ind1+ig1]-clustering2[ind1+ig1+1+ig2], ord=1)/nclus2
            print(E1, '\t', E2)
            d += np.abs(E1-E2)
            if ig2 == 1000:
                break
        break
    return d


def ARI_Hullermeier2012(gene_list, cluster_table1, cluster_table2):
    ngenes = len(gene_list)

    # create gene X cluster tables
    nclus1 = len(cluster_table1)
    clustering1 = pd.DataFrame(index=gene_list, columns=cluster_table1.index,
                               data=np.zeros((ngenes, nclus1)))
    for clus1 in cluster_table1.index:
        # read cluster gene list
        genes = cluster_table1.loc[clus1, 'genes'].replace("['", "").replace("']", "").split("', '")
        clustering1.loc[genes, clus1] = 1

    nclus2 = len(cluster_table2)
    clustering2 = pd.DataFrame(index=gene_list, columns=cluster_table2.index,
                               data=np.zeros((ngenes, nclus2)))
    for clus2 in cluster_table2.index:
        # read cluster gene list
        genes = cluster_table2.loc[clus2, 'genes'].replace("['", "").replace("']", "").split("', '")
        clustering2.loc[genes, clus2] = 1

    # divide to small tasks and run in parallel
    clustering1 = np.array(clustering1)
    clustering2 = np.array(clustering2)
    tic = time.perf_counter()
    with mp.Pool(mp.cpu_count()) as pool:
        d_list = pool.starmap(compute_sum_of_diff,
                              [(clustering1, clustering2, nclus1, nclus2, gene_list, inds[0], inds[1]) for inds in
                               [(0, 1000),
                                (1000, 2000),
                                (2000, 3000),
                                (3000, 4000),
                                (4000, 5000),
                                (5000, 6000),
                                (6000, 7000),
                                (7000, 8000),
                                (8000, 9000),
                                (9000, 10000)]])
    print(d_list)
    toc = time.perf_counter()
    print(f"{toc - tic:0.4f} seconds")
    return 1 - sum(d_list)/(ngenes*(ngenes-1)/2)


def purity(cluster_table1, cluster_table2):    
    similarity = pd.DataFrame(np.zeros((len(cluster_table1), len(cluster_table2))),
                             columns=cluster_table2.index, index=cluster_table1.index) 
    for clus1 in cluster_table1.index:
        # read cluster gene list
        genes_clus1 = cluster_table1.loc[clus1, 'genes'].replace("['", "").replace("']", "").split("', '")
        for clus2 in cluster_table2.index:
            # read cluster gene list
            genes_clus2 = cluster_table2.loc[clus2, 'genes'].replace("['", "").replace("']", "").split("', '")
            # compute similarity index
            similarity.loc[clus1, clus2] = len(set(genes_clus1).intersection(set(genes_clus2)))
    return sum(similarity.max(axis=1))/similarity.sum().sum(), sum(similarity.max())/similarity.sum().sum()


def create_heatmap(cluster_table1, cluster_table2, similarity_func, similarity_func_name, 
                   name_method1, name_method2, name_clus1, name_clus2, out_folder):
    sns.set(font_scale=1.4)
    
    similarity = np.empty(shape=(len(cluster_table1), len(cluster_table2)), dtype=float) 
    for clus1_ind, clus1 in enumerate(cluster_table1.index):
        # read cluster gene list
        genes_clus1 = cluster_table1.loc[clus1, 'genes'].replace("['", "").replace("']", "").split("', '")
        for clus2_ind, clus2 in enumerate(cluster_table2.index):
            # read cluster gene list
            genes_clus2 = cluster_table2.loc[clus2, 'genes'].replace("['", "").replace("']", "").split("', '")

            # compute similarity index
            similarity[clus1_ind][clus2_ind] = similarity_func(set(genes_clus1), set(genes_clus2))
        
    similarity_df = pd.DataFrame(similarity, index=cluster_table1.index, columns=cluster_table2.index)
    
    # save df
    similarity_df.to_excel(os.path.join(out_folder, similarity_func_name + '_similarities_' 
                                        + name_method1 + '_vs_' + name_method2 + '.xlsx'))
    
    if similarity.shape[0]<16 and similarity.shape[1]<16:
        annot = True
    else:
        annot = False
        
    # save heatmap
    plt.figure(figsize=(20, 10))
    b = sns.heatmap(similarity, annot=annot, fmt='.2f', 
                yticklabels=[name_clus1 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(cluster_table1.index, cluster_table1['#genes'])], 
                xticklabels=[name_clus2 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(cluster_table2.index, cluster_table2['#genes'])])
    if similarity.shape[1]>=16:
        b.set_xticklabels(labels=[name_clus2 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(cluster_table2.index, cluster_table2['#genes'])], size=12)
        
    plt.savefig(os.path.join(out_folder, similarity_func_name + '_heatmap_' 
                             + name_method1 + '_vs_' + name_method2 + '.jpg'), bbox_inches = "tight")
    plt.savefig(os.path.join(out_folder, similarity_func_name + '_heatmap_' 
                             + name_method1 + '_vs_' + name_method2 + '.eps'), bbox_inches = "tight", format='eps')
    plt.close()
    
    # pairwise ordering of the heatmap
    ordered_clus1 = []
    ordered_clus2 = []

    left_clus1 = list(cluster_table1.index)
    left_clus2 = list(cluster_table2.index)

    while left_clus1 and left_clus2:
        #maxx = np.max(similarity[np.c_[np.array(left_clus1)], np.array(left_clus2)])
        #inds = np.where(similarity == maxx)
        
        maxx = similarity_df.loc[left_clus1, left_clus2].max().max()
        inds = [np.flatnonzero((similarity_df.loc[left_clus1, left_clus2]==maxx).values)//similarity_df.loc[left_clus1, left_clus2].shape[1], 
                np.unique(np.flatnonzero((similarity_df.loc[left_clus1, left_clus2]==maxx).values)%similarity_df.loc[left_clus1, left_clus2].shape[1])]
        
        for ind1, ind2 in zip(inds[0], inds[1]):
            clus1 = list(left_clus1)[ind1]
            clus2 = list(left_clus2)[ind2]
            if clus1 in left_clus1 and clus2 in left_clus2:
                break
                
        ordered_clus1.append(clus1)
        ordered_clus2.append(clus2)

        left_clus1.remove(clus1)
        left_clus2.remove(clus2)

    if left_clus1:
        ordered_clus1 += left_clus1
    if left_clus2:
        ordered_clus2 += left_clus2

    #similarity_ordered = similarity[np.c_[np.array(ordered_clus1)], np.array(ordered_clus2)]
    similarity_df_ordered = similarity_df.loc[ordered_clus1, ordered_clus2]
    similarity_ordered = np.array(similarity_df_ordered)
    
    if similarity_ordered.shape[0]<16 and similarity_ordered.shape[1]<16:
        annot = True
    else:
        annot = False
    
    plt.figure(figsize=(20, 10))
    b = sns.heatmap(similarity_ordered, annot=annot, fmt='.2f',
                yticklabels=[name_clus1 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(ordered_clus1, cluster_table1.loc[ordered_clus1, "#genes"])], 
                xticklabels=[name_clus2 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(ordered_clus2, cluster_table2.loc[ordered_clus2, "#genes"])])
    
    if similarity_ordered.shape[1]>=16:
        b.set_xticklabels(labels=[name_clus2 + ' ' + str(cl) + ' (' + str(ngenes) + ' genes)' for cl, ngenes in 
                             zip(ordered_clus2, cluster_table2.loc[ordered_clus2, "#genes"])], size=12)

    plt.savefig(os.path.join(out_folder, similarity_func_name + '_heatmap_' 
                             + name_method1 + '_vs_' + name_method2 + '_ordered.jpg'), 
                bbox_inches = "tight")
    plt.savefig(os.path.join(out_folder, similarity_func_name + '_heatmap_' 
                             + name_method1 + '_vs_' + name_method2 + '_ordered.eps'), 
                bbox_inches = "tight", format='eps')
    plt.close()
    
    return similarity_df

if __name__ == '__main__':
    data_n = 7    #Splatter- simulated dataset
    out_folder = "comparison_of_all_methods_to_gound_truth_of_Splatter_dataset"
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    # load true cluster table
    true_file = 'D:\\OneDrive - Technion\\SPIRAL_simulated_synthetic_data\\true_cluster_table.csv'
    true_cluster_table = pd.read_csv(true_file, index_col=0)

    # Compare methods to ground truth and produce heatmaps
    data_path = "D:\\OneDrive - Technion\\SPIRAL\\static\\analysis\\data" + str(data_n)
    species = open(os.path.join(data_path, 'species.txt'), "r").read()
    impute_method = open(os.path.join(data_path, 'imputation_method.txt'), "r").read()

    #get genelist for this dataset
    '''
    genelist_file = os.path.join(data_path, impute_method + '_genetable.p')
    # Load genes_table
    with open(genelist_file, 'rb') as fp:
        genes_table = pickle.load(fp)
    gene_list = list(genes_table.index)
    '''
    raw_data = pd.read_csv(os.path.join(data_path, 'counts.txt'), index_col=0)
    gene_list = list(raw_data.index)

    # For each true cluster, what is its best Jaccard index
    best_jaccard = pd.DataFrame(index=range(len(true_cluster_table)))
    method_sizes = dict()

    purity1_dict = dict()
    purity2_dict = dict()

    rand_ind_dict = dict()

    for method in [
        'SPIRAL',
        #'gene_clustering',
        #'cell_clustering',
        'Seurat',
        'Slingshot+tradeSeq (nPoints=20)',
        'Slingshot+tradeSeq (nPoints=200)',
        "Hotspot (n_neighbors=30)",
        "Hotspot (n_neighbors=300)",
        #"Hotspot (model='danb', n_neighbors=30)",
        #"Hotspot (model='danb', n_neighbors=300)",
        #"Hotspot (model='bernoulli', n_neighbors=30)",
        #"Hotspot (model='bernoulli', n_neighbors=300)",
        'nsNMF'
                  ]:

        flag = False

        if method in [
            #'gene_clustering', 'cell_clustering',
            'nsNMF', 'Seurat',
            'Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)',
            'Hotspot', "Hotspot (n_neighbors=30)", "Hotspot (n_neighbors=300)",
            #"Hotspot (model='danb', n_neighbors=30)",
            #"Hotspot (model='danb', n_neighbors=300)",
            #"Hotspot (model='bernoulli', n_neighbors=30)",
            #"Hotspot (model='bernoulli', n_neighbors=300)",
                     ]:
            '''
            if method == 'gene_clustering':
                method_folder = 'comparison_to_gene_clustering_filter_genes_exp_and_coef/data' + str(data_n) + '/'
                name_clus2 = 'cluster'
            elif method == 'cell_clustering':
                method_folder = 'comparison_to_DE_genes/data' + str(data_n) + '/'
                name_clus2 = 'DE gene set'
            '''
            if method == 'nsNMF':
                method_folder = '../SPIRAL_for_paper/comparison_to_nsNMF/data' + str(data_n) + '/'
                name_clus2 = 'module'
            elif method == 'Seurat':
                method_folder = '../SPIRAL_for_paper/comparison_to_Seurat_DE_genes/data' + str(data_n) + '/'
                name_clus2 = 'DE gene set'
            elif method in ['Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)']:
                method_folder = '../SPIRAL_for_paper/comparison_to_Slingshot/data' + str(data_n) + '/'
                name_clus2 = 'gene cluster'
            elif 'Hotspot' in method:
                method_folder = '../SPIRAL_for_paper/comparison_to_Hotspot/data' + str(data_n) + '/'
                name_clus2 = 'gene module'

            if method == 'Slingshot+tradeSeq (nPoints=20)':
                file = 'cluster_table_20.csv'
                method_cluster_table_file = os.path.join(method_folder, file)
                flag = True
            elif method == 'Slingshot+tradeSeq (nPoints=200)':
                file = 'cluster_table_200.csv'
                method_cluster_table_file = os.path.join(method_folder, file)
                flag = True
            elif 'Hotspot' in method:
                if method in ["Hotspot (model='danb', n_neighbors=30)", "Hotspot (n_neighbors=30)"]:
                    file = 'cluster_table_danb_30.csv'
                elif method in ["Hotspot (model='danb', n_neighbors=300)", "Hotspot (n_neighbors=300)"]:
                    file = 'cluster_table_danb_300.csv'
                elif method=="Hotspot (model='bernoulli', n_neighbors=30)":
                    file = 'cluster_table_bernoulli_30.csv'
                elif method=="Hotspot (model='bernoulli', n_neighbors=300)":
                    file = 'cluster_table_bernoulli_300.csv'
                method_cluster_table_file = os.path.join(method_folder, file)
                flag = True
            else:
                # load method's cluster table
                for file in os.listdir(method_folder):
                    if file.startswith("cluster_table"):
                        print(file)
                        method_cluster_table_file = os.path.join(method_folder, file)
                        flag = True

            if flag:
                method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)

                if method in ['Seurat', 'Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)']:
                    method_cluster_table.loc[:, 'genes'] = [
                        str(genes.replace("['", "").replace("']", "").split(", ")) for genes in method_cluster_table.loc[:, 'genes']]

                print('\n\n', method, method_cluster_table_file)

        elif method == 'SPIRAL':
            flag = True
            folder = 'D:\\OneDrive - Technion\\SPIRAL\\static\\analysis\\data' + str(data_n)
            impute_method = 'agg_wald'
            name_clus2 = 'structure'

            struct_table_file = os.path.join(folder, impute_method + '_sigtable_filt_GO_vis.xlsx')
            ANALYSIS_FOLDER = 'D:\\OneDrive - Technion\\SPIRAL\\static\\analysis'
            data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
            genetable_file = os.path.join(data_path, impute_method + '_genetable.p')
            with open(genetable_file, 'rb') as fp:
                genes_table = pickle.load(fp)

            struct_table = pd.read_excel(struct_table_file, index_col=0)

            #mode = 'all_structs'
            mode = 'shorter_list'

            if mode == 'shorter_list':
                struct_list = [1, 2, 3, 9, 10, 43]
                struct_list_str = str(struct_list)

            elif mode == 'all_structs':
                struct_list = list(struct_table.index)
                struct_list_str = 'all'

            method_cluster_table = struct_table.loc[struct_list, ['num_genes_in_struct', 'genes']].rename(columns={'num_genes_in_struct': '#genes'})
            method_cluster_table.loc[:, 'genes'] = [str(a.split(',')) for a in method_cluster_table.loc[:, 'genes']]

        if flag:
            print(method, '\n', method_cluster_table.head())

            method_sizes[method] = list(method_cluster_table.loc[:, '#genes'])
            '''
            # create heatmap- Jaccard
            similarity_df_jaccard = create_heatmap(cluster_table1=true_cluster_table, cluster_table2=method_cluster_table,
                                   similarity_func=jaccard, similarity_func_name='jaccard',
                                   name_method1='true', name_method2=method, 
                                   name_clus1='true cluster', name_clus2=name_clus2,
                                   out_folder=out_folder)
            
            # create heatmap- Dice
            create_heatmap(cluster_table1=true_cluster_table, cluster_table2=method_cluster_table,
                           similarity_func=dice, similarity_func_name='Dice',
                           name_method1='true', name_method2=method, 
                           name_clus1='true cluster', name_clus2=name_clus2,
                           out_folder=out_folder)
    
            # create heatmap- szymkiewicz_simpson
            create_heatmap(cluster_table1=true_cluster_table, cluster_table2=method_cluster_table,
                           similarity_func=szymkiewicz_simpson, similarity_func_name='szymkiewicz_simpson',
                           name_method1='true', name_method2=method, 
                           name_clus1='true cluster', name_clus2=name_clus2,
                           out_folder=out_folder)
    
            # Fill in the best jaccard index for each true structure
            best_jaccard.loc[:, method] = list(similarity_df_jaccard.max(axis=1))
            '''
            # compute purity
            purity1, purity2 = purity(cluster_table1=true_cluster_table, cluster_table2=method_cluster_table)
            purity1_dict[method] = purity1
            purity2_dict[method] = purity2

            # compute our rand index
            rand_ind_dict[method] = ARI_Hullermeier2012(gene_list=gene_list,
                                                        cluster_table1=true_cluster_table,
                                                        cluster_table2=method_cluster_table)
            print("ari:", rand_ind_dict[method])

    # write dict to file
    json_object = json.dumps(rand_ind_dict, indent=4)
    with open("synthetic_rand_ind_dict.json", 'w') as file:
        file.write(json_object)

    json_object = json.dumps(purity1_dict, indent=4)
    with open("synthetic_purity1_dict.json", 'w') as file:
        file.write(json_object)

    json_object = json.dumps(purity2_dict, indent=4)
    with open("synthetic_purity2_dict.json", 'w') as file:
        file.write(json_object)

    #method_sizes['True clusters'] = list(true_cluster_table.loc[:, '#genes'])

