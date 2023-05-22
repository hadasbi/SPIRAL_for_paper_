#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import seaborn as sns; sns.set_theme()
import matplotlib.pyplot as plt
import pickle5 as pickle
import multiprocessing as mp
import json
import time
#from funcs import *

def method_in_2_lines(method):
    return method.replace(' (', '\n(')


def jaccard(set1, set2):
    return len(set1.intersection(set2))/len(set1.union(set2))


def compute_sum_of_diff(clustering1, clustering2, gene_list, ind1, ind2):
    d = 0
    nclus1 = clustering1.shape[1]
    nclus2 = clustering2.shape[1]

    # add score for every gene pair
    for ig1, g1 in enumerate(gene_list[ind1:ind2]):
        if (ig1 + 1) % 100 == 0:
            print(ig1 + 1, '/', ind2 - ind1)
        for ig2, g2 in enumerate(gene_list[ind1+ig1+1:]):
            #print(ind1+ig1, ind1+ig1+1+ig2)
            E1 = 1 - np.linalg.norm(clustering1[ind1+ig1]-clustering1[ind1+ig1+1+ig2], ord=1)/nclus1
            E2 = 1 - np.linalg.norm(clustering2[ind1+ig1]-clustering2[ind1+ig1+1+ig2], ord=1)/nclus2
            #print(E1, '\t', E2)
            d += np.abs(E1-E2)
    return d


#ANALYSIS_FOLDER = "../SPIRAL_webtool/static/analysis/"
ANALYSIS_FOLDER = "../SPIRAL_webtool/static/analysis/"
def filter_struct_lst_for_result_panel(data_n, new_Jaccard_thr_genes):
    # produces a list of structures, in which there is no pair of structures with Jaccard index larger than
    # Jaccard_thr_genes. The function relies on the fact that the structure numbers are
    # range(1, len(Jaccard_mat_genes) + 1), and that structure i corresponds to row i-1 in Jaccard_mat_genes.
    print('filter_struct_lst_for_result_panel!!!')

    data_path = os.path.join(ANALYSIS_FOLDER, 'data' + str(data_n))
    # get the Jaccard matrix of gene lists
    Jaccard_mat_genes = np.load(os.path.join(data_path, 'Jaccard_mat_genes.npy'))

    inds = list(range(1, len(Jaccard_mat_genes) + 1))
    thr_struct_lst = []

    while inds:
        i = inds.pop(0)
        thr_struct_lst.append(i)
        similar_to_i = []

        for j in inds:
            if Jaccard_mat_genes[i - 1][j - 1] >= new_Jaccard_thr_genes:
                similar_to_i.append(j)

        inds = [j for j in inds if j not in similar_to_i]

    return thr_struct_lst


def cluster_table_to_matrix(gene_list, cluster_table):
    ngenes = len(gene_list)
    
    # create gene X cluster table
    nclus = len(cluster_table)
    clustering = pd.DataFrame(index=gene_list, columns=cluster_table.index,
                               data=np.zeros((ngenes, nclus)))
    for clus in cluster_table.index:
        # read cluster gene list
        genes = cluster_table.loc[clus, 'genes'].replace("['", "").replace("']", "").split("', '")
        clustering.loc[genes, clus] = 1
        
    return clustering


def ARI_Hullermeier2012(gene_list, clustering1, clustering2):   
    # clustering1, clustering2 are clustering matrices (the output of cluster_table_to_matrix)
    ngenes = len(gene_list)
    nclus1 = clustering1.shape[1]
    nclus2 = clustering2.shape[1]

    # divide to small tasks and run in parallel
    clustering1 = np.array(clustering1)
    clustering2 = np.array(clustering2)
    tic = time.perf_counter()
    with mp.Pool(mp.cpu_count()) as pool:
        d_list = pool.starmap(compute_sum_of_diff,
                              [(clustering1, clustering2, gene_list, inds[0], inds[1]) for inds in
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
    print(toc - tic, "seconds")
    return 1 - sum(d_list) / (ngenes * (ngenes - 1) / 2)


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
                   name_method1, name_method2, name_clus1, name_clus2, out_folder, save_to_file=True):
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
    
    if save_to_file:
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
    '''
    # # Compute the true cluster table for each of the synthetic datasets
    for data_n in range(51, 61):
        print(data_n)
        folder = '../SPIRAL_simulated_synthetic_data/data' + str(data_n)
        #counts_file = os.path.join(folder, 'splatter_branching_path.txt')
        coldata_file = os.path.join(folder, 'splatter_branching_path_colData.txt')
        rowdata_file = os.path.join(folder, 'splatter_branching_path_rowData.txt')
        #data = pd.read_csv(counts_file, index_col=0)
        coldata = pd.read_csv(coldata_file, sep="\t", index_col=0)
        rowdata = pd.read_csv(rowdata_file, sep="\t", index_col=0)
        low_thr = 0.9
        high_thr = 1.1
        for path in [1, 2, 3]:
            rowdata.loc[:, 'Path' + str(path) + '_up'] = (rowdata.loc[:, 'DEFacPath' + str(path)] >= high_thr)
            rowdata.loc[:, 'Path' + str(path) + '_down'] = (rowdata.loc[:, 'DEFacPath' + str(path)] <= low_thr)
        true_cluster_table = pd.DataFrame(index=range(6), columns=['#genes', 'genes'])
        true_cluster_table.index.name = 'cluster'
        for path in [1, 2, 3]:
            path_up_genes = list(rowdata[rowdata.loc[:, 'Path' + str(path) + '_up']].index)
            true_cluster_table.loc[2*(path - 1), '#genes'] = len(path_up_genes)
            true_cluster_table.loc[2*(path - 1), 'genes'] = str(path_up_genes)
    
            path_down_genes = list(rowdata[rowdata.loc[:, 'Path' + str(path) + '_down']].index)
            true_cluster_table.loc[2*(path - 1) + 1, '#genes'] = len(path_down_genes)
            true_cluster_table.loc[2*(path - 1) + 1, 'genes'] = str(path_down_genes)
        true_cluster_table.to_csv(os.path.join(folder, 'true_cluster_table.csv'))
    '''

    # # Compute the adjusted rand index for each method, for each dataset
    methods = [
                'SPIRAL (Jaccard index=0.05)',
                'SPIRAL (Jaccard index=0.5)',
                'SPIRAL (Jaccard index=0.75)',
                'Seurat',
                'Slingshot+tradeSeq (nPoints=20)',
                'Slingshot+tradeSeq (nPoints=200)',
                "Hotspot (n_neighbors=30)",
                "Hotspot (n_neighbors=300)",
                'nsNMF',
                'SingleCellHaystack (k=3)',
                'SingleCellHaystack (k=5)',
                'SingleCellHaystack (k=6)'
    ]

    method_sizes = pd.DataFrame(index=range(51, 61), columns=['True clusters'] + methods)
    rand_ind = pd.DataFrame(index=range(51, 61), columns=methods)
    best_jaccard_51 = pd.DataFrame(index=range(6))
    jaccard_list = pd.DataFrame(columns=['method', 'best_Jaccard_index_for_cluster'])

    for data_n in range(51, 61):
        out_folder = "../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset/data" + str(data_n)
        if not os.path.exists("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset"):
            os.mkdir("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset")
        if not os.path.exists(out_folder):
            os.mkdir(out_folder)

        # get gene list
        folder = '../SPIRAL_simulated_synthetic_data/data' + str(data_n)
        rowdata_file = os.path.join(folder, 'splatter_branching_path_rowData.txt')
        rowdata = pd.read_csv(rowdata_file, sep="\t", index_col=0)
        gene_list = list(rowdata.index)

        # load true cluster table
        folder = '../SPIRAL_simulated_synthetic_data/data' + str(data_n)
        true_cluster_table = pd.read_csv(os.path.join(folder, 'true_cluster_table.csv'), index_col=0)
        method_sizes.loc[data_n, 'True clusters'] = str(list(true_cluster_table.loc[:, '#genes']))

        # compute true clustering matrix
        true_clustering = cluster_table_to_matrix(gene_list, true_cluster_table)

        for method in methods:
            flag = False
            # load method's cluster table
            if method == 'nsNMF':
                method_folder = '../SPIRAL_for_paper/nsNMF/data' + str(data_n) + '/'
                name_clus2 = 'module'
                for file in os.listdir(method_folder):
                    if file.startswith("cluster_table"):
                        #print(file)
                        method_cluster_table_file = os.path.join(method_folder, file)
                        flag = True
                method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)
                if data_n != 51:
                    method_cluster_table.loc[:, 'genes'] = [
                        str(a.split(', ')) for a in method_cluster_table.loc[:, 'genes']]
            elif method == 'Seurat':
                method_folder = '../SPIRAL_for_paper/Seurat_DE_genes/data' + str(data_n) + '/'
                name_clus2 = 'DE gene set'
                file = "cluster_table.csv"
                method_cluster_table_file = os.path.join(method_folder, file)
                flag = True
                method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)
                method_cluster_table.loc[:, 'genes'] = [
                        str(genes.replace("['", "").replace("']", "").split(", ")) for genes in method_cluster_table.loc[:, 'genes']]
            elif method in ['Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)']:
                method_folder = '../SPIRAL_for_paper/Slingshot/data' + str(data_n) + '/'
                name_clus2 = 'gene cluster'
                if method == 'Slingshot+tradeSeq (nPoints=20)':
                    file = 'cluster_table_20.csv'
                    method_cluster_table_file = os.path.join(method_folder, file)
                    flag = True
                elif method == 'Slingshot+tradeSeq (nPoints=200)':
                    file = 'cluster_table_200.csv'
                    method_cluster_table_file = os.path.join(method_folder, file)
                    flag = True
                method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)
                method_cluster_table.loc[:, 'genes'] = [
                        str(genes.replace("['", "").replace("']", "").split(", ")) for genes in method_cluster_table.loc[:, 'genes']]
            elif 'Hotspot' in method:
                method_folder = '../SPIRAL_for_paper/Hotspot/data' + str(data_n) + '/'
                name_clus2 = 'gene module'
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
                method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)
            elif 'SingleCellHaystack' in method:
                method_folder = '../SPIRAL_for_paper/singlecellhaystack/data' + str(data_n) + '/'
                name_clus2 = 'cluster'
                if method=='SingleCellHaystack (k=3)':
                    file = 'data' + str(data_n) + '.singleCellHaystack.k_3.csv'
                elif method=='SingleCellHaystack (k=5)':
                    file = 'data' + str(data_n) + '.singleCellHaystack.k_5.csv'
                elif method=='SingleCellHaystack (k=6)':
                    file = 'data' + str(data_n) + '.singleCellHaystack.k_6.csv'
                method_cluster_table_file = os.path.join(method_folder, file)
                flag = True
                method_cluster_table = pd.read_csv(method_cluster_table_file)
                method_cluster_table.columns = ['#genes', 'genes']
                method_cluster_table.loc[:, 'genes'] = [
                    str(a.split(' ')) for a in method_cluster_table.loc[:, 'genes']]
            elif method in ['SPIRAL (Jaccard index=0.05)',
                            'SPIRAL (Jaccard index=0.5)',
                            'SPIRAL (Jaccard index=0.75)']:
                flag = True
                #folder = "../SPIRAL_webtool/static/analysis/data" + str(data_n)
                #impute_method = 'agg_wald'
                name_clus2 = 'structure'

                # load SPIRAL cluster table
                folder = "../SPIRAL_webtool/static/analysis/data" + str(data_n)
                #species = open(os.path.join(folder, 'species.txt'), "r").read()
                impute_method = open(os.path.join(folder, 'imputation_method.txt'), "r").read()

                #get genelist for this dataset
                '''
                genelist_file = os.path.join(data_path, impute_method + '_genetable.p')
                # Load genes_table
                with open(genelist_file, 'rb') as fp:
                    genes_table = pickle.load(fp)
                gene_list = list(genes_table.index)
                '''
                #raw_data = pd.read_csv(os.path.join(data_path, 'counts.txt'), index_col=0)
                #gene_list = list(raw_data.index)

                struct_table_file = os.path.join(folder, impute_method + '_sigtable_filt_GO_vis.xlsx')
                #genetable_file = os.path.join(folder, impute_method + '_genetable.p')
                #with open(genetable_file, 'rb') as fp:
                #    genes_table = pickle.load(fp)

                struct_table = pd.read_excel(struct_table_file, index_col=0)

                if method=='SPIRAL (Jaccard index=0.05)':
                    struct_list = filter_struct_lst_for_result_panel(data_n=data_n,
                                                                     new_Jaccard_thr_genes=0.05)
                    struct_list_str = str(struct_list)

                elif method=='SPIRAL (Jaccard index=0.5)':
                    struct_list = filter_struct_lst_for_result_panel(data_n=data_n,
                                                                     new_Jaccard_thr_genes=0.5)
                    struct_list_str = str(struct_list)

                elif method=='SPIRAL (Jaccard index=0.75)':
                    struct_list = filter_struct_lst_for_result_panel(data_n=data_n,
                                                                     new_Jaccard_thr_genes=0.75)
                    struct_list_str = str(struct_list)

                method_cluster_table = struct_table.loc[struct_list, ['num_genes_in_struct', 'genes']].rename(
                    columns={'num_genes_in_struct': '#genes'})
                method_cluster_table.loc[:, 'genes'] = [str(a.split(',')) for a in method_cluster_table.loc[:, 'genes']]

            # compute indices
            if flag:
                print(data_n, method, '\n', method_cluster_table.head())

                method_sizes.loc[data_n, method] = str(list(method_cluster_table.loc[:, '#genes']))

                # Fill in the best jaccard index for each true structure
                # create heatmap- Jaccard
                similarity_df_jaccard = create_heatmap(cluster_table1=true_cluster_table,
                                                       cluster_table2=method_cluster_table,
                                                       similarity_func=jaccard, similarity_func_name='jaccard',
                                                       name_method1='true', name_method2=method,
                                                       name_clus1='true cluster', name_clus2=name_clus2,
                                                       out_folder=out_folder, save_to_file=True)

                if data_n==51:
                    best_jaccard_51.loc[:, method] = list(similarity_df_jaccard.max(axis=1))
                df = pd.DataFrame(columns=['method', 'best_Jaccard_index_for_cluster'])
                df.loc[:, 'best_Jaccard_index_for_cluster'] = list(similarity_df_jaccard.max(axis=1))
                df.loc[:, 'method'] = method
                jaccard_list = pd.concat([jaccard_list, df], ignore_index=True)

                '''    
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
                '''


                # compute the adjusted rand index (ARI)
                method_clustering = cluster_table_to_matrix(gene_list, method_cluster_table)
                rand_ind.loc[data_n, method] = ARI_Hullermeier2012(gene_list=gene_list, 
                                                                   clustering1=true_clustering,
                                                                   clustering2=method_clustering)

    rand_ind.to_csv("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset/rand_ind.csv")
    method_sizes.to_csv("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset/method_sizes.csv")
    jaccard_list.to_csv("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset/jaccard_list.csv")

    rand_ind2 = pd.DataFrame(columns=['ARI', 'method'])
    for method in methods:
        df = pd.DataFrame(columns=['ARI', 'method'])
        df.loc[:, 'ARI'] = rand_ind.loc[:, method]
        df.loc[:, 'method'] = method
        rand_ind2 = pd.concat([rand_ind2, df], ignore_index=True)
    rand_ind2.to_csv("../SPIRAL_for_paper/comparison_of_all_methods_to_ground_truth_of_Splatter_dataset/rand_ind2.csv")

