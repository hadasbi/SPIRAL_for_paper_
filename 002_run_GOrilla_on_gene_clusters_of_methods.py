#!/usr/bin/env python
# coding: utf-8

import os
import numpy as np
import pandas as pd
import warnings; warnings.simplefilter('ignore')
import pickle
import mechanicalsoup
from time import time, sleep
import re
from datetime import datetime
import json


def method_in_2_lines(method):
    return method.replace(' (', '\n(')


ANALYSIS_FOLDER = "../SPIRAL.web.tool/static/analysis/"
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


def get_GO_terms(gorilla_url, impute_method, species, background_list, geneslist,
                 pval=0.001):
    
    #p-value threshold is 0.001 since this is the highest allowed by GOrilla
    
    go_dict = dict()
    
    target_list = str(geneslist).replace("'", "").replace('[', '').replace(']', '').replace(', ', '\n')
    #print('target_list:', target_list)

    browser = mechanicalsoup.StatefulBrowser()
    browser.open(gorilla_url, verify=False)

    browser.select_form()

    browser["species"] = species
    # <option value="ARABIDOPSIS_THALIANA">Arabidopsis thaliana</option>
    # <option value="SACCHAROMYCES_CEREVISIAE">Saccharomyces cerevisiae</option>
    # <option value="CAENORHABDITIS_ELEGANS">Caenorhabditis elegans</option>
    # <option value="DROSOPHILA_MELANOGASTER">Drosophila melanogaster</option>
    # <option value="DANIO_RERIO">Danio rerio (Zebrafish)</option>
    # <option value="HOMO_SAPIENS" selected="selected">Homo sapiens</option>
    # <option value="MUS_MUSCULUS">Mus musculus</option>
    # <option value="RATTUS_NORVEGICUS">Rattus norvegicus</option>

    browser["run_mode"] = "hg"
    # "hg" for Two unranked lists of genes (target and background lists) 
    # "mhg" for Single ranked list of genes

    browser["target_set"] = target_list
    browser["background_set"] = background_list

    browser["db"] = "all"
    # "all" OR "proc" (process) OR "func" (function) OR "comp" (component)

    browser["pvalue_thresh"] = "0.001"
    browser["analysis_name"] = ""
    browser["user_email"] = ""

    browser["output_excel"] = True

    # Uncomment to launch a real web browser on the current page.
    # browser.launch_browser()

    # Uncomment to display a summary of the filled-in form
    # browser.form.print_summary()

    response = browser.submit_selected()
    # time.sleep(30)
    # browser.get_current_page()
    # print(response.text)
    new_link = browser.get_url()
    # print('new_link:', new_link)

    sleep(10)
    browser = mechanicalsoup.StatefulBrowser()
    browser.open(new_link, verify=False)
    new_link2 = browser.get_url()

    c = 0
    while ("GOResults" not in new_link2) and (c < 10):
        sleep(3)
        browser = mechanicalsoup.StatefulBrowser()
        browser.open(new_link, verify=False)
        new_link2 = browser.get_url()
        c += 1
    if ("GOResults" not in new_link2):
        print('PROBLEM: "GOResults" not in new_link2')

    # print('new_link2:', new_link2)

    ind = [m.start() for m in re.finditer('/', new_link2)][-1]
    for ontology_init, ontology_html in zip(['proc', 'func', 'comp'],
                                            ['GOResultsPROCESS.html', 'GOResultsFUNCTION.html',
                                             'GOResultsCOMPONENT.html']):
        new_link3 = new_link2[:ind + 1] + ontology_html
        print(new_link3)
        
        try:
            dfs = pd.read_html(new_link3)
            df = dfs[1]
            df.columns = df.iloc[0, :]
            df = df.drop([0])
            # print(df)
            df = df[df['P-value'].astype(float) <= pval]
            if len(df) > 0:
                go_dict[ontology_init] = [(a + ':' + b + ' (qval' + str(c) + ')') for a, b, c in
                     zip(df['GO term'], df['Description'], df['FDR q-value'])]
            else:
                go_dict[ontology_init] = []
        except ValueError as error:
            if str(error) == 'No tables found':
                go_dict[ontology_init] = []
            else:
                go_dict[ontology_init] = str(error)
    print('Gorilla_access_time:', datetime.now().strftime("%d/%m/%Y, %H:%M:%S"))

    return go_dict


# Build dataset_dict
dataset_dict = dict(zip([1, 3, 2, 5, 4, 6], 
                        ['lymphoblastoid cell line GM12891', 
                        'Zebrafish differentiation (7 time points)',
                        'sagittal-posterior section of a mouse brain',
                        'human prostate',
                        'human differentiation (days 0 and 7 − 21)',
                        'mouse B cells, treated with anti-IgM mAb (25 time points in the 6 hours post stimulation)'
                       ]))
dataset_type = dict(zip([1, 3, 2, 5, 4, 6], 
                        ['single cell', 'single cell', 
                        'spatial', 'spatial',
                        'bulk', 'bulk'
                       ]))

dataset_dict[20] = 'Slide-seqV2- diabetic kidney mouse'
dataset_type[20] = 'spatial'

dataset_dict[21] = 'Slide-seqV2- normal kidney mouse'
dataset_type[21] = 'spatial'

for data_n in range(61, 81):
    f = open("../SPIRAL.web.tool/static/analysis/data" + str(data_n) + "/dataset_name.txt", 'r')
    name = f.readlines()[0]
    #print(name)
    f.close()
    
    dataset_dict[data_n] = 'mouse ' + name.replace('__', ' ')
    dataset_type[data_n] = 'spatial'


# Build methods_dict
methods_dict = {'single cell': ['SPIRAL (Jaccard index=0.05)',
                                'SPIRAL (Jaccard index=0.5)',
                                'SPIRAL (Jaccard index=0.75)',
                                'Seurat',
                                'Slingshot+tradeSeq (nPoints=20)',
                                'Slingshot+tradeSeq (nPoints=200)',
                                "Hotspot (n_neighbors=30)",
                                "Hotspot (n_neighbors=300)",
                                'nsNMF',
                                'SingleCellHaystack (k=3)',
                                'SingleCellHaystack (k=5)'],
               'spatial':  ['SPIRAL (Jaccard index=0.05)',
                            'SPIRAL (Jaccard index=0.5)',
                            'SPIRAL (Jaccard index=0.75)',
                            'Seurat',
                            "Hotspot (n_neighbors=30)",
                            "Hotspot (n_neighbors=300)",
                            'nsNMF',
                            'SpatialDE (C=3)',
                            'SpatialDE (C=5)',
                            'SingleCellHaystack (k=3)',
                            'SingleCellHaystack (k=5)'],
               'bulk': ['SPIRAL (Jaccard index=0.05)',
                        'SPIRAL (Jaccard index=0.5)',
                        'SPIRAL (Jaccard index=0.75)',
                        'Seurat',
                        'nsNMF']}


out_folder = "../SPIRAL_for_paper_/comparison_of_all_methods_real_datasets/GO_terms/"
if not os.path.exists(out_folder):
    os.mkdir(out_folder)


# Find and save list of GO terms for every dataset and method
gorilla_url = 'http://cbl-gorilla.cs.technion.ac.il/'

for data_n in dataset_dict.keys():
    dataset = dataset_dict[data_n]
    typee = dataset_type[data_n]
    methods = methods_dict[typee]
    
    for method in methods:
        print(dataset, typee, method)
        outfile = os.path.join(out_folder, dataset + ' - ' + method + ' - GO_terms.json')
        if not os.path.exists(outfile):
            try:
                data_path = "../SPIRAL.web.tool/static/analysis/data" + str(data_n)
                species = open(os.path.join(data_path, 'species.txt'), "r").read()
                impute_method = open(os.path.join(data_path, 'imputation_method.txt'), "r").read()

                #get genelist for this dataset
                genelist_file = os.path.join(data_path, impute_method + '_genetable.p')
                # Load genes_table
                with open(genelist_file, 'rb') as fp:
                    genes_table = pickle.load(fp)
                background_list = str(list(genes_table.index)).replace("'", "").replace('[', '').replace(']', '').replace(', ','\n')

                flag = False
                if method == 'nsNMF':
                    method_folder = 'nsNMF/data' + str(data_n) + '/'
                    name_clus2 = 'module'
                    for file in os.listdir(method_folder):
                        if file.startswith("cluster_table"):
                            #print(file)
                            method_cluster_table_file = os.path.join(method_folder, file)
                            flag = True
                elif method == 'Seurat':
                    method_folder = 'Seurat_DE_genes/data' + str(data_n) + '/'
                    name_clus2 = 'DE gene set'
                    for file in os.listdir(method_folder):
                        if file.startswith("cluster_table"):
                            #print(file)
                            method_cluster_table_file = os.path.join(method_folder, file)
                            flag = True
                elif method in ['Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)']:
                    method_folder = 'Slingshot/data' + str(data_n) + '/'
                    name_clus2 = 'gene cluster'
                    if method == 'Slingshot+tradeSeq (nPoints=20)':
                        file = 'cluster_table_20.csv'
                        method_cluster_table_file = os.path.join(method_folder, file)
                        flag = True
                    elif method == 'Slingshot+tradeSeq (nPoints=200)':
                        file = 'cluster_table_200.csv'
                        method_cluster_table_file = os.path.join(method_folder, file)
                        flag = True
                elif 'Hotspot' in method:
                    method_folder = 'Hotspot/data' + str(data_n) + '/'
                    name_clus2 = 'gene module'
                    if method == "Hotspot (n_neighbors=30)":
                        file = 'cluster_table_danb_30.csv'     
                    elif method == "Hotspot (n_neighbors=300)":
                        file = 'cluster_table_danb_300.csv'
                    method_cluster_table_file = os.path.join(method_folder, file)
                    flag = True
                elif 'SpatialDE' in method:
                    method_folder = 'SpatialDE/data' + str(data_n) + '/'
                    name_clus2 = 'pattern'
                    if method=='SpatialDE (C=3)':
                        file = 'cluster_table_3.csv'
                    elif method=='SpatialDE (C=5)':
                        file = 'cluster_table_5.csv'
                    method_cluster_table_file = os.path.join(method_folder, file)
                    flag = True
                elif 'SingleCellHaystack' in method:
                    method_folder = 'singlecellhaystack/data' + str(data_n) + '/'
                    name_clus2 = 'cluster'
                    if method=='SingleCellHaystack (k=3)':
                        file = 'data' + str(data_n) + '.singleCellHaystack.k_3.csv'
                    elif method=='SingleCellHaystack (k=5)':
                        file = 'data' + str(data_n) + '.singleCellHaystack.k_5.csv'
                    method_cluster_table_file = os.path.join(method_folder, file)
                    flag = True


                if flag:
                    method_cluster_table = pd.read_csv(method_cluster_table_file, index_col=0)

                    if method in ['Seurat', 'Slingshot+tradeSeq (nPoints=20)', 'Slingshot+tradeSeq (nPoints=200)']:
                        method_cluster_table.loc[:, 'genes'] = [
                            str(genes.replace("['", "").replace("']", "").split(", ")) for genes in method_cluster_table.loc[:, 'genes']]
                    elif 'SingleCellHaystack' in method:
                        method_cluster_table = pd.read_csv(method_cluster_table_file)
                        method_cluster_table.columns = ['#genes', 'genes']
                        method_cluster_table.loc[:, 'genes'] = [
                            str(a.split(' ')) for a in method_cluster_table.loc[:, 'genes']]
                    print('\n\n', method, method_cluster_table_file)

                elif method in ['SPIRAL (Jaccard index=0.05)', 
                                'SPIRAL (Jaccard index=0.5)', 
                                'SPIRAL (Jaccard index=0.75)']:
                    flag = True

                    struct_table_file = os.path.join(data_path, impute_method + '_sigtable_filt_GO_vis.xlsx')

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

                    method_cluster_table = struct_table.loc[struct_list, ['num_genes_in_struct', 'genes']].rename(columns={'num_genes_in_struct': '#genes'})
                    method_cluster_table.loc[:, 'genes'] = [str(a.split(',')) for a in method_cluster_table.loc[:, 'genes']]

                if not flag:
                    print('\n\n', dataset, '+', method, flag, '!!!!!!!!!!!\n\n')
                if flag:
                    print(dataset, '+', method, '\n', method_cluster_table.head())

                    #module_GO_terms[dataset][method] = dict()
                    GO_terms = dict()

                    cluster_table = method_cluster_table
                    for clus_ind, clus in enumerate(cluster_table.index):
                        # read cluster gene list
                        if method == 'nsNMF' and data_n in range(20, 81):
                            geneslist = cluster_table.loc[clus, 'genes'].split(", ")
                        else:
                            geneslist = cluster_table.loc[clus, 'genes'].replace("['", "").replace("']", "").split("', '")
                        print(clus, geneslist)
                        GO_terms[clus] = get_GO_terms(
                            gorilla_url, impute_method, species, background_list, geneslist)

                    # write dict to file
                    json_object = json.dumps(GO_terms, indent=4)
                    with open(outfile, 'w') as file:
                         file.write(json_object)
            except:
                pass

