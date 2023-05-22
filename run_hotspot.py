#!/usr/bin/env python
# coding: utf-8

# https://hotspot.readthedocs.io/en/latest/CD4_Tutorial.html
# https://yoseflab.github.io/Hotspot/Spatial_Tutorial.html


import os
import sys
import hotspot
import anndata
import numpy as np
import pandas as pd
import hotspot
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns
import mplscience
from scipy.io import mmread
from scipy.sparse import csr_matrix, csc_matrix
import scanpy as sc
import warnings; warnings.simplefilter('ignore')

    
for data_n in [3]:
    try:
        print('\n\n\nData', data_n)
        folder = '../SPIRAL_webtool/static/analysis/data' + str(data_n)

        out_folder = 'Hotspot/data' + str(data_n) + '/'
        if not os.path.isdir('./Hotspot'):
            os.mkdir('./Hotspot')
        if not os.path.isdir(out_folder):   
            os.mkdir(out_folder)

        # load normalized and filtered file
        counts_file = os.path.join(folder, 'counts_norm_filt.csv')
        counts = pd.read_csv(counts_file, index_col=0, sep='\t')
        counts = counts.loc[counts.sum(axis=1) != 0]

        # Filter genes
        gene_counts = (counts > 0).sum(axis=1)
        valid_genes = gene_counts >= 50
        counts = counts.loc[valid_genes]

        if data_n in [2, 5, 20, 21, 24]+list(range(61, 81)): #spatial
            pos_file = os.path.join(folder, 'spatial_coors.csv')
            if data_n==5:
                pos = pd.read_csv(pos_file, index_col=0, sep='\t', header=None)
            else:
                pos = pd.read_csv(pos_file, index_col=0)
            pos.columns = ['X', 'Y']
            pos.index.name = None
            pos = pos.loc[list(counts), :]

        #######################################################################################################
        for n_neighbors in [30, 300]:
            model = 'danb'  # model='bernoulli' does not work
            if data_n in [1, 3, 7] + list(range(51, 61)): #single cell
                adata = anndata.AnnData(counts.T)
            elif data_n in [2, 5] + [20, 21, 24] + list(range(61, 81)): #spatial
                adata = anndata.AnnData(counts.T, obsm={"spatial": pos})
            adata.obs_names = list(counts)
            adata.obs = pd.DataFrame(index=list(counts), data=counts.sum(), columns=['total_counts'])
            adata.var_names = list(counts.index)

            adata.layers["counts"] = adata.X.copy()
            sc.pp.normalize_total(adata)
            sc.pp.log1p(adata)
            adata.layers["log_normalized"] = adata.X.copy()
            sc.pp.scale(adata)

            if data_n in [1, 3, 7] + list(range(51, 61)): #single cell
                sc.tl.pca(adata)

                with mplscience.style_context():
                    sc.pl.pca_variance_ratio(adata)
                    plt.show()

                # rerun with fewer components
                sc.tl.pca(adata, n_comps=10)

            # Create the Hotspot object and the neighborhood graph
            # hotspot works a lot faster with a csc matrix!
            adata.layers["counts_csc"] = csc_matrix(adata.layers["counts"])
            if data_n in [1, 3, 7] + list(range(51, 61)): #single cell
                hs = hotspot.Hotspot(
                    adata,
                    layer_key="counts_csc",
                    model=model,
                    latent_obsm_key="X_pca",
                    umi_counts_obs_key="total_counts"
                )
            elif data_n in [2, 5] + [20, 21, 24] + list(range(61, 81)): #spatial
                hs = hotspot.Hotspot(
                    adata,
                    layer_key="counts_csc",
                    model=model,
                    latent_obsm_key="spatial",
                    umi_counts_obs_key="total_counts"
                )

            hs.create_knn_graph(
                weighted_graph=False, n_neighbors=n_neighbors,
            )

            hs_results = hs.compute_autocorrelations(jobs=15)
            print(hs_results.head(15))

            # Select the genes with significant lineage autocorrelation
            hs_genes = hs_results.index[hs_results.FDR < 0.05]

            # Compute pair-wise local correlations between these genes
            lcz = hs.compute_local_correlations(hs_genes, jobs=15)

            modules = hs.create_modules(min_gene_threshold=15, core_only=True, fdr_threshold=0.05)
            print(modules.value_counts())

            module_list = list(set(modules[modules != -1]))
            method_cluster_table = pd.DataFrame(columns=['#genes', 'genes'], index=module_list)
            for m in module_list:
                genes = list(modules[modules == m].index)
                method_cluster_table.loc[m, "#genes"] = len(genes)
                method_cluster_table.loc[m, "genes"] = str(genes)
            method_cluster_table.to_csv(os.path.join(out_folder, 'cluster_table_' + model + '_' + str(n_neighbors) + '.csv'))

            hs.plot_local_correlations(vmin=-12, vmax=12)
        
    except Exception as e:
        print('\n\n\n!!!! Data', data_n, 'did not work!!!!!!!')
        print(e)
        print('\n')
            

