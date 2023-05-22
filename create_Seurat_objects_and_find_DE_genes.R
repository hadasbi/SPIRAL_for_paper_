library(Matrix)
library(tidyverse)
library(rliger)
library(dplyr)
library(patchwork)
library(Seurat)
library(tidyverse) 
library(scTenifoldNet)

setwd('D:/OneDrive - Technion/SPIRAL_for_paper')

# load counts
#data_n = 1    #Zhang2019
#data_n = 2    #Visium mouse brain
#data_n = 3    #Zebrafish Wagner 2018
###data_n = 4    #Bulk RNA-seq data of human differentiation- Mendel-Gutfreund lab
#data_n = 5    #Spatial transcriptomics data of a normal human prostate (FFPE) (10x Genomics- Visium)
#data_n = 6    #Bulk RNA-seq data of mouse B-cells that were treated with anti-IgM mAb- 25 time points in the 6 hours post stimulation (Chiang et al. 2020)
#data_n = 51-60    #Splatter- simulated datasets
#data_n = 20, 21, 24
#data_n = 61-80  

for (data_n in c(4)) {
  print(paste("data", as.character(data_n)))
  
  data_dir <- paste0('D:/OneDrive - Technion/SPIRAL/static/analysis/data', as.character(data_n))
  #data_dir <- paste0('../SPIRAL_webtool/static/analysis/data', as.character(data_n))
  
  seurat_object_file = paste0("./Seurat_objects/data", as.character(data_n), ".rds")
  out_folder_for_lists_of_DE_genes = paste0("./Seurat_DE_genes/data", as.character(data_n))
  if (!file.exists(out_folder_for_lists_of_DE_genes)){
    dir.create(out_folder_for_lists_of_DE_genes)
  }
  
  #######################################################################
  # construct a Seurat object (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html)
  
  if (data_n %in% c(1, 2, 3, 5, 7, 20, 21, 24, 
                    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                    71, 72, 73, 74, 75, 76, 77, 78, 79, 80)) {
    # Load the normalized and filtered dataset
    counts = read.table(paste0(data_dir, '/counts_norm_filt.csv'), row.names=1, header = TRUE)
  } else {
    # Load the raw dataset
    if (data_n==4) {
      counts = read.table(paste0(data_dir, '/counts.csv'), row.names=1, header = TRUE)
    } else {
      counts = read.table(paste0(data_dir, '/counts.txt'), row.names=1, header = TRUE)
    }
    counts = as.data.frame(lapply(counts, round))
  }
  
  # Initialize the Seurat object with the normalized and filtered data
  data <- CreateSeuratObject(counts = counts, 
                             project = paste0('data', as.character(data_n)))
  if (data_n %in% c(1, 2, 3, 5, 7, 20, 21, 24, 
                    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                    71, 72, 73, 74, 75, 76, 77, 78, 79, 80)) {
    data <- NormalizeData(data)
  } else {
    #data@scale.data <- as.matrix(log(counts + 1))
    data <- SetAssayData(object = data, slot = "scale.data", 
                         new.data = as.matrix(log(cpmNormalization(counts) + 1)))
  }
  
  # Identification of highly variable features (feature selection)
  data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  #top10 <- head(VariableFeatures(data), 10)
  
  # Scaling the data
  all.genes <- rownames(data)
  data <- ScaleData(data, features = all.genes)
  
  # Perform linear dimensional reduction
  data <- RunPCA(data, features = VariableFeatures(object = data), 
                 npcs=min(50, ncol(x = data)-1))
  
  #VizDimLoadings(data, dims = 1:2, reduction = "pca")
  #DimPlot(data, reduction = "pca")
  #DimHeatmap(data, dims = 1, cells = 500, balanced = TRUE)
  #DimHeatmap(data, dims = 1:15, cells = 500, balanced = TRUE)
  
  # Find the true dimensionality of the data
  #ElbowPlot(data)
  
  #if (data_n == 1) {
  #  dimensionality = 15
  #} else if (data_n == 2) {
  #  dimensionality = 15
  #} else if (data_n == 3) {
  #  dimensionality = 15
  #} else if (data_n == 5) {
  #  dimensionality = 15
  #} else if (data_n == 7) {
  #  dimensionality = 6
  #}
  
  dimensionality = 30
  
  # Cluster the cells
  data <- FindNeighbors(data, dims = 1:dimensionality)
  data <- FindClusters(data, resolution = 0.5)
  
  # UMAP
  data <- RunUMAP(data, dims = 1:dimensionality)
  DimPlot(data, reduction = "umap")
  
  # save the object to file
  saveRDS(data, file = seurat_object_file)
  
  ################################################################
  # find all markers of each cluster and save them to file
  if (data_n %in% c(1, 2, 3, 5, 7, 20, 21, 24, 
                    51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
                    61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 
                    71, 72, 73, 74, 75, 76, 77, 78, 79, 80)) {
    data.markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25)
  } else {
    data.markers <- FindAllMarkers(data, min.pct = 0.25, logfc.threshold = 0.25,
                                   test.use="DESeq2")
  }
  cluster_list = unique(data.markers$cluster)
  
  cluster_table = data.frame(matrix(nrow=0,ncol=6))
  colnames(cluster_table) <- c('cluster', 'cell_cluster', 'ncells', 
                               'up_or_down', 'ngenes', 'genes')
  
  ind = 1
  for (c in cluster_list) {
    d = data.markers[data.markers$cluster == c,]
    
    genes_up = d[d$avg_log2FC > 1,]$gene
    if (length(genes_up) != 0) {
      cluster_table[nrow(cluster_table) + 1,] = c(ind, c,
                                                  sum(data$seurat_clusters == c),
                                                  'up', length(genes_up),
                                                  toString(genes_up))
      ind = ind+1
    }
    
    genes_down = d[d$avg_log2FC < 1,]$gene
    if (length(genes_down) != 0) {
      cluster_table[nrow(cluster_table) + 1,] = c(ind, c,
                                                  sum(data$seurat_clusters == c),
                                                  'down', length(genes_down),
                                                  toString(genes_down))
      ind = ind+1
    }
  }
  
  colnames(cluster_table) <- c('cluster', 'cell_cluster', '#cells', 
                               'up_or_down', '#genes', 'genes')
  
  write.csv(cluster_table, 
            paste0(out_folder_for_lists_of_DE_genes,"/cluster_table.csv"), 
            row.names = FALSE)
}
