library(Matrix)
library(tidyverse)
#library(rliger)
library(dplyr)
#library(patchwork)
library(Seurat)
library(tidyverse) 
library(SingleCellExperiment)
library(slingshot, quietly = FALSE)
library(tradeSeq)
library(clusterExperiment)

setwd('D:/OneDrive - Technion/SPIRAL_for_paper')

# load counts
#data_n = 1    #Zhang2019
######data_n = 2    #Visium mouse brain
#data_n = 3    #Zebrafish Wagner 2018
######data_n = 4    #Bulk RNA-seq data of human differentiation- Mendel-Gutfreund lab
######data_n = 5    #Spatial transcriptomics data of a normal human prostate (FFPE) (10x Genomics- Visium)
######data_n = 6    #Bulk RNA-seq data of mouse B-cells that were treated with anti-IgM mAb- 25 time points in the 6 hours post stimulation (Chiang et al. 2020)
#data_n = 7    #Splatter- simulated dataset

#######################################################################
# Tutorial: https://nbisweden.github.io/workshop-archive/workshop-scRNAseq/2020-01-27/labs/compiled/slingshot/slingshot.html

# Define a color pallete to use
pal <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))

suppressPackageStartupMessages({
  library(Seurat)
  library(cowplot)
})

for (data_n in c(56, 57, 58, 59, 60)) {
  print(paste("data", as.character(data_n)))
  
  data_dir <- paste0('D:/OneDrive - Technion/SPIRAL/static/analysis/data', as.character(data_n))
  #data_dir <- paste0('../SPIRAL/static/analysis/data', as.character(data_n))
  
  #seurat_object_file = paste0("./Seurat_objects/data", as.character(data_n), ".rds")
  out_folder_for_gene_modules = paste0("./comparison_to_Slingshot/data", as.character(data_n))
  if (!file.exists(out_folder_for_gene_modules)){
    dir.create(out_folder_for_gene_modules)
  }
  
  counts = read.table(paste0(data_dir, '/counts_norm_filt.csv'), row.names=1, header = TRUE)
  data <- CreateSeuratObject(counts = counts, 
                             project = paste0('data', as.character(data_n)))
  data <- NormalizeData(data)
  data <- FindVariableFeatures(data, nfeatures = 2000)
  data <- ScaleData(data)
  data <- RunPCA(data)
  
  data <- FindNeighbors(data)
  data <- FindClusters(data, resolution = 1)
  
  data <- RunUMAP(data, n.neighbors = 10, dims = 1:50, spread = 2, min.dist = 0.3)
  
  # Plot the clusters
  DimPlot(data, group.by = "RNA_snn_res.1")
  
  # Save the objects as separate matrices for input in slingshot
  dimred <- data@reductions$umap@cell.embeddings
  clustering <- data$RNA_snn_res.1
  counts <- as.matrix(data@assays$RNA@counts[data@assays$RNA@var.features, ])
  
  suppressPackageStartupMessages({
    library(slingshot)
  })
  
  # Run default Slingshot lineage identification
  set.seed(1)
  lineages <- getLineages(data = dimred, clusterLabels = clustering)
  
  lineages
  
  # Plot the lineages
  par(mfrow = c(1, 2))
  plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
  for (i in levels(clustering)) {
    text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
  }
  plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
  #lines(lineages, lwd = 3, col = "black")
  
  curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
  curves
  
  plot(dimred, col = pal[clustering], asp = 1, pch = 16)
  #lines(curves, lwd = 3, col = "black")
  
  sce <- fitGAM(counts = as.matrix(counts), sds = curves)
  
  plotGeneCount(curves, counts, clusters = clustering, models = sce)
  
  # Define function to plot
  plot_differential_expression <- function(feature_id) {
    feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
    cowplot::plot_grid(plotGeneCount(curves, counts, gene = feature_id[1], clusters = clustering, models = sce) + ggplot2::theme(legend.position = "none"), 
                       plotSmoothers(sce, as.matrix(counts), gene = feature_id[1]))
  }
  
  patternTest_association <- patternTest(sce)
  patternTest_association$feature_id <- rownames(patternTest_association)
  gene_nona <- rownames(patternTest_association[!is.na(patternTest_association$df), ])
  
  
  # Tutorial: https://statomics.github.io/tradeSeq/articles/tradeSeq.html
  for (nPointsClus in c(20, 200)) {
    clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                         genes = gene_nona)
    clusterLabels <- primaryCluster(clusPat$rsec)
    
    cUniq <- unique(clusterLabels)
    cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
    
    cluster_table = data.frame(matrix(nrow=0,ncol=3))
    colnames(cluster_table) <- c('cluster', 'ngenes', 'genes')
    
    for (clus in cUniq) {
      print(sum(clusterLabels == clus))
      genes = gene_nona[clusterLabels == clus]
      if (length(genes)>=20) {
        cluster_table[nrow(cluster_table) + 1,] = c(clus,
                                                    length(genes),
                                                    toString(genes))
      }
    }
    
    colnames(cluster_table) <- c('cluster', '#genes', 'genes')
    
    write.csv(cluster_table, 
              paste0(out_folder_for_gene_modules,"/cluster_table_", 
                     toString(nPointsClus), ".csv"), 
              row.names = FALSE)
  }
}

