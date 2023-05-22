

library(singleCellHaystack)
library(ggplot2)
library(Seurat)

## Application to single-cell RNA-seq data
 #Preparing  data

  #This data should include the following objects:
  
  #dat.expression: a matrix object with the expression of genes (rows) in each cell (columns).
  #dat.pca: the output of PCA. This data.frame contains the first 20 pricipal components (PCs).
  #dat.tsne: a data.frame with t-SNE coordinates (2D) based on the first 20 PCs.


############################################################################################################################
# create data

dat.expression <- read.csv("counts.csv")

####################
  
  ## using Seurat to creat PCA and tsney
  seurat.data <- CreateSeuratObject(counts = dat.expression, project = data_name, min.cells = 3, min.features = 200)
  seurat.data <- NormalizeData(seurat.data, normalization.method = "LogNormalize", scale.factor = 10000)
  
  #Identification of highly variable genes 
  seurat.data <- FindVariableFeatures(seurat.data, selection.method = "vst", nfeatures = 2000)
  
  #Scaling the data
  all.genes <- rownames(dat.expression)
  seurat.data <- ScaleData(seurat.data, features = all.genes)
  
  #Perform linear dimensional reduction PCA
  seurat.data <- RunPCA(seurat.data, features = VariableFeatures(object = seurat.data),npcs = 20)
  dat.pca <- Embeddings(seurat.data, reduction = "pca")[, 1:20] 
  
  #Running haystack on multi-dimensional coordinates 
  set.seed(123)
  res.pc20 <- haystack(x = dat.pca, expression = dat.expression)
  
  # get the 1000 most significant DEGs, and cluster them by their distribution pattern in the 2D plot
  sorted.table <- show_result_haystack(res.haystack = res.pc20, n = 1000)
  gene.subset <- row.names(sorted.table)
  
  # cluster the genes by their expression pattern in the input space using hierarchical clustering
  res.hc <- hclust_haystack(dat.pca, dat.expression[gene.subset, ], grid.coordinates=res.pc20$info$grid.coordinates)
  
  ### cut into clusters using the cutree and output to table
  
  number.of.clustrs <- c(3,5,6)
  
  for (n in 1:length(number.of.clustrs)){
    res.hc.clusters <- cutree(res.hc, k=number.of.clustrs[n])
  
    clust.table <- data.frame(clust<-integer(),genes=character())
    clust.table <- matrix(nrow=0,ncol=2)
    
    tt <- table(res.hc.clusters)
    for (c in unique(res.hc.clusters)) {
      cluster.names <- paste(names(which(res.hc.clusters == c)),collapse = ' ')
      clust.table <- rbind(clust.table,c(tt[c],cluster.names))
      
    }
    
    colnames(clust.table)<- c("genes_number","genes_names")
    
    write.csv(clust.table, file ="output_file.csv",row.names=F,quote=F)
  
  }






