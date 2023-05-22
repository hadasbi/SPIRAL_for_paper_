

library(singleCellHaystack)
library(ggplot2)
library(Seurat)

## Application to  spatial transcriptomics data
  #The data should include the following objects:
  
  #dat.expression: a matrix object with the expression of genes (rows) in each cell (columns).
  #dat.coord: the spatial coordinates of the Visium spots. Please note that we are not using an embedding as input space here, but the actual 2D coordinates of 
            #spots inside the tissue. 

############################################################################################################################
# create data
data_name <- "data"
dat.expression <- read.csv("ounts.csv")
dat.coord <- read.csv("spatial_coors.csv")

############################################################################################################################

  ## running singleCellHaystack 
  
  set.seed(123)
  res <- haystack(dat.coord, dat.expression)
  
  
  # get the 1000 most significant DEGs, and cluster them by their distribution pattern in the 2D plot
  sorted.table <- show_result_haystack(res.haystack = res, n = 1000)
  gene.subset <- row.names(sorted.table)
  
  # cluster the genes by their expression pattern in the input space using hierarchical clustering
  res.hc <- hclust_haystack(dat.coord, dat.expression[gene.subset, ], grid.coordinates=res$info$grid.coordinates)
  # hclust_haystack returns as result a hclust tree, which we can cut into clusters using the cutree function. Here, we arbitrarily set the 
  
  
  ### cut into clusters using the cutree and output to table
  
  number.of.clustrs <- c(3,5)
  
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

