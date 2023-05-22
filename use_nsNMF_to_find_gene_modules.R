# Follow the procedure to find gene modules as in:
# Moncada, R., Barkley, D., Wagner, F., Chiodin, M., Devlin, J. C., Baron, M., ... & Yanai, I. (2020). Integrating microarray-based spatial transcriptomics and single-cell RNA-seq reveals tissue architecture in pancreatic ductal adenocarcinomas. Nature biotechnology, 38(3), 333-342.

#setwd('D:/OneDrive - Technion/SPIRAL_for_paper/')
setwd('/Data2/Hadas/SPIRAL_for_paper/')

library(NMF)

# load counts
#data_n = 1    #Zhang2019
#data_n = 2    #Visium mouse brain
#data_n = 3    #Zebrafish Wagner 2018
#data_n = 4    #Bulk RNA-seq data of human differentiation- Mendel-Gutfreund lab
#data_n = 5    #Spatial transcriptomics data of a normal human prostate (FFPE) (10x Genomics- Visium)
#data_n = 6    #Bulk RNA-seq data of mouse B-cells that were treated with anti-IgM mAb- 25 time points in the 6 hours post stimulation (Chiang et al. 2020)
#data_n = 7    #Splatter- simulated dataset


for (data_n in c(53, 54, 55, 56, 57, 58, 59, 60, 20, 21, 24, 61, 62, 63, 64, 65,
                 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80)) {
  print(paste("data", as.character(data_n)))
  
  #folder = paste0('../SPIRAL/static/analysis/data', as.character(data_n))
  folder = paste0('../SPIRAL_webtool/static/analysis/data', as.character(data_n))
  counts_file = file.path(folder, 'counts_norm_filt.csv') #already normalized to median count
  counts = read.csv(counts_file, sep='\t', row.names=1)
  # remove empty rows
  counts = counts[which(rowSums(counts != 0) > 0), ]
  
  # set out_folder
  out_folder = paste0('./comparison_to_nsNMF/data', as.character(data_n))
  if (!file.exists('./comparison_to_nsNMF')) {
    dir.create('./comparison_to_nsNMF')
  }
  if (!file.exists(out_folder)) {
    dir.create(out_folder)
  }
  
  # Normalize to 100000 per cell
  #size_factor = apply(counts, 2, function(x) sum(x))
  #counts = sweep(counts, 2, size_factor, "/")*100000
  
  # Normalize genes to z-scores
  mean_factor = apply(counts, 1, function(x) mean(x))
  counts = sweep(counts, 1, mean_factor, "-")
  std_factor = apply(counts, 1, function(x) sd(x))
  counts = sweep(counts, 1, std_factor, "/")
  
  # Set negative values to zero
  counts[counts<0] = 0
  
  # nsNMF
  rank = 20
  nout = nmf(counts, rank, "nsNMF", seed="nndsvd")
  
  cluster_table = data.frame(matrix(nrow=0,ncol=3))
  colnames(cluster_table) <- c('cluster', 'ngenes', 'genes')
  ind = 1
  
  # selecting the genes in the gene modules
  W = as.data.frame(basis(nout))
  gene_modules = list()
  for (col in colnames(W)) {
    W = W[order(W[,col], decreasing = TRUE),]
    gene_list = list()
    for (gene in rownames(W)) {
      if (max(W[gene, ]) == W[gene, col]) {
        gene_list = append(gene_list, gene)
      }
      else {
        break
      }
    }
    # keep modules with at least 20 genes
    if (length(gene_list) >= 20) {
      gene_modules[[length(gene_modules)+1]] = gene_list
      # write module to file
      fileConn = file(file.path(out_folder, paste0("module_", col, ".txt")))
      write(toString(gene_list), fileConn)
      close(fileConn)
      # add module to cluster table
      cluster_table[nrow(cluster_table) + 1,] = c(ind, length(gene_list),
                                                  toString(gene_list))
      ind = ind+1
    }
  }
  print(lengths(gene_modules))
  
  colnames(cluster_table) <- c('cluster', '#genes', 'genes')
  write.csv(cluster_table, paste0(out_folder,"/cluster_table.csv"), 
            row.names = FALSE)
}
