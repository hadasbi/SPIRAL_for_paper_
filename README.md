# SPIRAL_for_paper_
 This code was used to produce the synthetic data and figures for the paper:
 
 SPIRAL: Significant Process InfeRence ALgorithm for single cell RNA-sequencing and spatial transcriptomics\ Hadas Biran, Tamar Hashimshony, Tamar Lahav, Or Efrat, Yael Mandel-Gutfreund and Zohar Yakhini.

## Description of zip file
 - SPIRAL_simulated_synthetic_data.zip: includes:
   * run_splatter.R: code to create synthetic datasets.
   * 01_analysis_of_Splatter_data.ipynb: basic analysis of synthetic dataset, and creation of the true cluster table.
   * data51, data52, ..., data60: folders of synthetic datasets (10 repetitions).

## Description of folders
 - comparison_of_all_methods_real_datasets: related to Figures 3b, 3c, 4c, 4d, 5b, 6b, 7b.
 - comparison_of_all_methods_to_ground_truth_of_Splatter_dataset: related to Figures 2b, 2c.
 - Hotspot: Hotspot output for relevant datasets (see datasets numbering below).
 - nsNMF: nsNMF output for relevant datasets (see datasets numbering below).
 - Seurat_DE_genes: Seurat output for relevant datasets (see datasets numbering below).
 - singlecellhaystack: singlecellhaystack output for relevant datasets (see datasets numbering below).
 - Slingshot: Slingshot output for relevant datasets (see datasets numbering below).
 - SpatialDE: SpatialDE output for relevant datasets (see datasets numbering below).

## Description of files
- 001a_compare_all_methods_to_ground_truth_of_Splatter_dataset.py + 001b_compare_all_methods_to_ground_truth_of_Splatter_dataset.ipynb - code to run methods comparison for synthetic data (Figure 2)
- 002_run_GOrilla_on_gene_clusters_of_methods.py - find enriched GO terms for the gene modules of all benchmark methods and SPIRAL (for all real datasets)
- 003_compare_GO_terms - compare GO terms between methods. The code uses accessory files: fb.gpad, goa_human.gpad, mgi.gpad, zfin.gpad, go-basic.obo.
- 004_compare_gene_module_sizes - compare gene module sizes between methods
- 005_compare_overall_number_of_participating_genes - compare overall number of participating genes between methods (not in the paper)
- run_hotspot.py - code to run Hotspot
- run_SpatialDE.py - code to run SpatialDE
- create_Seurat_objects_and_find_DE_genes.R - code to tun Seurat
- run_singleCellHaystack_on_Single_cell_RNAseq.R - code to run singleCellHaystack on single cell data
- run_singleCellHaystack_on_spatial_data.R - code to run singleCellHaystack on spatial data
- run_slingshot_tradeseq.R - code to run slingshot+tradeseq
- use_nsNMF_to_find_gene_modules.R - code to run nsNMF

## Datasets numbering:
- data1: single-cell RNAseq dataset of lymphoblastoid cells (Figure 3)
- data2: spatial transcriptomics dataset (Visium) of a mouse brain (Figure 5)
- data3: single cell dataset of Zebrafish differentiation (Figure 4)
- data4: bulk RNA-seq dataset of human embryonic stem cells during differentiation (Supplementary Figure 2)
- data5: spatial transcriptomics dataset (Visium) of a normal human prostate (Figure 7)
- data6: bulk RNA-seq dataset of mouse B cells that were treated with anti-IgM mAb to mimic B cell receptor stimulation (Supplementary Figure 3)
- data51, data52, ... , data60 - 10 repetitions of synthetic dataset (Figure 2)
- data66: spatial transcriptomics dataset (ST) of a coronal section of an 18 month old mouse brain with Alzheimer disease (Figure 6)   
 
