#####BIOINFORMATICS-PROJECT-COURSE#####
R STUDIO CODE FOR BIOINFORMATICS PROJECT 
#Differential Gene Expression Analysis and Visualization#

###Project Overview###
This project performs a differential gene expression analysis to identify significant changes in treated and control samples. RNA sequencing data is analyzed to answer specific biological questions using bioinformatics tools. The project uses Galaxy for pre-processing and R for downstream analysis, applying statistical methods to identify and visualize differentially expressed genes.



## Features
1. **Preprocessing and Normalization**:
   - Utilized normalized RNA-seq data from Galaxy.
   - Filtered low-expressed genes and retained the top 500 most variable genes.

2. **Visualization Techniques**:
   - **Box Plot**: To assess normalization and overall distribution of expression levels.
   - **Volcano Plot**: To visualize differentially expressed genes (Log Fold Change > 1, adjusted p-value < 0.05).
   - **Violin Plot**: To show the distribution of expression levels for significant genes.
   - **Heatmap**: To cluster significant genes and samples based on expression profiles.
   - **PCA, t-SNE, and UMAP**: For dimensionality reduction and clustering analysis of treated vs. control samples.
   - **Scatter Plot**: To compare mean expression levels between conditions.

3. **Functional Analysis**:
   - **KEGG Pathway Enrichment**: Identified enriched biological pathways such as cancer-related pathways and cellular senescence.

 
###Prerequisites###
- R software
- RStudio software
- Required R packages:
         library(edgeR)
         library(tidyverse)
         library(affy)
         library(ggplot2)
         library(ggpubr)
         library(MASS)
         library(proxy)
         library(Rtsne)
         library(sva)
         library(gplots)
         library(clusterProfiler)
         library(enrichplot)
         library(org.Hs.eg.db) 
         library(pheatmap)
         library(reshape2)
         library(ggrepel)
         library(AnnotationDbi)
         library(umap)


