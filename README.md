**BIOINFORMATICS-PROJECT-COURSE**

R STUDIO CODE FOR BIOINFORMATICS PROJECT 

**Differential Gene Expression Analysis and Visualization**

**Project Overview**
This project performs a differential gene expression analysis to identify significant changes in treated and control samples. RNA sequencing data is analyzed to answer specific biological questions using bioinformatics tools. The project uses Galaxy for pre-processing and R for downstream analysis, applying statistical methods to identify and visualize differentially expressed genes.


**Pre-processing steps**

1.The SRR accesion number from NCBI (sample identifiers) are pasted onto GALAXY.EU with its appropiate CONTROL/TREATED TAG.

2.Reads are extracted in FASTQ format

3.The data im working with in this project is single end data

4.The data is then subjectd to FASTQ Quality control for both treated and control

5.Another tool MULTIQC is also used

6.The samples will be trimmed using fastp/cutadapt tool in galaxy

7.HISAT2/STAR tools are used as splice aware aligners to map data to genome

8.MultiQC is ran again for the data as a quality control

9.The data is then subjected to quantification

10.A tool called featureCounts is used to measure gene expression in RNA-Seq (Gene annotation)

11.edgeR is used to perform diffrential expression of count data. And the data is filtered with a minimum CPM of '1' and minimum samples of '3'

12.A tool called annotatemyID is also executed for ease of handling data in R



 **Visualization Techniques**:
   - **Scatter Plot**: To compare mean expression levels between conditions.
   - **PCA, t-SNE, and UMAP**: For dimensionality reduction and clustering analysis of treated and control samples.
   - **Box Plot**: To assess normalization and overall distribution of expression levels.
   - **Volcano Plot**: To visualize differentially expressed genes (Log Fold Change > 1, adjusted p-value < 0.05).
   - **Violin Plot**: To show the distribution of expression levels for significant genes.
   - **Heatmap**: To cluster significant genes and samples based on expression profiles.
     
 **Analysis**:
   - **KEGG Pathway analysis**: Identified enriched biological pathways such as cancer-related pathways and cellular senescence.

 
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
  - Annotated Data from Galaxy
    


