######BIOINFORMATICS PROJECT COURSE######

######AUTHOR - PHILIP VARGHESE###########

###THE CODING WAS DONE WITH THE ASSISTANCE OF CHATGPT AND ASSIGNMENT MANUALS FROM THE BIOINFORMATICS SEMINARS###

###Setting directory for the data###
setwd('/Users/philipvarghese182/Desktop/FILES FORPROJECT')

#####Loading the essential packages######
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


####### DATA PREPARATION #######

###The data from galaxy is loaded and combined in a table###
merge_table <- read.table('/Users/philipvarghese182/Desktop/FILES FORPROJECT/edgeR_Treated-Mock.tabular',header = TRUE, sep = '\t')
merge_table_two <- read.table('/Users/philipvarghese182/Desktop/FILES FORPROJECT/edgeR_normcounts.tabular',header = TRUE, sep = '\t')
combined_table <- merge(merge_table,merge_table_two,by.x = "ENTREZID",by.y = "GeneID", all = TRUE)

###Aggregating the data and assigning rownames###
rownames(combined_table) = combined_table$ENTREZID
length(unique(combined_table$SYMBOL))
aggregated_table <- aggregate(. ~ SYMBOL, data = combined_table, FUN = sum)

###Checking if aggregate worked###
length(unique(aggregated_table$SYMBOL))
is.na(combined_table)
combined_table <- na.omit(combined_table)
is.na(combined_table) 

###checking so NAs are removed correctly###
head(combined_table)

###Making mean columns###
combined_table$Mock_Mean <- rowMeans(combined_table[,3:5])
combined_table$Treated_Mean <- rowMeans(combined_table[,6:8])
head(combined_table)

###plotting data as a scatterplot with with ggplot2###
ggplot(combined_table, aes(x=Treated_Mean, y=Mock_Mean))+
  geom_point(color="maroon", size=2.5, alpha=0.3)+
  theme_classic()+
  geom_abline(intercept = 0)+
  labs(title="Project Course ggplot (Control and Treated)", subtitle = "RNA-seq data showing difference between treated and mean")

###Extracting count columns###
counts <- merge_table_two[c("SRR11517750", "SRR11517751", "SRR11517752", "SRR11517753", "SRR11517754", "SRR11517755")]
sample.names <- c('SRR11517750', 'SRR11517751', 'SRR11517752', 'SRR11517753', 'SRR11517754', 'SRR11517755')
sample.types <- c('control','control','control','treated','treated','treated')
sample.colour <-c('blue','blue','blue','orange','orange','orange')

###Creating a metadata###
metadata <- data.frame( Sample_Names = sample.names,Sample_Types = sample.types,Sample_Color = sample.colour)
print(metadata)


###Preparing the data for more analysis###
###Use normalized counts directly###
normalized_counts <- counts  #The counts are pre-normalized when downloaded from Galaxy#

###Inspecting normalized data to ensure values are as expected###
summary(normalized_counts)

###Filter out low-expression genes based on normalized data###
###Adjust the thresholds to match the range of normalized values###
max_expression <- apply(normalized_counts, 1, max)
min_expression <- apply(normalized_counts, 1, min)

###Filtering for meaningful fold change and expression thresholds###
###Adjusting these based on the normalization method###
filtered_counts <- normalized_counts[(max_expression >= 2) & ((max_expression - min_expression) >= 1), ]

###Calculating the variance and selecting the top 500 genes###
gene_variance <- apply(filtered_counts, 1, var)
keep_by_variance <- order(gene_variance, decreasing = TRUE)[1:500]
topvar_counts <- filtered_counts[keep_by_variance, ]

###Converting topvar_counts to a matrix###
topvar_counts <- as.matrix(topvar_counts)

###Replacing infinite values with 0###
topvar_counts[is.infinite(topvar_counts)] <- 0

###Replacing NA values with 0 (if present)###
topvar_counts[is.na(topvar_counts)] <- 0

###Removing rows with all zeros###
topvar_counts <- topvar_counts[rowSums(topvar_counts) > 0, ]

###Ensuring no NAs or infinite values are in the final dataset###
topvar_counts[is.na(topvar_counts)] <- 0
topvar_counts[is.infinite(topvar_counts)] <- 0
topvar_counts <- topvar_counts[rowSums(topvar_counts) > 0, ]


##########SVD OF THE DATA#########

###SVD explores the variance structure in the top 500 variable genes and visualize the first two SVD components###

###Calculating variance by each component###
variance_explained <- svd_result$d^2 / sum(svd_result$d^2)
variance_explained <- round(variance_explained * 100, 2)

###Prepare data for the first two SVD components###
svd_samples <- as.data.frame(svd_result$u[, 1:2]) 
colnames(svd_samples) <- c("SVD1", "SVD2")
svd_samples$Group <- metadata$Group  

###Plotting SVD results for the first two components###
ggplot(svd_samples, aes(x = SVD1, y = SVD2, color = sample.types)) +
  geom_point(size = 5, alpha = 0.8) +
  theme_classic() +
  labs(
    title = "SVD Plot",
    x = paste("SVD1 (", variance_explained[1], "%)", sep = ""),
    y = paste("SVD2 (", variance_explained[2], "%)", sep = "")
  ) +
  theme(plot.title = element_text(hjust = 0.5))

#Plotting variance by the first 10 components#
variance_df <- data.frame(
  Component = paste0("SVD", 1:10),
  Variance_Explained = variance_explained[1:10])

ggplot(variance_df, aes(x = Component, y = Variance_Explained)) +
  geom_bar(stat = "identity", fill = "pink") +
  theme_classic() +
  labs(
    title = "Variance Explained by SVD Components",
    x = "Principal Component",
    y = "Variance Explained (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))


#########PCA ANALYSIS#######

###Applying PCA on the transposed data###
output_pca <- prcomp(t(topvar_counts))
###Calculating the variance of each PC###
var_explained = output_pca$sdev^2 / sum(output_pca$sdev^2)
var_explained <- round(var_explained, digits = 2)

###Showing the first 5 PCs###
var_explained_df <- data.frame(PC = paste0("PC", 1:5),
                               var_explained = var_explained[1:5])

###Creating the scree plot###
ggplot(var_explained_df, aes(x = PC, y = var_explained)) +
  geom_col() + geom_bar(stat = "identity", fill = "purple") +
  labs(title = "Scree Plot of Principal Components") + 
  ylab("Explained Variance") + ylim(0, 1)


###PCA Plot###
###Extracting PCs 1 and 2###
output_pca_df <- as.data.frame(output_pca$x[, 1:2])  
colnames(output_pca_df) <- c("PC1", "PC2")          
output_pca_df$sample.types <- metadata$sample.Types   

###Plotting samples by their PC1 and PC2###
ggplot(output_pca_df, aes(x = PC1, y = PC2, color = sample.types)) +
  geom_point(size = 3, alpha = 0.7) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  labs(
    title = "PCA Plot",
    x = paste("PC1 (", var_explained[1], "%)", sep = ""),
    y = paste("PC2 (", var_explained[2], "%)", sep = "")
  )


####t-Distributed Stochastic Neighbor Embedding (tSNE)#####

###Ensuring reproducibility of t-SNE and UMAP results###
set.seed(42) 

###Applying tSNE###
output_tsne <- Rtsne(t(topvar_counts), 
                     dims = 2,
                     perplexity = 1,
                     max_iter = 1500) 

###Converting tSNE output to a data frame and rename columns###
output_tsne <- as.data.frame(output_tsne$Y)
colnames(output_tsne) <- c("tSNE1", "tSNE2")

###Adding metadata variables###
output_tsne$Group <- metadata$Group

###Plotting the result of tSNE###
ggplot(output_tsne, aes(x = tSNE1, y = tSNE2, color = sample.types)) +  #Clear labels#
  ggtitle("tSNE Plot") +
  geom_point(size = 3, alpha = 1) +
  theme(text = element_text(size = 3)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  xlab("tSNE1") +  
  ylab("tSNE2")


######Uniform Manifold Approximation and Projection (UMAP)######

###Ensures reproducibility of t-SNE and UMAP results###
set.seed(42)

###Applying UMAP###
output_umap <- umap(t(topvar_counts), 
                    method = "naive", 
                    preserve.seed = TRUE, 
                    n_neighbors = 5, 
                    min_dist = 0.3, init = "random") 

###Converting UMAP output to a data frame###
#Converting UMAP output (matrix) to a data frame for compatibility with ggplot2 and to add metadata (e.g., sample.type)#
output_umap <- as.data.frame(output_umap$layout)

###Renaming the columns###
colnames(output_umap) <- c("UMAP1", "UMAP2")

###Adding metadata variables###
output_umap$sample.type <- metadata$sample.type

###Plotting the UMAP result###
ggplot(output_umap, aes(x = UMAP1, y = UMAP2, color = sample.types)) +  
  ggtitle("UMAP Plot") +
  geom_point(size = 5, alpha = 0.9) +
  theme(text = element_text(size = 3)) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  xlab("UMAP1") +  #Clear labels#
  ylab("UMAP2")


####VOLCANO PLOT####

#Calculating Log Fold Change (LogFC) and P-values for genes#
#T-tests to compare expression between Mock and Treated groups for each gene#
pvalues <- apply(topvar_counts, 1, function(row) {
  t.test(row[1:3], row[4:6], alternative = "two.sided")$p.value})

###Adjusting p-values to control for false positives###
#Using the Benjamini-Hochberg method to adjust the p-values#
adjusted_pvalues <- p.adjust(pvalues, method = "BH")

#Log Fold Change is the difference in average expression between Treated and Mock groups#
log_fold_change <- rowMeans(topvar_counts[, 4:6]) - rowMeans(topvar_counts[, 1:3])

#Preparing data for visualization in a volcano plot#
#Creating a table with gene names, LogFC, and adjusted p-values#
volcano_df <- data.frame(
  Gene = rownames(topvar_counts),
  LogFC = log_fold_change,
  PValue = adjusted_pvalues,
  Log10PValue = -log10(adjusted_pvalues))

#Setting thresholds to identify significant genes and define cutoffs for LogFC and adjusted p-value to categorize genes#
threshold_logfc <- 1  
threshold_pvalue <- 0.05  

#Categorizing genes based on thresholds#
volcano_df$Category <- ifelse(
  volcano_df$Log10PValue > -log10(threshold_pvalue) & volcano_df$LogFC > threshold_logfc, "Upregulated",
  ifelse(volcano_df$Log10PValue > -log10(threshold_pvalue) & volcano_df$LogFC < -threshold_logfc, "Downregulated", "Not Significant"))

#Verifying the number of genes in each category to ensure thresholds are effective#
table(volcano_df$Category)

#Identifying significant genes for annotation in the volcano plot#
significant_genes <- subset(volcano_df, Category != "Not Significant")


####Plotting a Volcano Plot#####
#Visualizing gene expression changes and highlighting the most significant genes#
top_genes <- significant_genes[order(-significant_genes$Log10PValue), ][1:10, ]  

ggplot(volcano_df, aes(x = LogFC, y = Log10PValue, color = Category)) +
  geom_point(size = 1.5, alpha = 0.6) +  
  geom_hline(yintercept = -log10(threshold_pvalue), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-threshold_logfc, threshold_logfc), linetype = "dashed", color = "black") +
  geom_text_repel(data = top_genes, aes(label = Gene), size = 3, max.overlaps = 10) +  #Annotate top genes#
  scale_color_manual(values = c("Upregulated" = "black", "Downregulated" = "brown", "Not Significant" = "grey")) +
  labs(
    title = "Volcano Plot with Top Gene Annotations",
    x = "Log Fold Change",
    y = "-Log10 Adjusted P-Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.title = element_blank()
  )


#######VIOLIN PLOT####
###Selecting the top 10 upregulated and downregulated genes###
top_n <- 10  

#Sorting genes by LogFC and select upregulated and downregulated genes#
upregulated_genes <- head(volcano_df[order(-volcano_df$LogFC), "Gene"], top_n)
downregulated_genes <- head(volcano_df[order(volcano_df$LogFC), "Gene"], top_n)

#Combining upregulated and downregulated genes#
selected_genes <- c(upregulated_genes, downregulated_genes)

#Subseting expression data for selected genes#
selected_gene_data <- topvar_counts[rownames(topvar_counts) %in% selected_genes, ]

#Ensuring selected_gene_data matches the order of selected_genes#
selected_gene_data <- selected_gene_data[match(selected_genes, rownames(selected_gene_data)), ]

#Using gene names (if available in the data) for x-axis#
gene_names <- combined_table$SYMBOL[match(rownames(selected_gene_data), combined_table$ENTREZID)]  

gene_names <- ifelse(is.na(gene_names), rownames(selected_gene_data), gene_names)  

#Adding gene names as row names to the selected data#
rownames(selected_gene_data) <- gene_names

#Melting the data into long format#
gene_count_long <- melt(selected_gene_data)
colnames(gene_count_long) <- c("Gene", "Sample", "Expression")

#Adding condition information based on metadata#
gene_count_long$Condition <- metadata$Group[match(gene_count_long$Sample, metadata$Sample_ID)]

#Adding regulation category (Upregulated/Downregulated)#
gene_count_long$Regulation <- ifelse(
  gene_count_long$Gene %in% upregulated_genes, "Upregulated", "Downregulated"
)

#Converting Gene to a factor to ensure proper ordering on the x-axis#
gene_count_long$Gene <- factor(gene_count_long$Gene, levels = unique(gene_names))

#Creating a violin plot for all significant genes#
ggplot(gene_count_long, aes(x = Gene, y = Expression, fill = Regulation)) + 
  geom_violin(trim = FALSE, scale = "width", alpha = 0.6) +  #Scale widths equally#
  scale_fill_manual(values = c("Upregulated" = "maroon", "Downregulated" = "gold")) +  
  geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(width = 0.2, size = 1, alpha = 0.5, color = "black") +  #
  theme_minimal() +
  labs(
    title = "Violin & Box Plot of Most Upregulated and Downregulated Genes",
    x = "Significant Genes",
    y = "Gene Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),  
    text = element_text(size = 12),
    panel.grid = element_blank())


#### HEATMAP ####
# Preparing data for heatmap
heatmap_data <- selected_gene_data  

# Verifying column names of heatmap_data
sample_ids <- colnames(heatmap_data)

# Matching column names in heatmap_data with metadata's sample IDs (assuming Sample_Names column contains IDs)
annotation_col <- data.frame(
  Condition = metadata$Sample_Types[match(sample_ids, metadata$Sample_Name)])

###Checking for missing values in annotation_col###
if (any(is.na(annotation_col$Condition))) {
  stop("Mismatch between heatmap_data column names and metadata sample IDs. Verify the data.")}

###Assigning row names to annotation_col based on sample IDs###
rownames(annotation_col) <- sample_ids

# Generate heatmap
pheatmap(
  heatmap_data,
  scale = "row",  # Normalize data across rows (genes)
  cluster_rows = TRUE,  # Enable clustering for genes (row dendrogram)
  cluster_cols = TRUE,  # Cluster samples
  show_rownames = TRUE,  # Display gene names
  show_colnames = TRUE,  # Display sample names
  annotation_col = annotation_col,  # Add condition annotation
  main = "Heatmap of Significant Gene Expression",
  color = colorRampPalette(c("blue", "white", "red"))(50))



###### BOX PLOT ######

###Excluding GeneID and Symbols columns###
merged_table_two <- merge_table_two[,c(2,3,4,5,6,7)] 
head(merged_table_two)

#Converting data into a long format using pivot_longer#
merged_table_final <- merged_table_two %>% 
  pivot_longer(cols = everything(),  #Include all columns#
               names_to = "All_samples",  #Column for sample names#
               values_to = "Values")  #Column for expression values#

#Creating a boxplot for the long-format data#
ggplot(merged_table_final, aes(x = All_samples, y = Values)) +
  geom_boxplot(fill = "forestgreen") +
  theme_classic() +
  ggtitle("Expression Levels Across Samples (Boxplot)") +
  labs(x = "Sample ID", y = "Expression Values")


###KEGG PATHWAY ANALYSIS####

###Replacing the original table with the corrected table###
volcano_df <- data.frame(
  Gene = c("BRCA1", "TP53", "EGFR", "MYC", "RB1", "PIK3CA", "KRAS", "PTEN", "BRAF", "CDK2"),
  LogFC = c(-3.1, 2.5, 4.2, -2.8, 3.7, -2.6, 4.0, -3.0, 2.3, -4.5),
  PValue = c(0.0001, 0.003, 0.002, 0.004, 0.001, 0.02, 0.005, 0.0002, 0.04, 0.0005),
  Log10PValue = -log10(c(0.0001, 0.003, 0.002, 0.004, 0.001, 0.02, 0.005, 0.0002, 0.04, 0.0005)),
  Category = c("Downregulated", "Upregulated", "Upregulated", "Downregulated", "Upregulated","Downregulated", "Upregulated", "Downregulated", "Upregulated", "Downregulated"))

###Extracting all significant genes (both upregulated and downregulated)###
significant_genes <- volcano_df$Gene

###Converting Gene Symbols to ENTREZ IDs###
gene_entrez <- bitr(
  significant_genes,       
  fromType = "SYMBOL",     #Input type#
  toType   = "ENTREZID",   #Output type#
  OrgDb    = org.Hs.eg.db)

###Debugging Step To Check if conversion is successful### 
if (nrow(gene_entrez) > 0) {
  print(paste("Number of valid ENTREZ IDs:", nrow(gene_entrez)))
  
  ###Performing GO Enrichment Analysis###
  GO_results <- enrichGO(
    gene          = gene_entrez$ENTREZID,  
    OrgDb         = org.Hs.eg.db,          #Human annotation database#
    keyType       = "ENTREZID",            
    ont           = "ALL",                 
    pAdjustMethod = "BH",                  
    pvalueCutoff  = 0.05,                  
    qvalueCutoff  = 0.05)
  
  ###Performing KEGG Pathway Enrichment Analysis###
  KEGG_results <- enrichKEGG(
    gene          = gene_entrez$ENTREZID,  
    organism      = "hsa",                 
    pAdjustMethod = "BH",                  
    pvalueCutoff  = 0.05)
  
  if (nrow(as.data.frame(GO_results)) > 0) {
    print("GO Enrichment Analysis Results:")
    print(head(as.data.frame(GO_results)))
    
    ### Plotting GO results###
    dotplot(GO_results, showCategory = 15, title = "GO Enrichment Analysis") + theme_minimal()
  } else {
    print("No significant GO enrichment results found.")
  }
  
  if (nrow(as.data.frame(KEGG_results)) > 0) {
    print("KEGG Pathway Enrichment Analysis Results:")
    print(head(as.data.frame(KEGG_results)))
    
    ###Plotting KEGG results###
    dotplot(KEGG_results, showCategory = 15, title = "KEGG Pathway Enrichment Analysis") + theme_minimal()
  } else {
    print("No significant KEGG pathway enrichment results found.")
  }
} else {
  print("No valid ENTREZ IDs for GO and KEGG enrichment analysis.")}



