# Install and load necessary libraries
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

# BiocManager::install("fgsea")
# BiocManager::install("msigdbr")

library(fgsea)
library(msigdbr)

# Read your input file, and rank your genes based on logFC
df <- read.csv("/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/placental tissue RNA-seq data DGE Analysis Results Using Limma Voom and CPM Filter of 1 in _min_samples Samples_Adjusted for Sex RIN Batch_Filtered for gene symbols.csv", row.names = 1)

# Sort by logFC to get gene rankings
gene_ranking <- df$logFC
names(gene_ranking) <- rownames(df)
gene_ranking <- sort(gene_ranking, decreasing = TRUE) # High logFC at the top

# Load the gene sets (KEGG pathways in this case)
kegg_gene_sets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")

# Convert the DataFrame to the correct format for fgsea
kegg_gene_sets_list <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)

# Set a seed for reproducibility
set.seed(1234)

# Run fgsea with your ranked genes and the KEGG pathways
fgsea_res <- fgsea(pathways = kegg_gene_sets_list, stats = gene_ranking)


# Order by adjusted p-value
fgsea_res <- fgsea_res[order(fgsea_res$padj), ]

# View top significant pathways
head(fgsea_res)
#str(fgsea_res)


# Convert the list to a data frame
fgsea_res_df <- do.call(rbind, fgsea_res)
fgsea_res_df<-t(fgsea_res_df)
# Write the data frame to a .csv file
write.csv(fgsea_res_df, "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/Functional Annotation/fgsea/fgsea_results.csv", row.names = FALSE)
