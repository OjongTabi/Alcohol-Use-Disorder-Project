# Load required libraries
library(pathfindR)
library(dplyr)
library(clusterProfiler)
set.seed(123) # for reproducibility

# Load the data
df <- read.csv("/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/placental tissue RNA-seq data DGE Analysis Results Using Limma Voom and CPM Filter of 1 in _min_samples Samples_Adjusted for Sex RIN Batch_Filtered for gene symbols.csv", row.names = 1)
head(df)

# Filter for significant genes
df <- dplyr::filter(df, P.Value < 0.05 & abs(logFC) >= 0.58496)
head(df)

# Get the gene list
gene_list <- rownames(df)
gene_list

# Create data frame with randomly generated fold-changes and adjusted p-values less than 0.05
df <- data.frame(Gene_Symbol = gene_list,
                 FC = 2^(df$logFC), # random normal variables for fold-changes
                 adj_p_val = df$P.Value) # random uniform variables for adjusted p-values

# Run pathfindR with the significant genes data frame
tryCatch({
  results <- run_pathfindR(df)
}, error = function(e) {
  print("There was an error:")
  print(e)
})
#str(results)

# Cluster the enriched terms
enrichment_results <- cluster_enriched_terms(results)


