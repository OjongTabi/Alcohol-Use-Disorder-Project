library(clusterProfiler)
library(org.Hs.eg.db)

# Load the data
df <- read.csv("/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/placental tissue RNA-seq data DGE Analysis Results Using Limma Voom and CPM Filter of 1 in _min_samples Samples_Adjusted for Sex RIN Batch_Filtered for gene symbols.csv", row.names = 1)

# Filter for significant genes
df <- dplyr::filter(df, P.Value < 0.05 & abs(logFC) >= 0.58496)

# Get the gene list (gene symbols)
gene_list_symbols <- rownames(df)

# Convert gene symbols to Entrez IDs
gene_list <- mapIds(org.Hs.eg.db, 
                    keys = gene_list_symbols,
                    column = "ENTREZID", 
                    keytype = "SYMBOL", 
                    multiVals = "first")

# Remove NA values (if any) after conversion
gene_list <- gene_list[!is.na(gene_list)]

# Perform KEGG enrichment analysis
res <- enrichKEGG(gene = gene_list,
                  organism = 'hsa',
                  keyType = 'kegg',
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)

# Check if there are any enriched terms
if(nrow(res) > 0){
  # Visualize the result
  dotplot(res)
} else {
  print("No significantly enriched pathways found.")
}
