#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")

#BiocManager::install("EnsDb.Hsapiens.v79")

library(enrichR)
library(readxl)
library(EnsDb.Hsapiens.v79)

# Read your input file, and rank your genes based on logFC
df <- read.csv("/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/placental tissue RNA-seq data DGE Analysis Results Using Limma Voom and CPM Filter of 1 in _min_samples Samples_Adjusted for Sex RIN Batch_Filtered for gene symbols.csv", row.names = 1)

# Upregulated genes
DEGs_upregulated = df[df$P.Value < 0.05 & df$logFC >= 0.58496, ]

# Downregulated genes
DEGs_downregulated = df[df$P.Value < 0.05 & df$logFC <= -0.58496, ]

# Include KEGG database for enrichment
dbs = c("GO_Molecular_Function_2023", "GO_Cellular_Component_2023", "GO_Biological_Process_2023", "KEGG_2023")

# Perform enrichment for upregulated genes
Y = rownames(DEGs_upregulated)
enriched = enrichr(Y, dbs)

# Save results to Excel files
library(writexl)
write_xlsx(data.frame(enriched[["GO_Molecular_Function_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_UP_MF_3.xlsx")
write_xlsx(data.frame(enriched[["GO_Cellular_Component_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_UP_CC_3.xlsx")
write_xlsx(data.frame(enriched[["GO_Biological_Process_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_UP_BP_3.xlsx")
write_xlsx(data.frame(enriched[["KEGG_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_UP_KEGG_3.xlsx")

# Perform enrichment for downregulated genes
Y = rownames(DEGs_downregulated)
enriched = enrichr(Y, dbs)

# Save results to Excel files
write_xlsx(data.frame(enriched[["GO_Molecular_Function_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_DOWN_MF.xlsx")
write_xlsx(data.frame(enriched[["GO_Cellular_Component_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_DOWN_CC.xlsx")
write_xlsx(data.frame(enriched[["GO_Biological_Process_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_DOWN_BP.xlsx")
write_xlsx(data.frame(enriched[["KEGG_2023"]]), "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/enrichr_DOWN_KEGG.xlsx")
