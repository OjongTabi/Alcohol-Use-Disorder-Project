# Load necessary packages
if (!require(goseq)) install.packages("goseq")
library(goseq)

# Define the path
data_path <- "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/LimmaVoom_CPM1_Covariates_2023.07.17/placental tissue RNA-seq data DGE Analysis Results Using Limma Voom and CPM Filter of 1 in _min_samples Samples_Adjusted for Sex RIN Batch_Filtered for gene symbols.csv"

# Check if the file exists
if (!file.exists(data_path)) {
  stop("File does not exist. Please check the path.")
}

# Read the data
df <- read.csv(data_path, row.names = 1)

# Check if data is loaded correctly
if (is.null(df)) {
  stop("Failed to load data. Please check the file.")
}

# Extract gene names from the dataframe
gene.names <- row.names(df)

# Check gene ID's
if (length(gene.names) == 0) {
  stop("No gene ID's found. Please check the data.")
}

# Determine differentially expressed genes
#DEgenes <- abs(df$logFC) >= 0.58496
DEgenes <- df$P.Value < 0.05 & abs(df$logFC) >= 0.58496

# Generate a named binary vector indicating differential expression
names(DEgenes) <- gene.names

# Get gene lengths data
gene.lengths <- getlength(gene.names, "hg19", "geneSymbol")

# Check gene lengths
if (is.null(gene.lengths)) {
  stop("Could not obtain gene lengths. Please check the gene names and database.")
}

# Run GOseq analysis
pwf <- nullp(DEgenes, bias.data = gene.lengths)

# Check if the pwf object is created
if (is.null(pwf)) {
  stop("Failed to create the pwf object. Please check the input to nullp().")
}

# Correct for multiple testing
GO.wall <- goseq(pwf, "hg19", "geneSymbol")

# Check if the GO.wall object is created
if (is.null(GO.wall)) {
  stop("Failed to create the GO.wall object. Please check the input to goseq().")
}

# Set your output directory
output_dir <- "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/Functional Annotation/GSEA/goseq/"

# Write GO.wall to a csv file in the specified directory
write.csv(GO.wall, file = paste0(output_dir, "GO_wall_results.csv"))

