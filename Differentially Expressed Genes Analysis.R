#-------------------------------------------------------------------------------
#Description: R script to merge the counts matrices, perform DE to get significant genes by brain region, and generate visualizations
#-------------------------------------------------------------------------------

# Load the necessary libraries
# tidyverse includes dplyr, ggplot2 and stringr
# variancePartition: used for gene expression data analysis
# limma and edgeR: used for differential expression analysis
# rtracklayer: used for genomic annotations
# EnhancedVolcano: used for plotting volcano plots
# magrittr: provides the pipe operator %>% used for chaining operations
library(tidyverse)
library(variancePartition)
library(limma)
library(edgeR)
library(rtracklayer)
library(EnhancedVolcano)
library(magrittr)

# Define the working directory where the data files are located
setwd("/projectnb/zhangh/Ojong_Tabi/DCC_Run_CAU/")

# Define the notin operator
`%notin%` <- negate(`%in%`)

# Import the GTF file (contains gene and transcript annotations) and filter the entries
finalgtf <- rtracklayer::import("/projectnb/zhangh/JP_HumBrain/reference/tst/gencode_v39.lncipedia_v52.annotationALLCOMBINED.gtf", format="gtf")
finalgtf <- data.frame(finalgtf)
finalgtf <- subset(finalgtf, finalgtf$type=="gene", select=c("gene_id", "gene_type", "gene_name"))

# Load the counts matrix (contains gene expression data)
counts_matrix <- read.table('/projectnb/zhangh/Ojong_Tabi/DCC_Run_Amy/Amy_merged_file.tsv', header = TRUE, sep = '\t')
# Generate a unique identifier for each gene by concatenating the chromosome, start, and end coordinates
counts_matrix$Coordinates <- paste(counts_matrix$Chr, counts_matrix$Start, counts_matrix$End, sep="_")
# Remove the original coordinate columns
counts_matrix <- counts_matrix[, -c(1:3)]
# Reorder the columns to move "Coordinates" to the first column
counts_matrix <- counts_matrix[, c("Coordinates", setdiff(names(counts_matrix), "Coordinates"))]
# Filter out rows where all counts are zero
counts_matrix_ready <- counts_matrix[apply(counts_matrix[,c(2:25)], 1, function(x) !all(x==0)), ] 

# Load the metadata file
meta <- read.csv("/projectnb/zhangh/JP_HumBrain/mRNAseq/PCA_DE/Set1_2_Phenotype_192samples_metadata.csv", sep = ",", header = TRUE)
# Remove unnecessary columns
meta <- meta[,-c(1,5,8,10,17,20)]
# Replace missing values in the "Alcohol_Intake" column with 0
meta$Alcohol_Intake[is.na(meta$Alcohol_Intake)] <- 0
# Convert the "LeftRightBrain" column into a factor and relabel it to "Left2_Right1_Brain"
meta$LeftRightBrain <- factor(meta$LeftRightBrain, levels=c("Right", "Left"), labels=c(1,2)) ; colnames(meta)[9] <- "Left2_Right1_Brain"
# Convert the "Left2_Right1_Brain" column into integer
meta$Left2_Right1_Brain <- as.integer(meta$Left2_Right1_Brain)
# Convert the "Classification" column into a factor with specified levels
meta$Classification <- factor(meta$Classification, levels = c("Control", "AlcoholUseDisorder"), labels=c("Control", "AlcoholUseDisorder"))
# Create a new column "LiverClass2" based on the "LiverClass" column
meta$LiverClass2 <- ifelse(meta$LiverClass=="Normal", 1, 2)
# Convert the "LiverClass2" column into integer
meta$LiverClass2 <- as.integer(factor(meta$LiverClass2, levels=c(1,2)))
# Remove the original "LiverClass" column
meta <- subset(meta, select=-c(LiverClass))
# Filter the data to retain only the samples that start with "AMY"
meta <- meta %>%
  filter(grepl("^AMY", SampleID))

# Replace "HIPPO_" with "HIP_" in the "SampleID" column
final_meta <- meta
final_meta$SampleID <- str_replace(final_meta$SampleID, pattern = "HIPPO_", replacement = "HIP_")
# Move the "SampleID" column to rownames
rownames(final_meta) <- NULL
final_meta <- final_meta %>% column_to_rownames(var="SampleID")

# Check the data structure of the final metadata
str(meta)

# Perform variance partitioning following the documentation
# First, convert categorical variables into factors
catCols <- c("Set", "BrainRegion", "M1_F2", "Left2_Right1_Brain", "Smoking2_NoSmoking1", "LiverClass2")
final_meta2 <- final_meta %<>% mutate_at(catCols, factor)

# Change rownames to start with "AMY_" if they don't already
rownames(final_meta2) <- ifelse(grepl("^AMY_", rownames(final_meta2)), 
                                rownames(final_meta2), 
                                paste("AMY", gsub("^X", "", rownames(final_meta2)), sep = "_"))

# Define a function to extract the sample ID from the column names
get_sample_id <- function(col_name) {
  # Look for "AMY_###" at the beginning of the column name
  if (grepl("^AMY_\\d+", col_name)) {
    return(sub("^(AMY_\\d+).*", "\\1", col_name))
  } 
  # Look for "X###_AMY" within the column name
  else if (grepl("X\\d+_AMY", col_name)) {
    parts <- strsplit(col_name, "_")[[1]]
    return(paste("AMY", gsub("^X", "", parts[1]), sep = "_"))
  }
  return(NA)
}

# Apply this function to the column names of 'counts_matrix_ready'
sample_ids <- sapply(colnames(counts_matrix_ready), get_sample_id)

# Retain only the columns whose extracted sample ids exist in 'final_meta2'
counts_matrix_ready_amy <- counts_matrix_ready[, sample_ids %in% rownames(final_meta2)]
# Renaming the columns
names(counts_matrix_ready_amy) <- sample_ids[names(counts_matrix_ready_amy)]                                        # First, create a separate metadata frame for AMY region
final_meta2_amy <- final_meta2[grepl("^AMY_", rownames(final_meta2)), ]

# Factorize categorical variables for variancePartition
catCols <- c("Set", "BrainRegion", "M1_F2", "Left2_Right1_Brain", "Smoking2_NoSmoking1", "LiverClass2")
final_meta2_amy <- final_meta2_amy %>% mutate_at(catCols, factor)

# Store gene identifiers separately
gene_ids <- counts_matrix_ready_amy$Coordinates
# Create the counts matrix without the gene identifiers column
counts_matrix_ready_numeric_amy <- counts_matrix_ready_amy[, -1]

# Create a DGEList object, which is suitable for downstream differential expression analysis
gExpr <- DGEList(counts = counts_matrix_ready_amy)
# Add the gene IDs back to the DGEList object
gExpr$genes <- gene_ids

# Filter out the lowly expressed genes
fltr <- filterByExpr(gExpr, robust=TRUE)
gExpr <- gExpr[fltr,]

# Normalize the library sizes
gExpr <- calcNormFactors(gExpr)

# Create the design matrix for the differential expression analysis
design <- model.matrix( ~ Classification, final_meta2_amy)

# Apply the voom transformation, which is necessary for limma's linear model fitting
vobjGenes <- voom(gExpr, design)

# Fit the variance partitioning model
form <- ~ Age + PMI + Brain_pH + (1|Classification) + RIN + Alcohol_Intake + (1|Set) + (1|M1_F2) + BrainWeight + (1|LiverClass2) + (1|Left2_Right1_Brain) + (1|Smoking2_NoSmoking1)
varPart <- fitExtractVarPartModel(vobjGenes, form, final_meta2_amy)

# Plot the variance partitioning results
plotVarPart(sortCols(varPart))
# Load necessary libraries
library(edgeR)
library(limma)

# Assign unique row names to counts_matrix_ready_amy
# Create a unique identifier by combining "Coordinates" and a sequence
unique_id <- paste(final_meta2_amy$Coordinates, seq_len(nrow(final_meta2_amy)), sep="_")
rownames(final_meta2_amy) <- unique_id

# Confirm that the sample labels in metadata (row names) are the same as in count data (column names)
identical(rownames(final_meta2_amy), colnames(counts_matrix_ready))

# Construct the design matrix
design <- model.matrix(~0 + Set + M1_F2 + RIN + PMI + Age, final_meta2_amy)

# Generate DGEList object
# Assuming counts_matrix_ready's first column is geneid
# Append row number to ensure uniqueness
rownames(counts_matrix_ready) <- paste(counts_matrix_ready[, 1], seq_len(nrow(counts_matrix_ready)), sep="_")
counts_matrix_ready <- counts_matrix_ready[, -1]  # remove the first column

# Convert the dataframe to a matrix
counts_matrix_ready <- as.matrix(counts_matrix_ready)

# Now generate the DGEList object
dge <- DGEList(counts=counts_matrix_ready)



# Filter out low expressed genes
keep <- filterByExpr(dge, design, robust=TRUE)
dge <- dge[keep,,keep.lib.sizes=FALSE]
dge <- calcNormFactors(dge, method = "TMM")

# voom transformation
v <- voom(dge, design)

# Determine correlation between repeated measurements from the same patient
patient <- sapply(final_meta2_amy$SubjectID, function(x) strsplit(x, "_")[[1]], USE.NAMES = F)[2,]
patient.block = as.factor(patient)
dupcor <- duplicateCorrelation(v, design, block = patient.block)

# Linear model fit
fit <- lmFit(v, design, block = patient.block, correlation = dupcor$consensus.correlation)

# Create contrasts of interest
contrast_levels <- levels(final_meta2_amy$Set)
cont.matrix <- makeContrasts(paste("Diff=", contrast_levels[1], "-", contrast_levels[2], sep=""), levels = design)

# Fit contrasts matrix
fit2 <- contrasts.fit(fit, cont.matrix)

# Apply empirical Bayes smoothing
fit2 <- eBayes(fit2, robust=TRUE)

# Extract top differentially expressed genes
amyT <- topTable(fit2, n=Inf, coef = 'Diff=1-2', adjust="BH")
head(amyT)
# Extract top differentially expressed genes
#amyT <- topTable(fit2, n=Inf, adjust="BH")




# Write output to file
#write.csv(amyT, file="AmyT_diff_expr_results.csv")

# Load necessary libraries
library(ggplot2)

# Generate a column for -log10 adjusted p-values
#amyT$PValueAdjLog <- -log10(amyT$adj.P.Val)

# Generate a column for log2 fold changes
#amyT$logFC <- log2(amyT$logFC)

# Generate the Volcano Plot
ggplot(data = amyT, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(alpha = 0.5, size = 1.5) +
  theme_minimal() +
  theme(text = element_text(size=20)) +
  ggtitle("Volcano Plot of Differentially Expressed circRNAs") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-Value") +
  geom_hline(yintercept = -log10(0.05), color = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), color = "blue", linetype = "dashed")

# Save the plot as a PNG image
#ggsave(filename = "VolcanoPlot_DGE.png", dpi = 300)

#  filter  raw counts matrix to only include genes from the topTable output
filtered_counts_matrix <- counts_matrix_raw[rownames(counts_matrix_raw) %in% rownames(amyT), ]


