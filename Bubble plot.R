# Load the necessary package
library(readxl)
library(ggplot2)
# Define the file path
file_path <- "/projectnb/zhangh/Ojong_Tabi/HumanPlacenta_RNAseq/Functional Annotation/DAVID/2023.7.22_DAVID Functional Annotation_244 genes_P0.05.xlsx"

# Read the third sheet of the Excel file
data <- read_excel(file_path, "Chart3")

# Print the data
print(data)

# Create the bubble plot
ggplot(data, aes(x = `[-log10(Pvalue)]`, y = reorder(Term, `[-log10(Pvalue)]`), size = Count)) +
  geom_point(alpha = 0.6) +
  scale_size(range = c(1, 8)) +
  theme(
    axis.text.y = element_text(size = 6, angle = 0, hjust = 1),
    title = element_text(size = 12),
    axis.title = element_text(size = 10)
  ) +
  labs(
    title = "KEGG Pathways, Gene Counts, and Significance",
    x = "-log10(P-value)",
    y = "KEGG Pathway",
    size = "Number of Genes"
  )



