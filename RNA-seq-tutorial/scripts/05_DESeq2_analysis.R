# ==============================================================================
# RNA-seq Analysis Pipeline - Step 5: Differential Expression Analysis
# ==============================================================================
# Author: Alessandro Mieli
# Purpose: Identify differentially expressed genes using DESeq2
# ==============================================================================

# Install required packages (run once)
# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("pheatmap")
# install.packages("ggplot2")

# Load libraries
library(DESeq2)
library(pheatmap)
library(ggplot2)

# ==============================================================================
# 1. Import count data
# ==============================================================================

# Set working directory to folder containing count files
setwd("/path/to/your/data")

# Import gene counts (adjust filename as needed)
# Assumes file has GeneID in first column and counts in remaining columns
counts <- read.delim("gene_counts.csv", header = TRUE, row.names = 1, sep = ',')

# Filter low-count genes (optional but recommended)
counts <- counts[which(rowSums(counts) > 50), ]

# ==============================================================================
# 2. Define experimental conditions
# ==============================================================================

# Example: 4 treatment groups with 3 biological replicates each
# Adjust according to your experimental design
condition <- factor(c("Control", "Treatment1", "Treatment2",
                      "Control", "Treatment1", "Treatment2",
                      "Control", "Treatment1", "Treatment2",))

# Create sample metadata
coldata <- data.frame(row.names = colnames(counts), condition)

# ==============================================================================
# 3. Create DESeq2 dataset and run analysis
# ==============================================================================

dds <- DESeqDataSetFromMatrix(countData = counts, 
                               colData = coldata, 
                               design = ~condition)

# Run differential expression analysis
dds <- DESeq(dds)

# ==============================================================================
# 4. Quality control: PCA plot
# ==============================================================================

# Variance stabilizing transformation for visualization
vsdata <- vst(dds, blind = FALSE)

# PCA plot to check sample clustering
plotPCA(vsdata, intgroup = "condition")

# Dispersion plot
plotDispEsts(dds)

# ==============================================================================
# 5. Extract differentially expressed genes
# ==============================================================================

# Pairwise comparison: Treatment1 vs Control
res <- results(dds, contrast = c("condition", "Treatment1", "Control"))

# Remove genes with missing values
sigs <- na.omit(res)

# Filter by significance (FDR < 0.05)
sigs <- sigs[sigs$padj < 0.05, ]

# View top differentially expressed genes
head(sigs[order(sigs$padj), ])

# ==============================================================================
# 6. Generate heatmap of top DEGs
# ==============================================================================

# Get top 50 most significant genes
top_genes <- rownames(sigs[order(sigs$padj), ])[1:50]

# Extract counts for these genes
top_counts <- counts[top_genes, ]

# Generate heatmap
pheatmap(top_counts, 
         scale = "row", 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         main = "Heatmap of Top 50 DEGs")

# ==============================================================================
# 7. Create volcano plot
# ==============================================================================

# Prepare data for plotting
volcano_data <- as.data.frame(res)
volcano_data$diffexpressed <- "NO"
volcano_data$diffexpressed[volcano_data$log2FoldChange > 0.6 & volcano_data$pvalue < 0.05] <- "UP"
volcano_data$diffexpressed[volcano_data$log2FoldChange < -0.6 & volcano_data$pvalue < 0.05] <- "DOWN"

# Generate volcano plot
ggplot(data = volcano_data, aes(x = log2FoldChange, y = -log10(pvalue), col = diffexpressed)) +
  geom_point() +
  scale_color_manual(values = c("NO" = "black", "UP" = "red", "DOWN" = "blue")) +
  theme_minimal() +
  geom_vline(xintercept = c(-0.6, 0.6), col = "red", linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  labs(title = "Volcano Plot: Treatment1 vs Control",
       x = "Log2 Fold Change",
       y = "-Log10 p-value")

# ==============================================================================
# 8. Export results
# ==============================================================================

# Save significant DEGs to CSV file
write.csv(as.data.frame(sigs), file = "DEGs_Treatment1_vs_Control.csv")

# Save normalized counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, file = "normalized_counts.csv")

print("Analysis complete!")
