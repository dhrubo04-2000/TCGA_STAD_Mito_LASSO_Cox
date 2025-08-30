# ==============================
# Differential Expression Analysis of TCGA-STAD (Mitochondrial Genes)
# Author: Your Name
# Date: YYYY-MM-DD
# ==============================

# --- Load packages ---
library(readr)
library(tidyverse)
library(DESeq2)
library(EnhancedVolcano)

# --- Set working directory ---
setwd("D:/Gastric_GSE")

# --- Load data ---
STAD_mito_normal <- read_csv("tcga_stad/STAD_mito_normal.csv")
STAD_mito_tumor  <- read_csv("tcga_stad/STAD_mito_tumor_filtered_os.csv")

# --- Preprocessing ---
# Set rownames and remove gene annotation columns
STAD_mito_normal <- STAD_mito_normal %>% column_to_rownames("Genes") %>% select(-any_of(c("SYMBOL", "Symbol")))
STAD_mito_tumor  <- STAD_mito_tumor  %>% column_to_rownames("Genes") %>% select(-any_of(c("SYMBOL", "Symbol")))

# Filter low-expression genes (sum > 10)
STAD_mito_normal <- STAD_mito_normal[rowSums(STAD_mito_normal) > 10, ]
STAD_mito_tumor  <- STAD_mito_tumor[rowSums(STAD_mito_tumor) > 10, ]

# Keep only common genes between normal and tumor
common_genes <- intersect(rownames(STAD_mito_normal), rownames(STAD_mito_tumor))
STAD_mito_normal <- STAD_mito_normal[common_genes, ]
STAD_mito_tumor  <- STAD_mito_tumor[common_genes, ]

# Combine into one matrix
all_counts_unpaired <- cbind(STAD_mito_normal, STAD_mito_tumor)
write.csv(all_counts_unpaired, "tcga_stad/STAD_mito_all_counts_unpaired_os.csv", row.names = TRUE)

# --- Metadata preparation (colData) ---
# Extract batch info from barcodes (6th field)
extract_field <- function(barcodes, index) {
  sapply(strsplit(barcodes, "-"), function(x) x[index])
}

barcodes <- colnames(all_counts_unpaired)
condition <- factor(c(rep("Normal", ncol(STAD_mito_normal)),
                      rep("Tumor", ncol(STAD_mito_tumor))))

coldata_batch <- data.frame(
  row.names = barcodes,
  condition = condition,
  batch     = factor(extract_field(barcodes, 6)),
  TSS       = extract_field(barcodes, 2),
  plate     = extract_field(barcodes, 5),
  vial      = extract_field(barcodes, 7)
)

# --- Differential expression analysis ---
dds_batch <- DESeqDataSetFromMatrix(
  countData = all_counts_unpaired,
  colData   = coldata_batch,
  design    = ~ batch + condition
)

dds_batch <- DESeq(dds_batch)
res_batch <- results(dds_batch, contrast = c("condition", "Tumor", "Normal"))

# Significant DEGs
resOrdered <- res_batch[order(res_batch$padj), ]
resSig     <- subset(resOrdered, padj < 0.05 & abs(log2FoldChange) > 1)

# Save results
write.csv(as.data.frame(resSig), "tcga_stad/STAD_mito_DE_os.csv", row.names = TRUE)

# --- Visualization ---
EnhancedVolcano(
  resOrdered,
  lab = rownames(resOrdered),
  x = "log2FoldChange",
  y = "padj",
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "DEGs Volcano Plot",
  subtitle = "TCGA STAD Mitochondrial Genes",
  pointSize = 3.0,
  labSize = 3.0,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  colConnectors = "black",
  colAlpha = 0.5
)

# --- PCA plots for batch effect check ---
vsd <- vst(dds_batch, blind = FALSE)
plotPCA(vsd, intgroup = "condition")
plotPCA(vsd, intgroup = "batch")
plotPCA(vsd, intgroup = "TSS")
plotPCA(vsd, intgroup = "plate")

# ==============================
# End of Script
# ==============================



