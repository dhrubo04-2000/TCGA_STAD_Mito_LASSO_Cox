# ==============================
# Mitochondrial Genes and Immune Cell Correlation Analysis
# TCGA-STAD
# Author: Dhrubo Chowdhury
# ==============================

# --- Load libraries ---
library(readr)
library(tidyverse)
library(DESeq2)
library(pheatmap)
library(reshape2)

# --- Set working directory ---
setwd("D:/Gastric_GSE")

# --- Load tumor mitochondrial gene expression ---
STAD_mito_tumor <- read_csv("tcga_stad/STAD_mito_tumor.csv") %>%
  select(-Symbol, -SYMBOL)

# --- Load DEGs ---
STAD_mito_DE_os <- read_csv("tcga_stad/STAD_mito_DE_os.csv")
colnames(STAD_mito_DE_os)[1] <- "Genes"
rownames(STAD_mito_DE_os) <- STAD_mito_DE_os$Genes
rownames(STAD_mito_tumor) <- STAD_mito_tumor$Genes

# Filter tumor matrix for DEGs
matching_genes <- intersect(rownames(STAD_mito_tumor), rownames(STAD_mito_DE_os))
STAD_mito_tumor_filtered <- STAD_mito_tumor[matching_genes, ]
STAD_mito_tumor_filtered <- STAD_mito_tumor_filtered %>%
  remove_rownames() %>% column_to_rownames("Genes")
write.csv(STAD_mito_tumor_filtered, "tcga_stad/STAD_mito_matrix_DEGs.csv", row.names = TRUE)

# --- VST normalization ---
count_matrix <- as.matrix(STAD_mito_tumor_filtered)
mode(count_matrix) <- "numeric"
sample_info <- data.frame(row.names = colnames(count_matrix),
                          condition = rep("Tumor", ncol(count_matrix)))
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~1)
vst_data <- varianceStabilizingTransformation(dds, blind = TRUE)
vst_expr_matrix <- assay(vst_data)
write.csv(vst_expr_matrix, "tcga_stad/STAD_mito_vstDEG_normalized.csv", row.names = TRUE)

# --- Load CIBERSORT output ---
cibersort_output <- read_csv("CIBERSORT_stad/cibersort_output_STAD_final.csv")
colnames(cibersort_output)[1] <- "sample_id"

# Keep only samples with P-value < 0.05
cibersort_significant <- cibersort_output %>% filter(`P-value` < 0.05)
rownames(cibersort_significant) <- cibersort_significant$sample_id
cibersort_significant <- cibersort_significant %>% select(-sample_id)

# --- Match samples ---
common_samples <- intersect(colnames(vst_expr_matrix), rownames(cibersort_significant))
expr_mat_common <- vst_expr_matrix[, common_samples]
ciber_common <- cibersort_significant[common_samples, ]

# --- Spearman correlation ---
cor_results <- data.frame()
for (gene in rownames(expr_mat_common)) {
  for (immune_cell in colnames(ciber_common)[1:22]) {
    cor_test <- cor.test(as.numeric(expr_mat_common[gene, ]),
                         as.numeric(ciber_common[, immune_cell]),
                         method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      Immune_Cell = immune_cell,
      Spearman_r = cor_test$estimate,
      p_value = cor_test$p.value
    ))
  }
}
cor_results$FDR <- p.adjust(cor_results$p_value, method = "BH")
sig_cor <- cor_results %>% filter(FDR < 0.05)

# Save results
write.csv(cor_results, "tcga_stad/STAD_mito_immune_correlation_results.csv", row.names = FALSE)
write.csv(sig_cor, "tcga_stad/STAD_mito_immune_correlation_significant.csv", row.names = FALSE)

# --- Heatmap of significant correlations ---
heatmap_matrix <- reshape2::acast(sig_cor, Gene ~ Immune_Cell, value.var = "Spearman_r")
heatmap_matrix[is.na(heatmap_matrix)] <- 0

pheatmap(heatmap_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         clustering_method = "complete",
         scale = "none",
         main = "Significant Spearman Correlation: mt Genes vs Immune Cells")




















