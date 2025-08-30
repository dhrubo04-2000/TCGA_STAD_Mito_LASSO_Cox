# ==============================
# Random Survival Forest (RSF) for TCGA-STAD Survival Analysis
# Author: Dhrubo Chowdhury
# ==============================

# --- Load libraries ---
BiocManager::install("randomForestSRC")
library(randomForestSRC)
library(survival)
library(survminer)
library(timeROC)
library(dplyr)

# --- Prepare data ---
# Remove previous risk columns
cox_clean <- cox_filtered_clean %>% select(-risk_group, -risk_score)

# Define survival object
surv_obj <- Surv(cox_clean$OS_time, cox_clean$OS_event)

# Use only gene expression columns
gene_cols <- colnames(cox_clean)[4:ncol(cox_clean)]
gene_data <- cox_clean[, gene_cols]

# Combine for RSF input
rsf_input <- data.frame(OS_time = cox_clean$OS_time,
                        OS_event = cox_clean$OS_event,
                        gene_data)

# --- Fit RSF model ---
set.seed(42)
rsf_model <- rfsrc(Surv(OS_time, OS_event) ~ ., 
                   data = rsf_input,
                   ntree = 1000,
                   importance = TRUE)

# Model summary and error plot
print(rsf_model)
plot(rsf_model)

# --- Variable importance ---
rsf_importance <- sort(rsf_model$importance, decreasing = TRUE)
top_genes_df <- data.frame(Gene = names(rsf_importance),
                           Importance = rsf_importance)
top_20_genes <- head(top_genes_df, 20)
write.csv(top_20_genes, "RSF/top_20_genes_rsf.csv", row.names = FALSE)

# --- RSF-based risk score ---
cox_clean$rsf_risk <- rsf_model$predicted
cox_clean$rsf_risk_group <- ifelse(cox_clean$rsf_risk > median(cox_clean$rsf_risk),
                                   "High Risk", "Low Risk")

# Kaplan-Meier plot
fit_rsf <- survfit(Surv(OS_time, OS_event) ~ rsf_risk_group, data = cox_clean)
ggsurvplot(fit_rsf, data = cox_clean, pval = TRUE, risk.table = TRUE,
           title = "RSF-based Risk Group Survival",
           palette = c("red", "blue"))

# --- OOB Concordance Index ---
print(paste("RSF OOB Concordance:", round(rsf_model$concordance, 3)))

# --- Time-dependent ROC ---
roc_rsf <- timeROC(T = cox_clean$OS_time,
                   delta = cox_clean$OS_event,
                   marker = cox_clean$rsf_risk,
                   cause = 1,
                   times = c(365, 1095, 1825),
                   ROC = TRUE)
plot(roc_rsf, time = 365, col = "red", title = "RSF Time-dependent ROC")
plot(roc_rsf, time = 1095, col = "blue", add = TRUE)
plot(roc_rsf, time = 1825, col = "green", add = TRUE)
legend("bottomright", legend = c("1yr", "3yr", "5yr"),
       col = c("red", "blue", "green"), lwd = 2)

# --- RSF with top 15 genes ---
top15_genes <- c("PUS1","ABCA9","BIK","HPDL","POLQ","PDE2A","CYP11A1",
                 "MSRB3","TOMM40","ARMCX1","ATAD3A","PDK4","BNIP3",
                 "FKBP10","MICU3")
rsf_input_top15 <- data.frame(OS_time = cox_clean$OS_time,
                              OS_event = cox_clean$OS_event,
                              cox_clean[, top15_genes])
rsf_model_top15 <- rfsrc(Surv(OS_time, OS_event) ~ ., data = rsf_input_top15,
                         ntree = 1000, importance = TRUE)
cox_clean$rsf_top15_risk <- rsf_model_top15$predicted
cox_clean$rsf_top15_group <- ifelse(cox_clean$rsf_top15_risk > median(cox_clean$rsf_top15_risk),
                                    "High Risk", "Low Risk")

# KM plot for top 15 genes
fit_top15 <- survfit(Surv(OS_time, OS_event) ~ rsf_top15_group, data = cox_clean)
ggsurvplot(fit_top15, data = cox_clean, pval = TRUE, risk.table = TRUE,
           title = "RSF (Top 15 Genes) Risk Stratification",
           palette = c("red", "blue"))

# C-index for top 15 gene model
cindex_rsf <- concordance.index(x = cox_clean$rsf_top15_risk,
                                surv.time = cox_clean$OS_time,
                                surv.event = cox_clean$OS_event,
                                method = "noether")
print(paste("C-index:", round(cindex_rsf$c.index, 3)))
print(paste("Standard Error:", round(cindex_rsf$se, 3)))









