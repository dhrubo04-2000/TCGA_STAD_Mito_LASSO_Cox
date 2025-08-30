# ==============================
# LASSO-Cox Regression for TCGA-STAD Survival Analysis
# Author: Your Name
# Date: YYYY-MM-DD
# ==============================

# --- Load libraries ---
library(readr)
library(tidyverse)
library(glmnet)
library(survival)
library(survminer)
library(survcomp)
library(timeROC)

# --- Set working directory and load data ---
setwd("D:/Gastric_GSE")
cox_data <- read_csv("cox_ready/cox_filtered_normalized_DE_os.csv")

# --- Preprocessing ---
# Remove patients with OS_time = 0
cox_data <- cox_data %>% filter(OS_time != 0)

# Remove irrelevant clinical columns
cox_clean <- cox_data %>% 
  select(-c(cases.submitter_id,
            demographic.days_to_death,
            diagnoses.days_to_last_follow_up,
            demographic.vital_status,
            demographic.age_at_index,
            demographic.gender,
            diagnoses.ajcc_pathologic_stage,
            diagnoses.tumor_grade))

# --- Prepare predictors and survival outcome ---
gene_cols <- colnames(cox_clean)[4:ncol(cox_clean)]
x <- as.matrix(cox_clean[, gene_cols])
y <- Surv(cox_clean$OS_time, cox_clean$OS_event)

# --- LASSO-Cox regression ---
set.seed(123)
cv_lasso <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10)

# Extract best lambda and coefficients
best_lambda <- cv_lasso$lambda.min
lasso_coef <- coef(cv_lasso, s = best_lambda)
selected_genes <- lasso_coef[lasso_coef != 0]
selected_genes_df <- data.frame(Gene = rownames(selected_genes), Coefficient = as.numeric(selected_genes))
write.csv(selected_genes_df, "Lasso_cox/selected_genes_lasso_cox.csv", row.names = FALSE)

# --- Compute risk scores ---
x_selected <- x[, selected_genes_df$Gene]
cox_clean$risk_score <- as.numeric(x_selected %*% selected_genes_df$Coefficient)

# --- Model performance ---
cindex <- concordance.index(
  x = cox_clean$risk_score,
  surv.time = cox_clean$OS_time,
  surv.event = cox_clean$OS_event,
  method = "noether"
)
print(paste("C-index:", round(cindex$c.index, 3)))

# --- Kaplan-Meier curve ---
cox_clean$risk_group <- ifelse(cox_clean$risk_score > median(cox_clean$risk_score), "High Risk", "Low Risk")
fit <- survfit(Surv(OS_time, OS_event) ~ risk_group, data = cox_clean)
ggsurvplot(fit, data = cox_clean, pval = TRUE, risk.table = TRUE,
           title = "Survival by LASSO-based Risk Score",
           palette = c("red", "blue"))

# --- Time-dependent ROC ---
roc_result <- timeROC(
  T = cox_clean$OS_time,
  delta = cox_clean$OS_event,
  marker = cox_clean$risk_score,
  cause = 1,
  times = c(365, 1095, 1825),
  ROC = TRUE
)
plot(roc_result, time = 365, col = "red", title = "Time-dependent ROC")
plot(roc_result, time = 1095, col = "blue", add = TRUE)
plot(roc_result, time = 1825, col = "green", add = TRUE)
legend("bottomright", legend = c("1-year", "3-year", "5-year"),
       col = c("red", "blue", "green"), lwd = 2)













