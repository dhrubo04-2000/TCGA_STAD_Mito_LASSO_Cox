# ==============================
# Multivariate Cox Regression Analysis
# TCGA-STAD: Gene Expression + Clinical Variables
# Author: Your Name
# Date: YYYY-MM-DD
# ==============================

# --- Load libraries ---
library(readr)
library(tidyverse)
library(survival)
library(gtsummary)
library(flextable)
library(officer)
library(broom)

# --- Set working directory ---
setwd("D:/Gastric_GSE")

# --- Load processed Cox dataset ---
cox_data <- read_csv("cox_ready/cox_filtered_normalized_DE_os_model_ready.csv")

# --- Data preprocessing ---
cox_data <- cox_data %>%
  filter(OS_time > 0) %>%                     # remove OS_time = 0
  filter(!is.na(diagnoses.ajcc_pathologic_stage) & !is.na(diagnoses.tumor_grade)) %>%
  mutate(
    age = suppressWarnings(as.numeric(demographic.age_at_index)),
    gender = as.factor(demographic.gender),
    stage = as.factor(diagnoses.ajcc_pathologic_stage),
    grade = as.factor(diagnoses.tumor_grade),
    FKBP10 = as.numeric(FKBP10)              # convert genes to numeric
  )

# --- Fit univariate Cox for a single gene (example: MICU3) ---
cox_model <- coxph(Surv(OS_time, OS_event) ~ MICU3 + age + gender + stage + grade,
                   data = cox_data)
summary(cox_model)

# --- Generate publication-ready table using gtsummary ---
cox_table <- tbl_regression(
  cox_model,
  exponentiate = TRUE,
  label = list(
    MICU3 ~ "MICU3 Expression",
    age ~ "Age at Diagnosis",
    gender ~ "Gender",
    stage ~ "Pathologic Stage",
    grade ~ "Tumor Grade"
  )
) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_significance_stars() %>%
  bold_p(t = 0.05)

# --- Export table to Word ---
ft <- as_flex_table(cox_table)
doc <- read_docx() %>%
  body_add_flextable(ft)
print(doc, target = "cox_model_table.docx")

# --- Multivariate Cox for multiple genes + clinical variables ---
gene_list <- c("POLQ", "PDK4", "FKBP10", "CYP11A1", "SPHKAP", 
               "PDE2A", "ACADL", "DMGDH", "BIK", "MICU3", 
               "ARMCX1", "BNIP3", "ABCA9", "MSRB3")

# Convert gene columns to numeric
cox_data[gene_list] <- lapply(cox_data[gene_list], as.numeric)

# Create formula string dynamically
formula_string <- paste("Surv(OS_time, OS_event) ~",
                        paste(gene_list, collapse = " + "),
                        "+ age + stage + grade + gender")

# Fit multivariate Cox model
cox_model_multi <- coxph(as.formula(formula_string), data = cox_data)
summary(cox_model_multi)

# Generate table for multiple genes
tbl_regression(cox_model_multi, exponentiate = TRUE) %>%
  bold_labels() %>%
  italicize_levels() %>%
  add_significance_stars() %>%
  bold_p(t = 0.05)







