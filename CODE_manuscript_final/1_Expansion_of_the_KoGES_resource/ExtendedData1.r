#!/usr/bin/env Rscript
# ExtendedData1.r
#
# Purpose:
#   Baseline phenotype summary tables for KoGES (Extended Data 1):
#   categorical case/control counts and continuous descriptive statistics,
#   for (1) all participants and (2) genotyped participants (.fam IDs).
#
# Inputs:
#   Categorical:
#     DATA_manuscript_final/Extended_Data/raw/Merge_cate_cuminci.csv
#     (B_ID, B_Cohort; disease vars coded 1=control, 2=case)
#   Continuous (already QC'd: NC lung-function set to NA; LDL <= 0 set to NA):
#     DATA_manuscript_final/Extended_Data/raw/Merge_conti.csv
#   Genotyped IDs:
#     .../KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated.fam
#
# Outputs (analysis intermediates):
#   DATA_manuscript_final/ExtendedData/ExtendedData1_categorical.csv
#   DATA_manuscript_final/ExtendedData/ExtendedData1_continuous.csv
#
# Publication sheet (caption + category/unit layout):
#   After this script, run:
#     python3 CODE_manuscript_final/ExtendedData_export/export_ExtendedData_publication.py
#   → RESULT_manuscript_final/ExtendedData/Extended Data1.xlsx
#
# Usage:
#   Rscript ExtendedData1.r

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
cate_file <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/raw/Merge_cate_cuminci.csv"
conti_file <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/raw/Merge_conti.csv"
fam_file <- "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated.fam"
out_dir <- "/home2/sjim/openKoGES/DATA_manuscript_final/ExtendedData"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Load
# -----------------------------------------------------------------------------
cate_data <- read.csv(cate_file, stringsAsFactors = FALSE)
conti_data <- read.csv(conti_file, stringsAsFactors = FALSE)

fam_data <- read.table(fam_file, stringsAsFactors = FALSE, sep = "\t", header = FALSE)
fam_ids <- unique(fam_data[[1]])
cat("Loaded", length(fam_ids), "IDs from .fam\n")

# -----------------------------------------------------------------------------
# Target variables
# -----------------------------------------------------------------------------
target_vars_cate <- c(
  "B_SEX",
  "B_HTN2", "B_DM2", "B_LIP", "B_THY", "B_GT", "B_OSTE", "B_BPH", "B_MI",
  "B_ANG", "B_CAD", "B_CVA", "B_STR", "B_CLD", "B_TB", "B_COPD2", "B_BRON",
  "B_ASTH", "B_ULCER", "B_LIV", "B_FLIV", "B_STOMUL", "B_GASTRO", "B_POL",
  "B_GB", "B_CKD", "B_KD", "B_UB", "B_UT", "B_THYCA", "B_UTCA", "B_GCA",
  "B_COLCA", "B_LCA", "B_BRCA", "B_NOI", "B_FRAC", "B_GLAU", "B_ALLER",
  "B_CATA", "B_ARTH", "B_PER"
)

target_vars_conti <- c(
  "B_AGE",
  "B_RBC_B", "B_WBC_B", "B_MCV", "B_MCH", "B_MCHC", "B_PLAT", "B_HCT", "B_HB",
  "B_FEV1best", "B_FEV1pp", "B_FFRbest", "B_FVCbest", "B_FVCpp",
  "B_HDL", "B_LDL", "B_TG", "B_TCHL", "B_APOA", "B_APOB",
  "B_VITD", "B_IRON", "B_FOLATE", "B_ALBUMIN", "B_URICACID", "B_HOMOCYST",
  "B_HIP", "B_BMI", "B_WHR", "B_WAIST", "B_WEIGHT", "B_HEIGHT",
  "B_AST", "B_ALT", "B_ALP", "B_R_GTP", "B_T_BIL",
  "B_BUN", "B_eGFR", "B_CYSTATINC", "B_PHOSPHATE", "B_CREATINE",
  "B_HBA1C", "B_INSUL", "B_GLU0", "B_CPEPTIDE",
  "B_SBP", "B_DBP", "B_PULSE"
)

available_vars_cate <- intersect(target_vars_cate, names(cate_data))
available_vars_conti <- intersect(target_vars_conti, names(conti_data))

# -----------------------------------------------------------------------------
# Summary helpers (1=control, 2=case for categorical)
# -----------------------------------------------------------------------------
calculate_cate_summary <- function(data, var_list, group_label = "All") {
  rows <- lapply(var_list, function(var) {
    if (!var %in% names(data)) return(NULL)
    x <- data[[var]]
    control_n <- sum(x == 1, na.rm = TRUE)
    case_n <- sum(x == 2, na.rm = TRUE)
    data.frame(
      Var_name = var,
      Total_n = sum(!is.na(x)),
      Case_n = case_n,
      Control_n = control_n,
      Case_Control_Ratio = if (control_n > 0) round(case_n / control_n, 4) else NA_real_,
      NA_n = sum(is.na(x)),
      Group = group_label,
      stringsAsFactors = FALSE
    )
  })
  rbindlist(rows, fill = TRUE)
}

calculate_conti_summary <- function(data, var_list, group_label = "All") {
  rows <- lapply(var_list, function(var) {
    if (!var %in% names(data)) return(NULL)
    x <- data[[var]]
    x_ok <- x[!is.na(x)]
    data.frame(
      Var_name = var,
      Total_n = length(x_ok),
      Mean = if (length(x_ok)) round(mean(x_ok), 4) else NA_real_,
      SD = if (length(x_ok)) round(sd(x_ok), 4) else NA_real_,
      Min = if (length(x_ok)) round(min(x_ok), 4) else NA_real_,
      Max = if (length(x_ok)) round(max(x_ok), 4) else NA_real_,
      NA_n = sum(is.na(x)),
      Group = group_label,
      stringsAsFactors = FALSE
    )
  })
  rbindlist(rows, fill = TRUE)
}

# -----------------------------------------------------------------------------
# Summaries: All + genotyped (.fam)
# -----------------------------------------------------------------------------
cate_all <- calculate_cate_summary(cate_data, available_vars_cate, "All")
conti_all <- calculate_conti_summary(conti_data, available_vars_conti, "All")

cate_fam <- calculate_cate_summary(
  cate_data[cate_data$B_ID %in% fam_ids, ],
  available_vars_cate,
  "FAM_IDs"
)
conti_fam <- calculate_conti_summary(
  conti_data[conti_data$B_ID %in% fam_ids, ],
  available_vars_conti,
  "FAM_IDs"
)

cate_summary <- rbind(cate_all, cate_fam)
conti_summary <- rbind(conti_all, conti_fam)

setnames(
  cate_summary,
  c("Var_name", "Total(n)", "Case(n)", "Control(n)", "Case-Control Ratio", "NA(n)", "Group")
)
setnames(
  conti_summary,
  c("Var_name", "Total(n)", "Mean", "SD", "Min", "Max", "NA(n)", "Group")
)

# Drop SEX/AGE rows from published tables
cate_summary <- cate_summary[Group %in% c("All", "FAM_IDs") & Var_name != "B_SEX"]
conti_summary <- conti_summary[Group %in% c("All", "FAM_IDs") & Var_name != "B_AGE"]

# -----------------------------------------------------------------------------
# Phenotype display labels
# -----------------------------------------------------------------------------
pheno_map_cate <- data.table(
  Var_name = c(
    "B_HTN2", "B_DM2", "B_LIP", "B_THY", "B_GT", "B_OSTE", "B_BPH", "B_MI",
    "B_ANG", "B_CAD", "B_CVA", "B_STR", "B_CLD", "B_TB", "B_COPD2", "B_BRON",
    "B_ASTH", "B_ULCER", "B_LIV", "B_FLIV", "B_STOMUL", "B_GASTRO", "B_POL",
    "B_GB", "B_CKD", "B_KD", "B_UB", "B_UT", "B_THYCA", "B_UTCA", "B_GCA",
    "B_COLCA", "B_LCA", "B_BRCA", "B_NOI", "B_FRAC", "B_GLAU", "B_ALLER",
    "B_CATA", "B_ARTH", "B_PER"
  ),
  Phenotype_label = c(
    "Hypertension", "Diabetes mellitus", "Hyperlipidemia", "Thyroid disease",
    "Gout", "Osteoporosis", "Benign prostatic hyperplasia", "Myocardial infarction",
    "Angina pectoris", "Coronary artery disease", "Cerebrovascular disease", "Stroke",
    "Chronic lung disease", "Pulmonary tuberculosis",
    "Chronic obstructive pulmonary disease (COPD)", "Chronic bronchitis", "Asthma",
    "Peptic ulcer", "Acute liver disease", "Fatty liver",
    "Chronic gastritis / Peptic ulcer", "Chronic gastritis", "Intestinal polyp",
    "Gallbladder disease", "Chronic kidney disease (CKD)", "Kidney disease",
    "Cystitis", "Urinary tract infection", "Thyroid cancer", "Uterine cancer",
    "Gastric cancer", "Colorectal cancer", "Lung cancer", "Breast cancer",
    "Neurosis", "Fracture", "Glaucoma", "Allergic disease", "Cataract",
    "Arthritis", "Periodontal Disease"
  )
)

pheno_map_conti <- data.table(
  Var_name = c(
    "B_RBC_B", "B_WBC_B", "B_MCV", "B_MCH", "B_MCHC", "B_PLAT", "B_HCT", "B_HB",
    "B_FEV1best", "B_FEV1pp", "B_FFRbest", "B_FVCbest", "B_FVCpp",
    "B_HDL", "B_LDL", "B_TG", "B_TCHL", "B_APOA", "B_APOB",
    "B_VITD", "B_IRON", "B_FOLATE", "B_ALBUMIN", "B_URICACID", "B_HOMOCYST",
    "B_HIP", "B_BMI", "B_WHR", "B_WAIST", "B_WEIGHT", "B_HEIGHT",
    "B_AST", "B_ALT", "B_ALP", "B_R_GTP", "B_T_BIL",
    "B_BUN", "B_eGFR", "B_CYSTATINC", "B_PHOSPHATE", "B_CREATINE",
    "B_HBA1C", "B_INSUL", "B_GLU0", "B_CPEPTIDE",
    "B_SBP", "B_DBP", "B_PULSE"
  ),
  Phenotype_label = c(
    "Red blood cell count (RBC)", "White blood cell count (WBC)",
    "Mean corpuscular volume (MCV)", "Mean corpuscular hemoglobin (MCH)",
    "Mean corpuscular hemoglobin concentration (MCHC)",
    "Platelet count (Plt)", "Hematocrit (HCT)", "Hemoglobin (Hb)",
    "FEV1",
    "FEV1_%pred (Percentage of predicted forced expiratory volume in 1 second)",
    "FEV1/FVC", "FVC",
    "FVC_%pred (Percentage of predicted forced vital capacity)",
    "HDL-cholesterol (HDL-C)", "LDL-cholesterol (LDL-C)", "Triglycerides (TG)",
    "Total cholesterol (TC)", "Apolipoprotein A", "Apolipoprotein B",
    "Vit.D", "Iron", "Folate", "Albumin (Alb)", "Uric acids (UA)", "Homocysteine",
    "Hip", "Body mass index (BMI)", "Waist-to-hip ratio (WHR)", "Waist",
    "Body weight (BW)", "Body height",
    "Aspartate transaminase (AST)", "Alanine aminotransferase (ALT)",
    "Alkaline phosphatase (ALP)", "r_Gtp (Gamma glutamyl transferase)",
    "Total Bilirubin (TBil)", "Blood urea nitrogen (BUN)", "eGFR (CKD-EPI)",
    "Cystatin C", "Phosphate", "Serum creatinine (sCr)",
    "HbA1C", "Insulin", "Glucose", "C-peptide",
    "Systolic blood pressure (SBP)", "Diastolic blood pressure (DBP)", "Pulse"
  )
)

final_cate <- left_join(pheno_map_cate, as.data.frame(cate_summary), by = "Var_name") %>%
  arrange(Group)
final_conti <- left_join(pheno_map_conti, as.data.frame(conti_summary), by = "Var_name") %>%
  arrange(Group)

# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
out_cate <- file.path(out_dir, "ExtendedData1_categorical.csv")
out_conti <- file.path(out_dir, "ExtendedData1_continuous.csv")
write.csv(final_cate, out_cate, row.names = FALSE)
write.csv(final_conti, out_conti, row.names = FALSE)

cat("[DONE]", out_cate, "\n")
cat("[DONE]", out_conti, "\n")
