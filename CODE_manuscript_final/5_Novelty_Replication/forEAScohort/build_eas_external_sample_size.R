#!/usr/bin/env Rscript
################################################################################
# Build phenotype-level external N + EAS panel phenotype eligibility
#
# Source
#   Phenotype × cohort mapping table (KBA160K VarNM + BBJ2021 / TWB / TPMI N).
#   N values are hardcoded below from authoritative public / sumstat sources:
#     - BBJ: PheWeb.jp downloads (Sakaue & Kanai Nat Genet 2021); Cont sumstat
#       has no N column. MCH=128028, MCV=129832, GGT=133471.
#     - TWB: GCST*.tsv.gz-meta.yaml sample_size (= constant n in sumstat).
#       Most continuous = 92615; FEV1/FVC (GCST90278627) = 62901.
#     - TPMI: unique N column in /data6/TPMIsumstat/{trait}.saige.txt.gz
#       (constant across variants). HR -> B_PULSE N=141961.
#
# Outputs
#   1) eas_external_phenotype_N.tsv
#        output_stem-level analytic N for BBJ / TWB / TPMI
#   2) eas_panel_phenotype_eligibility.tsv
#        koges_pheno_code-level flag for whether panel meta is allowed
#
# Notes
#   - Do NOT use koges_n / NCHROBS_koges.
#   - output_stem follows EAS Cont file names:
#       B_GLU0 -> B_GLUO ; B_FFRbest -> B_FFR_best
#   - KoGES Pulse (heart rate) is NOT matched to BBJ/TPMI PP (pulse pressure).
#   - TPMI heart rate uses trait HR -> output_stem B_PULSE (N=141961).
################################################################################

suppressPackageStartupMessages(library(data.table))

out_n_path <- Sys.getenv(
  "EAS_EXTERNAL_N_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eas_external_phenotype_N.tsv"
)
out_elig_path <- Sys.getenv(
  "EAS_PANEL_ELIGIBILITY_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eas_panel_phenotype_eligibility.tsv"
)

# ------------------------------------------------------------------------------
# 1. Phenotype N table
# ------------------------------------------------------------------------------
pheno_n <- fread(text = "
clinical\tvarnm\tbbj_n\ttwb_n\ttpmi_n
Alanine aminotransferase (ALT)\tB_ALT\t150545\t92615\t203766
Albumin (Alb)\tB_ALBUMIN\t120539\t92615\tNA
Alkaline phosphatase (ALP)\tB_ALP\t118886\tNA\tNA
Apolipoprotein A\tB_APOA\tNA\tNA\tNA
Apolipoprotein B\tB_APOB\tNA\tNA\tNA
Aspartate transaminase (AST)\tB_AST\t150068\t92615\t152196
Blood urea nitrogen (BUN)\tB_BUN\t148767\t92615\t140894
Body height\tB_HEIGHT\t165056\t92615\t213895
Body mass index (BMI)\tB_BMI\t163835\t92615\t191458
Body weight (BW)\tB_WEIGHT\t165419\t92615\tNA
C-peptide\tB_CPEPTIDE\tNA\tNA\tNA
Cystatin C\tB_CYSTATINC\tNA\tNA\tNA
Diastolic blood pressure (DBP)\tB_DBP\t145515\t92615\t196215
eGFR\tB_eGFR\t154633\tNA\tNA
FEV1_best\tB_FEV1best\tNA\tNA\tNA
FEV1_pred\tB_FEV1pred\tNA\tNA\tNA
FFRbest (FEV1/FVC_BEST)\tB_FFRbest\tNA\t62901\tNA
FFRpred (FEV1/FVC_PRED)\tB_FFRpred\tNA\tNA\tNA
Folate\tB_FOLATE\tNA\tNA\tNA
FVC_best\tB_FVCbest\tNA\tNA\tNA
FVC_pred\tB_FVCpred\tNA\tNA\tNA
Glucose\tB_GLU0\t133336\t92615\t111635
HbA1C\tB_HBA1C\t71221\t92615\t90509
HDL-cholesterol (HDL-C)\tB_HDL\t74970\t92615\t142676
Hematocrit (HCT)\tB_HCT\t153015\t92615\t177395
Hemoglobin (Hb)\tB_HB\t152447\t92615\t194120
Hip\tB_HIP\tNA\t92615\tNA
Homocysteine\tB_HOMOCYST\tNA\tNA\tNA
Insulin\tB_INSUL\tNA\tNA\tNA
Iron\tB_IRON\tNA\tNA\tNA
LDL-cholesterol (LDL-C)\tB_LDL\t72866\t92615\t165202
Mean corpuscular hemoglobin (MCH)\tB_MCH\t128028\tNA\tNA
Mean corpuscular hemoglobin concentration (MCHC)\tB_MCHC\t135482\tNA\tNA
Mean corpuscular volume (MCV)\tB_MCV\t129832\tNA\tNA
Phosphatae\tB_PHOSPHATE\t46848\tNA\tNA
Platelet count (Plt)\tB_PLAT\t148623\t92615\t168102
Pulse (heart rate)\tB_PULSE\tNA\tNA\t141961
r_Gtp (Gamma glutamyl transferase)\tB_R_GTP\t133471\t92615\tNA
Red blood cell count (RBC)\tB_RBC_B\t153512\t92615\t173813
Serum creatinine (sCr)\tB_CREATINE\t150266\t92615\t217954
Sytolic blood pressure (SBP)\tB_SBP\t145505\t92615\t196273
Total Bilirubin (TBil)\tB_T_BIL\t124341\t92615\t91506
Total cholesterol (TC)\tB_TCHL\t135808\t92615\t161535
Triglycerides (TG)\tB_TG\t111667\t92615\t161062
Uric acids (UA)\tB_URICACID\t129405\t92615\t116933
Vit.D\tB_VITD\tNA\tNA\tNA
Waist\tB_WAIST\tNA\t92615\tNA
Waist-to-hip ratio (WHR)\tB_WHR\tNA\t92615\tNA
White blood cell count (WBC)\tB_WBC_B\t154355\t92615\t180365
", na.strings = c("NA", ""))

stem_map <- c(
  B_GLU0 = "B_GLUO",
  B_FFRbest = "B_FFR_best"
)

pheno_n[, output_stem := fifelse(
  varnm %chin% names(stem_map),
  unname(stem_map[varnm]),
  varnm
)]
pheno_n[, `:=`(
  bbj_n = as.numeric(bbj_n),
  twb_n = as.numeric(twb_n),
  tpmi_n = as.numeric(tpmi_n)
)]
pheno_n[, `:=`(
  bbj_n_source = fifelse(
    is.finite(bbj_n), "phenotype_map_table:BBJ2021", NA_character_
  ),
  twb_n_source = fifelse(
    is.finite(twb_n), "phenotype_map_table:TWB", NA_character_
  ),
  tpmi_n_source = fifelse(
    is.finite(tpmi_n), "phenotype_map_table:TPMI", NA_character_
  )
)]

if (anyDuplicated(pheno_n$output_stem)) {
  stop("Duplicated output_stem in phenotype N map table.")
}

out_n <- pheno_n[, .(
  output_stem,
  clinical,
  varnm,
  bbj_n,
  bbj_n_source,
  twb_n,
  twb_n_source,
  tpmi_n,
  tpmi_n_source
)]
setorder(out_n, output_stem)

# ------------------------------------------------------------------------------
# 2. koges_pheno_code eligibility for EAS panel meta
#    (bolt code -> external stems used in novelty scripts)
# ------------------------------------------------------------------------------
elig <- fread(text = "
koges_pheno_code\tbbj_output_stem\ttwb_output_stem\ttpmi_output_stem\teas_panel_meta_eligible\tmeta_exclusion_reason
B_FFRbest_noNC_BOLT_result\tNA\tB_FFR_best\tNA\tTRUE\t
RIN_B_ALBUMIN_BOLT_result\tB_ALBUMIN\tB_ALBUMIN\tNA\tTRUE\t
RIN_B_ALP_BOLT_result\tB_ALP\tNA\tNA\tTRUE\t
RIN_B_ALT_BOLT_result\tB_ALT\tB_ALT\tB_ALT\tTRUE\t
RIN_B_AST_BOLT_result\tB_AST\tB_AST\tB_AST\tTRUE\t
RIN_B_BMI_BOLT_result\tB_BMI\tB_BMI\tB_BMI\tTRUE\t
RIN_B_BUN_BOLT_result\tB_BUN\tB_BUN\tB_BUN\tTRUE\t
RIN_B_CREATINE_BOLT_result\tB_CREATINE\tB_CREATINE\tB_CREATINE\tTRUE\t
RIN_B_DBP_BOLT_result\tB_DBP\tB_DBP\tB_DBP\tTRUE\t
RIN_B_eGFR_BOLT_result\tB_eGFR\tNA\tNA\tTRUE\t
RIN_B_GLU0_BOLT_result\tB_GLUO\tB_GLUO\tB_GLUO\tTRUE\t
RIN_B_HB_BOLT_result\tB_HB\tB_HB\tB_HB\tTRUE\t
RIN_B_HBA1C_BOLT_result\tB_HBA1C\tB_HBA1C\tB_HBA1C\tTRUE\t
RIN_B_HCT_BOLT_result\tB_HCT\tB_HCT\tB_HCT\tTRUE\t
RIN_B_HDL_BOLT_result\tB_HDL\tB_HDL\tB_HDL\tTRUE\t
RIN_B_HEIGHT_BOLT_result\tB_HEIGHT\tB_HEIGHT\tB_HEIGHT\tTRUE\t
RIN_B_HIP_BOLT_result\tNA\tB_HIP\tNA\tTRUE\t
RIN_B_LDLclean_BOLT_result\tB_LDL\tB_LDL\tB_LDL\tTRUE\t
RIN_B_MCH_BOLT_result\tB_MCH\tNA\tNA\tTRUE\t
RIN_B_MCHC_BOLT_result\tB_MCHC\tNA\tNA\tTRUE\t
RIN_B_MCV_BOLT_result\tB_MCV\tNA\tNA\tTRUE\t
RIN_B_PHOSPHATE_BOLT_result\tB_PHOSPHATE\tNA\tNA\tTRUE\t
RIN_B_PLAT_BOLT_result\tB_PLAT\tB_PLAT\tB_PLAT\tTRUE\t
RIN_B_PULSE_BOLT_result\tNA\tNA\tB_PULSE\tTRUE\t
RIN_B_R_GTP_BOLT_result\tB_R_GTP\tB_R_GTP\tNA\tTRUE\t
RIN_B_RBC_B_BOLT_result\tB_RBC_B\tB_RBC_B\tB_RBC_B\tTRUE\t
RIN_B_SBP_BOLT_result\tB_SBP\tB_SBP\tB_SBP\tTRUE\t
RIN_B_T_BIL_BOLT_result\tB_T_BIL\tB_T_BIL\tB_T_BIL\tTRUE\t
RIN_B_TCHL_BOLT_result\tB_TCHL\tB_TCHL\tB_TCHL\tTRUE\t
RIN_B_TG_BOLT_result\tB_TG\tB_TG\tB_TG\tTRUE\t
RIN_B_URICACID_BOLT_result\tB_URICACID\tB_URICACID\tB_URICACID\tTRUE\t
RIN_B_WAIST_BOLT_result\tNA\tB_WAIST\tNA\tTRUE\t
RIN_B_WBC_B_BOLT_result\tB_WBC_B\tB_WBC_B\tB_WBC_B\tTRUE\t
RIN_B_WEIGHT_BOLT_result\tB_WEIGHT\tB_WEIGHT\tNA\tTRUE\t
RIN_B_WHR_BOLT_result\tNA\tB_WHR\tNA\tTRUE\t
", na.strings = c("NA", ""))

elig[, eas_panel_meta_eligible := as.logical(eas_panel_meta_eligible)]
if (anyDuplicated(elig$koges_pheno_code)) {
  stop("Duplicated koges_pheno_code in eligibility table.")
}

# Validate stems referenced by eligibility exist in N table (when non-NA)
for (col in c("bbj_output_stem", "twb_output_stem", "tpmi_output_stem")) {
  stems <- unique(na.omit(elig[[col]]))
  miss <- setdiff(stems, out_n$output_stem)
  if (length(miss) > 0L) {
    stop("Eligibility stem(s) missing from N table (", col, "): ",
         paste(miss, collapse = ", "))
  }
}

dir.create(dirname(out_n_path), recursive = TRUE, showWarnings = FALSE)
fwrite(out_n, out_n_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(elig, out_elig_path, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:", out_n_path, "\n")
cat("Wrote:", out_elig_path, "\n")
cat("N phenotypes:", nrow(out_n),
    "| eligible koges codes:", sum(elig$eas_panel_meta_eligible),
    "| excluded:", sum(!elig$eas_panel_meta_eligible), "\n")
if (any(!elig$eas_panel_meta_eligible)) {
  cat("Excluded phenotypes:\n")
  print(elig[eas_panel_meta_eligible == FALSE,
             .(koges_pheno_code, meta_exclusion_reason)])
}
