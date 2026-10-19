#!/usr/bin/env Rscript
# ============================================================
# count_gws_variants.R
#
# Purpose:
#   Count genome-wide significant variants (P < 5e-8) per phenotype
#   for continuous (BOLT-LMM) and binary (SAIGE) GWAS sumstats.
#   This counts all SNPs below threshold (not LD-clumped leads).
#   Results are printed to the console only (no output files).
#
# Inputs (both under the same sumstat directory):
#   Continuous:
#     .../Result/GWAS_sumstat/*_BOLT_result
#     P column: P_BOLT_LMM if present, else P_BOLT_LMM_INF
#   Binary:
#     .../Result/GWAS_sumstat/*_SAIGE_cumulative_inci.SAIGE.SEfixed.txt
#     P column: p.value
#
# Usage:
#   Rscript count_gws_variants.R
#   NCORES=8 Rscript count_gws_variants.R
#   P_THR=5e-8 GROUP=cont Rscript count_gws_variants.R
#   GROUP=bin Rscript count_gws_variants.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
})

###############################################################################
# Settings
###############################################################################

P_THR <- as.numeric(Sys.getenv("P_THR", "5e-8"))
if (!is.finite(P_THR) || P_THR <= 0 || P_THR > 1) {
  stop("Invalid P_THR: ", Sys.getenv("P_THR"), call. = FALSE)
}

GROUP <- tolower(trimws(Sys.getenv("GROUP", "all")))
if (!GROUP %in% c("all", "cont", "continuous", "bin", "binary", "cate")) {
  stop(
    "Invalid GROUP: ", GROUP,
    ". Use all | cont | bin",
    call. = FALSE
  )
}

n_cores <- suppressWarnings(as.integer(Sys.getenv("NCORES", unset = NA)))
if (is.na(n_cores) || n_cores < 1L) {
  n_cores <- max(1L, parallel::detectCores() - 1L)
}

SUMSTAT_DIR <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"

###############################################################################
# Phenotype display labels (same manuscript set as clump scripts)
###############################################################################

cont_label <- data.table(
  file_stem = c(
    "B_FEV1best_noNC", "B_FEV1pp_noNC", "B_FFRbest_noNC",
    "B_FVCbest_noNC", "B_FVCpp_noNC",
    "RIN_B_ALBUMIN", "RIN_B_ALP", "RIN_B_ALT", "RIN_B_APOA", "RIN_B_APOB",
    "RIN_B_AST", "RIN_B_BMI", "RIN_B_BUN", "RIN_B_CPEPTIDE", "RIN_B_CREATINE",
    "RIN_B_CYSTATINC", "RIN_B_DBP", "RIN_B_eGFR", "RIN_B_FOLATE", "RIN_B_GLU0",
    "RIN_B_HB", "RIN_B_HBA1C", "RIN_B_HCT", "RIN_B_HDL", "RIN_B_HEIGHT",
    "RIN_B_HIP", "RIN_B_HOMOCYST", "RIN_B_INSUL", "RIN_B_IRON",
    "RIN_B_LDLclean", "RIN_B_MCH", "RIN_B_MCHC", "RIN_B_MCV",
    "RIN_B_PHOSPHATE", "RIN_B_PLAT", "RIN_B_PULSE", "RIN_B_R_GTP",
    "RIN_B_RBC_B", "RIN_B_SBP", "RIN_B_T_BIL", "RIN_B_TCHL", "RIN_B_TG",
    "RIN_B_URICACID", "RIN_B_VITD", "RIN_B_WAIST", "RIN_B_WBC_B",
    "RIN_B_WEIGHT", "RIN_B_WHR",
    "RIN_B_FEV1best_noNC", "RIN_B_FEV1pp_noNC", "RIN_B_FFRbest_noNC",
    "RIN_B_FVCbest_noNC", "RIN_B_FVCpp_noNC"
  ),
  Phenotype_label = c(
    "FEV1", "FEV1 %pred", "FEV1/FVC", "FVC", "FVC %pred",
    "Albumin (Alb)", "Alkaline phosphatase (ALP)",
    "Alanine aminotransferase (ALT)", "Apolipoprotein A", "Apolipoprotein B",
    "Aspartate transaminase (AST)", "Body mass index (BMI)",
    "Blood urea nitrogen (BUN)", "C-peptide", "Serum creatinine (sCr)",
    "Cystatin C", "Diastolic blood pressure (DBP)", "eGFR (CKD-EPI)",
    "Folate", "Glucose", "Hemoglobin (Hb)", "HbA1C", "Hematocrit (HCT)",
    "HDL-cholesterol (HDL-C)", "Body height", "Hip", "Homocysteine",
    "Insulin", "Iron", "LDL-cholesterol (LDL-C)",
    "Mean corpuscular hemoglobin (MCH)",
    "Mean corpuscular hemoglobin concentration (MCHC)",
    "Mean corpuscular volume (MCV)", "Phosphate", "Platelet count (Plt)",
    "Pulse", "r_Gtp (Gamma glutamyl transferase)",
    "Red blood cell count (RBC)", "Systolic blood pressure (SBP)",
    "Total Bilirubin (TBil)", "Total cholesterol (TC)", "Triglycerides (TG)",
    "Uric acids (UA)", "Vit.D", "Waist", "White blood cell count (WBC)",
    "Body weight (BW)", "Waist-to-hip ratio (WHR)",
    "FEV1", "FEV1 %pred", "FEV1/FVC", "FVC", "FVC %pred"
  )
)

bin_label <- data.table(
  Pheno_code = c(
    "B_ALLER", "B_BRCA", "B_CKD", "B_CAD", "B_DM2", "B_HTN2",
    "B_LIP", "B_GT", "B_FLIV", "B_THY", "B_THYCA", "B_ULCER",
    "B_FRAC", "B_GCA", "B_NOI", "B_OSTE", "B_STOMUL", "B_CATA",
    "B_CLD", "B_CVA", "B_GASTRO", "B_GB", "B_PER", "B_UTCA",
    "B_ANG", "B_ARTH", "B_ASTH", "B_BPH", "B_BRON", "B_COLCA",
    "B_COPD2", "B_GLAU", "B_KD", "B_LCA", "B_LIV", "B_MI",
    "B_POL", "B_STR", "B_TB", "B_UB", "B_UT"
  ),
  Phenotype_label = c(
    "Allergic disease", "Breast cancer", "Chronic kidney disease (CKD)",
    "Coronary artery disease", "Diabetes mellitus", "Hypertension",
    "Hyperlipidemia", "Gout", "Fatty liver", "Thyroid disease",
    "Thyroid cancer", "Peptic ulcer", "Fracture", "Gastric cancer",
    "Neurosis", "Osteoporosis", "Chronic gastritis / Peptic ulcer",
    "Cataract", "Chronic lung disease", "Cerebrovascular disease",
    "Chronic gastritis", "Gallbladder disease", "Periodontal Disease",
    "Uterine cancer",
    "Angina", "Arthritis", "Asthma", "Benign prostatic hyperplasia",
    "Bronchitis", "Colorectal cancer", "COPD", "Glaucoma",
    "Kidney disease", "Lung cancer", "Liver disease", "Myocardial infarction",
    "Polyp", "Stroke", "Tuberculosis", "Urinary bladder disease",
    "Urinary tract disease"
  )
)

###############################################################################
# Helpers
###############################################################################

read_header <- function(path) {
  line <- readLines(path, n = 1L, warn = FALSE)
  if (length(line) == 0L || !nzchar(line)) {
    stop("Could not read header from: ", path, call. = FALSE)
  }
  strsplit(line, "\t", fixed = TRUE)[[1]]
}

resolve_bolt_p_field <- function(cols, phenotype = NULL) {
  cols <- as.character(cols)
  if ("P_BOLT_LMM" %in% cols) return("P_BOLT_LMM")
  if ("P_BOLT_LMM_INF" %in% cols) return("P_BOLT_LMM_INF")
  stop(
    "No usable BOLT P column found",
    if (!is.null(phenotype)) paste0(" in ", phenotype) else "",
    ". Expected P_BOLT_LMM or P_BOLT_LMM_INF.",
    call. = FALSE
  )
}

count_p_lt <- function(path, p_field, p_thr) {
  dt <- fread(path, select = p_field, showProgress = FALSE)
  p <- suppressWarnings(as.numeric(dt[[1L]]))
  list(
    n_total = length(p),
    n_non_missing_p = sum(is.finite(p)),
    n_gws = sum(is.finite(p) & p < p_thr)
  )
}

run_one <- function(job, p_thr) {
  t0 <- proc.time()[["elapsed"]]
  out <- tryCatch(
    {
      cnt <- count_p_lt(job$path, job$p_field, p_thr)
      data.table(
        group = job$group,
        Phenotype = job$phenotype,
        Phenotype_label = job$label,
        file = basename(job$path),
        p_field = job$p_field,
        n_total_variants = as.integer(cnt$n_total),
        n_non_missing_p = as.integer(cnt$n_non_missing_p),
        n_gws_variants = as.integer(cnt$n_gws),
        status = "ok",
        message = NA_character_,
        elapsed_sec = round(proc.time()[["elapsed"]] - t0, 2)
      )
    },
    error = function(e) {
      data.table(
        group = job$group,
        Phenotype = job$phenotype,
        Phenotype_label = job$label,
        file = basename(job$path),
        p_field = job$p_field,
        n_total_variants = NA_integer_,
        n_non_missing_p = NA_integer_,
        n_gws_variants = NA_integer_,
        status = "error",
        message = conditionMessage(e),
        elapsed_sec = round(proc.time()[["elapsed"]] - t0, 2)
      )
    }
  )
  cat(
    sprintf(
      "[%s] %s  n_gws=%s  n_total=%s  p=%s  status=%s\n",
      job$group,
      job$phenotype,
      ifelse(is.na(out$n_gws_variants), "NA", format(out$n_gws_variants, big.mark = ",")),
      ifelse(is.na(out$n_total_variants), "NA", format(out$n_total_variants, big.mark = ",")),
      out$p_field,
      out$status
    )
  )
  out
}

build_cont_jobs <- function() {
  files <- list.files(SUMSTAT_DIR, pattern = "_BOLT_result$", full.names = TRUE)
  if (length(files) == 0L) {
    stop("No *_BOLT_result files found in: ", SUMSTAT_DIR, call. = FALSE)
  }

  rbindlist(lapply(files, function(path) {
    fname <- basename(path)
    stem <- sub("_BOLT_result$", "", fname)
    header <- read_header(path)
    p_field <- resolve_bolt_p_field(header, fname)
    label <- cont_label[file_stem == stem, Phenotype_label]
    if (length(label) == 0L) label <- NA_character_ else label <- label[1L]
    data.table(
      group = "continuous",
      phenotype = stem,
      label = label,
      path = path,
      p_field = p_field
    )
  }))
}

build_bin_jobs <- function() {
  files <- list.files(
    SUMSTAT_DIR,
    pattern = "_SAIGE_cumulative_inci\\.SAIGE\\.SEfixed\\.txt$",
    full.names = TRUE
  )
  if (length(files) == 0L) {
    stop(
      "No *_SAIGE_cumulative_inci.SAIGE.SEfixed.txt files found in: ",
      SUMSTAT_DIR,
      call. = FALSE
    )
  }

  rbindlist(lapply(files, function(path) {
    fname <- basename(path)
    code <- sub("_SAIGE_cumulative_inci\\.SAIGE\\.SEfixed\\.txt$", "", fname)
    header <- read_header(path)
    if (!"p.value" %in% header) {
      stop("Missing p.value column in: ", fname, call. = FALSE)
    }
    label <- bin_label[Pheno_code == code, Phenotype_label]
    if (length(label) == 0L) label <- NA_character_ else label <- label[1L]
    data.table(
      group = "binary",
      phenotype = code,
      label = label,
      path = path,
      p_field = "p.value"
    )
  }))
}

log_group_summary <- function(res, group_name) {
  cat("\n===== ", group_name, " summary =====\n", sep = "")
  print(res[, .(
    Phenotype,
    Phenotype_label,
    n_gws_variants,
    n_total_variants,
    p_field,
    status
  )])
  cat(
    group_name, " phenotypes: ", nrow(res),
    " | with GWS: ",
    sum(res$status == "ok" & res$n_gws_variants > 0L, na.rm = TRUE),
    " | sum n_gws: ",
    format(sum(res$n_gws_variants, na.rm = TRUE), big.mark = ","),
    "\n",
    sep = ""
  )
}

###############################################################################
# Main
###############################################################################

cat("P threshold:", format(P_THR, scientific = TRUE), "\n")
cat("Cores:", n_cores, "\n")
cat("Sumstat dir:", SUMSTAT_DIR, "\n")
cat("No output files will be written (console log only).\n")

run_cont <- GROUP %in% c("all", "cont", "continuous")
run_bin  <- GROUP %in% c("all", "bin", "binary", "cate")

if (run_cont) {
  jobs <- build_cont_jobs()
  setorder(jobs, phenotype)
  cat("Continuous phenotypes:", nrow(jobs), "\n")
  res_list <- mclapply(
    seq_len(nrow(jobs)),
    function(i) run_one(jobs[i], P_THR),
    mc.cores = n_cores
  )
  res_cont <- rbindlist(res_list, fill = TRUE)
  setorder(res_cont, -n_gws_variants, Phenotype)
  log_group_summary(res_cont, "continuous")
}

if (run_bin) {
  jobs <- build_bin_jobs()
  setorder(jobs, phenotype)
  cat("Binary phenotypes:", nrow(jobs), "\n")
  res_list <- mclapply(
    seq_len(nrow(jobs)),
    function(i) run_one(jobs[i], P_THR),
    mc.cores = n_cores
  )
  res_bin <- rbindlist(res_list, fill = TRUE)
  setorder(res_bin, -n_gws_variants, Phenotype)
  log_group_summary(res_bin, "binary")
}

cat("Done.\n")
