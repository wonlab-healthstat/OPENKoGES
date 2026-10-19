#!/usr/bin/env Rscript
# SAIGE.R
#
# Purpose:
#   Genome-wide association analysis of cumulative-incidence binary traits
#   in KoGES using the legacy SAIGE package:
#   (1) Step1  fitNULLGLMM   — null model + variance ratio
#   (2) Step2  SPAGMMATtest  — single-variant association
#   (3) SE fix               — correct SE for SPA-applied variants
#                              SE = |BETA| / |qnorm(p.value / 2)|
#
#   Covariate sets (auto by phenotype, or --mode=...):
#     default       : AGE, AGE^2, AGE×SEX, AGE^2×SEX, SEX, COHORT, PC1–20
#     unique_sex    : drop SEX terms (B_BPH, B_BRCA, B_UTCA)
#     unique_cohort : drop COHORT     (B_PER, B_STR, B_TB, B_ANG)
#
# Inputs:
#   Phenotype + covariates:
#     .../Phenotype/DATA/SAIGE_cumulative_incidence/
#       SAIGE_withcovar_cumulative_inci_<PHENO>.txt
#   Genotypes (PLINK):
#     .../Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated.{bed,bim,fam}
#   Sparse GRM:
#     .../Result/SAIGE/SAIGE_grm_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx[.sampleIDs.txt]
#
# Outputs:
#   Step1/2:
#     .../Result/SAIGE_cumulative_incidence/<PHENO>_SAIGE_cumulative_inci{.rda,.varianceRatio.txt,.SAIGE.txt}
#   SE-fixed sumstats:
#     .../Result/GWAS_sumstat/<PHENO>_SAIGE_cumulative_inci.SAIGE.SEfixed.txt
#     .../Result/GWAS_sumstat/SEfixed_QC/
#


library(xts)
library(data.table)
library(SAIGE)

cat("SAIGE ", as.character(packageVersion("SAIGE")), "\n", sep = "")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
name.grm   <- "/data6/OpenGWAS/KBA160K/Result/SAIGE/SAIGE_grm"
name.bfile <- "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated"
pheno.dir  <- "/data6/OpenGWAS/KBA160K_longitudinal/Phenotype/DATA/SAIGE_cumulative_incidence"
saige.dir  <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/SAIGE_cumulative_incidence"
sumstat.dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"
qc.dir      <- file.path(sumstat.dir, "SEfixed_QC")

dir.create(saige.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(sumstat.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc.dir, recursive = TRUE, showWarnings = FALSE)

phenoCol <- "B_disease_binary"

# -----------------------------------------------------------------------------
# Covariate modes
# -----------------------------------------------------------------------------
pheno_unique_sex    <- c("B_BPH", "B_BRCA", "B_UTCA")
pheno_unique_cohort <- c("B_PER", "B_STR", "B_TB", "B_ANG")

covar_sets <- list(
  default = list(
    covarColList = c(
      "B_AGE", "B_AGE_SQ", "B_AGE_SEX", "B_AGE_SQ_SEX",
      "B_SEX_BINARY", "COHORT", paste0("PC", 1:20)
    ),
    qCovarList = c("B_SEX_BINARY", "COHORT")
  ),
  unique_sex = list(
    covarColList = c("B_AGE", "B_AGE_SQ", "COHORT", paste0("PC", 1:20)),
    qCovarList = c("COHORT")
  ),
  unique_cohort = list(
    covarColList = c(
      "B_AGE", "B_AGE_SQ", "B_AGE_SEX", "B_AGE_SQ_SEX",
      "B_SEX_BINARY", paste0("PC", 1:20)
    ),
    qCovarList = c("B_SEX_BINARY")
  )
)

infer_mode <- function(pheno_name) {
  if (pheno_name %in% pheno_unique_sex) return("unique_sex")
  if (pheno_name %in% pheno_unique_cohort) return("unique_cohort")
  "default"
}

get_covar_for_pheno <- function(pheno_name, mode_override = "auto") {
  mode <- if (identical(mode_override, "auto")) infer_mode(pheno_name) else mode_override
  if (!mode %in% names(covar_sets)) {
    stop("Unknown --mode: ", mode)
  }
  list(mode = mode, covars = covar_sets[[mode]])
}

pheno_list_binary_default <- c("B_DM2")

# -----------------------------------------------------------------------------
# CLI: --stage --njob --mode --force ; remaining args = phenotypes
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
njob <- 1L
mode_override <- "auto"
stage <- "all"
force_overwrite <- FALSE
pheno_list_binary <- pheno_list_binary_default

parse_flag_eq <- function(args, name) {
  idx <- grep(paste0("^--", name, "="), args)
  if (length(idx) == 0) return(list(args = args, value = NULL))
  list(args = args[-idx[1]], value = sub(paste0("^--", name, "="), "", args[idx[1]]))
}

parse_flag_space <- function(args, name) {
  idx <- which(args == paste0("--", name))
  if (length(idx) == 0) return(list(args = args, value = NULL))
  if (idx[1] == length(args)) stop("--", name, " requires a value")
  list(args = args[-c(idx[1], idx[1] + 1)], value = args[idx[1] + 1])
}

if (length(args) > 0) {
  force_idx <- which(args %in% c("--force", "-f"))
  if (length(force_idx) > 0) {
    force_overwrite <- TRUE
    args <- args[-force_idx]
  }

  p <- parse_flag_eq(args, "njob"); args <- p$args
  if (is.null(p$value)) {
    p <- parse_flag_space(args, "njob"); args <- p$args
  }
  if (!is.null(p$value)) {
    njob_val <- suppressWarnings(as.integer(p$value))
    if (is.na(njob_val) || njob_val < 1) stop("Invalid --njob: ", p$value)
    njob <- njob_val
  }

  p <- parse_flag_eq(args, "mode"); args <- p$args
  if (is.null(p$value)) {
    p <- parse_flag_space(args, "mode"); args <- p$args
  }
  if (!is.null(p$value)) mode_override <- p$value

  p <- parse_flag_eq(args, "stage"); args <- p$args
  if (is.null(p$value)) {
    p <- parse_flag_space(args, "stage"); args <- p$args
  }
  if (!is.null(p$value)) stage <- p$value

  if (length(args) > 0) {
    pheno_list_binary <- unique(trimws(args))
    pheno_list_binary <- pheno_list_binary[nzchar(pheno_list_binary)]
  }
}

if (!mode_override %in% c("auto", names(covar_sets))) {
  stop("Invalid --mode: ", mode_override)
}
if (!stage %in% c("step1", "step2", "sefix", "all")) {
  stop("Invalid --stage: ", stage)
}

mode_plan <- vapply(
  pheno_list_binary,
  function(p) get_covar_for_pheno(p, mode_override)$mode,
  character(1)
)
cat("stage:", stage, " njob:", njob, " mode:", mode_override, "\n")
cat("phenotypes:", paste(pheno_list_binary, collapse = ", "), "\n")
for (i in seq_along(pheno_list_binary)) {
  cat("  ", pheno_list_binary[i], " -> ", mode_plan[i], "\n", sep = "")
}

run_parallel <- function(x, FUN) {
  if (njob <= 1L) {
    lapply(x, FUN)
  } else {
    parallel::mclapply(x, FUN, mc.cores = njob, mc.preschedule = FALSE)
  }
}

# -----------------------------------------------------------------------------
# Step1: fitNULLGLMM
# -----------------------------------------------------------------------------
run_step1 <- function(pheno.iter) {
  cov <- get_covar_for_pheno(pheno.iter, mode_override)
  mode <- cov$mode
  covarColList <- cov$covars$covarColList
  qCovarList   <- cov$covars$qCovarList

  cat("\n[Step1]", pheno.iter, "mode=", mode, "\n")

  pheno.file <- file.path(
    pheno.dir,
    paste0("SAIGE_withcovar_cumulative_inci_", pheno.iter, ".txt")
  )
  name.out <- file.path(saige.dir, paste0(pheno.iter, "_SAIGE_cumulative_inci"))

  if (!file.exists(pheno.file)) {
    cat("[SKIP] missing phenotype file:", pheno.file, "\n")
    return(FALSE)
  }

  tryCatch({
    SAIGE::fitNULLGLMM(
      plinkFile = name.bfile,
      traitType = "binary",
      phenoFile = pheno.file,
      phenoCol = phenoCol,
      covarColList = covarColList,
      qCovarCol = qCovarList,
      sampleIDColinphenoFile = "IID",
      LOCO = FALSE,
      IsOverwriteVarianceRatioFile = TRUE,
      useSparseGRMtoFitNULL = TRUE,
      sparseGRMFile = paste0(
        name.grm,
        "_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx"
      ),
      sparseGRMSampleIDFile = paste0(
        name.grm,
        "_relatednessCutoff_0.125_2000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
      ),
      outputPrefix = name.out,
      nThreads = 1
    )
    cat("[DONE] Step1", pheno.iter, "\n")
    TRUE
  }, error = function(e) {
    cat("[ERROR] Step1", pheno.iter, ":", conditionMessage(e), "\n")
    FALSE
  })
}

# -----------------------------------------------------------------------------
# Step2: SPAGMMATtest
# -----------------------------------------------------------------------------
run_step2 <- function(pheno_name) {
  cat("\n[Step2]", pheno_name, "\n")

  name.out   <- file.path(saige.dir, paste0(pheno_name, "_SAIGE_cumulative_inci"))
  model.file <- paste0(name.out, ".rda")
  var.file   <- paste0(name.out, ".varianceRatio.txt")
  out.file   <- paste0(name.out, ".SAIGE.txt")

  if (!file.exists(model.file)) {
    cat("[SKIP] missing model:", model.file, "\n")
    return(NULL)
  }
  if (!file.exists(var.file)) {
    cat("[SKIP] missing variance ratio:", var.file, "\n")
    return(NULL)
  }

  tryCatch({
    SAIGE::SPAGMMATtest(
      bedFile = paste0(name.bfile, ".bed"),
      bimFile = paste0(name.bfile, ".bim"),
      famFile = paste0(name.bfile, ".fam"),
      LOCO = FALSE,
      is_Firth_beta = TRUE,
      pCutoffforFirth = 0.05,
      is_output_moreDetails = TRUE,
      GMMATmodelFile = model.file,
      varianceRatioFile = var.file,
      SAIGEOutputFile = out.file
    )
    cat("[DONE] Step2", pheno_name, "->", out.file, "\n")
    out.file
  }, error = function(e) {
    cat("[ERROR] Step2", pheno_name, ":", conditionMessage(e), "\n")
    NULL
  })
}

# -----------------------------------------------------------------------------
# SE fix (SPA variants only)
# -----------------------------------------------------------------------------
required_cols <- c("BETA", "SE", "p.value", "Is.SPA")
qc_cols <- c(
  "CHR", "POS", "MarkerID", "BETA", "SE_original", "SE_corrected",
  "SE_ratio", "p.value", "p.value.NA", "Is.SPA"
)

is_complete_text_file <- function(path) {
  info <- file.info(path)
  if (is.na(info$size) || info$size == 0) return(FALSE)
  con <- file(path, "rb")
  on.exit(close(con), add = TRUE)
  seek(con, where = info$size - 1L, origin = "start")
  identical(readBin(con, what = "raw", n = 1L), as.raw(10L))
}

empty_summary <- function(pheno, infile, outfile, status, note = "") {
  data.table(
    pheno = pheno, status = status, note = note,
    infile = infile, outfile = outfile,
    n_total = NA_integer_, n_spa = NA_integer_, n_corrected = NA_integer_,
    n_spa_not_corrected = NA_integer_,
    n_diff_gt_1pct = NA_integer_, n_diff_gt_10pct = NA_integer_,
    median_SE_ratio = NA_real_, mean_SE_ratio = NA_real_,
    max_SE_ratio = NA_real_, min_SE_ratio = NA_real_
  )
}

saige_txt_path <- function(pheno) {
  file.path(saige.dir, paste0(pheno, "_SAIGE_cumulative_inci.SAIGE.txt"))
}

fix_one_pheno <- function(pheno) {
  infile  <- saige_txt_path(pheno)
  outfile <- file.path(
    sumstat.dir,
    paste0(pheno, "_SAIGE_cumulative_inci.SAIGE.SEfixed.txt")
  )
  qcfile  <- file.path(qc.dir, paste0(pheno, "_SEfixed_QC.txt"))
  topfile <- file.path(qc.dir, paste0(pheno, "_SEfixed_top20.txt"))

  cat("\n[SE fix]", pheno, "\n")

  if (!file.exists(infile)) {
    return(empty_summary(pheno, infile, outfile, "SKIP", "input not found"))
  }
  if (!force_overwrite && file.exists(outfile)) {
    return(empty_summary(pheno, infile, outfile, "SKIP", "outfile exists"))
  }
  if (!is_complete_text_file(infile)) {
    return(empty_summary(pheno, infile, outfile, "SKIP", "incomplete input"))
  }

  tryCatch({
    dt <- fread(infile, showProgress = FALSE)

    missing_cols <- setdiff(required_cols, names(dt))
    if (length(missing_cols) > 0) {
      stop("missing columns: ", paste(missing_cols, collapse = ", "))
    }

    dt[, SE_original := SE]
    dt[, SE_corrected := SE]

    spa <- tolower(as.character(dt$Is.SPA)) == "true"
    idx <- spa &
      !is.na(dt$p.value) & dt$p.value > 0 & dt$p.value < 1 &
      is.finite(dt$p.value) & !is.na(dt$BETA) & dt$BETA != 0

    dt[idx, SE_corrected := abs(BETA / qnorm(p.value / 2))]

    finite_ok <- idx & is.finite(dt$SE_corrected) & dt$SE_corrected > 0
    dt[finite_ok, SE := SE_corrected]
    dt[, SE_ratio := SE_original / SE_corrected]

    n_total <- nrow(dt)
    n_spa <- as.integer(sum(spa, na.rm = TRUE))
    n_corrected <- as.integer(sum(finite_ok, na.rm = TRUE))
    n_spa_not_corrected <- as.integer(sum(spa & !finite_ok, na.rm = TRUE))

    ratio_ok <- finite_ok & is.finite(dt$SE_ratio) & dt$SE_ratio > 0
    n_diff_gt_1pct  <- as.integer(sum(ratio_ok & abs(dt$SE_ratio - 1) > 0.01, na.rm = TRUE))
    n_diff_gt_10pct <- as.integer(sum(ratio_ok & abs(dt$SE_ratio - 1) > 0.10, na.rm = TRUE))

    se_ratio_vals <- dt$SE_ratio[ratio_ok]
    med_ratio  <- if (length(se_ratio_vals) > 0) as.numeric(median(se_ratio_vals)) else NA_real_
    mean_ratio <- if (length(se_ratio_vals) > 0) as.numeric(mean(se_ratio_vals)) else NA_real_
    max_ratio  <- if (length(se_ratio_vals) > 0) as.numeric(max(se_ratio_vals)) else NA_real_
    min_ratio  <- if (length(se_ratio_vals) > 0) as.numeric(min(se_ratio_vals)) else NA_real_

    spa_tab <- dt[, .N, by = .(Is.SPA = as.character(Is.SPA))]
    qc_con <- file(qcfile, open = "wt")
    sink(qc_con)
    cat("pheno:", pheno, "\n")
    cat("infile:", infile, "\n")
    cat("outfile:", outfile, "\n")
    cat("time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
    cat("table(Is.SPA):\n"); print(spa_tab); cat("\n")
    cat("n_total:", n_total, "\n")
    cat("n_spa:", n_spa, "\n")
    cat("n_corrected:", n_corrected, "\n")
    cat("n_spa_not_corrected:", n_spa_not_corrected, "\n")
    cat("n_|SE_ratio-1|>0.01:", n_diff_gt_1pct, "\n")
    cat("n_|SE_ratio-1|>0.10:", n_diff_gt_10pct, "\n")
    cat("median_SE_ratio:", med_ratio, "\n")
    cat("mean_SE_ratio:", mean_ratio, "\n")
    cat("min_SE_ratio:", min_ratio, "\n")
    cat("max_SE_ratio:", max_ratio, "\n")
    sink(); close(qc_con)

    keep_qc_cols <- intersect(qc_cols, names(dt))
    top20 <- dt[ratio_ok][order(-abs(log(SE_ratio)))]
    if (nrow(top20) > 20L) top20 <- top20[1:20]
    if (nrow(top20) > 0) {
      fwrite(top20[, ..keep_qc_cols], topfile, sep = "\t", quote = FALSE, na = "NA")
    } else {
      fwrite(data.table(note = "no SPA variants eligible"), topfile, sep = "\t", quote = FALSE)
    }

    dt[, c("SE_corrected", "SE_ratio") := NULL]
    fwrite(dt, outfile, sep = "\t", quote = FALSE, na = "NA", showProgress = FALSE)
    cat("[DONE] SE fix", pheno, "SPA=", n_spa, "corrected=", n_corrected, "\n")

    data.table(
      pheno = pheno, status = "OK", note = "",
      infile = infile, outfile = outfile,
      n_total = n_total, n_spa = n_spa, n_corrected = n_corrected,
      n_spa_not_corrected = n_spa_not_corrected,
      n_diff_gt_1pct = n_diff_gt_1pct, n_diff_gt_10pct = n_diff_gt_10pct,
      median_SE_ratio = med_ratio, mean_SE_ratio = mean_ratio,
      max_SE_ratio = max_ratio, min_SE_ratio = min_ratio
    )
  }, error = function(e) {
    cat("[ERROR] SE fix", pheno, ":", conditionMessage(e), "\n")
    empty_summary(pheno, infile, outfile, "ERROR", conditionMessage(e))
  })
}

run_sefix <- function(phenos) {
  if (njob > 1L) setDTthreads(1L) else setDTthreads(0L)

  summaries <- run_parallel(phenos, fix_one_pheno)
  summary_dt <- rbindlist(summaries, fill = TRUE)

  summary_file <- file.path(qc.dir, "SEfixed_QC_summary.tsv")
  if (file.exists(summary_file)) {
    old_summary <- fread(summary_file)
    old_summary <- old_summary[!pheno %in% summary_dt$pheno]
    summary_dt <- rbindlist(list(old_summary, summary_dt), fill = TRUE)
  }
  setorder(summary_dt, pheno)
  fwrite(summary_dt, summary_file, sep = "\t", quote = FALSE, na = "NA")
  cat("[DONE] SE fix summary ->", summary_file, "\n")
  invisible(summary_dt)
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
if (stage %in% c("step1", "all")) run_parallel(pheno_list_binary, run_step1)
if (stage %in% c("step2", "all")) run_parallel(pheno_list_binary, run_step2)
if (stage %in% c("sefix", "all")) run_sefix(pheno_list_binary)

cat("[DONE] stage=", stage, "\n", sep = "")
