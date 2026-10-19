#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   In-memory builders for Cont Supplementary Data 3 Summary and cohort-detail
#   tables. Returns data.tables only; does NOT write Excel.
#
# Inputs (callers supply fread() results)
#   - EAS Cont novelty TSVs: BBJ / TWB / TPMI
#   - EUR Cont novelty TSVs: UKB / FinnGen / MVP
#   - Single-cohort novelty TSV for build_cohort_detail_cont()
#
# Outputs (returned, not written)
#   - build_summary_eas_cont()  -> Summary_EAS-shaped data.table
#   - build_summary_eur_cont()  -> Summary_EUR-shaped data.table
#   - build_cohort_detail_cont() -> one cohort detail sheet data.table
#
# Usage
#   source("ExtendedData2_builders_cont.R")
#   Called by ExtendedData2_cont.R after external_validation + panel meta.
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("package 'data.table' is required")
}
suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Shared helpers (Cont summary labels)
# ------------------------------------------------------------------------------

supple3_cont_fmt_p <- function(p) {
  ifelse(is.na(p), "—", formatC(p, format = "e", digits = 2))
}

supple3_cont_get_novelty_label <- function(novel_flag, mapping_status, cohort_prefix) {
  no_match_val <- paste0("no_matching_", cohort_prefix, "_phenotype")
  ifelse(
    !is.na(mapping_status) & mapping_status == no_match_val,
    "—",
    ifelse(
      is.na(novel_flag),
      "—",
      ifelse(novel_flag == TRUE, "novel", "known")
    )
  )
}

supple3_cont_get_replication_label <- function(
    mapping_status_col,
    replication_testable_col,
    not_testable_reason_col,
    replicated_col,
    direction_concordant_col,
    cohort_prefix
) {
  no_match_val <- paste0("no_matching_", cohort_prefix, "_phenotype")
  absent_pattern <- paste0("lead_variant_absent_after_", cohort_prefix, "_basic_qc")
  result <- character(length(mapping_status_col))
  for (i in seq_along(result)) {
    ms  <- mapping_status_col[i]
    rt  <- replication_testable_col[i]
    ntr <- not_testable_reason_col[i]
    rp  <- replicated_col[i]
    dc  <- direction_concordant_col[i]
    if (!is.na(ms) && ms == no_match_val) {
      result[i] <- "—"
    } else if (!is.na(rt) && rt == FALSE) {
      if (!is.na(ntr) && grepl(absent_pattern, ntr)) {
        result[i] <- "No variant"
      } else {
        result[i] <- "Not analyzed"
      }
    } else if (!is.na(rt) && rt == TRUE) {
      if (!is.na(rp) && rp == TRUE && !is.na(dc) && dc == TRUE) {
        result[i] <- "Replicated"
      } else {
        result[i] <- "Not replicated"
      }
    } else {
      result[i] <- "—"
    }
  }
  result
}

supple3_cont_get_combined_novelty <- function(...) {
  labels <- list(...)
  n <- length(labels[[1L]])
  result <- character(n)
  for (i in seq_len(n)) {
    vals <- vapply(labels, function(x) x[[i]], character(1L))
    vals <- vals[vals != "—"]
    if (length(vals) == 0L) {
      result[i] <- "—"
    } else if (any(vals == "known")) {
      result[i] <- "known"
    } else {
      result[i] <- "novel"
    }
  }
  result
}

supple3_cont_get_any_cohort_replication <- function(...) {
  labels <- list(...)
  n <- length(labels[[1L]])
  result <- character(n)
  for (i in seq_len(n)) {
    vals_all <- vapply(labels, function(x) x[[i]], character(1L))
    if (all(vals_all == "—")) {
      result[i] <- "—"
    } else if (any(vals_all == "Replicated")) {
      result[i] <- "Replicated"
    } else {
      vals <- vals_all[vals_all != "—"]
      if (all(vals == "No variant")) {
        result[i] <- "No variant"
      } else if (all(vals %in% c("No variant", "Not analyzed"))) {
        result[i] <- "Not analyzed"
      } else {
        result[i] <- "Not replicated"
      }
    }
  }
  result
}

supple3_cont_as_logical_flag <- function(x) {
  !is.na(x) & (x == TRUE | x == 1L | x == "TRUE" | x == "true")
}

supple3_cont_dir_aware_raw_p <- function(testable, concordant, exact_p) {
  is_t <- supple3_cont_as_logical_flag(testable)
  out <- rep(NA_real_, length(testable))
  if (!any(is_t)) return(out)
  out[is_t] <- 1
  ok <- is_t & supple3_cont_as_logical_flag(concordant) & is.finite(exact_p)
  out[ok] <- exact_p[ok]
  out
}

supple3_cont_get_integrated_replication <- function(
    testable_list,
    concordant_list,
    exact_p_list,
    pheno_code,
    alpha = 0.05
) {
  n <- length(pheno_code)
  n_cohort <- length(testable_list)
  if (length(concordant_list) != n_cohort || length(exact_p_list) != n_cohort) {
    stop("testable_list / concordant_list / exact_p_list length mismatch")
  }
  testable_mat <- do.call(cbind, lapply(testable_list, supple3_cont_as_logical_flag))
  k <- as.integer(rowSums(testable_mat))
  concordant_mat <- do.call(cbind, lapply(concordant_list, supple3_cont_as_logical_flag))
  n_concordant <- as.integer(rowSums(testable_mat & concordant_mat))
  p_mat <- do.call(
    cbind,
    Map(supple3_cont_dir_aware_raw_p, testable_list, concordant_list, exact_p_list)
  )
  min_p <- suppressWarnings(apply(p_mat, 1L, function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0L) NA_real_ else min(x)
  }))
  min_p[k == 0L] <- NA_real_
  p_integrated <- ifelse(
    k >= 1L & is.finite(min_p),
    pmin(1, k * min_p),
    NA_real_
  )
  M_global <- as.integer(sum(k >= 1L))
  M <- rep(M_global, n)
  thr_scalar <- if (M_global > 0L) alpha / M_global else NA_real_
  threshold <- rep(thr_scalar, n)
  label <- rep("—", n)
  idx <- which(k >= 1L)
  label[idx] <- ifelse(
    is.finite(p_integrated[idx]) &
      is.finite(threshold[idx]) &
      p_integrated[idx] <= threshold[idx],
    "Replicated",
    "Not replicated"
  )
  list(
    k = k,
    n_concordant = n_concordant,
    min_raw_p = min_p,
    p_integrated = p_integrated,
    M = M,
    M_global = M_global,
    threshold = threshold,
    label = label
  )
}

supple3_cont_yes_no <- function(x) {
  x <- as.logical(x)
  fifelse(is.na(x), NA_character_, fifelse(x, "Yes", "No"))
}

# ------------------------------------------------------------------------------
# build_summary_eas_cont(bbj, twb, tpmi)
# ------------------------------------------------------------------------------
build_summary_eas_cont <- function(bbj, twb, tpmi) {
  bbj <- as.data.table(bbj)
  twb <- as.data.table(twb)
  tpmi <- as.data.table(tpmi)

  if (!identical(bbj$lead_row_id, twb$lead_row_id) ||
      !identical(bbj$lead_row_id, tpmi$lead_row_id)) {
    stop("BBJ / TWB / TPMI lead_row_id order differs.")
  }
  if (!identical(bbj$koges_lead_variant, twb$koges_lead_variant) ||
      !identical(bbj$koges_lead_variant, tpmi$koges_lead_variant)) {
    stop("BBJ / TWB / TPMI koges_lead_variant order differs.")
  }

  required_bonf_cols <- c(
    "replicated_global_testable",
    "replication_threshold_global_testable",
    "replication_testable",
    "direction_concordant"
  )
  required_exact_p <- c(bbj = "bbj_exact_p", twb = "twb_exact_p", tpmi = "tpmi_exact_p")
  for (nm in c("bbj", "twb", "tpmi")) {
    dt <- get(nm)
    missing <- setdiff(c(required_bonf_cols, required_exact_p[[nm]]), names(dt))
    if (length(missing) > 0L) {
      stop(nm, " input is missing required column(s): ", paste(missing, collapse = ", "))
    }
  }

  bbj_thr <- unique(na.omit(bbj$replication_threshold_global_testable))
  twb_thr <- unique(na.omit(twb$replication_threshold_global_testable))
  tpmi_thr <- unique(na.omit(tpmi$replication_threshold_global_testable))
  if (length(bbj_thr) == 0L || length(twb_thr) == 0L || length(tpmi_thr) == 0L) {
    stop("replication_threshold_global_testable is empty.")
  }
  message(sprintf(
    "EAS Cont cohort-specific global_testable thresholds | BBJ: %s | TWB: %s | TPMI: %s",
    paste(formatC(bbj_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(twb_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(tpmi_thr, format = "e", digits = 6), collapse = ", ")
  ))

  bbj_novelty_label <- supple3_cont_get_novelty_label(
    bbj$novel_relative_to_bbj, bbj$bbj_mapping_status, "bbj"
  )
  twb_novelty_label <- supple3_cont_get_novelty_label(
    twb$novel_relative_to_twb, twb$twb_mapping_status, "twb"
  )
  tpmi_novelty_label <- supple3_cont_get_novelty_label(
    tpmi$novel_relative_to_tpmi, tpmi$tpmi_mapping_status, "tpmi"
  )
  combined_novelty_label <- supple3_cont_get_combined_novelty(
    bbj_novelty_label, twb_novelty_label, tpmi_novelty_label
  )

  bbj_rep_nominal <- supple3_cont_get_replication_label(
    bbj$bbj_mapping_status, bbj$replication_testable,
    bbj$replication_not_testable_reason, bbj$replicated_nominal,
    bbj$direction_concordant, "bbj"
  )
  twb_rep_nominal <- supple3_cont_get_replication_label(
    twb$twb_mapping_status, twb$replication_testable,
    twb$replication_not_testable_reason, twb$replicated_nominal,
    twb$direction_concordant, "twb"
  )
  tpmi_rep_nominal <- supple3_cont_get_replication_label(
    tpmi$tpmi_mapping_status, tpmi$replication_testable,
    tpmi$replication_not_testable_reason, tpmi$replicated_nominal,
    tpmi$direction_concordant, "tpmi"
  )
  any_rep_nominal <- supple3_cont_get_any_cohort_replication(
    bbj_rep_nominal, twb_rep_nominal, tpmi_rep_nominal
  )

  bbj_rep_bonf <- supple3_cont_get_replication_label(
    bbj$bbj_mapping_status, bbj$replication_testable,
    bbj$replication_not_testable_reason, bbj$replicated_global_testable,
    bbj$direction_concordant, "bbj"
  )
  twb_rep_bonf <- supple3_cont_get_replication_label(
    twb$twb_mapping_status, twb$replication_testable,
    twb$replication_not_testable_reason, twb$replicated_global_testable,
    twb$direction_concordant, "twb"
  )
  tpmi_rep_bonf <- supple3_cont_get_replication_label(
    tpmi$tpmi_mapping_status, tpmi$replication_testable,
    tpmi$replication_not_testable_reason, tpmi$replicated_global_testable,
    tpmi$direction_concordant, "tpmi"
  )
  any_rep_bonf <- supple3_cont_get_any_cohort_replication(
    bbj_rep_bonf, twb_rep_bonf, tpmi_rep_bonf
  )

  integ <- supple3_cont_get_integrated_replication(
    testable_list = list(
      bbj$replication_testable, twb$replication_testable, tpmi$replication_testable
    ),
    concordant_list = list(
      bbj$direction_concordant, twb$direction_concordant, tpmi$direction_concordant
    ),
    exact_p_list = list(bbj$bbj_exact_p, twb$twb_exact_p, tpmi$tpmi_exact_p),
    pheno_code = bbj$koges_pheno_code,
    alpha = 0.05
  )
  # Integrated label is computed in the source script for QC messaging only;
  # Summary_EAS sheet columns do not include it.
  message(sprintf(
    "Integrated EAS Cont replication | M_global=%d | threshold=%s",
    integ$M_global,
    format(if (integ$M_global > 0L) 0.05 / integ$M_global else NA_real_, digits = 6)
  ))

  bbj_min_p <- ifelse(
    bbj$bbj_mapping_status == "no_matching_bbj_phenotype", NA_real_, bbj$bbj_top_window_p
  )
  twb_min_p <- ifelse(
    twb$twb_mapping_status == "no_matching_twb_phenotype", NA_real_, twb$twb_top_window_p
  )
  tpmi_min_p <- ifelse(
    tpmi$tpmi_mapping_status == "no_matching_tpmi_phenotype", NA_real_, tpmi$tpmi_top_window_p
  )

  out <- data.table(
    Phenotype = bbj$koges_phenotype,
    Chromosome = bbj$koges_chr,
    `Base pair position` = bbj$koges_bp,
    rsID = bbj$koges_lead_variant,
    `P-value in KoGES` = supple3_cont_fmt_p(bbj$koges_p),
    `Minimum P-value within ±1 Mb in BBJ` = supple3_cont_fmt_p(bbj_min_p),
    `Minimum P-value within ±1 Mb in TWB` = supple3_cont_fmt_p(twb_min_p),
    `Minimum P-value within ±1 Mb in TPMI` = supple3_cont_fmt_p(tpmi_min_p),
    `Novel relative to assessed EastAsian-ancestry cohorts` = combined_novelty_label,
    `Variant-level replication - BBJ` = bbj_rep_nominal,
    `Variant-level replication - TWB` = twb_rep_nominal,
    `Variant-level replication - TPMI` = tpmi_rep_nominal,
    `Nominal replication in ≥1 cohort` = any_rep_nominal,
    `Bonferroni replication - BBJ` = bbj_rep_bonf,
    `Bonferroni replication - TWB` = twb_rep_bonf,
    `Bonferroni replication - TPMI` = tpmi_rep_bonf,
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort` = any_rep_bonf
  )
  setorderv(out, c("Phenotype", "Chromosome", "Base pair position"))
  out
}

# ------------------------------------------------------------------------------
# build_summary_eur_cont(ukb, finngen, mvp)
# ------------------------------------------------------------------------------
build_summary_eur_cont <- function(ukb, finngen, mvp) {
  ukb <- as.data.table(ukb)
  finngen <- as.data.table(finngen)
  mvp <- as.data.table(mvp)

  if (!identical(ukb$lead_row_id, finngen$lead_row_id) ||
      !identical(ukb$lead_row_id, mvp$lead_row_id)) {
    stop("UKB / FinnGen / MVP lead_row_id order differs.")
  }
  if (!identical(ukb$koges_lead_variant, finngen$koges_lead_variant) ||
      !identical(ukb$koges_lead_variant, mvp$koges_lead_variant)) {
    stop("UKB / FinnGen / MVP koges_lead_variant order differs.")
  }

  required_bonf_cols <- c(
    "replicated_global_testable",
    "replication_threshold_global_testable",
    "replication_testable",
    "direction_concordant"
  )
  required_exact_p <- c(ukb = "ukb_exact_p", finngen = "finngen_exact_p", mvp = "mvp_exact_p")
  for (nm in c("ukb", "finngen", "mvp")) {
    dt <- get(nm)
    missing <- setdiff(c(required_bonf_cols, required_exact_p[[nm]]), names(dt))
    if (length(missing) > 0L) {
      stop(nm, " input is missing required column(s): ", paste(missing, collapse = ", "))
    }
  }

  ukb_thr <- unique(na.omit(ukb$replication_threshold_global_testable))
  finngen_thr <- unique(na.omit(finngen$replication_threshold_global_testable))
  mvp_thr <- unique(na.omit(mvp$replication_threshold_global_testable))
  if (length(ukb_thr) == 0L || length(finngen_thr) == 0L || length(mvp_thr) == 0L) {
    stop("replication_threshold_global_testable is empty.")
  }
  message(sprintf(
    "EUR Cont cohort-specific global_testable thresholds | UKB: %s | FinnGen: %s | MVP: %s",
    paste(formatC(ukb_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(finngen_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(mvp_thr, format = "e", digits = 6), collapse = ", ")
  ))

  ukb_novelty_label <- supple3_cont_get_novelty_label(
    ukb$novel_relative_to_ukb, ukb$ukb_mapping_status, "ukb"
  )
  finngen_novelty_label <- supple3_cont_get_novelty_label(
    finngen$novel_relative_to_finngen, finngen$finngen_mapping_status, "finngen"
  )
  mvp_novelty_label <- supple3_cont_get_novelty_label(
    mvp$novel_relative_to_mvp, mvp$mvp_mapping_status, "mvp"
  )
  combined_novelty_label <- supple3_cont_get_combined_novelty(
    ukb_novelty_label, finngen_novelty_label, mvp_novelty_label
  )

  ukb_rep_nominal <- supple3_cont_get_replication_label(
    ukb$ukb_mapping_status, ukb$replication_testable,
    ukb$replication_not_testable_reason, ukb$replicated_nominal,
    ukb$direction_concordant, "ukb"
  )
  finngen_rep_nominal <- supple3_cont_get_replication_label(
    finngen$finngen_mapping_status, finngen$replication_testable,
    finngen$replication_not_testable_reason, finngen$replicated_nominal,
    finngen$direction_concordant, "finngen"
  )
  mvp_rep_nominal <- supple3_cont_get_replication_label(
    mvp$mvp_mapping_status, mvp$replication_testable,
    mvp$replication_not_testable_reason, mvp$replicated_nominal,
    mvp$direction_concordant, "mvp"
  )
  any_rep_nominal <- supple3_cont_get_any_cohort_replication(
    ukb_rep_nominal, finngen_rep_nominal, mvp_rep_nominal
  )

  ukb_rep_bonf <- supple3_cont_get_replication_label(
    ukb$ukb_mapping_status, ukb$replication_testable,
    ukb$replication_not_testable_reason, ukb$replicated_global_testable,
    ukb$direction_concordant, "ukb"
  )
  finngen_rep_bonf <- supple3_cont_get_replication_label(
    finngen$finngen_mapping_status, finngen$replication_testable,
    finngen$replication_not_testable_reason, finngen$replicated_global_testable,
    finngen$direction_concordant, "finngen"
  )
  mvp_rep_bonf <- supple3_cont_get_replication_label(
    mvp$mvp_mapping_status, mvp$replication_testable,
    mvp$replication_not_testable_reason, mvp$replicated_global_testable,
    mvp$direction_concordant, "mvp"
  )
  any_rep_bonf <- supple3_cont_get_any_cohort_replication(
    ukb_rep_bonf, finngen_rep_bonf, mvp_rep_bonf
  )

  integ <- supple3_cont_get_integrated_replication(
    testable_list = list(
      ukb$replication_testable, finngen$replication_testable, mvp$replication_testable
    ),
    concordant_list = list(
      ukb$direction_concordant, finngen$direction_concordant, mvp$direction_concordant
    ),
    exact_p_list = list(ukb$ukb_exact_p, finngen$finngen_exact_p, mvp$mvp_exact_p),
    pheno_code = ukb$koges_pheno_code,
    alpha = 0.05
  )
  # Integrated label is computed in the source script for QC messaging only;
  # Summary_EUR sheet columns do not include it.
  message(sprintf(
    "Integrated EUR Cont replication | M_global=%d | threshold=%s",
    integ$M_global,
    format(if (integ$M_global > 0L) 0.05 / integ$M_global else NA_real_, digits = 6)
  ))

  ukb_min_p <- ifelse(
    ukb$ukb_mapping_status == "no_matching_ukb_phenotype", NA_real_, ukb$ukb_top_window_p
  )
  finngen_min_p <- ifelse(
    finngen$finngen_mapping_status == "no_matching_finngen_phenotype",
    NA_real_, finngen$finngen_top_window_p
  )
  mvp_min_p <- ifelse(
    mvp$mvp_mapping_status == "no_matching_mvp_phenotype", NA_real_, mvp$mvp_top_window_p
  )

  out <- data.table(
    Phenotype = ukb$koges_phenotype,
    Chromosome = ukb$koges_chr,
    `Base pair position` = ukb$koges_bp,
    rsID = ukb$koges_lead_variant,
    `P-value in KoGES` = supple3_cont_fmt_p(ukb$koges_p),
    `Minimum P-value within ±1 Mb in UKB` = supple3_cont_fmt_p(ukb_min_p),
    `Minimum P-value within ±1 Mb in FinnGen` = supple3_cont_fmt_p(finngen_min_p),
    `Minimum P-value within ±1 Mb in MVP` = supple3_cont_fmt_p(mvp_min_p),
    `Novel relative to assessed European-ancestry cohorts` = combined_novelty_label,
    `Variant-level replication - UKB` = ukb_rep_nominal,
    `Variant-level replication - FinnGen` = finngen_rep_nominal,
    `Variant-level replication - MVP` = mvp_rep_nominal,
    `Nominal replication in ≥1 cohort` = any_rep_nominal,
    `Bonferroni replication - UKB` = ukb_rep_bonf,
    `Bonferroni replication - FinnGen` = finngen_rep_bonf,
    `Bonferroni replication - MVP` = mvp_rep_bonf,
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort` = any_rep_bonf
  )
  setorderv(out, c("Phenotype", "Chromosome", "Base pair position"))
  out
}

# ------------------------------------------------------------------------------
# build_cohort_detail_cont(novelty_dt, cohort_prefix)
#   cohort_prefix: bbj | twb | tpmi | ukb | finngen | mvp
# ------------------------------------------------------------------------------
build_cohort_detail_cont <- function(novelty_dt, cohort_prefix) {
  dt <- as.data.table(novelty_dt)
  prefix <- tolower(as.character(cohort_prefix))
  required_common <- c(
    "lead_row_id", "koges_phenotype", "koges_pheno_code",
    "koges_lead_variant", "koges_variant_id", "koges_chr", "koges_bp",
    "koges_nearest_gene", "koges_ea", "koges_nea", "koges_eaf", "koges_beta",
    "koges_se", "koges_p", "koges_beta_aligned_a1",
    "replication_testable", "direction_concordant", "replicated_nominal",
    "replicated_global_testable", "replication_threshold_global_testable",
    "replication_status"
  )
  source_cols <- paste0(prefix, c(
    "_source_trait", "_mapping_status", "_n_variants_in_window",
    "_n_gws_hits_in_window", "_top_window_snp", "_top_window_rsid_koges",
    "_top_window_chr", "_top_window_pos", "_top_window_p",
    "_top_window_distance_bp", "_novelty_status",
    "_exact_chrpos_present", "_exact_snp", "_exact_rsid_koges",
    "_exact_a1_koges", "_exact_a2_koges", "_exact_beta_aligned_a1",
    "_exact_eaf_aligned_a1", "_exact_se", "_exact_p"
  ))
  missing <- setdiff(c(required_common, source_cols), names(dt))
  if (length(missing) > 0L) {
    stop(
      "Cont cohort detail (", prefix, ") missing required column(s): ",
      paste(missing, collapse = ", ")
    )
  }
  C <- function(suffix) paste0(prefix, suffix)
  data.table(
    `KoGES Phenotype` = dt$koges_phenotype,
    `KoGES Lead Variant (rsID)` = dt$koges_lead_variant,
    `KoGES Variant ID` = dt$koges_variant_id,
    `KoGES Chromosome` = dt$koges_chr,
    `KoGES Position (bp)` = dt$koges_bp,
    `Nearest Gene` = dt$koges_nearest_gene,
    `KoGES Effect Allele` = dt$koges_ea,
    `KoGES Other Allele` = dt$koges_nea,
    `KoGES EAF` = dt$koges_eaf,
    `KoGES Beta` = dt$koges_beta,
    `KoGES Standard Error` = dt$koges_se,
    `KoGES P Value` = dt$koges_p,
    `Lowest-P External Variant rsID Within ±1 Mb` = dt[[C("_top_window_rsid_koges")]],
    `Lowest-P External Variant Chromosome` = dt[[C("_top_window_chr")]],
    `Lowest-P External Variant Position (bp)` = dt[[C("_top_window_pos")]],
    `Minimum External P Value Within ±1 Mb` = dt[[C("_top_window_p")]],
    `Locus Novelty Assessment` = dt[[C("_novelty_status")]],
    `Exact-Position External Variant Available` = supple3_cont_yes_no(dt[[C("_exact_chrpos_present")]]),
    `Exact-Position External Variant rsID` = dt[[C("_exact_rsid_koges")]],
    `Aligned Effect Allele (A1)` = dt[[C("_exact_a1_koges")]],
    `Aligned Other Allele (A2)` = dt[[C("_exact_a2_koges")]],
    `KoGES Beta Aligned to A1` = dt$koges_beta_aligned_a1,
    `External Beta Aligned to A1` = dt[[C("_exact_beta_aligned_a1")]],
    `External EAF Aligned to A1` = dt[[C("_exact_eaf_aligned_a1")]],
    `External Standard Error` = dt[[C("_exact_se")]],
    `External P Value at Exact Position` = dt[[C("_exact_p")]],
    `Exact-Variant Replication Testable` = supple3_cont_yes_no(dt$replication_testable),
    `Effect Direction Concordant` = supple3_cont_yes_no(dt$direction_concordant),
    `Nominal Replication` = supple3_cont_yes_no(dt$replicated_nominal),
    `Primary Bonferroni Threshold` = dt$replication_threshold_global_testable,
    `Bonferroni-Significant Replication` = supple3_cont_yes_no(dt$replicated_global_testable),
    `Replication Assessment` = dt$replication_status
  )
}
