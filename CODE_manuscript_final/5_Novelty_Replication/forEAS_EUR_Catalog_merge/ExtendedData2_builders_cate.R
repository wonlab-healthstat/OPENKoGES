#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   In-memory builders for Cate (binary/SAIGE) Supplementary Data 3 Summary and
#   cohort-detail tables. Returns data.tables only; does NOT write Excel.
#
# Inputs (callers supply fread() results)
#   - EAS Cate novelty + per_match: BBJ, TPMI
#   - EUR Cate novelty + per_match: UKB, FinnGen, MVP
#   - EAF lookup (Cate_lead_variants_GWAScatalog_novel.tsv) for cohort detail when koges_eaf absent
#
# Outputs (returned, not written)
#   - build_summary_eas_cate() / build_summary_eur_cate()
#   - build_cohort_detail_cate()
#
# Usage
#   source("ExtendedData2_builders_cate.R")
#   Called by ExtendedData2_cate.R after external_validation + panel meta.
#
# Notes
#   - EAS = BBJ + TPMI only (no TWB).
#   - EUR validation TSVs use cohort prefixes (ukb_/finngen_/mvp_); rename_tpmi kept for legacy.
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("package 'data.table' is required")
}
suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# Shared helpers
# ------------------------------------------------------------------------------
supple3_cate_fmt_p <- function(p) {
  ifelse(is.na(p), "—", formatC(p, format = "e", digits = 2))
}

supple3_cate_get_novelty_label <- function(novel_flag, mapping_status, cohort_prefix) {
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

supple3_cate_get_replication_label <- function(
    mapping_status_col,
    replication_testable_col,
    not_testable_reason_col,
    replicated_col,
    direction_concordant_col,
    cohort_prefix
) {
  no_match_val <- paste0("no_matching_", cohort_prefix, "_phenotype")
  absent_pattern <- paste0(
    "lead_variant_absent_after_", cohort_prefix, "_basic_qc"
  )
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

supple3_cate_get_combined_novelty <- function(...) {
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

supple3_cate_get_any_cohort_replication <- function(...) {
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

supple3_cate_yes_dash <- function(flag) {
  ifelse(flag %in% TRUE, "Yes", "—")
}

supple3_cate_expand_rep_with_additional <- function(direct_label, add_flag) {
  ifelse(add_flag %in% TRUE, "Replicated", direct_label)
}

supple3_cate_as_logical_flag <- function(x) {
  !is.na(x) & (x == TRUE | x == 1L | x == "TRUE" | x == "true")
}

supple3_cate_as_flag <- function(x) {
  if (is.logical(x)) return(x)
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %chin% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %chin% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

supple3_cate_yes_no <- function(x) {
  x <- as.logical(x)
  fifelse(is.na(x), NA_character_, fifelse(x, "Yes", "No"))
}

supple3_cate_dir_aware_raw_p <- function(testable, concordant, exact_p) {
  is_t <- supple3_cate_as_logical_flag(testable)
  out <- rep(NA_real_, length(testable))
  if (!any(is_t)) return(out)
  out[is_t] <- 1
  ok <- is_t & supple3_cate_as_logical_flag(concordant) & is.finite(exact_p)
  out[ok] <- exact_p[ok]
  out
}

supple3_cate_get_integrated_replication <- function(
    testable_list,
    concordant_list,
    exact_p_list,
    alpha = 0.05
) {
  n_cohort <- length(testable_list)
  testable_mat <- do.call(cbind, lapply(testable_list, supple3_cate_as_logical_flag))
  k <- as.integer(rowSums(testable_mat))
  concordant_mat <- do.call(cbind, lapply(concordant_list, supple3_cate_as_logical_flag))
  n_concordant <- as.integer(rowSums(testable_mat & concordant_mat))
  p_mat <- do.call(
    cbind,
    Map(supple3_cate_dir_aware_raw_p, testable_list, concordant_list, exact_p_list)
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
  n <- length(k)
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
    M_global = M_global,
    threshold = threshold,
    label = label
  )
}

supple3_cate_pick_any_match <- function(direct, other) {
  fifelse(supple3_cate_as_logical_flag(direct), direct, other)
}

supple3_cate_pick_any_p <- function(testable_direct, p_direct, testable_other, p_other) {
  fifelse(supple3_cate_as_logical_flag(testable_direct), p_direct, p_other)
}

supple3_cate_min_window_p <- function(pm, p_col, mapping_col, scope) {
  x <- copy(as.data.table(pm))
  need_pm <- c("lead_row_id", "is_direct_match", p_col, mapping_col)
  miss_pm <- setdiff(need_pm, names(x))
  if (length(miss_pm)) {
    stop("per-match missing column(s): ", paste(miss_pm, collapse = ", "))
  }
  if (identical(scope, "direct")) {
    x <- x[is_direct_match %in% TRUE]
  } else if (identical(scope, "other")) {
    x <- x[is_direct_match %in% FALSE]
  }
  x <- x[get(mapping_col) == "mapped_file_found" & is.finite(get(p_col))]
  if (!nrow(x)) {
    return(data.table(lead_row_id = integer(), min_p = numeric()))
  }
  x[, .(min_p = min(get(p_col))), by = lead_row_id]
}

supple3_cate_rename_tpmi_prefix <- function(dt, cohort_prefix) {
  dt <- as.data.table(dt)
  nms <- names(dt)
  hit <- grepl("^tpmi_", nms)
  if (any(hit)) {
    setnames(dt, old = nms[hit], new = sub("^tpmi_", paste0(cohort_prefix, "_"), nms[hit]))
  }
  dt
}

supple3_cate_first_present <- function(dt, candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0L) return(NULL)
  dt[[hit[[1L]]]]
}

supple3_cate_pick_direct_one <- function(pm, prefix, label) {
  stem_col <- paste0(prefix, "_output_stem")
  map_col <- paste0(prefix, "_mapping_status")
  p_col <- paste0(prefix, "_exact_p")
  need <- c(
    "lead_row_id", "is_direct_match", "replication_testable",
    stem_col, map_col,
    paste0(prefix, "_exact_chrpos_present"),
    paste0(prefix, "_exact_snp"),
    paste0(prefix, "_exact_rsid_koges"),
    paste0(prefix, "_exact_a1_koges"),
    paste0(prefix, "_exact_a2_koges"),
    paste0(prefix, "_exact_beta_aligned_a1"),
    paste0(prefix, "_exact_eaf_aligned_a1"),
    paste0(prefix, "_exact_se"),
    p_col,
    paste0(prefix, "_top_window_rsid_koges"),
    paste0(prefix, "_top_window_chr"),
    paste0(prefix, "_top_window_pos"),
    paste0(prefix, "_top_window_p"),
    paste0(prefix, "_source_trait"),
    paste0(prefix, "_n_variants_in_window"),
    paste0(prefix, "_n_gws_hits_in_window"),
    paste0(prefix, "_novelty_status")
  )
  soft <- paste0(prefix, "_novelty_status")
  need_hard <- setdiff(need, soft)
  miss <- setdiff(need_hard, names(pm))
  if (length(miss) > 0L) {
    stop(label, " per_match missing required column(s): ", paste(miss, collapse = ", "))
  }
  x <- copy(as.data.table(pm))
  if ("koges_pheno_code.x" %chin% names(x) && !"koges_pheno_code" %chin% names(x)) {
    setnames(x, "koges_pheno_code.x", "koges_pheno_code")
  }
  x <- x[is_direct_match %in% TRUE]
  if (!nrow(x)) {
    return(x[0L])
  }
  x[, rank_testable := as.integer(!(supple3_cate_as_flag(replication_testable) %in% TRUE))]
  x[, rank_p := fifelse(is.finite(as.numeric(get(p_col))), as.numeric(get(p_col)), Inf)]
  setorderv(x, c("lead_row_id", "rank_testable", "rank_p", stem_col))
  x[, .SD[1L], by = lead_row_id]
}

supple3_cate_resolve_replication_fields <- function(nov) {
  pick_col <- function(cands) {
    hit <- cands[cands %in% names(nov)]
    if (length(hit) == 0L) return(rep(NA, nrow(nov)))
    nov[[hit[[1L]]]]
  }
  list(
    replication_testable = pick_col(c(
      "replication_testable_direct", "replication_testable"
    )),
    direction_concordant = pick_col(c(
      "direction_concordant_direct", "direction_concordant"
    )),
    replicated_nominal = pick_col(c(
      "replicated_nominal_direct", "replicated_nominal"
    )),
    replicated_primary = pick_col(c(
      "replicated_primary_direct", "replicated_primary",
      "replicated_global_testable"
    )),
    replication_threshold = pick_col(c(
      "replication_threshold_primary_direct",
      "replication_threshold_primary",
      "replication_threshold_global_testable_direct",
      "replication_threshold_global_testable"
    )),
    replication_status = pick_col(c(
      "replication_status_direct", "replication_status"
    ))
  )
}

supple3_cate_load_eaf_lookup <- function(path) {
  eaf_lookup_raw <- fread(
    path, sep = "\t", na.strings = c("NA", ""),
    select = c("row_id", "EAF")
  )
  miss_eaf <- setdiff(c("row_id", "EAF"), names(eaf_lookup_raw))
  if (length(miss_eaf) > 0L) {
    stop("EAF lookup missing column(s): ", paste(miss_eaf, collapse = ", "))
  }
  eaf_lookup <- eaf_lookup_raw[, .(
    lead_row_id = as.integer(row_id),
    koges_eaf_lookup = as.numeric(EAF)
  )]
  if (anyDuplicated(eaf_lookup$lead_row_id)) {
    stop("Duplicate row_id in EAF lookup: ", path)
  }
  if (anyNA(eaf_lookup$koges_eaf_lookup)) {
    stop("NA EAF values in lookup: ", path)
  }
  eaf_lookup
}

# ------------------------------------------------------------------------------
# build_summary_eas_cate(bbj, tpmi, bbj_pm, tpmi_pm)
# ------------------------------------------------------------------------------
build_summary_eas_cate <- function(bbj, tpmi, bbj_pm, tpmi_pm) {
  bbj <- as.data.table(bbj)
  tpmi <- as.data.table(tpmi)
  bbj_pm <- as.data.table(bbj_pm)
  tpmi_pm <- as.data.table(tpmi_pm)

  if (!identical(bbj$lead_row_id, tpmi$lead_row_id)) {
    stop("BBJ / TPMI lead_row_id order differs.")
  }
  if (!identical(bbj$koges_lead_variant, tpmi$koges_lead_variant) ||
      !identical(bbj$koges_chr, tpmi$koges_chr) ||
      !identical(bbj$koges_bp, tpmi$koges_bp)) {
    stop("BBJ / TPMI KoGES lead coordinates differ.")
  }

  need_lead <- c(
    "lead_row_id", "koges_phenotype", "koges_pheno_code",
    "koges_chr", "koges_bp", "koges_lead_variant", "koges_p",
    "mapping_status_direct", "mapping_status_any",
    "additionally_known_from_other",
    "replication_testable_direct", "replication_testable_other",
    "not_testable_reason_direct",
    "replicated_nominal_direct", "replicated_primary_direct",
    "direction_concordant_direct", "direction_concordant_other",
    "best_testable_p_direct", "best_testable_p_other",
    "additionally_replicated_nominal_from_other",
    "additionally_replicated_primary_from_other",
    "replication_threshold_primary_direct"
  )
  need_bbj <- c(need_lead, "novel_relative_to_bbj_direct", "novel_relative_to_bbj_any")
  need_tpmi <- c(need_lead, "novel_relative_to_tpmi_direct", "novel_relative_to_tpmi_any")
  miss_bbj <- setdiff(need_bbj, names(bbj))
  miss_tpmi <- setdiff(need_tpmi, names(tpmi))
  if (length(miss_bbj)) stop("BBJ lead TSV missing: ", paste(miss_bbj, collapse = ", "))
  if (length(miss_tpmi)) stop("TPMI lead TSV missing: ", paste(miss_tpmi, collapse = ", "))

  bbj_p_direct <- supple3_cate_min_window_p(
    bbj_pm, "bbj_top_window_p", "bbj_mapping_status", "direct"
  )
  bbj_p_other <- supple3_cate_min_window_p(
    bbj_pm, "bbj_top_window_p", "bbj_mapping_status", "other"
  )
  tpmi_p_direct <- supple3_cate_min_window_p(
    tpmi_pm, "tpmi_top_window_p", "tpmi_mapping_status", "direct"
  )
  tpmi_p_other <- supple3_cate_min_window_p(
    tpmi_pm, "tpmi_top_window_p", "tpmi_mapping_status", "other"
  )
  bbj[bbj_p_direct, min_window_p_direct := i.min_p, on = "lead_row_id"]
  bbj[bbj_p_other, min_window_p_other := i.min_p, on = "lead_row_id"]
  tpmi[tpmi_p_direct, min_window_p_direct := i.min_p, on = "lead_row_id"]
  tpmi[tpmi_p_other, min_window_p_other := i.min_p, on = "lead_row_id"]

  bbj_thr <- unique(na.omit(bbj$replication_threshold_primary_direct))
  tpmi_thr <- unique(na.omit(tpmi$replication_threshold_primary_direct))
  message(sprintf(
    "EAS Cate direct-match Bonferroni thresholds | BBJ: %s | TPMI: %s",
    paste(formatC(bbj_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(tpmi_thr, format = "e", digits = 6), collapse = ", ")
  ))

  bbj_nov_direct <- supple3_cate_get_novelty_label(
    bbj$novel_relative_to_bbj_direct, bbj$mapping_status_direct, "bbj"
  )
  tpmi_nov_direct <- supple3_cate_get_novelty_label(
    tpmi$novel_relative_to_tpmi_direct, tpmi$mapping_status_direct, "tpmi"
  )
  bbj_nov_any <- supple3_cate_get_novelty_label(
    bbj$novel_relative_to_bbj_any, bbj$mapping_status_any, "bbj"
  )
  tpmi_nov_any <- supple3_cate_get_novelty_label(
    tpmi$novel_relative_to_tpmi_any, tpmi$mapping_status_any, "tpmi"
  )
  combined_nov_direct <- supple3_cate_get_combined_novelty(bbj_nov_direct, tpmi_nov_direct)
  combined_nov_any <- supple3_cate_get_combined_novelty(bbj_nov_any, tpmi_nov_any)

  bbj_rep_nom_direct <- supple3_cate_get_replication_label(
    bbj$mapping_status_direct, bbj$replication_testable_direct,
    bbj$not_testable_reason_direct, bbj$replicated_nominal_direct,
    bbj$direction_concordant_direct, "bbj"
  )
  tpmi_rep_nom_direct <- supple3_cate_get_replication_label(
    tpmi$mapping_status_direct, tpmi$replication_testable_direct,
    tpmi$not_testable_reason_direct, tpmi$replicated_nominal_direct,
    tpmi$direction_concordant_direct, "tpmi"
  )
  any_rep_nom_direct <- supple3_cate_get_any_cohort_replication(
    bbj_rep_nom_direct, tpmi_rep_nom_direct
  )
  bbj_rep_nom_any <- supple3_cate_expand_rep_with_additional(
    bbj_rep_nom_direct, bbj$additionally_replicated_nominal_from_other
  )
  tpmi_rep_nom_any <- supple3_cate_expand_rep_with_additional(
    tpmi_rep_nom_direct, tpmi$additionally_replicated_nominal_from_other
  )
  any_rep_nom_any <- supple3_cate_get_any_cohort_replication(
    bbj_rep_nom_any, tpmi_rep_nom_any
  )

  bbj_rep_bonf_direct <- supple3_cate_get_replication_label(
    bbj$mapping_status_direct, bbj$replication_testable_direct,
    bbj$not_testable_reason_direct, bbj$replicated_primary_direct,
    bbj$direction_concordant_direct, "bbj"
  )
  tpmi_rep_bonf_direct <- supple3_cate_get_replication_label(
    tpmi$mapping_status_direct, tpmi$replication_testable_direct,
    tpmi$not_testable_reason_direct, tpmi$replicated_primary_direct,
    tpmi$direction_concordant_direct, "tpmi"
  )
  any_rep_bonf_direct <- supple3_cate_get_any_cohort_replication(
    bbj_rep_bonf_direct, tpmi_rep_bonf_direct
  )
  bbj_rep_bonf_any <- supple3_cate_expand_rep_with_additional(
    bbj_rep_bonf_direct, bbj$additionally_replicated_primary_from_other
  )
  tpmi_rep_bonf_any <- supple3_cate_expand_rep_with_additional(
    tpmi_rep_bonf_direct, tpmi$additionally_replicated_primary_from_other
  )
  any_rep_bonf_any <- supple3_cate_get_any_cohort_replication(
    bbj_rep_bonf_any, tpmi_rep_bonf_any
  )

  integ_direct <- supple3_cate_get_integrated_replication(
    testable_list = list(bbj$replication_testable_direct, tpmi$replication_testable_direct),
    concordant_list = list(bbj$direction_concordant_direct, tpmi$direction_concordant_direct),
    exact_p_list = list(bbj$best_testable_p_direct, tpmi$best_testable_p_direct)
  )
  integ_direct_label <- integ_direct$label
  idx_k0 <- which(integ_direct$k == 0L)
  if (length(idx_k0) > 0L) integ_direct_label[idx_k0] <- any_rep_nom_direct[idx_k0]

  integ_any <- supple3_cate_get_integrated_replication(
    testable_list = list(
      supple3_cate_pick_any_match(bbj$replication_testable_direct, bbj$replication_testable_other),
      supple3_cate_pick_any_match(tpmi$replication_testable_direct, tpmi$replication_testable_other)
    ),
    concordant_list = list(
      supple3_cate_pick_any_match(bbj$direction_concordant_direct, bbj$direction_concordant_other),
      supple3_cate_pick_any_match(tpmi$direction_concordant_direct, tpmi$direction_concordant_other)
    ),
    exact_p_list = list(
      supple3_cate_pick_any_p(
        bbj$replication_testable_direct, bbj$best_testable_p_direct,
        bbj$replication_testable_other, bbj$best_testable_p_other
      ),
      supple3_cate_pick_any_p(
        tpmi$replication_testable_direct, tpmi$best_testable_p_direct,
        tpmi$replication_testable_other, tpmi$best_testable_p_other
      )
    )
  )
  integ_any_label <- integ_any$label
  idx_k0_any <- which(integ_any$k == 0L)
  if (length(idx_k0_any) > 0L) integ_any_label[idx_k0_any] <- any_rep_nom_any[idx_k0_any]

  message(sprintf(
    "Integrated EAS Cate (direct) | M_global=%d | 0.05/M=%s",
    integ_direct$M_global,
    format(if (integ_direct$M_global > 0L) 0.05 / integ_direct$M_global else NA_real_, digits = 6)
  ))
  message(sprintf(
    "Integrated EAS Cate (direct+other) | M_global=%d | 0.05/M=%s",
    integ_any$M_global,
    format(if (integ_any$M_global > 0L) 0.05 / integ_any$M_global else NA_real_, digits = 6)
  ))

  out <- data.table(
    Phenotype = bbj$koges_phenotype,
    Chromosome = bbj$koges_chr,
    `Base pair position` = bbj$koges_bp,
    rsID = bbj$koges_lead_variant,
    `P-value in KoGES` = supple3_cate_fmt_p(bbj$koges_p),
    `Minimum P-value within ±1 Mb in BBJ (직접 매칭)` = supple3_cate_fmt_p(bbj$min_window_p_direct),
    `Minimum P-value within ±1 Mb in BBJ (other matching)` = supple3_cate_fmt_p(bbj$min_window_p_other),
    `Minimum P-value within ±1 Mb in TPMI (직접 매칭)` = supple3_cate_fmt_p(tpmi$min_window_p_direct),
    `Minimum P-value within ±1 Mb in TPMI (other matching)` = supple3_cate_fmt_p(tpmi$min_window_p_other),
    `Novel relative to BBJ (직접 매칭)` = bbj_nov_direct,
    `Novel relative to TPMI (직접 매칭)` = tpmi_nov_direct,
    `Additionally known from other matching - BBJ` = supple3_cate_yes_dash(bbj$additionally_known_from_other),
    `Additionally known from other matching - TPMI` = supple3_cate_yes_dash(tpmi$additionally_known_from_other),
    `Novel relative to assessed EastAsian-ancestry cohorts (직접 매칭)` = combined_nov_direct,
    `Novel relative to assessed EastAsian-ancestry cohorts (직접 매칭 + other)` = combined_nov_any,
    `Variant-level replication - BBJ (직접 매칭)` = bbj_rep_nom_direct,
    `Variant-level replication - TPMI (직접 매칭)` = tpmi_rep_nom_direct,
    `Nominal replication in ≥1 cohort (직접 매칭)` = any_rep_nom_direct,
    `Additionally replicated nominally from other matching - BBJ` = supple3_cate_yes_dash(
      bbj$additionally_replicated_nominal_from_other
    ),
    `Additionally replicated nominally from other matching - TPMI` = supple3_cate_yes_dash(
      tpmi$additionally_replicated_nominal_from_other
    ),
    `Nominal replication in ≥1 cohort (직접 매칭 + other)` = any_rep_nom_any,
    `Bonferroni replication - BBJ (직접 매칭)` = bbj_rep_bonf_direct,
    `Bonferroni replication - TPMI (직접 매칭)` = tpmi_rep_bonf_direct,
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort (직접 매칭)` =
      any_rep_bonf_direct,
    `Additionally Bonferroni-replicated from other matching - BBJ` = supple3_cate_yes_dash(
      bbj$additionally_replicated_primary_from_other
    ),
    `Additionally Bonferroni-replicated from other matching - TPMI` = supple3_cate_yes_dash(
      tpmi$additionally_replicated_primary_from_other
    ),
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort (직접 매칭 + other)` =
      any_rep_bonf_any,
    `Integrated EAS replication (직접 매칭)` = integ_direct_label,
    `Integrated EAS replication (직접 매칭 + other)` = integ_any_label
  )
  setorderv(out, c("Phenotype", "Chromosome", "Base pair position"))
  out
}

# ------------------------------------------------------------------------------
# build_summary_eur_cate(ukb, finngen, mvp, ukb_pm, finngen_pm, mvp_pm)
#   EUR Cate novelty/per_match use cohort prefixes (ukb_/finngen_/mvp_).
# ------------------------------------------------------------------------------
build_summary_eur_cate <- function(ukb, finngen, mvp, ukb_pm, finngen_pm, mvp_pm) {
  ukb <- as.data.table(ukb)
  finngen <- as.data.table(finngen)
  mvp <- as.data.table(mvp)
  ukb_pm <- as.data.table(ukb_pm)
  finngen_pm <- as.data.table(finngen_pm)
  mvp_pm <- as.data.table(mvp_pm)

  if (!identical(ukb$lead_row_id, finngen$lead_row_id) ||
      !identical(ukb$lead_row_id, mvp$lead_row_id)) {
    stop("UKB / FinnGen / MVP lead_row_id order differs.")
  }
  if (!identical(ukb$koges_lead_variant, finngen$koges_lead_variant) ||
      !identical(ukb$koges_lead_variant, mvp$koges_lead_variant) ||
      !identical(ukb$koges_chr, finngen$koges_chr) ||
      !identical(ukb$koges_chr, mvp$koges_chr) ||
      !identical(ukb$koges_bp, finngen$koges_bp) ||
      !identical(ukb$koges_bp, mvp$koges_bp)) {
    stop("UKB / FinnGen / MVP KoGES lead coordinates differ.")
  }

  need_common <- c(
    "lead_row_id", "koges_phenotype", "koges_pheno_code",
    "koges_chr", "koges_bp", "koges_lead_variant", "koges_p",
    "mapping_status_direct", "mapping_status_any",
    "additionally_known_from_other",
    "replication_testable_direct", "replication_testable_other",
    "not_testable_reason_direct",
    "replicated_nominal_direct", "replicated_primary_direct",
    "direction_concordant_direct", "direction_concordant_other",
    "best_testable_p_direct", "best_testable_p_other",
    "additionally_replicated_nominal_from_other",
    "additionally_replicated_primary_from_other",
    "replication_threshold_primary_direct"
  )
  for (nm in c("ukb", "finngen", "mvp")) {
    pfx <- nm
    need_lead <- c(
      need_common,
      paste0("novel_relative_to_", pfx, "_direct"),
      paste0("novel_relative_to_", pfx, "_any")
    )
    miss <- setdiff(need_lead, names(get(nm)))
    if (length(miss)) stop(toupper(nm), " lead TSV missing: ", paste(miss, collapse = ", "))
  }

  attach_window_p <- function(lead_dt, pm, label, pfx) {
    p_direct <- supple3_cate_min_window_p(
      pm, paste0(pfx, "_top_window_p"), paste0(pfx, "_mapping_status"), "direct"
    )
    p_other <- supple3_cate_min_window_p(
      pm, paste0(pfx, "_top_window_p"), paste0(pfx, "_mapping_status"), "other"
    )
    lead_dt[p_direct, min_window_p_direct := i.min_p, on = "lead_row_id"]
    lead_dt[p_other, min_window_p_other := i.min_p, on = "lead_row_id"]
    message(label, " direct window-P leads: ", nrow(p_direct), " / other: ", nrow(p_other))
    lead_dt
  }
  ukb <- attach_window_p(ukb, ukb_pm, "UKB", "ukb")
  finngen <- attach_window_p(finngen, finngen_pm, "FinnGen", "finngen")
  mvp <- attach_window_p(mvp, mvp_pm, "MVP", "mvp")

  ukb_thr <- unique(na.omit(ukb$replication_threshold_primary_direct))
  finngen_thr <- unique(na.omit(finngen$replication_threshold_primary_direct))
  mvp_thr <- unique(na.omit(mvp$replication_threshold_primary_direct))
  message(sprintf(
    "EUR Cate direct-match Bonferroni thresholds | UKB: %s | FinnGen: %s | MVP: %s",
    paste(formatC(ukb_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(finngen_thr, format = "e", digits = 6), collapse = ", "),
    paste(formatC(mvp_thr, format = "e", digits = 6), collapse = ", ")
  ))

  ukb_nov_direct <- supple3_cate_get_novelty_label(
    ukb$novel_relative_to_ukb_direct, ukb$mapping_status_direct, "ukb"
  )
  finngen_nov_direct <- supple3_cate_get_novelty_label(
    finngen$novel_relative_to_finngen_direct, finngen$mapping_status_direct, "finngen"
  )
  mvp_nov_direct <- supple3_cate_get_novelty_label(
    mvp$novel_relative_to_mvp_direct, mvp$mapping_status_direct, "mvp"
  )
  ukb_nov_any <- supple3_cate_get_novelty_label(
    ukb$novel_relative_to_ukb_any, ukb$mapping_status_any, "ukb"
  )
  finngen_nov_any <- supple3_cate_get_novelty_label(
    finngen$novel_relative_to_finngen_any, finngen$mapping_status_any, "finngen"
  )
  mvp_nov_any <- supple3_cate_get_novelty_label(
    mvp$novel_relative_to_mvp_any, mvp$mapping_status_any, "mvp"
  )
  combined_nov_direct <- supple3_cate_get_combined_novelty(
    ukb_nov_direct, finngen_nov_direct, mvp_nov_direct
  )
  combined_nov_any <- supple3_cate_get_combined_novelty(
    ukb_nov_any, finngen_nov_any, mvp_nov_any
  )

  ukb_rep_nom_direct <- supple3_cate_get_replication_label(
    ukb$mapping_status_direct, ukb$replication_testable_direct,
    ukb$not_testable_reason_direct, ukb$replicated_nominal_direct,
    ukb$direction_concordant_direct, "ukb"
  )
  finngen_rep_nom_direct <- supple3_cate_get_replication_label(
    finngen$mapping_status_direct, finngen$replication_testable_direct,
    finngen$not_testable_reason_direct, finngen$replicated_nominal_direct,
    finngen$direction_concordant_direct, "finngen"
  )
  mvp_rep_nom_direct <- supple3_cate_get_replication_label(
    mvp$mapping_status_direct, mvp$replication_testable_direct,
    mvp$not_testable_reason_direct, mvp$replicated_nominal_direct,
    mvp$direction_concordant_direct, "mvp"
  )
  any_rep_nom_direct <- supple3_cate_get_any_cohort_replication(
    ukb_rep_nom_direct, finngen_rep_nom_direct, mvp_rep_nom_direct
  )
  ukb_rep_nom_any <- supple3_cate_expand_rep_with_additional(
    ukb_rep_nom_direct, ukb$additionally_replicated_nominal_from_other
  )
  finngen_rep_nom_any <- supple3_cate_expand_rep_with_additional(
    finngen_rep_nom_direct, finngen$additionally_replicated_nominal_from_other
  )
  mvp_rep_nom_any <- supple3_cate_expand_rep_with_additional(
    mvp_rep_nom_direct, mvp$additionally_replicated_nominal_from_other
  )
  any_rep_nom_any <- supple3_cate_get_any_cohort_replication(
    ukb_rep_nom_any, finngen_rep_nom_any, mvp_rep_nom_any
  )

  ukb_rep_bonf_direct <- supple3_cate_get_replication_label(
    ukb$mapping_status_direct, ukb$replication_testable_direct,
    ukb$not_testable_reason_direct, ukb$replicated_primary_direct,
    ukb$direction_concordant_direct, "ukb"
  )
  finngen_rep_bonf_direct <- supple3_cate_get_replication_label(
    finngen$mapping_status_direct, finngen$replication_testable_direct,
    finngen$not_testable_reason_direct, finngen$replicated_primary_direct,
    finngen$direction_concordant_direct, "finngen"
  )
  mvp_rep_bonf_direct <- supple3_cate_get_replication_label(
    mvp$mapping_status_direct, mvp$replication_testable_direct,
    mvp$not_testable_reason_direct, mvp$replicated_primary_direct,
    mvp$direction_concordant_direct, "mvp"
  )
  any_rep_bonf_direct <- supple3_cate_get_any_cohort_replication(
    ukb_rep_bonf_direct, finngen_rep_bonf_direct, mvp_rep_bonf_direct
  )
  ukb_rep_bonf_any <- supple3_cate_expand_rep_with_additional(
    ukb_rep_bonf_direct, ukb$additionally_replicated_primary_from_other
  )
  finngen_rep_bonf_any <- supple3_cate_expand_rep_with_additional(
    finngen_rep_bonf_direct, finngen$additionally_replicated_primary_from_other
  )
  mvp_rep_bonf_any <- supple3_cate_expand_rep_with_additional(
    mvp_rep_bonf_direct, mvp$additionally_replicated_primary_from_other
  )
  any_rep_bonf_any <- supple3_cate_get_any_cohort_replication(
    ukb_rep_bonf_any, finngen_rep_bonf_any, mvp_rep_bonf_any
  )

  integ_direct <- supple3_cate_get_integrated_replication(
    testable_list = list(
      ukb$replication_testable_direct, finngen$replication_testable_direct,
      mvp$replication_testable_direct
    ),
    concordant_list = list(
      ukb$direction_concordant_direct, finngen$direction_concordant_direct,
      mvp$direction_concordant_direct
    ),
    exact_p_list = list(
      ukb$best_testable_p_direct, finngen$best_testable_p_direct,
      mvp$best_testable_p_direct
    )
  )
  integ_direct_label <- integ_direct$label
  idx_k0 <- which(integ_direct$k == 0L)
  if (length(idx_k0) > 0L) integ_direct_label[idx_k0] <- any_rep_nom_direct[idx_k0]

  integ_any <- supple3_cate_get_integrated_replication(
    testable_list = list(
      supple3_cate_pick_any_match(ukb$replication_testable_direct, ukb$replication_testable_other),
      supple3_cate_pick_any_match(finngen$replication_testable_direct, finngen$replication_testable_other),
      supple3_cate_pick_any_match(mvp$replication_testable_direct, mvp$replication_testable_other)
    ),
    concordant_list = list(
      supple3_cate_pick_any_match(ukb$direction_concordant_direct, ukb$direction_concordant_other),
      supple3_cate_pick_any_match(finngen$direction_concordant_direct, finngen$direction_concordant_other),
      supple3_cate_pick_any_match(mvp$direction_concordant_direct, mvp$direction_concordant_other)
    ),
    exact_p_list = list(
      supple3_cate_pick_any_p(
        ukb$replication_testable_direct, ukb$best_testable_p_direct,
        ukb$replication_testable_other, ukb$best_testable_p_other
      ),
      supple3_cate_pick_any_p(
        finngen$replication_testable_direct, finngen$best_testable_p_direct,
        finngen$replication_testable_other, finngen$best_testable_p_other
      ),
      supple3_cate_pick_any_p(
        mvp$replication_testable_direct, mvp$best_testable_p_direct,
        mvp$replication_testable_other, mvp$best_testable_p_other
      )
    )
  )
  integ_any_label <- integ_any$label
  idx_k0_any <- which(integ_any$k == 0L)
  if (length(idx_k0_any) > 0L) integ_any_label[idx_k0_any] <- any_rep_nom_any[idx_k0_any]

  message(sprintf(
    "Integrated EUR Cate (direct) | M_global=%d | 0.05/M=%s",
    integ_direct$M_global,
    format(if (integ_direct$M_global > 0L) 0.05 / integ_direct$M_global else NA_real_, digits = 6)
  ))
  message(sprintf(
    "Integrated EUR Cate (direct+other) | M_global=%d | 0.05/M=%s",
    integ_any$M_global,
    format(if (integ_any$M_global > 0L) 0.05 / integ_any$M_global else NA_real_, digits = 6)
  ))

  out <- data.table(
    Phenotype = ukb$koges_phenotype,
    Chromosome = ukb$koges_chr,
    `Base pair position` = ukb$koges_bp,
    rsID = ukb$koges_lead_variant,
    `P-value in KoGES` = supple3_cate_fmt_p(ukb$koges_p),
    `Minimum P-value within ±1 Mb in UKB (직접 매칭)` = supple3_cate_fmt_p(ukb$min_window_p_direct),
    `Minimum P-value within ±1 Mb in UKB (other matching)` = supple3_cate_fmt_p(ukb$min_window_p_other),
    `Minimum P-value within ±1 Mb in FinnGen (직접 매칭)` = supple3_cate_fmt_p(finngen$min_window_p_direct),
    `Minimum P-value within ±1 Mb in FinnGen (other matching)` = supple3_cate_fmt_p(finngen$min_window_p_other),
    `Minimum P-value within ±1 Mb in MVP (직접 매칭)` = supple3_cate_fmt_p(mvp$min_window_p_direct),
    `Minimum P-value within ±1 Mb in MVP (other matching)` = supple3_cate_fmt_p(mvp$min_window_p_other),
    `Novel relative to UKB (직접 매칭)` = ukb_nov_direct,
    `Novel relative to FinnGen (직접 매칭)` = finngen_nov_direct,
    `Novel relative to MVP (직접 매칭)` = mvp_nov_direct,
    `Additionally known from other matching - UKB` = supple3_cate_yes_dash(ukb$additionally_known_from_other),
    `Additionally known from other matching - FinnGen` = supple3_cate_yes_dash(finngen$additionally_known_from_other),
    `Additionally known from other matching - MVP` = supple3_cate_yes_dash(mvp$additionally_known_from_other),
    `Novel relative to assessed European-ancestry cohorts (직접 매칭)` = combined_nov_direct,
    `Novel relative to assessed European-ancestry cohorts (직접 매칭 + other)` = combined_nov_any,
    `Variant-level replication - UKB (직접 매칭)` = ukb_rep_nom_direct,
    `Variant-level replication - FinnGen (직접 매칭)` = finngen_rep_nom_direct,
    `Variant-level replication - MVP (직접 매칭)` = mvp_rep_nom_direct,
    `Nominal replication in ≥1 cohort (직접 매칭)` = any_rep_nom_direct,
    `Additionally replicated nominally from other matching - UKB` = supple3_cate_yes_dash(
      ukb$additionally_replicated_nominal_from_other
    ),
    `Additionally replicated nominally from other matching - FinnGen` = supple3_cate_yes_dash(
      finngen$additionally_replicated_nominal_from_other
    ),
    `Additionally replicated nominally from other matching - MVP` = supple3_cate_yes_dash(
      mvp$additionally_replicated_nominal_from_other
    ),
    `Nominal replication in ≥1 cohort (직접 매칭 + other)` = any_rep_nom_any,
    `Bonferroni replication - UKB (직접 매칭)` = ukb_rep_bonf_direct,
    `Bonferroni replication - FinnGen (직접 매칭)` = finngen_rep_bonf_direct,
    `Bonferroni replication - MVP (직접 매칭)` = mvp_rep_bonf_direct,
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort (직접 매칭)` = any_rep_bonf_direct,
    `Additionally Bonferroni-replicated from other matching - UKB` = supple3_cate_yes_dash(
      ukb$additionally_replicated_primary_from_other
    ),
    `Additionally Bonferroni-replicated from other matching - FinnGen` = supple3_cate_yes_dash(
      finngen$additionally_replicated_primary_from_other
    ),
    `Additionally Bonferroni-replicated from other matching - MVP` = supple3_cate_yes_dash(
      mvp$additionally_replicated_primary_from_other
    ),
    `Cohort-specific Bonferroni-significant replication in ≥1 cohort (직접 매칭 + other)` = any_rep_bonf_any,
    `Integrated EUR replication (직접 매칭)` = integ_direct_label,
    `Integrated EUR replication (직접 매칭 + other)` = integ_any_label
  )
  setorderv(out, c("Phenotype", "Chromosome", "Base pair position"))
  out
}

# ------------------------------------------------------------------------------
# build_cohort_detail_cate(nov, pm, cohort, prefix, eaf_lookup = NULL)
#   eaf_lookup: data.table(lead_row_id, koges_eaf_lookup) or NULL if koges_eaf present
#   For EUR cohorts, pass per_match already renamed OR set rename_tpmi=TRUE.
# ------------------------------------------------------------------------------
build_cohort_detail_cate <- function(
    nov, pm, cohort, prefix,
    eaf_lookup = NULL,
    rename_tpmi = FALSE
) {
  nov <- as.data.table(nov)
  pm <- as.data.table(pm)
  prefix <- tolower(as.character(prefix))
  if (isTRUE(rename_tpmi)) {
    pm <- supple3_cate_rename_tpmi_prefix(pm, prefix)
  }
  pm_one <- supple3_cate_pick_direct_one(pm, prefix, cohort)

  C <- function(suffix) paste0(prefix, suffix)
  join_cols <- intersect(
    c(
      "lead_row_id",
      "koges_beta_aligned_a1",
      C("_source_trait"),
      C("_mapping_status"),
      C("_n_variants_in_window"),
      C("_n_gws_hits_in_window"),
      C("_top_window_snp"),
      C("_top_window_rsid_koges"),
      C("_top_window_chr"),
      C("_top_window_pos"),
      C("_top_window_p"),
      C("_top_window_distance_bp"),
      C("_exact_chrpos_present"),
      C("_exact_snp"),
      C("_exact_rsid_koges"),
      C("_exact_a1_koges"),
      C("_exact_a2_koges"),
      C("_exact_beta_aligned_a1"),
      C("_exact_eaf_aligned_a1"),
      C("_exact_se"),
      C("_exact_p")
    ),
    names(pm_one)
  )
  dt <- copy(nov)
  if (length(join_cols) > 1L && nrow(pm_one) > 0L) {
    drop_overlap <- setdiff(join_cols, "lead_row_id")
    drop_overlap <- drop_overlap[drop_overlap %in% names(dt)]
    if (length(drop_overlap) > 0L) dt[, (drop_overlap) := NULL]
    dt <- merge(dt, pm_one[, ..join_cols], by = "lead_row_id", all.x = TRUE, sort = FALSE)
  }

  identity_req <- c(
    "lead_row_id", "koges_phenotype", "koges_pheno_code",
    "koges_lead_variant", "koges_variant_id", "koges_chr", "koges_bp",
    "koges_nearest_gene", "koges_ea", "koges_nea", "koges_beta",
    "koges_se", "koges_p"
  )
  miss_id <- setdiff(identity_req, names(dt))
  if (length(miss_id) > 0L) {
    stop(cohort, " novelty missing identity column(s): ", paste(miss_id, collapse = ", "))
  }

  novelty_status <- supple3_cate_first_present(dt, c(
    C("_novelty_status_direct"),
    C("_novelty_status"),
    "novelty_status_direct",
    "novelty_status"
  ))
  if (is.null(novelty_status)) novelty_status <- rep(NA_character_, nrow(dt))

  repf <- supple3_cate_resolve_replication_fields(dt)
  if ("koges_eaf" %in% names(dt) && !all(is.na(dt$koges_eaf))) {
    koges_eaf <- as.numeric(dt$koges_eaf)
  } else {
    if (is.null(eaf_lookup)) {
      stop(cohort, ": koges_eaf absent and eaf_lookup not provided.")
    }
    eaf_lookup <- as.data.table(eaf_lookup)
    idx <- match(as.integer(dt$lead_row_id), eaf_lookup$lead_row_id)
    if (anyNA(idx)) {
      stop(cohort, ": ", sum(is.na(idx)), " lead_row_id(s) missing from EAF lookup.")
    }
    koges_eaf <- eaf_lookup$koges_eaf_lookup[idx]
  }
  koges_beta_a1 <- if ("koges_beta_aligned_a1" %in% names(dt)) {
    dt$koges_beta_aligned_a1
  } else {
    rep(NA_real_, nrow(dt))
  }
  get_pm <- function(suffix) {
    col <- C(suffix)
    if (col %in% names(dt)) dt[[col]] else rep(NA, nrow(dt))
  }

  data.table(
    `KoGES Phenotype` = dt$koges_phenotype,
    `KoGES Lead Variant (rsID)` = dt$koges_lead_variant,
    `KoGES Variant ID` = dt$koges_variant_id,
    `KoGES Chromosome` = dt$koges_chr,
    `KoGES Position (bp)` = dt$koges_bp,
    `Nearest Gene` = dt$koges_nearest_gene,
    `KoGES Effect Allele` = dt$koges_ea,
    `KoGES Other Allele` = dt$koges_nea,
    `KoGES EAF` = koges_eaf,
    `KoGES Beta` = dt$koges_beta,
    `KoGES Standard Error` = dt$koges_se,
    `KoGES P Value` = dt$koges_p,
    `Lowest-P External Variant rsID Within ±1 Mb` = get_pm("_top_window_rsid_koges"),
    `Lowest-P External Variant Chromosome` = get_pm("_top_window_chr"),
    `Lowest-P External Variant Position (bp)` = get_pm("_top_window_pos"),
    `Minimum External P Value Within ±1 Mb` = get_pm("_top_window_p"),
    `Locus Novelty Assessment` = novelty_status,
    `Exact-Position External Variant Available` = supple3_cate_yes_no(get_pm("_exact_chrpos_present")),
    `Exact-Position External Variant rsID` = get_pm("_exact_rsid_koges"),
    `Aligned Effect Allele (A1)` = get_pm("_exact_a1_koges"),
    `Aligned Other Allele (A2)` = get_pm("_exact_a2_koges"),
    `KoGES Beta Aligned to A1` = koges_beta_a1,
    `External Beta Aligned to A1` = get_pm("_exact_beta_aligned_a1"),
    `External EAF Aligned to A1` = get_pm("_exact_eaf_aligned_a1"),
    `External Standard Error` = get_pm("_exact_se"),
    `External P Value at Exact Position` = get_pm("_exact_p"),
    `Exact-Variant Replication Testable` = supple3_cate_yes_no(repf$replication_testable),
    `Effect Direction Concordant` = supple3_cate_yes_no(repf$direction_concordant),
    `Nominal Replication` = supple3_cate_yes_no(repf$replicated_nominal),
    `Primary Bonferroni Threshold` = repf$replication_threshold,
    `Bonferroni-Significant Replication` = supple3_cate_yes_no(repf$replicated_primary),
    `Replication Assessment` = repf$replication_status
  )
}
