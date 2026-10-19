# ==============================================================================
# Build EAS/EUR panel-replication tables from panel-meta TSVs (in memory).
#
# Does NOT re-run meta-analysis. Validates stored replication flags and returns
# lead-level tables used by ExtendedData2_{cont,cate}.R
# (Replication_EAS_meta / Replication_EUR_meta sheets).
#
# Usage
#   source(".../panel_replication_from_tsv.R")
#   tabs <- build_panel_replication_tables(eas_path, eur_path)
#   # tabs$eas, tabs$eur  -> lead-level data.tables
# ==============================================================================

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("package 'data.table' is required")
}
suppressPackageStartupMessages(library(data.table))

.panel_repl_as_flag <- function(x) {
  if (is.logical(x)) return(x)
  y <- toupper(trimws(as.character(x)))
  out <- rep(NA, length(y))
  out[y %chin% c("TRUE", "T", "1", "YES", "Y")] <- TRUE
  out[y %chin% c("FALSE", "F", "0", "NO", "N")] <- FALSE
  out
}

.panel_repl_fmt_npct <- function(n, d) {
  if (!is.finite(d) || d <= 0) return(NA_character_)
  sprintf("%d (%.1f%%)", as.integer(n), 100 * n / d)
}

.panel_repl_check_cols <- function(dt, cols, label) {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) {
    stop(label, " missing required columns: ", paste(miss, collapse = ", "))
  }
}

.panel_repl_panel_info <- function(dt, panel) {
  req <- c(
    paste0(panel, "_panel_analyzed"),
    paste0(panel, "_n_contributing_cohorts"),
    paste0(panel, "_panel_P"),
    paste0(panel, "_panel_koges_beta_aligned_a1"),
    paste0(panel, "_panel_direction_concordant"),
    paste0(panel, "_panel_replication_nominal"),
    paste0(panel, "_panel_replication_Bonferroni"),
    paste0(panel, "_panel_M_global"),
    paste0(panel, "_panel_alpha_Bonferroni"),
    paste0(panel, "_sensitivity_signedZ_Z"),
    paste0(panel, "_sensitivity_signedZ_P")
  )
  .panel_repl_check_cols(dt, req, panel)

  analyzed <- .panel_repl_as_flag(dt[[paste0(panel, "_panel_analyzed")]])
  direction <- .panel_repl_as_flag(dt[[paste0(panel, "_panel_direction_concordant")]])
  rep_nom <- .panel_repl_as_flag(dt[[paste0(panel, "_panel_replication_nominal")]])
  rep_bonf <- .panel_repl_as_flag(dt[[paste0(panel, "_panel_replication_Bonferroni")]])
  k <- as.integer(dt[[paste0(panel, "_n_contributing_cohorts")]])

  Mv <- unique(na.omit(as.integer(dt[[paste0(panel, "_panel_M_global")]])))
  av <- unique(na.omit(as.numeric(dt[[paste0(panel, "_panel_alpha_Bonferroni")]])))
  if (length(Mv) != 1L || length(av) != 1L) {
    stop(panel, ": M_global/alpha must each be unique.")
  }
  M <- Mv[[1L]]
  alpha <- av[[1L]]
  if (sum(analyzed %in% TRUE) != M) {
    stop(panel, ": M_global != number of analyzed lead rows.")
  }
  if (abs(alpha - 0.05 / M) > 1e-12) {
    stop(panel, ": stored Bonferroni alpha != 0.05/M_global.")
  }

  pp <- as.numeric(dt[[paste0(panel, "_panel_P")]])
  expected_nom <- rep(NA, nrow(dt))
  expected_bonf <- rep(NA, nrow(dt))
  idx <- analyzed %in% TRUE
  expected_nom[idx] <- (direction[idx] %in% TRUE) & is.finite(pp[idx]) & pp[idx] <= 0.05
  expected_bonf[idx] <- (direction[idx] %in% TRUE) & is.finite(pp[idx]) & pp[idx] <= alpha
  if (any(expected_nom[idx] != rep_nom[idx], na.rm = TRUE)) {
    stop(panel, ": nominal replication flag mismatch.")
  }
  if (any(expected_bonf[idx] != rep_bonf[idx], na.rm = TRUE)) {
    stop(panel, ": Bonferroni replication flag mismatch.")
  }

  z <- as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_Z")]])
  zp <- as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_P")]])
  kb <- as.numeric(dt[[paste0(panel, "_panel_koges_beta_aligned_a1")]])
  zdir <- rep(NA, nrow(dt))
  znom <- rep(NA, nrow(dt))
  zbonf <- rep(NA, nrow(dt))
  iz <- idx & is.finite(z) & is.finite(zp) & is.finite(kb) & z != 0 & kb != 0
  zdir[iz] <- sign(z[iz]) == sign(kb[iz])
  znom[iz] <- zdir[iz] & zp[iz] <= 0.05
  zbonf[iz] <- zdir[iz] & zp[iz] <= alpha

  list(
    analyzed = analyzed, direction = direction, rep_nom = rep_nom, rep_bonf = rep_bonf,
    k = k, M = M, alpha = alpha, signedZ_dir = zdir,
    signedZ_nom = znom, signedZ_bonf = zbonf
  )
}

.panel_repl_make_extended <- function(dt, panel, info, other_analyzed) {
  req <- c(
    "meta_exclusion_reason",
    paste0(panel, "_contributing_cohorts"),
    paste0(panel, "_panel_analyzed"),
    paste0(panel, "_panel_koges_beta_aligned_a1"),
    paste0(panel, "_panel_beta"), paste0(panel, "_panel_SE"), paste0(panel, "_panel_P"),
    paste0(panel, "_panel_direction_concordant"),
    paste0(panel, "_panel_replication_nominal"),
    paste0(panel, "_panel_replication_Bonferroni"),
    paste0(panel, "_sensitivity_signedZ_Z"), paste0(panel, "_sensitivity_signedZ_P"),
    paste0(panel, "_source_direction_pattern"),
    paste0(panel, "_panel_Q_P"), paste0(panel, "_panel_I2")
  )
  .panel_repl_check_cols(dt, req, panel)

  data.table(
    lead_row_id = dt$lead_row_id,
    koges_pheno_code = dt$koges_pheno_code,
    koges_phenotype = dt$koges_phenotype,
    koges_lead_variant = dt$koges_lead_variant,
    koges_chr = dt$koges_chr,
    koges_bp = dt$koges_bp,
    meta_exclusion_reason = dt$meta_exclusion_reason,
    panel_analyzed = info$analyzed,
    EAS_EUR_both_testable = (info$analyzed %in% TRUE) & (other_analyzed %in% TRUE),
    contributing_cohorts = dt[[paste0(panel, "_contributing_cohorts")]],
    koges_beta_aligned_a1 = as.numeric(dt[[paste0(panel, "_panel_koges_beta_aligned_a1")]]),
    external_panel_beta = as.numeric(dt[[paste0(panel, "_panel_beta")]]),
    external_panel_SE = as.numeric(dt[[paste0(panel, "_panel_SE")]]),
    external_panel_P = as.numeric(dt[[paste0(panel, "_panel_P")]]),
    direction_concordant_with_KoGES = info$direction,
    replication_nominal = info$rep_nom,
    replication_Bonferroni = info$rep_bonf,
    signedZ_Z = as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_Z")]]),
    signedZ_P = as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_P")]]),
    signedZ_replication_Bonferroni = info$signedZ_bonf,
    external_cohort_direction_pattern = dt[[paste0(panel, "_source_direction_pattern")]],
    external_cohort_Q_P = as.numeric(dt[[paste0(panel, "_panel_Q_P")]]),
    external_cohort_I2 = as.numeric(dt[[paste0(panel, "_panel_I2")]])
  )
}

#' Build validated EAS/EUR panel-replication tables from panel-meta TSVs.
#'
#' @param eas_path Path to EAS panel-meta TSV
#' @param eur_path Path to EUR panel-meta TSV
#' @return Named list with eas, eur, main_* blocks, N_common, alpha_common, ei, ui
build_panel_replication_tables <- function(eas_path, eur_path) {
  if (!file.exists(eas_path)) stop("EAS panel-meta TSV not found: ", eas_path)
  if (!file.exists(eur_path)) stop("EUR panel-meta TSV not found: ", eur_path)

  E <- fread(eas_path, sep = "\t", na.strings = c("NA", ""))
  U <- fread(eur_path, sep = "\t", na.strings = c("NA", ""))

  id_cols <- c(
    "lead_row_id", "koges_pheno_code", "koges_phenotype",
    "koges_lead_variant", "koges_chr", "koges_bp"
  )
  .panel_repl_check_cols(E, id_cols, "EAS")
  .panel_repl_check_cols(U, id_cols, "EUR")
  if (anyDuplicated(E$lead_row_id)) stop("EAS lead_row_id is not unique.")
  if (anyDuplicated(U$lead_row_id)) stop("EUR lead_row_id is not unique.")
  if (!setequal(E$lead_row_id, U$lead_row_id)) stop("EAS/EUR lead_row_id sets differ.")

  setkey(E, lead_row_id)
  setkey(U, lead_row_id)
  for (z in id_cols[-1]) {
    if (!identical(as.character(E[[z]]), as.character(U[[z]]))) {
      stop("Identity mismatch: ", z)
    }
  }

  ei <- .panel_repl_panel_info(E, "EAS")
  ui <- .panel_repl_panel_info(U, "EUR")
  common_analyzed <- (ei$analyzed %in% TRUE) & (ui$analyzed %in% TRUE)
  N_common <- sum(common_analyzed)
  alpha_common <- 0.05 / N_common

  primary_row <- function(panel, info) {
    n_sz <- sum(!is.na(info$signedZ_bonf))
    data.table(
      Panel = panel,
      N_analyzed = info$M,
      N_multi_cohort_meta_k_ge_2 = sum((info$analyzed %in% TRUE) & info$k >= 2L, na.rm = TRUE),
      Direction_concordance = .panel_repl_fmt_npct(sum(info$direction %in% TRUE, na.rm = TRUE), info$M),
      Nominal_replication = .panel_repl_fmt_npct(sum(info$rep_nom %in% TRUE, na.rm = TRUE), info$M),
      Bonferroni_replication = .panel_repl_fmt_npct(sum(info$rep_bonf %in% TRUE, na.rm = TRUE), info$M),
      SignedZ_Bonferroni_sensitivity = .panel_repl_fmt_npct(
        sum(info$signedZ_bonf %in% TRUE, na.rm = TRUE), n_sz
      )
    )
  }

  main_primary <- rbindlist(list(
    primary_row("EAS", ei),
    primary_row("EUR", ui)
  ))

  coverage_row <- function(panel, dt, info) {
    elig_col <- paste0(tolower(panel), "_panel_meta_eligible")
    .panel_repl_check_cols(dt, elig_col, panel)
    elig <- .panel_repl_as_flag(dt[[elig_col]])
    excluded <- elig %in% FALSE
    k0 <- info$k == 0L
    analyzed <- info$analyzed %in% TRUE
    not_available <- (elig %in% TRUE) & k0
    allele_blocked <- (elig %in% TRUE) & !k0 & !analyzed
    data.table(
      Panel = panel,
      N_leads = nrow(dt),
      N_phenotype_excluded_leads = sum(excluded, na.rm = TRUE),
      N_phenotypes_excluded = uniqueN(dt$koges_pheno_code[excluded]),
      N_not_available = sum(not_available, na.rm = TRUE),
      N_allele_anchor_blocked = sum(allele_blocked, na.rm = TRUE),
      N_analyzed = info$M,
      N_single_cohort_k1 = sum(analyzed & info$k == 1L, na.rm = TRUE),
      N_multi_cohort_k_ge_2 = sum(analyzed & info$k >= 2L, na.rm = TRUE)
    )
  }

  main_coverage <- rbindlist(list(
    coverage_row("EAS", E, ei),
    coverage_row("EUR", U, ui)
  ))

  common_row <- function(panel, dt, info) {
    pp <- as.numeric(dt[[paste0(panel, "_panel_P")]])
    zp <- as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_P")]])
    kb <- as.numeric(dt[[paste0(panel, "_panel_koges_beta_aligned_a1")]])
    z <- as.numeric(dt[[paste0(panel, "_sensitivity_signedZ_Z")]])
    d <- info$direction
    zdir <- is.finite(z) & is.finite(kb) & z != 0 & kb != 0 & sign(z) == sign(kb)
    idx <- common_analyzed
    n_sz <- sum(idx & is.finite(zp) & is.finite(z) & is.finite(kb) & z != 0 & kb != 0, na.rm = TRUE)
    data.table(
      Panel = panel,
      N_common_testable = N_common,
      Direction_concordance = .panel_repl_fmt_npct(sum(idx & (d %in% TRUE), na.rm = TRUE), N_common),
      Nominal_replication = .panel_repl_fmt_npct(
        sum(idx & (d %in% TRUE) & is.finite(pp) & pp <= 0.05, na.rm = TRUE), N_common
      ),
      Bonferroni_replication_shared_alpha = .panel_repl_fmt_npct(
        sum(idx & (d %in% TRUE) & is.finite(pp) & pp <= alpha_common, na.rm = TRUE), N_common
      ),
      SignedZ_Bonferroni_shared_alpha = .panel_repl_fmt_npct(
        sum(idx & zdir & is.finite(zp) & zp <= alpha_common, na.rm = TRUE), n_sz
      )
    )
  }

  main_common <- rbindlist(list(
    common_row("EAS", E, ei),
    common_row("EUR", U, ui)
  ))

  eas_bonf_primary <- ei$rep_bonf %in% TRUE
  eur_bonf_primary <- ui$rep_bonf %in% TRUE
  eas_p <- as.numeric(E$EAS_panel_P)
  eur_p <- as.numeric(U$EUR_panel_P)
  eas_bonf_shared <- common_analyzed & (ei$direction %in% TRUE) & is.finite(eas_p) & eas_p <= alpha_common
  eur_bonf_shared <- common_analyzed & (ui$direction %in% TRUE) & is.finite(eur_p) & eur_p <= alpha_common

  cross_block <- function(label, eas_ok, eur_ok, denom) {
    data.table(
      Definition = label,
      Denominator = denom,
      N_both = sum(eas_ok & eur_ok, na.rm = TRUE),
      N_EAS_only = sum(eas_ok & !eur_ok, na.rm = TRUE),
      N_EUR_only = sum(!eas_ok & eur_ok, na.rm = TRUE),
      N_neither = denom - sum(eas_ok | eur_ok, na.rm = TRUE),
      Both = .panel_repl_fmt_npct(sum(eas_ok & eur_ok, na.rm = TRUE), denom),
      EAS_only = .panel_repl_fmt_npct(sum(eas_ok & !eur_ok, na.rm = TRUE), denom),
      EUR_only = .panel_repl_fmt_npct(sum(!eas_ok & eur_ok, na.rm = TRUE), denom)
    )
  }

  main_cross <- rbindlist(list(
    cross_block("Primary panel-specific Bonferroni", eas_bonf_primary, eur_bonf_primary, nrow(E)),
    cross_block("Shared-alpha among common-testable", eas_bonf_shared, eur_bonf_shared, N_common)
  ))

  het_row <- function(panel, dt, info) {
    meta <- (info$analyzed %in% TRUE) & info$k >= 2L
    n_meta <- sum(meta, na.rm = TRUE)
    i2 <- as.numeric(dt[[paste0(panel, "_panel_I2")]])
    qp <- as.numeric(dt[[paste0(panel, "_panel_Q_P")]])
    pat <- dt[[paste0(panel, "_source_direction_pattern")]]
    analyzed <- info$analyzed %in% TRUE
    data.table(
      Panel = panel,
      N_multi_cohort_meta = n_meta,
      I2_ge_50 = .panel_repl_fmt_npct(sum(meta & is.finite(i2) & i2 >= 50, na.rm = TRUE), n_meta),
      I2_ge_75 = .panel_repl_fmt_npct(sum(meta & is.finite(i2) & i2 >= 75, na.rm = TRUE), n_meta),
      Q_P_lt_0.05 = .panel_repl_fmt_npct(sum(meta & is.finite(qp) & qp < 0.05, na.rm = TRUE), n_meta),
      Q_P_lt_0.001 = .panel_repl_fmt_npct(sum(meta & is.finite(qp) & qp < 0.001, na.rm = TRUE), n_meta),
      Mixed_source_direction_among_analyzed = .panel_repl_fmt_npct(
        sum(analyzed & pat == "Mixed", na.rm = TRUE), info$M
      )
    )
  }

  main_het <- rbindlist(list(
    het_row("EAS", E, ei),
    het_row("EUR", U, ui)
  ))

  E_ext <- .panel_repl_make_extended(E, "EAS", ei, ui$analyzed)
  U_ext <- .panel_repl_make_extended(U, "EUR", ui, ei$analyzed)

  list(
    eas = E_ext,
    eur = U_ext,
    main_primary = main_primary,
    main_coverage = main_coverage,
    main_common = main_common,
    main_cross = main_cross,
    main_het = main_het,
    N_common = N_common,
    alpha_common = alpha_common,
    ei = ei,
    ui = ui
  )
}
