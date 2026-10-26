#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   Build Supplementary Table 7 manuscript summaries of OpenKoGES external
#   validation from Extended Data 2 workbooks (Cont + Cate).
#
#   Extracted from the former summary-CSV blocks in ExtendedData2_cont.R and
#   ExtendedData2_cate.R so Excel merge and table export can be run separately.
#
# Run order
#   1. ExtendedData2_cont.R / ExtendedData2_cate.R (write ExtendedData2_*.xlsx)
#   2. Rscript SuppleTable7.R
#
# Inputs
#   - DATA_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx
#   - DATA_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx
#   Required sheets (each workbook):
#     Total_summary_Novelty_Repl, Replication_EAS_meta, Replication_EUR_meta
#
# Outputs
#   - RESULT_manuscript_final/SuppleTable7/SuppleTable7_cont.csv
#   - RESULT_manuscript_final/SuppleTable7/SuppleTable7_cate.csv
#
# Statistical definitions
#   1) Primary panel replication:
#      external-only IVW meta-analysis, KoGES-concordant direction, and
#      panel-specific Bonferroni significance.
#   2) Common-testable sensitivity:
#      associations analyzable in BOTH EAS and EUR panels, using one shared
#      Bonferroni threshold = 0.05 / N_common.
#   3) Strict novelty:
#      EAS locus novelty == Novel AND EUR locus novelty == Novel AND
#      GWAS Catalog novelty == Novel.
#   4) Heterogeneity:
#      I2 among external cohorts within each ancestry panel (not KoGES-vs-ext).
#   5) IVW vs signed-Z agreement of Bonferroni replication classifications.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

merge_dir <- "/home2/sjim/openKoGES/CODE_manuscript_final/5_Novelty_Replication/forEAS_EUR_Catalog_merge"
common_r <- file.path(merge_dir, "ExtendedData2_common.R")
if (!file.exists(common_r)) stop("Common helpers not found: ", common_r)
source(common_r)

xlsx_cont <- "/home2/sjim/openKoGES/DATA_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx"
xlsx_cate <- "/home2/sjim/openKoGES/DATA_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx"
out_dir <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable7"
csv_cont <- file.path(out_dir, "SuppleTable7_cont.csv")
csv_cate <- file.path(out_dir, "SuppleTable7_cate.csv")

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
as_num <- function(x) suppressWarnings(as.numeric(x))

n_contributing <- function(x) {
  y <- trimws(as.character(x))
  out <- integer(length(y))
  ok <- !is.na(y) & nzchar(y)
  if (any(ok)) {
    out[ok] <- lengths(strsplit(y[ok], "\\+"))
  }
  out
}

pct <- function(n, d) {
  if (is.na(d) || d == 0L) return(NA_real_)
  100 * n / d
}

fmt_count <- function(n) {
  format(as.integer(n), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_fraction <- function(n, d, show_pct = TRUE) {
  if (is.na(d) || d == 0L) return(NA_character_)
  if (show_pct) {
    sprintf(
      "%s/%s (%.1f%%)",
      fmt_count(n),
      fmt_count(d),
      pct(n, d)
    )
  } else {
    sprintf("%s/%s", fmt_count(n), fmt_count(d))
  }
}

make_row <- function(section, metric,
                     eas_n = NA_integer_, eas_d = NA_integer_,
                     eur_n = NA_integer_, eur_d = NA_integer_,
                     eas_count_only = FALSE, eur_count_only = FALSE,
                     eas_show_pct = TRUE, eur_show_pct = TRUE,
                     notes = NA_character_) {
  eas_display <- if (eas_count_only) {
    fmt_count(eas_n)
  } else if (!is.na(eas_d)) {
    fmt_fraction(eas_n, eas_d, eas_show_pct)
  } else {
    NA_character_
  }
  eur_display <- if (eur_count_only) {
    fmt_count(eur_n)
  } else if (!is.na(eur_d)) {
    fmt_fraction(eur_n, eur_d, eur_show_pct)
  } else {
    NA_character_
  }
  data.table(
    Section = section,
    Metric = metric,
    EAS_n = eas_n,
    EAS_denominator = eas_d,
    EAS_percent = if (!is.na(eas_d)) pct(eas_n, eas_d) else NA_real_,
    EAS = eas_display,
    EUR_n = eur_n,
    EUR_denominator = eur_d,
    EUR_percent = if (!is.na(eur_d)) pct(eur_n, eur_d) else NA_real_,
    EUR = eur_display,
    Notes = notes
  )
}

#' Build SuppleTable7 summary rows from an Extended Data 2 workbook.
#'
#' @param input_xlsx Path to ExtendedData2_continuous.xlsx or _categorical.xlsx
#' @param label Short label for messages ("Cont" / "Cate")
#' @param expected Optional named list of expected counts (Cont only); warnings if mismatch
#' @return data.table ready to fwrite
build_supple_table7 <- function(input_xlsx, label, expected = NULL) {
  if (!file.exists(input_xlsx)) {
    stop("[", label, "] Input file not found: ", input_xlsx)
  }

  total <- as.data.table(
    read.xlsx(
      input_xlsx,
      sheet = "Total_summary_Novelty_Repl",
      check.names = FALSE,
      sep.names = " "
    )
  )
  eas <- as.data.table(
    read.xlsx(
      input_xlsx,
      sheet = "Replication_EAS_meta",
      check.names = FALSE,
      sep.names = " "
    )
  )
  eur <- as.data.table(
    read.xlsx(
      input_xlsx,
      sheet = "Replication_EUR_meta",
      check.names = FALSE,
      sep.names = " "
    )
  )

  if (!("EAS_EUR_both_testable" %in% names(eas)) &&
      ("common_testable_both_panels" %in% names(eas))) {
    setnames(eas, "common_testable_both_panels", "EAS_EUR_both_testable")
  }
  if (!("EAS_EUR_both_testable" %in% names(eur)) &&
      ("common_testable_both_panels" %in% names(eur))) {
    setnames(eur, "common_testable_both_panels", "EAS_EUR_both_testable")
  }

  total_req <- c(
    "lead_row_id",
    "EAS locus novelty",
    "EUR locus novelty",
    "GWAS Catalog novelty"
  )
  meta_req <- c(
    "lead_row_id",
    "panel_analyzed",
    "EAS_EUR_both_testable",
    "contributing_cohorts",
    "external_panel_P",
    "direction_concordant_with_KoGES",
    "replication_nominal",
    "replication_Bonferroni",
    "signedZ_replication_Bonferroni",
    "external_cohort_direction_pattern",
    "external_cohort_I2"
  )
  require_cols(total, total_req, paste0(label, " Total_summary_Novelty_Repl"))
  require_cols(eas, meta_req, paste0(label, " Replication_EAS_meta"))
  require_cols(eur, meta_req, paste0(label, " Replication_EUR_meta"))

  total[, lead_row_id := as.integer(lead_row_id)]
  eas[, lead_row_id := as.integer(lead_row_id)]
  eur[, lead_row_id := as.integer(lead_row_id)]

  if (anyDuplicated(total$lead_row_id)) {
    stop("[", label, "] Duplicate lead_row_id in total sheet.")
  }
  if (anyDuplicated(eas$lead_row_id)) {
    stop("[", label, "] Duplicate lead_row_id in EAS sheet.")
  }
  if (anyDuplicated(eur$lead_row_id)) {
    stop("[", label, "] Duplicate lead_row_id in EUR sheet.")
  }
  if (!setequal(total$lead_row_id, eas$lead_row_id) ||
      !setequal(total$lead_row_id, eur$lead_row_id)) {
    stop("[", label, "] lead_row_id sets differ across sheets.")
  }

  eas <- eas[match(total$lead_row_id, eas$lead_row_id)]
  eur <- eur[match(total$lead_row_id, eur$lead_row_id)]
  if (!identical(total$lead_row_id, eas$lead_row_id) ||
      !identical(total$lead_row_id, eur$lead_row_id)) {
    stop("[", label, "] Failed to align EAS/EUR rows by lead_row_id.")
  }

  n_total <- nrow(total)
  eas_analyzed <- as_flag(eas$panel_analyzed) %in% TRUE
  eur_analyzed <- as_flag(eur$panel_analyzed) %in% TRUE
  eas_direction <- as_flag(eas$direction_concordant_with_KoGES) %in% TRUE
  eur_direction <- as_flag(eur$direction_concordant_with_KoGES) %in% TRUE
  eas_nominal <- as_flag(eas$replication_nominal) %in% TRUE
  eur_nominal <- as_flag(eur$replication_nominal) %in% TRUE
  eas_bonf <- as_flag(eas$replication_Bonferroni) %in% TRUE
  eur_bonf <- as_flag(eur$replication_Bonferroni) %in% TRUE
  eas_signed_bonf_raw <- as_flag(eas$signedZ_replication_Bonferroni)
  eur_signed_bonf_raw <- as_flag(eur$signedZ_replication_Bonferroni)
  eas_k <- n_contributing(eas$contributing_cohorts)
  eur_k <- n_contributing(eur$contributing_cohorts)
  eas_meta_k2 <- eas_analyzed & eas_k >= 2L
  eur_meta_k2 <- eur_analyzed & eur_k >= 2L

  strict_novel <-
    trimws(as.character(total[["EAS locus novelty"]])) == "Novel" &
    trimws(as.character(total[["EUR locus novelty"]])) == "Novel" &
    trimws(as.character(total[["GWAS Catalog novelty"]])) == "Novel"
  n_strict <- sum(strict_novel, na.rm = TRUE)

  common_eas <- as_flag(eas$EAS_EUR_both_testable)
  common_eur <- as_flag(eur$EAS_EUR_both_testable)
  common_equal <-
    (is.na(common_eas) & is.na(common_eur)) |
    (!is.na(common_eas) & !is.na(common_eur) & common_eas == common_eur)
  if (!all(common_equal)) {
    stop("[", label, "] EAS_EUR_both_testable differs between EAS and EUR sheets.")
  }

  common_testable <- common_eas %in% TRUE
  n_common <- sum(common_testable)
  if (n_common == 0L) stop("[", label, "] No common-testable associations.")
  alpha_common <- 0.05 / n_common

  eas_common_bonf <-
    common_testable &
    eas_analyzed &
    eas_direction &
    is.finite(as_num(eas$external_panel_P)) &
    as_num(eas$external_panel_P) <= alpha_common
  eur_common_bonf <-
    common_testable &
    eur_analyzed &
    eur_direction &
    is.finite(as_num(eur$external_panel_P)) &
    as_num(eur$external_panel_P) <= alpha_common

  strict_common <- strict_novel & common_testable
  n_strict_common <- sum(strict_common)
  eas_strict_common_bonf <- strict_common & eas_common_bonf
  eur_strict_common_bonf <- strict_common & eur_common_bonf

  eas_pattern <- trimws(as.character(eas$external_cohort_direction_pattern))
  eur_pattern <- trimws(as.character(eur$external_cohort_direction_pattern))
  eas_all_concordant <- !is.na(eas_pattern) & eas_pattern == "All concordant"
  eur_all_concordant <- !is.na(eur_pattern) & eur_pattern == "All concordant"
  eas_mixed <- !is.na(eas_pattern) & eas_pattern == "Mixed"
  eur_mixed <- !is.na(eur_pattern) & eur_pattern == "Mixed"
  eas_i2 <- as_num(eas$external_cohort_I2)
  eur_i2 <- as_num(eur$external_cohort_I2)
  eas_high_i2 <- eas_meta_k2 & is.finite(eas_i2) & eas_i2 >= 50
  eur_high_i2 <- eur_meta_k2 & is.finite(eur_i2) & eur_i2 >= 50

  eas_ivw_raw <- as_flag(eas$replication_Bonferroni)
  eur_ivw_raw <- as_flag(eur$replication_Bonferroni)
  eas_agreement_available <-
    eas_analyzed & !is.na(eas_ivw_raw) & !is.na(eas_signed_bonf_raw)
  eur_agreement_available <-
    eur_analyzed & !is.na(eur_ivw_raw) & !is.na(eur_signed_bonf_raw)
  eas_agree <- eas_agreement_available & (eas_ivw_raw == eas_signed_bonf_raw)
  eur_agree <- eur_agreement_available & (eur_ivw_raw == eur_signed_bonf_raw)

  rows <- list(
    make_row(
      "A. Overall external replication",
      "Analyzed phenotype-lead associations",
      eas_n = sum(eas_analyzed),
      eur_n = sum(eur_analyzed),
      eas_count_only = TRUE,
      eur_count_only = TRUE,
      notes = paste0(
        "Total OpenKoGES phenotype-lead associations in Extended Data: ",
        fmt_count(n_total), "."
      )
    ),
    make_row(
      "A. Overall external replication",
      "Direction concordant with KoGES",
      eas_n = sum(eas_analyzed & eas_direction),
      eas_d = sum(eas_analyzed),
      eur_n = sum(eur_analyzed & eur_direction),
      eur_d = sum(eur_analyzed)
    ),
    make_row(
      "A. Overall external replication",
      "Nominal replication",
      eas_n = sum(eas_analyzed & eas_nominal),
      eas_d = sum(eas_analyzed),
      eur_n = sum(eur_analyzed & eur_nominal),
      eur_d = sum(eur_analyzed)
    ),
    make_row(
      "A. Overall external replication",
      "Bonferroni-corrected replication",
      eas_n = sum(eas_analyzed & eas_bonf),
      eas_d = sum(eas_analyzed),
      eur_n = sum(eur_analyzed & eur_bonf),
      eur_d = sum(eur_analyzed),
      notes = "Primary panel-level external-only IVW replication classification."
    ),
    make_row(
      "A. Overall external replication",
      paste0(
        "Common-testable sensitivity (N=",
        fmt_count(n_common),
        "; shared alpha=0.05/",
        fmt_count(n_common),
        ")"
      ),
      eas_n = sum(eas_common_bonf),
      eas_d = n_common,
      eur_n = sum(eur_common_bonf),
      eur_d = n_common,
      notes = sprintf("Shared Bonferroni threshold = %.6g.", alpha_common)
    ),
    make_row(
      "B. Replication among strictly novel associations",
      "Strictly novel associations analyzable in panel",
      eas_n = sum(strict_novel & eas_analyzed),
      eas_d = n_strict,
      eur_n = sum(strict_novel & eur_analyzed),
      eur_d = n_strict,
      eas_show_pct = FALSE,
      eur_show_pct = FALSE,
      notes = paste0(
        "Strict novelty: Novel in EAS, EUR, and GWAS Catalog; N=",
        fmt_count(n_strict), "."
      )
    ),
    make_row(
      "B. Replication among strictly novel associations",
      "Nominal replication among analyzable strict-novel associations",
      eas_n = sum(strict_novel & eas_analyzed & eas_nominal),
      eas_d = sum(strict_novel & eas_analyzed),
      eur_n = sum(strict_novel & eur_analyzed & eur_nominal),
      eur_d = sum(strict_novel & eur_analyzed)
    ),
    make_row(
      "B. Replication among strictly novel associations",
      "Bonferroni replication among analyzable strict-novel associations",
      eas_n = sum(strict_novel & eas_analyzed & eas_bonf),
      eas_d = sum(strict_novel & eas_analyzed),
      eur_n = sum(strict_novel & eur_analyzed & eur_bonf),
      eur_d = sum(strict_novel & eur_analyzed),
      notes = "Uses each panel's primary Bonferroni replication classification."
    ),
    make_row(
      "B. Replication among strictly novel associations",
      paste0(
        "Common-testable strict-novel subset (N=",
        fmt_count(n_strict_common),
        "): Bonferroni replication"
      ),
      eas_n = sum(eas_strict_common_bonf),
      eas_d = n_strict_common,
      eur_n = sum(eur_strict_common_bonf),
      eur_d = n_strict_common,
      notes = paste0(
        "Restricted to strict-novel associations testable in both panels; ",
        "uses the SAME shared threshold as the full common-testable analysis ",
        "(0.05/N_common), not 0.05/N_strict_common."
      )
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "Associations with >=2 contributing external cohorts",
      eas_n = sum(eas_meta_k2),
      eas_d = sum(eas_analyzed),
      eur_n = sum(eur_meta_k2),
      eur_d = sum(eur_analyzed)
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "All contributing cohorts directionally concordant with KoGES (k>=2)",
      eas_n = sum(eas_meta_k2 & eas_all_concordant),
      eas_d = sum(eas_meta_k2),
      eur_n = sum(eur_meta_k2 & eur_all_concordant),
      eur_d = sum(eur_meta_k2)
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "Mixed external-cohort directions (k>=2)",
      eas_n = sum(eas_meta_k2 & eas_mixed),
      eas_d = sum(eas_meta_k2),
      eur_n = sum(eur_meta_k2 & eur_mixed),
      eur_d = sum(eur_meta_k2)
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "High between-cohort heterogeneity (I2 >=50%; k>=2)",
      eas_n = sum(eas_high_i2),
      eas_d = sum(eas_meta_k2),
      eur_n = sum(eur_high_i2),
      eur_d = sum(eur_meta_k2),
      notes = paste0(
        "I2 reflects heterogeneity among external cohorts within the panel, ",
        "not KoGES-vs-external heterogeneity."
      )
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "Panel direction concordant with KoGES among I2 >=50%",
      eas_n = sum(eas_high_i2 & eas_direction),
      eas_d = sum(eas_high_i2),
      eur_n = sum(eur_high_i2 & eur_direction),
      eur_d = sum(eur_high_i2)
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "All contributing cohorts concordant with KoGES among I2 >=50%",
      eas_n = sum(eas_high_i2 & eas_all_concordant),
      eas_d = sum(eas_high_i2),
      eur_n = sum(eur_high_i2 & eur_all_concordant),
      eur_d = sum(eur_high_i2)
    ),
    make_row(
      "C. Meta-analysis diagnostics and robustness",
      "IVW vs signed-Z Bonferroni replication classification agreement",
      eas_n = sum(eas_agree),
      eas_d = sum(eas_agreement_available),
      eur_n = sum(eur_agree),
      eur_d = sum(eur_agreement_available),
      notes = paste0(
        "Near-identical classification supports robustness of the replication ",
        "conclusion to the evidence-combination method; it does not imply that ",
        "a homogeneous common effect exists when I2 is high."
      )
    )
  )

  result <- rbindlist(rows, use.names = TRUE, fill = TRUE)

  if (!is.null(expected)) {
    if (!is.null(expected$n_total) && n_total != expected$n_total) {
      warning("[", label, "] Expected ", expected$n_total, " total rows; observed ", n_total, ".")
    }
    if (!is.null(expected$n_strict) && n_strict != expected$n_strict) {
      warning("[", label, "] Expected ", expected$n_strict, " strict-novel; observed ", n_strict, ".")
    }
    if (!is.null(expected$n_common) && n_common != expected$n_common) {
      warning("[", label, "] Expected ", expected$n_common, " common-testable; observed ", n_common, ".")
    }
    if (!is.null(expected$eas_analyzed) && sum(eas_analyzed) != expected$eas_analyzed) {
      warning("[", label, "] Expected ", expected$eas_analyzed, " EAS-analyzed; observed ", sum(eas_analyzed), ".")
    }
    if (!is.null(expected$eur_analyzed) && sum(eur_analyzed) != expected$eur_analyzed) {
      warning("[", label, "] Expected ", expected$eur_analyzed, " EUR-analyzed; observed ", sum(eur_analyzed), ".")
    }
  }

  message(
    "[", label, "] N_total=", n_total,
    "; N_strict=", n_strict,
    "; N_common=", n_common,
    "; EAS_analyzed=", sum(eas_analyzed),
    "; EUR_analyzed=", sum(eur_analyzed)
  )
  message(
    "[", label, "] Primary Bonferroni: EAS ",
    sum(eas_analyzed & eas_bonf), "/", sum(eas_analyzed),
    "; EUR ", sum(eur_analyzed & eur_bonf), "/", sum(eur_analyzed)
  )
  message(
    "[", label, "] Strict-novel Bonferroni: EAS ",
    sum(strict_novel & eas_analyzed & eas_bonf), "/",
    sum(strict_novel & eas_analyzed),
    "; EUR ",
    sum(strict_novel & eur_analyzed & eur_bonf), "/",
    sum(strict_novel & eur_analyzed)
  )
  message(
    "[", label, "] IVW vs signed-Z agreement: EAS ",
    sum(eas_agree), "/", sum(eas_agreement_available),
    "; EUR ",
    sum(eur_agree), "/", sum(eur_agreement_available)
  )
  message(
    "[", label, "] Shared common-testable Bonferroni alpha = ",
    format(alpha_common, digits = 8)
  )

  result
}

# ------------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cont_expected <- list(
  n_total = 8005L,
  n_strict = 77L,
  n_common = 4611L,
  eas_analyzed = 6170L,
  eur_analyzed = 5028L
)

result_cont <- build_supple_table7(xlsx_cont, "Cont", expected = cont_expected)
fwrite(result_cont, file = csv_cont, sep = ",", quote = TRUE, na = "")
message("Wrote: ", normalizePath(csv_cont))

result_cate <- build_supple_table7(xlsx_cate, "Cate", expected = NULL)
fwrite(result_cate, file = csv_cate, sep = ",", quote = TRUE, na = "")
message("Wrote: ", normalizePath(csv_cate))

message("SuppleTable7 done.")
