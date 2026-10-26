#!/usr/bin/env Rscript

# Threshold-based winner's-curse correction and external replication calibration
# for all OpenKoGES phenotype-lead associations.
#
# Input:  DATA_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx
# Outputs:
#   DATA_manuscript_final/ExtendedData/ExtendedData3.csv
#     (winner's-curse analysis dump for continuous leads; NOT Extended Data 3)
#   RESULT_manuscript_final/SuppleFigure9/SuppleFigure9.png
#   RESULT_manuscript_final/SuppleTable10/SuppleTable10.csv
#
# Publication sheets (run after this script):
#   python3 CODE_manuscript_final/ExtendedData_export/export_ExtendedData_publication.py
#   → RESULT_manuscript_final/ExtendedData/Extended Data6A.xlsx  (from ExtendedData3.csv)
#   → RESULT_manuscript_final/ExtendedData/Extended Data6B.xlsx
# Note: manuscript Extended Data 3 is the 172 COJO table, not this CSV.
#
# WC correction conditions on the GWS threshold applied to P_BOLT_LMM_INF
# (not a BETA/SE Wald Z):
#   Z_BOLT = sign(beta) * Phi^{-1}(1 - P_BOLT_LMM_INF / 2)
#   beta_WC = beta_BOLT * |Z_WC| / |Z_BOLT|
# SuppleTable10 uses Bonferroni replication (source-P classification) with
# cohort-specific expected power from raw and WC-corrected KoGES betas.

###############################################################################
# 0. Settings
###############################################################################

DATA_DIR <- "/home2/sjim/openKoGES/DATA_manuscript_final/ExtendedData"
INPUT_XLSX <- file.path(DATA_DIR, "ExtendedData2_continuous.xlsx")
RESULT_MANUSCRIPT_DIR <- "/home2/sjim/openKoGES/RESULT_manuscript_final"

DISCOVERY_SHEET <- "Total_summary_Novelty_Repl"
DISCOVERY_SELECTION_P <- 5e-8
SELECTION_IS_INCLUSIVE <- TRUE
COHORTS <- c("BBJ", "TWB", "TPMI", "UKB", "FinnGen", "MVP")

###############################################################################
# 1. Packages and helpers
###############################################################################

required_packages <- c("readxl", "dplyr", "ggplot2", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                               logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(
    "Install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

assert_columns <- function(data, required, sheet_name) {
  absent <- setdiff(required, names(data))
  if (length(absent) > 0) {
    stop(
      "Missing required column(s) in ", sheet_name, ": ",
      paste(absent, collapse = ", "),
      call. = FALSE
    )
  }
}

as_numeric <- function(x) suppressWarnings(as.numeric(x))

is_yes <- function(x) {
  toupper(trimws(as.character(x))) %in% c("YES", "TRUE", "1", "REPLICATED")
}

make_key <- function(phenotype, variant_id) {
  paste(as.character(phenotype), as.character(variant_id), sep = "||")
}

fmt_int <- function(x) {
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_1 <- function(x) {
  format(round(as.numeric(x), 1), nsmall = 1, big.mark = ",", scientific = FALSE, trim = TRUE)
}

fmt_n_pct <- function(n, rate, digits_n = 1L) {
  n_txt <- if (identical(digits_n, 0L)) fmt_int(n) else fmt_1(n)
  paste0(n_txt, " (", sprintf("%.1f%%", 100 * as.numeric(rate)), ")")
}

# Conditional-likelihood WC MLE on BOLT selection Z.
wc_correct_one <- function(beta_hat, p_bolt_lmm_inf, z_cutoff) {
  if (!is.finite(beta_hat) || beta_hat == 0 ||
      !is.finite(p_bolt_lmm_inf) || p_bolt_lmm_inf < 0 || p_bolt_lmm_inf > 1) {
    return(c(
      wc_beta = NA_real_, wc_z = NA_real_, raw_z = NA_real_,
      shrinkage_abs_ratio = NA_real_, boundary_solution = NA
    ))
  }
  p_bolt_lmm_inf <- max(p_bolt_lmm_inf, .Machine$double.xmin)

  z_bolt_signed <- sign(beta_hat) *
    stats::qnorm(p_bolt_lmm_inf / 2, lower.tail = FALSE)
  raw_z <- abs(z_bolt_signed)
  beta_sign <- ifelse(z_bolt_signed >= 0, 1, -1)

  conditional_nll <- function(theta) {
    p_selected <- pnorm(-z_cutoff - theta) + pnorm(theta - z_cutoff)
    p_selected <- max(p_selected, .Machine$double.xmin)
    -(dnorm(raw_z, mean = theta, sd = 1, log = TRUE) - log(p_selected))
  }

  upper <- max(raw_z + 12, z_cutoff + 12)
  fit <- stats::optimize(conditional_nll, interval = c(0, upper))
  theta_wc <- fit$minimum
  shrink_factor <- ifelse(raw_z > 0, theta_wc / raw_z, NA_real_)

  c(
    wc_beta = beta_hat * shrink_factor,
    wc_z = beta_sign * theta_wc,
    raw_z = z_bolt_signed,
    shrinkage_abs_ratio = shrink_factor,
    boundary_solution = theta_wc < 1e-7
  )
}

wc_correct_vector <- function(beta, p_bolt, z_cutoff) {
  n <- length(beta)
  ans <- vector("list", n)
  for (i in seq_len(n)) {
    if (i %% 500 == 0 || i == n) {
      message("Winner's-curse correction: ", i, " / ", n)
    }
    ans[[i]] <- wc_correct_one(beta[i], p_bolt[i], z_cutoff)
  }
  as.data.frame(do.call(rbind, ans), stringsAsFactors = FALSE)
}

# Direction-concordant two-sided replication power at alpha.
directional_replication_power <- function(beta_true, external_se, alpha_two_sided) {
  mu <- abs(beta_true) / external_se
  z_critical <- stats::qnorm(1 - alpha_two_sided / 2)
  stats::pnorm(mu - z_critical)
}

###############################################################################
# 2. Discovery leads → ExtendedData3.csv
###############################################################################

if (!file.exists(INPUT_XLSX)) {
  stop("Input workbook not found: ", INPUT_XLSX, call. = FALSE)
}

discovery_raw <- readxl::read_excel(INPUT_XLSX, sheet = DISCOVERY_SHEET)
assert_columns(
  discovery_raw,
  c(
    "lead_row_id", "Phenotype", "Variant ID", "rsID", "KoGES beta",
    "KoGES standard error", "KoGES P value", "KoGES N",
    "Integrated novelty status",
    "Strict novelty across EAS, EUR, and GWAS Catalog",
    "Novelty evidence pattern"
  ),
  DISCOVERY_SHEET
)

discovery <- discovery_raw %>%
  transmute(
    lead_row_id = as.character(.data[["lead_row_id"]]),
    phenotype = as.character(.data[["Phenotype"]]),
    variant_id = as.character(.data[["Variant ID"]]),
    rsid = as.character(.data[["rsID"]]),
    koges_beta = as_numeric(.data[["KoGES beta"]]),
    koges_se = as_numeric(.data[["KoGES standard error"]]),
    koges_p = as_numeric(.data[["KoGES P value"]]),
    koges_n = as_numeric(.data[["KoGES N"]]),
    integrated_novelty_status = as.character(.data[["Integrated novelty status"]]),
    strict_novelty_status = as.character(
      .data[["Strict novelty across EAS, EUR, and GWAS Catalog"]]
    ),
    novelty_evidence_pattern = as.character(.data[["Novelty evidence pattern"]])
  ) %>%
  mutate(
    join_key = make_key(phenotype, variant_id),
    koges_p_underflow_to_zero = is.finite(koges_p) & koges_p == 0,
    koges_p_raw = koges_p,
    koges_p = ifelse(koges_p_underflow_to_zero, .Machine$double.xmin, koges_p),
    valid_koges = is.finite(koges_beta) & koges_beta != 0 &
      is.finite(koges_p) & koges_p > 0 & koges_p <= 1 &
      is.finite(koges_se) & koges_se > 0,
    selected_by_p = if (SELECTION_IS_INCLUSIVE) {
      koges_p <= DISCOVERY_SELECTION_P
    } else {
      koges_p < DISCOVERY_SELECTION_P
    }
  )

if (anyDuplicated(discovery$join_key) > 0) {
  stop("Discovery join keys are not unique (Phenotype + Variant ID).", call. = FALSE)
}
if (!all(discovery$valid_koges)) {
  stop(
    "Invalid KoGES beta / P_BOLT_LMM_INF / SE in ",
    sum(!discovery$valid_koges), " discovery row(s).",
    call. = FALSE
  )
}

n_p_underflow <- sum(discovery$koges_p_underflow_to_zero, na.rm = TRUE)
if (n_p_underflow > 0) {
  message(
    "Floored ", n_p_underflow,
    " underflowed P_BOLT_LMM_INF values (0) to .Machine$double.xmin."
  )
}

z_selection <- stats::qnorm(1 - DISCOVERY_SELECTION_P / 2)
wc_matrix <- wc_correct_vector(discovery$koges_beta, discovery$koges_p, z_selection)

discovery_corrected <- bind_cols(discovery, tibble::as_tibble(wc_matrix)) %>%
  mutate(
    boundary_p_flag = abs(koges_p - DISCOVERY_SELECTION_P) < 1e-15,
    p_below_or_equal_selection_threshold = selected_by_p,
    wald_z = koges_beta / koges_se,
    bolt_vs_wald_z_diff = raw_z - wald_z
  )

extended_data3_csv <- file.path(DATA_DIR, "ExtendedData3.csv")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(discovery_corrected, extended_data3_csv, row.names = FALSE)
message("Wrote: ", normalizePath(extended_data3_csv, winslash = "/"))
message(
  "Discovery rows: ", nrow(discovery_corrected),
  "; selected: ", sum(discovery_corrected$selected_by_p, na.rm = TRUE),
  "; exactly at threshold: ", sum(discovery_corrected$boundary_p_flag, na.rm = TRUE)
)

###############################################################################
# 3. External cohorts → Bonferroni expected/observed (SuppleTable10)
###############################################################################

cohort_required_columns <- c(
  "KoGES Phenotype", "KoGES Variant ID", "Exact-Variant Replication Testable",
  "External Standard Error",
  "Bonferroni-Significant Replication", "Primary Bonferroni Threshold"
)

read_and_score_cohort <- function(cohort) {
  sheet <- paste0("Cohort_detail_", cohort)
  ext_raw <- readxl::read_excel(INPUT_XLSX, sheet = sheet)
  assert_columns(ext_raw, cohort_required_columns, sheet)

  ext <- ext_raw %>%
    transmute(
      cohort = cohort,
      phenotype = as.character(.data[["KoGES Phenotype"]]),
      variant_id = as.character(.data[["KoGES Variant ID"]]),
      testable = is_yes(.data[["Exact-Variant Replication Testable"]]),
      external_se = as_numeric(.data[["External Standard Error"]]),
      observed_bonferroni = is_yes(.data[["Bonferroni-Significant Replication"]]),
      bonferroni_alpha = as_numeric(.data[["Primary Bonferroni Threshold"]])
    ) %>%
    mutate(join_key = make_key(phenotype, variant_id)) %>%
    left_join(
      discovery_corrected %>%
        select(join_key, lead_row_id, koges_beta, wc_beta),
      by = "join_key"
    )

  if (any(is.na(ext$lead_row_id))) {
    stop(
      "Could not match ", sum(is.na(ext$lead_row_id)), " row(s) from ", sheet,
      " to discovery using Phenotype + Variant ID.",
      call. = FALSE
    )
  }

  ext %>%
    mutate(
      valid_for_power = testable &
        is.finite(wc_beta) & is.finite(koges_beta) &
        is.finite(external_se) & external_se > 0 &
        is.finite(bonferroni_alpha) & bonferroni_alpha > 0 & bonferroni_alpha < 1,
      expected_bonferroni_raw = ifelse(
        valid_for_power,
        directional_replication_power(koges_beta, external_se, bonferroni_alpha),
        NA_real_
      ),
      expected_bonferroni_wc = ifelse(
        valid_for_power,
        directional_replication_power(wc_beta, external_se, bonferroni_alpha),
        NA_real_
      )
    )
}

external_lead_level <- bind_rows(lapply(COHORTS, read_and_score_cohort))

supple_table9_core <- external_lead_level %>%
  filter(
    valid_for_power,
    is.finite(expected_bonferroni_raw),
    is.finite(expected_bonferroni_wc),
    !is.na(observed_bonferroni)
  ) %>%
  group_by(cohort) %>%
  summarise(
    n_testable = dplyr::n(),
    observed = sum(as.integer(observed_bonferroni)),
    expected_raw = sum(pmin(pmax(expected_bonferroni_raw, 0), 1)),
    expected_wc = sum(pmin(pmax(expected_bonferroni_wc, 0), 1)),
    .groups = "drop"
  ) %>%
  mutate(
    cohort = factor(cohort, levels = COHORTS),
    observed_rate = observed / n_testable,
    expected_raw_rate = expected_raw / n_testable,
    expected_wc_rate = expected_wc / n_testable,
    diff_pp = 100 * (observed_rate - expected_wc_rate),
    ratio = observed / expected_wc
  ) %>%
  arrange(cohort)

total_n <- sum(supple_table9_core$n_testable)
total_obs <- sum(supple_table9_core$observed)
total_raw <- sum(supple_table9_core$expected_raw)
total_wc <- sum(supple_table9_core$expected_wc)

supple_table9 <- bind_rows(
  mutate(supple_table9_core, cohort = as.character(cohort)),
  tibble::tibble(
    cohort = "All cohorts combined*",
    n_testable = total_n,
    observed = total_obs,
    expected_raw = total_raw,
    expected_wc = total_wc,
    observed_rate = total_obs / total_n,
    expected_raw_rate = total_raw / total_n,
    expected_wc_rate = total_wc / total_n,
    diff_pp = 100 * ((total_obs / total_n) - (total_wc / total_n)),
    ratio = total_obs / total_wc
  )
) %>%
  transmute(
    Cohort = cohort,
    `Testable, n` = fmt_int(n_testable),
    `Raw expected, n (%)` = fmt_n_pct(expected_raw, expected_raw_rate, digits_n = 1L),
    `WC expected, n (%)` = fmt_n_pct(expected_wc, expected_wc_rate, digits_n = 1L),
    `Observed, n (%)` = fmt_n_pct(observed, observed_rate, digits_n = 0L),
    `Observed − WC expected, %p` = fmt_1(diff_pp),
    `Observed / WC expected` = sprintf("%.3f", ratio)
  )

supple_table10_csv <- file.path(RESULT_MANUSCRIPT_DIR, "SuppleTable10", "SuppleTable10.csv")
dir.create(dirname(supple_table10_csv), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(supple_table10, supple_table10_csv, row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote: ", normalizePath(supple_table10_csv, winslash = "/"))

###############################################################################
# 4. Shrinkage attenuation figure → SuppleFigure9.png
###############################################################################

shrink_df <- discovery_corrected %>%
  filter(is.finite(shrinkage_abs_ratio), is.finite(raw_z)) %>%
  mutate(
    shrink = as.numeric(shrinkage_abs_ratio),
    shrink_plot = pmin(pmax(shrink, 0), 1),
    abs_z = abs(as.numeric(raw_z))
  )

n_leads <- nrow(shrink_df)
med_shrink <- stats::median(shrink_df$shrink, na.rm = TRUE)
iqr_shrink <- stats::quantile(shrink_df$shrink, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
n_lt_025 <- sum(shrink_df$shrink < 0.25, na.rm = TRUE)
pct_lt_025 <- 100 * n_lt_025 / n_leads

annot_a <- sprintf(
  "Median = %.3f\nIQR = %.3f–%.3f\n<0.25: %s / %s (%.1f%%)",
  med_shrink,
  iqr_shrink[1], iqr_shrink[2],
  format(n_lt_025, big.mark = ",", scientific = FALSE),
  format(n_leads, big.mark = ",", scientific = FALSE),
  pct_lt_025
)

panel_a <- ggplot(shrink_df, aes(x = shrink_plot)) +
  annotate(
    "rect",
    xmin = iqr_shrink[1], xmax = iqr_shrink[2],
    ymin = -Inf, ymax = Inf,
    fill = "#B8D4E8", alpha = 0.45
  ) +
  geom_histogram(
    breaks = seq(0, 1, by = 0.025),
    fill = "#1F4E79",
    colour = "#1F4E79",
    closed = "right"
  ) +
  geom_vline(xintercept = med_shrink, colour = "#C00000", linetype = "dashed", linewidth = 0.7) +
  annotate(
    "label",
    x = 0.02, y = Inf,
    label = annot_a,
    hjust = 0, vjust = 1.15,
    size = 3.2,
    label.size = 0.3,
    fill = "white",
    colour = "black",
    label.padding = grid::unit(0.35, "lines")
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2), expand = c(0.01, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.08))) +
  labs(
    title = "A  Distribution of winner's-curse attenuation",
    x = expression(paste("Shrinkage ratio, |", beta[WC] / beta[raw], "|")),
    y = "Number of lead associations"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    axis.title = element_text(size = 10),
    axis.text = element_text(colour = "black"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

z_breaks <- c(5.45, 5.6, 5.8, 6.0, 6.5, 7.0, 8.0, 10, 15, Inf)
z_labels <- c(
  "5.45–5.6", "5.6–5.8", "5.8–6.0", "6.0–6.5", "6.5–7.0",
  "7.0–8.0", "8.0–10", "10–15", "≥15"
)

shrink_by_z <- shrink_df %>%
  mutate(
    z_bin = cut(
      abs_z,
      breaks = z_breaks,
      labels = z_labels,
      right = FALSE,
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(z_bin)) %>%
  group_by(z_bin) %>%
  summarise(
    med = stats::median(shrink, na.rm = TRUE),
    q25 = stats::quantile(shrink, 0.25, na.rm = TRUE, names = FALSE),
    q75 = stats::quantile(shrink, 0.75, na.rm = TRUE, names = FALSE),
    .groups = "drop"
  )

panel_b <- ggplot(shrink_by_z, aes(x = z_bin, y = med, group = 1)) +
  geom_hline(yintercept = 1, colour = "#C00000", linetype = "dashed", linewidth = 0.65) +
  geom_errorbar(aes(ymin = q25, ymax = q75), width = 0.18, colour = "#1F4E79", linewidth = 0.55) +
  geom_line(colour = "#1F4E79", linewidth = 0.7) +
  geom_point(size = 2.4, colour = "#1F4E79") +
  annotate(
    "text",
    x = length(z_labels),
    y = 0.08,
    label = "— No attenuation",
    hjust = 1,
    vjust = 0.5,
    size = 3.1,
    colour = "#C00000"
  ) +
  scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, by = 0.2), expand = c(0, 0)) +
  labs(
    title = "B  Attenuation is strongest near the discovery threshold",
    x = expression(paste("Observed discovery |", Z[raw], "| bins")),
    y = expression(paste("Shrinkage ratio, |", beta[WC] / beta[raw], "|"))
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 11.5, hjust = 0),
    axis.title = element_text(size = 10),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

figure_wc_attenuation <- (panel_a | panel_b) +
  plot_annotation(
    title = "Threshold-based winner's-curse attenuation of KoGES discovery effects",
    caption = paste0(
      "All ", format(n_leads, big.mark = ",", scientific = FALSE),
      " phenotype–lead associations. Panel B shows the median (points) and ",
      "interquartile range (whiskers) of lead-specific shrinkage within discovery |Z| bins; ",
      "effect directions are retained."
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13, hjust = 0),
      plot.caption = element_text(size = 8.5, hjust = 0, colour = "grey20", margin = margin(t = 8))
    )
  )

figure_png <- file.path(RESULT_MANUSCRIPT_DIR, "SuppleFigure9", "SuppleFigure9.png")
dir.create(dirname(figure_png), recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = figure_png,
  plot = figure_wc_attenuation,
  width = 12.5,
  height = 5.6,
  units = "in",
  dpi = 400,
  bg = "white"
)
message("Wrote: ", normalizePath(figure_png, winslash = "/"))

message("\nCompleted.")
message("  ", normalizePath(extended_data3_csv, winslash = "/"))
message("  ", normalizePath(supple_table9_csv, winslash = "/"))
message("  ", normalizePath(figure_png, winslash = "/"))
