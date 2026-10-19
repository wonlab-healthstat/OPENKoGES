#!/usr/bin/env Rscript

# SuppleTable7: per-phenotype lead with the smallest KoGES P value
# (categorical traits from ExtendedData2_categorical.xlsx).
#
# Input:
#   RESULT_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx
#     - Total_summary_Novelty_Repl  (locus annotation, min-P, novelty)
#     - Ancestry_summary_EAS / EUR  (variant-level replication ≥1 cohort)
#
# Output:
#   RESULT_manuscript_final/SuppleTable7/SuppleTable7.csv
#
# Novel leads (Integrated novelty status == "Novel") are marked with "*"
# after the phenotype display name.
# EAS cohorts here are BBJ + TPMI only (no TWB in the categorical workbook).

###############################################################################
# 0. Settings
###############################################################################

INPUT_XLSX <- "/home2/sjim/openKoGES/RESULT_manuscript_final/ExtendedData/ExtendedData2_categorical.xlsx"
OUT_CSV <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable7/SuppleTable7.csv"

# Source Phenotype name → display name (and row order).
PHENOTYPE_DISPLAY <- c(
  "Allergic disease" = "Allergic disease",
  "Breast cancer" = "Breast cancer",
  "Cataract" = "Cataract",
  "Cerebrovascular disease" = "Cerebrovascular disease",
  "Chronic gastritis" = "Chronic gastritis",
  "Chronic gastritis / Peptic ulcer" = "Chronic gastritis / Peptic ulcer",
  "Chronic kidney disease (CKD)" = "Chronic kidney disease (CKD)",
  "Chronic lung disease" = "Chronic lung disease",
  "Coronary artery disease" = "Coronary artery disease",
  "Diabetes mellitus" = "Diabetes mellitus",
  "Fatty liver" = "Fatty liver",
  "Fracture" = "Fracture",
  "Gallbladder disease" = "Gallbladder disease",
  "Gastric cancer" = "Gastric cancer",
  "Gout" = "Gout",
  "Hyperlipidemia" = "Hyperlipidemia",
  "Hypertension" = "Hypertension",
  "Neurosis" = "Neurosis",
  "Osteoporosis" = "Osteoporosis",
  "Peptic ulcer" = "Peptic ulcer",
  "Periodontal Disease" = "Periodontal Disease",
  "Thyroid cancer" = "Thyroid cancer",
  "Thyroid disease" = "Thyroid disease",
  "Uterine cancer" = "Uterine cancer"
)

###############################################################################
# 1. Packages and helpers
###############################################################################

required_packages <- c("readxl", "dplyr")
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
})

as_numeric <- function(x) suppressWarnings(as.numeric(x))

fmt_p <- function(p) {
  p <- as_numeric(p)
  vapply(p, function(x) {
    if (!is.finite(x)) return("—")
    if (x <= 0 || x < 1e-300) return("<1.00E-300")
    sprintf("%.2E", x)
  }, character(1), USE.NAMES = FALSE)
}

fmt_dash <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == "" | toupper(trimws(x)) %in% c("NA", "NAN")] <- "—"
  x
}

pick_rep_col <- function(nms, ancestry_label) {
  candidates <- c(
    "Nominal replication in ≥1 cohort (직접 매칭)",
    "Nominal replication in ≥1 cohort"
  )
  hit <- intersect(candidates, nms)
  if (length(hit) == 0) {
    stop(
      "Missing variant-level replication column in ", ancestry_label, ". ",
      "Expected one of: ", paste(candidates, collapse = " / "),
      call. = FALSE
    )
  }
  hit[[1]]
}

###############################################################################
# 2. Read ExtendedData2 sheets
###############################################################################

if (!file.exists(INPUT_XLSX)) {
  stop("Input workbook not found: ", INPUT_XLSX, call. = FALSE)
}

total <- readxl::read_excel(INPUT_XLSX, sheet = "Total_summary_Novelty_Repl")
eas <- readxl::read_excel(INPUT_XLSX, sheet = "Ancestry_summary_EAS")
eur <- readxl::read_excel(INPUT_XLSX, sheet = "Ancestry_summary_EUR")

required_total <- c(
  "Phenotype", "Chromosome", "Base pair position (GRCh37)", "rsID",
  "Nearest gene", "KoGES P value",
  "Minimum P value within ±1 Mb in BBJ",
  "Minimum P value within ±1 Mb in TPMI",
  "Minimum P value within ±1 Mb in UKB",
  "Minimum P value within ±1 Mb in FinnGen",
  "Minimum P value within ±1 Mb in MVP",
  "Integrated novelty status"
)
absent_total <- setdiff(required_total, names(total))
if (length(absent_total) > 0) {
  stop(
    "Missing columns in Total_summary_Novelty_Repl: ",
    paste(absent_total, collapse = ", "),
    call. = FALSE
  )
}

eas_rep_col <- pick_rep_col(names(eas), "Ancestry_summary_EAS")
eur_rep_col <- pick_rep_col(names(eur), "Ancestry_summary_EUR")
message("EAS replication column: ", eas_rep_col)
message("EUR replication column: ", eur_rep_col)

###############################################################################
# 3. One lead per phenotype (minimum KoGES P)
###############################################################################

total_prep <- total %>%
  mutate(
    koges_p_num = as_numeric(.data[["KoGES P value"]]),
    chr_num = as_numeric(.data[["Chromosome"]]),
    bp_num = as_numeric(.data[["Base pair position (GRCh37)"]])
  ) %>%
  filter(is.finite(koges_p_num))

known_phenos <- intersect(names(PHENOTYPE_DISPLAY), unique(total_prep$Phenotype))
extra_phenos <- setdiff(unique(total_prep$Phenotype), names(PHENOTYPE_DISPLAY))
if (length(extra_phenos) > 0) {
  message(
    "Phenotypes present in ExtendedData2 but not in PHENOTYPE_DISPLAY ",
    "(appended at end): ", paste(extra_phenos, collapse = ", ")
  )
}
pheno_levels <- c(known_phenos, sort(extra_phenos))

best_leads <- total_prep %>%
  group_by(Phenotype) %>%
  slice_min(order_by = koges_p_num, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Phenotype = factor(Phenotype, levels = pheno_levels)) %>%
  arrange(Phenotype)

eas_join <- eas %>%
  transmute(
    Phenotype = as.character(.data[["Phenotype"]]),
    rsID = as.character(.data[["rsID"]]),
    eas_variant_level_replication = fmt_dash(.data[[eas_rep_col]])
  )

eur_join <- eur %>%
  transmute(
    Phenotype = as.character(.data[["Phenotype"]]),
    rsID = as.character(.data[["rsID"]]),
    eur_variant_level_replication = fmt_dash(.data[[eur_rep_col]])
  )

supple_table5 <- best_leads %>%
  mutate(
    Phenotype_src = as.character(Phenotype),
    display_base = ifelse(
      Phenotype_src %in% names(PHENOTYPE_DISPLAY),
      unname(PHENOTYPE_DISPLAY[Phenotype_src]),
      Phenotype_src
    ),
    is_novel = .data[["Integrated novelty status"]] %in% "Novel",
    Phenotype_out = ifelse(is_novel, paste0(display_base, "*"), display_base),
    Chromosomal_position = paste0(as.integer(chr_num), ":", as.integer(bp_num))
  ) %>%
  left_join(eas_join, by = c("Phenotype_src" = "Phenotype", "rsID" = "rsID")) %>%
  left_join(eur_join, by = c("Phenotype_src" = "Phenotype", "rsID" = "rsID")) %>%
  transmute(
    Phenotype = Phenotype_out,
    `Chromosomal position` = Chromosomal_position,
    rsID = as.character(rsID),
    `Nearest gene` = as.character(.data[["Nearest gene"]]),
    `P-value in KoGES` = fmt_p(koges_p_num),
    `Minimum P-value within ±1 Mb in BBJ` =
      fmt_p(.data[["Minimum P value within ±1 Mb in BBJ"]]),
    `Minimum P-value within ±1 Mb in TPMI` =
      fmt_p(.data[["Minimum P value within ±1 Mb in TPMI"]]),
    `EAS variant-level replication` = coalesce(eas_variant_level_replication, "—"),
    `Minimum P-value within ±1 Mb in UKB` =
      fmt_p(.data[["Minimum P value within ±1 Mb in UKB"]]),
    `Minimum P-value within ±1 Mb in FinnGen` =
      fmt_p(.data[["Minimum P value within ±1 Mb in FinnGen"]]),
    `Minimum P-value within ±1 Mb in MVP` =
      fmt_p(.data[["Minimum P value within ±1 Mb in MVP"]]),
    `EUR variant-level replication` = coalesce(eur_variant_level_replication, "—")
  )

n_novel <- sum(best_leads[["Integrated novelty status"]] %in% "Novel", na.rm = TRUE)
message(
  "SuppleTable5 rows: ", nrow(supple_table5),
  " phenotypes; novel (marked *): ", n_novel
)

###############################################################################
# 4. Write output
###############################################################################

dir.create(dirname(OUT_CSV), recursive = TRUE, showWarnings = FALSE)
utils::write.csv(supple_table5, OUT_CSV, row.names = FALSE, fileEncoding = "UTF-8")
message("Wrote: ", normalizePath(OUT_CSV, winslash = "/"))
