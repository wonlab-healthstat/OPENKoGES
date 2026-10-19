#!/usr/bin/env Rscript
################################################################################
# Build phenotype-level external N + EUR panel phenotype eligibility
#
# Source
#   KoGES_EURmatching_sjim.xlsx Continuous sheet
#   (UKB_N / FinnGen_N / MVP_N + cohort path flags via eur_koges_trait_map.R)
#
# Outputs
#   1) eur_external_phenotype_N.tsv
#        output_stem-level analytic N for UKB / FinnGen / MVP
#   2) eur_panel_phenotype_eligibility.tsv
#        koges_pheno_code-level flag for whether panel meta is allowed
#
# Notes
#   - Do NOT use KoGES N.
#   - A cohort stem is authorized only when the Excel path is present AND
#     that cohort's analytic N is finite and > 0.
#   - Phenotypes with no authorized EUR cohort remain in eligibility with
#     eur_panel_meta_eligible=FALSE and a filled meta_exclusion_reason.
################################################################################

suppressPackageStartupMessages(library(data.table))

script_dir <- tryCatch(
  dirname(sys.frame(1)$ofile),
  error = function(e) NA_character_
)
if (is.na(script_dir) || !nzchar(script_dir)) {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1L) {
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg)))
  } else {
    script_dir <- getwd()
  }
}
source("/home2/sjim/openKoGES/CODE_GWASsumstatQC/eur_koges_trait_map.R")

meta_excel <- Sys.getenv("EUR_MATCHING_EXCEL", META_EXCEL_DEFAULT)
out_n_path <- Sys.getenv(
  "EUR_EXTERNAL_N_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eur_external_phenotype_N.tsv"
)
out_elig_path <- Sys.getenv(
  "EUR_PANEL_ELIGIBILITY_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eur_panel_phenotype_eligibility.tsv"
)

meta <- read_eur_matching_meta(meta_excel, META_SHEET_DEFAULT)
need_n <- c("UKB_N", "FinnGen_N", "MVP_N")
miss_n <- setdiff(need_n, names(meta))
if (length(miss_n) > 0L) {
  stop("META_EXCEL missing N column(s): ", paste(miss_n, collapse = ", "))
}

meta[, `:=`(
  ukb_n = as.numeric(UKB_N),
  finngen_n = as.numeric(FinnGen_N),
  mvp_n = as.numeric(MVP_N)
)]

# Keep Continuous phenotypes with a KoGES BOLT stem
pheno <- meta[!is.na(VarNM) & nzchar(VarNM) & !is.na(output_stem) & nzchar(output_stem)]
if (anyDuplicated(pheno$output_stem)) {
  stop(
    "Duplicated output_stem in EUR matching Continuous sheet: ",
    paste(pheno$output_stem[duplicated(pheno$output_stem)], collapse = ", ")
  )
}
if (anyDuplicated(pheno$VarNM)) {
  stop(
    "Duplicated VarNM in EUR matching Continuous sheet: ",
    paste(pheno$VarNM[duplicated(pheno$VarNM)], collapse = ", ")
  )
}

# Authorize cohort only with path + finite analytic N
pheno[, `:=`(
  ukb_ok = has_ukb & is.finite(ukb_n) & ukb_n > 0,
  finngen_ok = has_finngen & is.finite(finngen_n) & finngen_n > 0,
  mvp_ok = has_mvp & is.finite(mvp_n) & mvp_n > 0
)]

out_n <- pheno[, .(
  output_stem,
  clinical = KoGES,
  varnm = trimws(as.character(var_nm)),
  koges_pheno_code = VarNM,
  ukb_n = fifelse(ukb_ok, ukb_n, NA_real_),
  ukb_n_source = fifelse(ukb_ok, "KoGES_EURmatching_sjim.xlsx:UKB_N", NA_character_),
  finngen_n = fifelse(finngen_ok, finngen_n, NA_real_),
  finngen_n_source = fifelse(
    finngen_ok, "KoGES_EURmatching_sjim.xlsx:FinnGen_N", NA_character_
  ),
  mvp_n = fifelse(mvp_ok, mvp_n, NA_real_),
  mvp_n_source = fifelse(mvp_ok, "KoGES_EURmatching_sjim.xlsx:MVP_N", NA_character_)
)]
setorder(out_n, output_stem)

elig <- pheno[, .(
  koges_pheno_code = VarNM,
  ukb_output_stem = fifelse(ukb_ok, output_stem, NA_character_),
  finngen_output_stem = fifelse(finngen_ok, output_stem, NA_character_),
  mvp_output_stem = fifelse(mvp_ok, output_stem, NA_character_),
  n_authorized = as.integer(ukb_ok) + as.integer(finngen_ok) + as.integer(mvp_ok)
)]
excl_reason <- "No scientifically harmonized quantitative phenotype in the EUR panel"
elig[, `:=`(
  eur_panel_meta_eligible = n_authorized >= 1L,
  meta_exclusion_reason = fifelse(
    n_authorized >= 1L,
    NA_character_,
    excl_reason
  )
)]
elig[, n_authorized := NULL]
# Keep eligible + explicitly excluded phenotypes (reason must never be blank)
elig[, meta_exclusion_reason := as.character(meta_exclusion_reason)]
setorder(elig, -eur_panel_meta_eligible, koges_pheno_code)

# Validate stems referenced by eligibility exist in N table with finite N
# (only for authorized / eligible stems)
elig_ok <- elig[eur_panel_meta_eligible == TRUE]
n_long <- rbindlist(list(
  out_n[is.finite(ukb_n) & ukb_n > 0, .(cohort = "UKB", output_stem, external_n = ukb_n)],
  out_n[is.finite(finngen_n) & finngen_n > 0,
        .(cohort = "FinnGen", output_stem, external_n = finngen_n)],
  out_n[is.finite(mvp_n) & mvp_n > 0, .(cohort = "MVP", output_stem, external_n = mvp_n)]
))
assert_stem_n <- function(stem_col, cohort_label) {
  stems <- unique(na.omit(elig_ok[[stem_col]]))
  stems <- stems[nzchar(stems)]
  miss <- setdiff(stems, n_long[cohort == cohort_label, output_stem])
  if (length(miss) > 0L) {
    stop("Eligibility ", stem_col, " missing finite ", cohort_label, " N: ",
         paste(miss, collapse = ", "))
  }
}
assert_stem_n("ukb_output_stem", "UKB")
assert_stem_n("finngen_output_stem", "FinnGen")
assert_stem_n("mvp_output_stem", "MVP")

dir.create(dirname(out_n_path), recursive = TRUE, showWarnings = FALSE)
fwrite(out_n, out_n_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(elig, out_elig_path, sep = "\t", quote = FALSE, na = "NA")

cat("Wrote:", out_n_path, "\n")
cat("Wrote:", out_elig_path, "\n")
cat(
  "N phenotypes:", nrow(out_n),
  "| eligible koges codes:", sum(elig$eur_panel_meta_eligible),
  "| excluded:", sum(!elig$eur_panel_meta_eligible),
  "| UKB stems:", sum(!is.na(elig$ukb_output_stem)),
  "| FinnGen stems:", sum(!is.na(elig$finngen_output_stem)),
  "| MVP stems:", sum(!is.na(elig$mvp_output_stem)),
  "\n"
)
if (any(!elig$eur_panel_meta_eligible)) {
  cat("Excluded phenotypes:\n")
  print(elig[eur_panel_meta_eligible == FALSE,
             .(koges_pheno_code, meta_exclusion_reason)])
}
