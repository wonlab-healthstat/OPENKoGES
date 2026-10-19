#!/usr/bin/env Rscript

# ==============================================================================
# Step3d: compare primary 1000G-clump-only and rescue-augmented EAS COJO runs
# ------------------------------------------------------------------------------
# The primary integrated table is the master denominator.  Sensitivity-only
# rescue targets are reported explicitly but never promoted to a primary claim.
# Before comparing results, this script proves the intended nesting relation:
# for every primary-ready target, all primary condition SNPs must be present in
# the rescue-augmented condition set.
# ============================================================================

suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv(
  "BASE_DIR",
  "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
)
PRIMARY_GCTA_DIR <- Sys.getenv(
  "PRIMARY_GCTA_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V5_EAS_primary_clumponly")
)
SENSITIVITY_GCTA_DIR <- Sys.getenv(
  "SENSITIVITY_GCTA_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V5_EAS_sensitivity_rescue")
)
PRIMARY_STEP2_DIR <- Sys.getenv(
  "PRIMARY_STEP2_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "ancestry_clumped_1000G_V5_EAS_primary_clumponly")
)
SENSITIVITY_STEP2_DIR <- Sys.getenv(
  "SENSITIVITY_STEP2_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "ancestry_clumped_1000G_V5_EAS_sensitivity_rescue")
)

PRIMARY_RESULTS <- Sys.getenv(
  "PRIMARY_RESULTS",
  file.path(PRIMARY_GCTA_DIR, "joint_COJO_conditional_results_integrated.tsv")
)
SENSITIVITY_RESULTS <- Sys.getenv(
  "SENSITIVITY_RESULTS",
  file.path(SENSITIVITY_GCTA_DIR, "joint_COJO_conditional_results_integrated.tsv")
)
PRIMARY_STEP2_MANIFEST <- Sys.getenv(
  "PRIMARY_STEP2_MANIFEST",
  file.path(PRIMARY_STEP2_DIR, "joint_COJO_condition_manifest_ancestry_clumped.tsv")
)
SENSITIVITY_STEP2_MANIFEST <- Sys.getenv(
  "SENSITIVITY_STEP2_MANIFEST",
  file.path(SENSITIVITY_STEP2_DIR, "joint_COJO_condition_manifest_ancestry_clumped.tsv")
)
OUT_FILE <- Sys.getenv(
  "OUT_FILE",
  file.path(BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V5_EAS_dualmode_comparison.tsv")
)
OUT_SUMMARY <- Sys.getenv(
  "OUT_SUMMARY",
  sub("\\.tsv$", "_summary.tsv", OUT_FILE)
)

STATUS_EXACT <- "known_exact_variant"
STATUS_INDEPENDENT <- "conditionally_independent_given_EAS_prior_signals"
STATUS_NOT_GWS <- "not_conditionally_GWS_given_EAS_prior_signals"
STATUS_HIGH_LD <- "not_independent_high_LD_with_prior_signal"
STATUS_UNASSESSABLE <- c(
  "conditionally_unassessable", "conditionally_unassessable_collinearity"
)

required_files <- c(
  PRIMARY_RESULTS, SENSITIVITY_RESULTS,
  PRIMARY_STEP2_MANIFEST, SENSITIVITY_STEP2_MANIFEST
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0L) {
  stop("Missing required dual-mode input file(s):\n", paste(missing_files, collapse = "\n"))
}

split_snps <- function(x) {
  if (is.na(x) || !nzchar(trimws(as.character(x)))) return(character())
  unique(trimws(unlist(strsplit(as.character(x), ";", fixed = TRUE))))
}

as_logical_flag <- function(x) {
  if (is.logical(x)) return(x)
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x_chr))
  out[x_chr %in% c("1", "true", "t", "yes", "y")] <- TRUE
  out[x_chr %in% c("0", "false", "f", "no", "n")] <- FALSE
  out
}

require_unique <- function(dt, file_label) {
  if (anyDuplicated(dt$target_id)) {
    stop(file_label, " contains duplicate target_id rows.")
  }
}

primary <- fread(PRIMARY_RESULTS, showProgress = FALSE)
sensitivity <- fread(SENSITIVITY_RESULTS, showProgress = FALSE)
for (nm in c("target_id", "EAS_run_mode", "final_cojo_status")) {
  if (!nm %in% names(primary)) stop("Primary results missing column: ", nm)
  if (!nm %in% names(sensitivity)) stop("Sensitivity results missing column: ", nm)
}
primary[, target_id := as.character(target_id)]
sensitivity[, target_id := as.character(target_id)]
require_unique(primary, "Primary results")
require_unique(sensitivity, "Sensitivity results")
if (!setequal(unique(primary$EAS_run_mode), "primary")) {
  stop("PRIMARY_RESULTS must contain only EAS_run_mode=primary.")
}
if (!setequal(unique(sensitivity$EAS_run_mode), "sensitivity")) {
  stop("SENSITIVITY_RESULTS must contain only EAS_run_mode=sensitivity.")
}

manifest_columns <- c(
  "target_id", "EAS_run_mode", "EAS_primary_ready",
  "condition_snp_ids_available", "primary_gate_failure_reason"
)
primary_m <- fread(PRIMARY_STEP2_MANIFEST, select = manifest_columns, showProgress = FALSE)
sensitivity_m <- fread(SENSITIVITY_STEP2_MANIFEST, select = manifest_columns, showProgress = FALSE)
for (dt in list(primary_m, sensitivity_m)) {
  if (anyDuplicated(dt$target_id)) stop("Step2 manifest contains duplicate target_id rows.")
}
primary_m[, target_id := as.character(target_id)]
sensitivity_m[, target_id := as.character(target_id)]
primary_m[, EAS_primary_ready := as_logical_flag(EAS_primary_ready)]
sensitivity_m[, EAS_primary_ready := as_logical_flag(EAS_primary_ready)]
if (!setequal(unique(primary_m$EAS_run_mode), "primary")) {
  stop("PRIMARY_STEP2_MANIFEST must contain only EAS_run_mode=primary.")
}
if (!setequal(unique(sensitivity_m$EAS_run_mode), "sensitivity")) {
  stop("SENSITIVITY_STEP2_MANIFEST must contain only EAS_run_mode=sensitivity.")
}

primary_ready <- primary_m[EAS_primary_ready %in% TRUE]
sensitivity_ready_ids <- sensitivity_m[EAS_primary_ready %in% TRUE, target_id]
missing_sensitivity_ready <- setdiff(primary_ready$target_id, sensitivity_ready_ids)
if (length(missing_sensitivity_ready) > 0L) {
  stop(
    "Primary-ready targets are absent from the sensitivity ready set (n=",
    length(missing_sensitivity_ready), "). Examples: ",
    paste(head(missing_sensitivity_ready, 20L), collapse = ", ")
  )
}

sensitivity_ready <- sensitivity_m[target_id %chin% primary_ready$target_id]
setkey(primary_ready, target_id)
setkey(sensitivity_ready, target_id)
nesting <- sensitivity_ready[primary_ready, nomatch = 0L]
nesting_failures <- nesting[vapply(seq_len(.N), function(i) {
  primary_snps <- split_snps(i.condition_snp_ids_available[[i]])
  sensitivity_snps <- split_snps(condition_snp_ids_available[[i]])
  !all(primary_snps %chin% sensitivity_snps)
}, logical(1L))]
if (nrow(nesting_failures) > 0L) {
  stop(
    "Primary condition SNP set is not a subset of the sensitivity set for ",
    nrow(nesting_failures), " target(s). Examples: ",
    paste(head(nesting_failures$target_id, 20L), collapse = ", ")
  )
}

sens_keep <- c(
  "target_id", "final_cojo_status", "conditional_result_code",
  "conditional_evidence_grade", "conditional_p",
  "condition_set_effectively_complete_GCTA", "high_ld_evidence_valid",
  "step2_primary_gate_failure_reason"
)
sens_keep <- intersect(sens_keep, names(sensitivity))
sensitivity_compare <- copy(sensitivity[, ..sens_keep])
setnames(
  sensitivity_compare,
  setdiff(names(sensitivity_compare), "target_id"),
  paste0("sensitivity_", setdiff(names(sensitivity_compare), "target_id"))
)
comparison <- merge(primary, sensitivity_compare, by = "target_id", all.x = TRUE, sort = FALSE)

sens_status <- comparison$sensitivity_final_cojo_status
comparison[, dualmode_robustness_status := fcase(
  final_cojo_status == STATUS_EXACT,
  "not_applicable_exact_known_variant",
  final_cojo_status == STATUS_INDEPENDENT & sens_status == STATUS_INDEPENDENT,
  "robust_independent_after_rescue_sensitivity",
  final_cojo_status == STATUS_INDEPENDENT &
    sens_status %chin% c(STATUS_NOT_GWS, STATUS_HIGH_LD),
  "independence_sensitive_to_direct_KoGES_rescue",
  final_cojo_status == STATUS_INDEPENDENT & sens_status %chin% STATUS_UNASSESSABLE,
  "sensitivity_unassessable",
  final_cojo_status == STATUS_INDEPENDENT & is.na(sens_status),
  "sensitivity_result_missing",
  final_cojo_status %chin% c(STATUS_NOT_GWS, STATUS_HIGH_LD),
  "primary_not_independent",
  step2_primary_gate_failure_reason == "no_EAS_1000G_clump_defined_signal_for_primary" &
    !is.na(sens_status) & !(sens_status %chin% STATUS_UNASSESSABLE),
  "sensitivity_only_rescue_target_no_primary_claim",
  final_cojo_status %chin% STATUS_UNASSESSABLE,
  "primary_unassessable",
  default = "not_primary_or_unassessable"
)]
comparison[, robust_additional_independent_association :=
  dualmode_robustness_status == "robust_independent_after_rescue_sensitivity"]

setorder(comparison, lead_row_id, target_id)
dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
fwrite(comparison, OUT_FILE, sep = "\t", quote = FALSE, na = "NA")
summary_dt <- comparison[, .N, by = dualmode_robustness_status][order(-N)]
setnames(summary_dt, "N", "n_targets")
fwrite(summary_dt, OUT_SUMMARY, sep = "\t", quote = FALSE, na = "NA")

cat("Primary results: ", PRIMARY_RESULTS, "\n", sep = "")
cat("Sensitivity results: ", SENSITIVITY_RESULTS, "\n", sep = "")
cat("Primary-ready subset checks: ", nrow(primary_ready), " passed\n", sep = "")
cat("Dual-mode comparison: ", OUT_FILE, "\n", sep = "")
cat("Comparison summary: ", OUT_SUMMARY, "\n", sep = "")
