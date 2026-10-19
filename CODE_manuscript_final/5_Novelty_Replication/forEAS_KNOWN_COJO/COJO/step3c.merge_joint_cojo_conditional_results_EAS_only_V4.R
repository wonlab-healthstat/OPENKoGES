#!/usr/bin/env Rscript

# ==============================================================================
# EAS-only Step3c: integrate Step1 + EAS-primary GCTA-COJO (full denominator)
# ------------------------------------------------------------------------------
# Builds one row per EAS-only Step1 target.  EUR_supporting is not required and
# is not used to change the EAS-primary verdict.
#
# Joint engine (do not point this OUT_FILE at a live joint V4 merge):
#   ../joint_EAS_EUR/step3c.merge_joint_cojo_conditional_results_FINAL_V4.R
#
# Step1 vs Step2 (EAS-only):
#   - Every non-exact Step1 target enters Step2 regardless of direct KoGES rsID
#     mapping. Step2 clumps exact 1000G EAS matches and, when
#     DIRECT_KOGES_RESCUE_IN_CORE=TRUE, inserts 1000G-unmatched KoGES-exact
#     rescues into the EAS-primary condition set.
#   - EAS_primary_ready is EAS external GWS prior evidence plus a usable KoGES
#     condition set (1000G clump representatives and/or direct rescue). A 1000G
#     clump is not required. EUR coverage is never a primary gate here.
#     Unresolved raw hits (partial_unresolved_prior) do not fail the gate.
#     Do not interpret primary-ready as "every reported EAS GWS row was a
#     covariate"; it means every representable EAS prior signal defined by
#     the conditioning workflow was accounted for.
#   - Same-coordinate SNP-indel and shared-allele encodings remain
#     representation-discordant unless exact identity or an LD proxy is shown.
#
# Shard mode (BATCH_RESULTS_GLOB):
#   - merges matching Step3b shard tables
#   - hard-stops on conflicting duplicate target_id rows
#   - requires every Step2 EAS-primary ready target to be present in Step3b
#
# Final paper-facing status (final_cojo_status):
#   known_exact_variant
#   conditionally_independent_given_EAS_prior_signals              (primary)
#   not_conditionally_GWS_given_EAS_prior_signals                  (primary)
#   not_independent_high_LD_with_prior_signal
#   conditionally_unassessable[_collinearity]
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv(
  "BASE_DIR",
  "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
)

STEP1_MANIFEST <- Sys.getenv(
  "STEP1_MANIFEST",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "step1",
    "joint_COJO_condition_manifest_preLD.tsv"
  )
)

STEP2_MANIFEST <- Sys.getenv(
  "STEP2_MANIFEST",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "ancestry_clumped_1000G_V4",
    "joint_COJO_condition_manifest_ancestry_clumped.tsv"
  )
)

BATCH_RESULTS <- Sys.getenv(
  "BATCH_RESULTS",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V4",
    "joint_COJO_conditional_results.EAS_primary.tsv"
  )
)

BATCH_RESULTS_GLOB <- Sys.getenv("BATCH_RESULTS_GLOB", "")

# Optional explicitly named EUR-supporting result table.  This avoids relying
# on a broad glob that can accidentally absorb output from another run.
EUR_BATCH_RESULTS <- Sys.getenv("EUR_BATCH_RESULTS", "")

READY_TARGETS_FILE <- Sys.getenv(
  "READY_TARGETS_FILE",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "ancestry_clumped_1000G_V4",
    "joint_COJO_EAS_primary_targets.tsv"
  )
)

# Default TRUE in all modes so incomplete Step3b merges fail loudly instead of
# silently absorbing missing ready targets as unassessable.
REQUIRE_READY_TARGETS_QC <- tolower(Sys.getenv(
  "REQUIRE_READY_TARGETS_QC",
  "true"
)) %in% c("1", "true", "t", "yes", "y")

# EAS-only: EUR supporting is not part of this analysis.
REQUIRE_EUR_SUPPORTING_QC <- tolower(Sys.getenv(
  "REQUIRE_EUR_SUPPORTING_QC",
  "false"
)) %in% c("1", "true", "t", "yes", "y")

OUT_FILE <- Sys.getenv(
  "OUT_FILE",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V4",
    "joint_COJO_conditional_results_integrated.tsv"
  )
)

OUT_SUMMARY <- Sys.getenv(
  "OUT_SUMMARY",
  file.path(
    BASE_DIR, "COJO", "EAS-only", "GCTA_conditional_V4",
    "joint_COJO_conditional_status_summary.tsv"
  )
)

OUT_PAPER <- Sys.getenv(
  "OUT_PAPER",
  file.path(
    dirname(OUT_FILE),
    "joint_COJO_conditional_results_paper.tsv"
  )
)

STATUS_EXACT <- "known_exact_variant"
STATUS_PRIMARY <- "preLD_candidates_available"
STATUS_INDEPENDENT <- "conditionally_independent_at_known_locus"
STATUS_NOT_GWS <- "not_conditionally_GWS_given_prior_signals"
STATUS_UNASSESSABLE <- "conditionally_unassessable"
STATUS_UNASSESSABLE_COLLINEARITY <- "conditionally_unassessable_collinearity"
STATUS_HIGH_LD_KNOWN <- "not_independent_high_LD_with_prior_signal"
STATUS_EUR_ONLY_NO_EAS_PRIMARY <- "EUR_only_prior_no_EAS_primary_test"

# The primary claim is EAS prior signals in KoGES (1000G clumps and/or
# direct KoGES-exact rescue). EUR supporting results are not produced here.
STATUS_INDEPENDENT_PRIMARY <- "conditionally_independent_given_EAS_prior_signals"
STATUS_NOT_GWS_PRIMARY <- "not_conditionally_GWS_given_EAS_prior_signals"

step3b_cols <- c(
  "analysis_id", "analysis_set", "target_id", "phenotype", "lead_SNP",
  "n_condition_SNP_requested", "n_condition_SNP_in_ma",
  "n_condition_SNP_retained_by_GCTA",
  "n_condition_SNP_missing", "n_condition_SNP_excluded_by_gcta",
  "condition_set_complete", "condition_set_complete_GCTA",
  "condition_set_effectively_complete_GCTA", "condition_set_status_GCTA",
  "condition_SNP_dropped_by_GCTA", "condition_SNP_retained_by_GCTA_ids",
  "n_condition_SNP_dropped_redundant_highLD",
  "condition_SNP_dropped_redundant_highLD",
  "condition_SNP_unresolved_after_GCTA",
  "max_dropped_condition_r2_to_retained", "condition_drop_LD_status",
  "max_lead_condition_r2_KoGES", "lead_highLD_condition_snps",
  "lead_condition_LD_status", "collinearity_suspected",
  "conditional_reason_code",
  "gwas_N_median", "ld_reference_N", "ld_reference_to_GWAS_N_ratio",
  "ld_reference_N_warning",
  "n_condition_SNP_step2", "n_condition_SNP_used",
  "marginal_beta", "marginal_se", "marginal_p",
  "conditional_beta", "conditional_se", "conditional_p",
  "conditional_status", "GCTA_run_status"
)

conflict_cols <- setdiff(step3b_cols, "analysis_id")

same_value <- function(x) {
  x_chr <- ifelse(is.na(x), "__NA__", as.character(x))
  length(unique(x_chr)) == 1L
}

as_logical_flag <- function(x) {
  if (is.logical(x)) {
    return(x)
  }
  x_chr <- tolower(trimws(as.character(x)))
  out <- rep(NA, length(x_chr))
  out[x_chr %in% c("1", "true", "t", "yes", "y")] <- TRUE
  out[x_chr %in% c("0", "false", "f", "no", "n")] <- FALSE
  out
}

resolve_step3b_files <- function() {
  if (nzchar(BATCH_RESULTS_GLOB)) {
    files <- Sys.glob(BATCH_RESULTS_GLOB)
    files <- files[file.exists(files)]
    if (length(files) == 0L) {
      stop("No files matched BATCH_RESULTS_GLOB: ", BATCH_RESULTS_GLOB)
    }
    return(sort(unique(files)))
  }
  if (!file.exists(BATCH_RESULTS)) {
    stop("Missing Step3b results file: ", BATCH_RESULTS)
  }
  BATCH_RESULTS
}

if (!file.exists(STEP1_MANIFEST)) {
  stop("Missing Step1 manifest: ", STEP1_MANIFEST)
}
if (!file.exists(STEP2_MANIFEST)) {
  stop("Missing Step2 manifest: ", STEP2_MANIFEST)
}

input_files <- resolve_step3b_files()
if (nzchar(EUR_BATCH_RESULTS)) {
  if (!file.exists(EUR_BATCH_RESULTS)) {
    stop("Missing EUR-supporting Step3b results file: ", EUR_BATCH_RESULTS)
  }
  input_files <- unique(c(input_files, EUR_BATCH_RESULTS))
}
shard_mode <- nzchar(BATCH_RESULTS_GLOB) || length(input_files) > 1L
cat("Step1 manifest:", STEP1_MANIFEST, "\n")
cat("Step2 manifest:", STEP2_MANIFEST, "\n")
cat("Merging", length(input_files), "Step3b result file(s)\n")
for (f in input_files) cat("  -", f, "\n")

parts <- lapply(input_files, function(path) {
  dt <- fread(path, showProgress = FALSE)
  missing <- setdiff(step3b_cols, names(dt))
  if (length(missing) > 0L) {
    stop(path, " missing columns: ", paste(missing, collapse = ", "))
  }
  dt[, .SD, .SDcols = step3b_cols]
})

step3b <- rbindlist(parts, use.names = TRUE, fill = TRUE)
step3b[, analysis_id := as.character(analysis_id)]
step3b[, analysis_set := as.character(analysis_set)]
step3b[, target_id := as.character(target_id)]

if (anyDuplicated(step3b$analysis_id)) {
  dup_ids <- unique(step3b$analysis_id[duplicated(step3b$analysis_id)])
  conflict_ids <- character()
  for (analysis_id_i in dup_ids) {
    sub <- step3b[analysis_id == analysis_id_i]
    cols_ok <- vapply(
      conflict_cols,
      function(nm) same_value(sub[[nm]]),
      logical(1L)
    )
    if (!all(cols_ok)) {
      conflict_ids <- c(conflict_ids, analysis_id_i)
    }
  }
  if (length(conflict_ids) > 0L) {
    stop(
      "Conflicting duplicate analysis_id across Step3b inputs (",
      length(conflict_ids), " analyses). Refusing to silently pick one. ",
      "Examples: ",
      paste(head(conflict_ids, 20L), collapse = ", "),
      if (length(conflict_ids) > 20L) " ..." else ""
    )
  }
  step3b <- unique(step3b, by = "analysis_id")
  cat(
    "Duplicate analysis_id rows were identical; kept one copy each (",
    length(dup_ids), " analyses).\n",
    sep = ""
  )
}

step2_wanted <- c(
  "target_id", "cojo_ready", "cojo_condition_status",
  "EAS_primary_ready", "EUR_supporting_ready",
  "EAS_prior_present", "EUR_prior_present", "EUR_only_prior",
  "EAS_clump_signal_set_complete", "EUR_clump_signal_set_complete",
  "EAS_clump_component_complete",
  "all_source_clump_signal_set_complete",
  "EUR_supporting_gate_failure_reason",
  "prior_coverage_status", "prior_set_coverage_status", "prior_identity_status",
  "raw_prior_reference_or_exact_KoGES_fraction",
  "proxy_transferability_status", "prior_signal_set_complete",
  "prior_signal_selection_scope",
  "n_raw_prior_hits_input", "n_raw_prior_hits_invalid",
  "n_raw_prior_hits_represented_1000G",
  "n_raw_prior_hits_unclumpable_1000G",
  "n_raw_prior_hits_direct_KoGES_rescued",
  "n_raw_prior_hits_unresolved_after_direct_KoGES_rescue",
  "n_external_clump_signals", "n_clump_representatives_unrepresented",
  "EAS_n_external_clump_signals", "EAS_n_clump_representatives_unrepresented",
  "EUR_n_external_clump_signals", "EUR_n_clump_representatives_unrepresented",
  "n_clump_representatives_exact_variant", "n_clump_representatives_ld_proxy",
  "n_clump_representatives_panel_wide_proxy",
  "n_clump_representatives_eur_only_ld_proxy",
  "strict_primary_ready", "primary_gate_failure_reason",
  "sensitivity_ready",
  "prior_input_failure_reasons", "reference_coverage_reasons",
  "unclumpable_prior_reasons",
  "unrepresented_clump_signal_reasons",
  "n_condition_snps_available",
  "n_condition_snps_ancestry_clumped",
  "n_condition_snps_selected",
  "n_condition_snps_from_1000G_clumps",
  "n_condition_snps_direct_KoGES_rescue",
  "n_direct_KoGES_rescued_hits_selected",
  "EAS_n_condition_snps_direct_KoGES_rescue",
  "n_EAS_selected", "n_EUR_selected"
)
step2_header <- names(fread(STEP2_MANIFEST, nrows = 0L, showProgress = FALSE))
step2_select <- intersect(step2_wanted, step2_header)
missing_step2_core <- setdiff(
  c("target_id", "cojo_ready", "cojo_condition_status"),
  step2_select
)
if (length(missing_step2_core) > 0L) {
  stop(
    "Step2 manifest missing required column(s): ",
    paste(missing_step2_core, collapse = ", ")
  )
}
step2 <- fread(STEP2_MANIFEST, select = step2_select, showProgress = FALSE)
if (!"prior_coverage_status" %in% names(step2)) {
  step2[, prior_coverage_status := NA_character_]
}
if (!"prior_set_coverage_status" %in% names(step2)) {
  step2[, prior_set_coverage_status := NA_character_]
}
if (!"prior_identity_status" %in% names(step2)) {
  step2[, prior_identity_status := NA_character_]
}
if (!"proxy_transferability_status" %in% names(step2)) {
  step2[, proxy_transferability_status := NA_character_]
}
if (!"n_raw_prior_hits_represented_1000G" %in% names(step2)) {
  step2[, n_raw_prior_hits_represented_1000G := NA_integer_]
}
if (!"reference_coverage_reasons" %in% names(step2)) {
  step2[, reference_coverage_reasons := NA_character_]
}
step2[, target_id := as.character(target_id)]
if (anyDuplicated(step2$target_id)) {
  stop("Step2 manifest has duplicate target_id.")
}
step2[, cojo_ready := as_logical_flag(cojo_ready)]
if (anyNA(step2$cojo_ready)) {
  stop("Step2 manifest cojo_ready contains non-TRUE/FALSE values.")
}
rename_from <- c(
  "cojo_condition_status",
  "EAS_primary_ready", "EUR_supporting_ready",
  "EAS_prior_present", "EUR_prior_present", "EUR_only_prior",
  "EAS_clump_signal_set_complete", "EUR_clump_signal_set_complete",
  "EAS_clump_component_complete",
  "all_source_clump_signal_set_complete",
  "EUR_supporting_gate_failure_reason",
  "prior_coverage_status",
  "prior_set_coverage_status",
  "prior_identity_status",
  "raw_prior_reference_or_exact_KoGES_fraction",
  "proxy_transferability_status",
  "prior_signal_set_complete",
  "prior_signal_selection_scope",
  "n_raw_prior_hits_input",
  "n_raw_prior_hits_invalid",
  "n_raw_prior_hits_represented_1000G",
  "n_raw_prior_hits_unclumpable_1000G",
  "n_raw_prior_hits_direct_KoGES_rescued",
  "n_raw_prior_hits_unresolved_after_direct_KoGES_rescue",
  "n_external_clump_signals",
  "n_clump_representatives_unrepresented",
  "EAS_n_external_clump_signals",
  "EAS_n_clump_representatives_unrepresented",
  "EUR_n_external_clump_signals",
  "EUR_n_clump_representatives_unrepresented",
  "n_clump_representatives_exact_variant",
  "n_clump_representatives_ld_proxy",
  "n_clump_representatives_panel_wide_proxy",
  "n_clump_representatives_eur_only_ld_proxy",
  "strict_primary_ready",
  "primary_gate_failure_reason",
  "sensitivity_ready",
  "prior_input_failure_reasons",
  "reference_coverage_reasons",
  "unclumpable_prior_reasons",
  "unrepresented_clump_signal_reasons",
  "n_condition_snps_available",
  "n_condition_snps_ancestry_clumped",
  "n_condition_snps_selected",
  "n_condition_snps_from_1000G_clumps",
  "n_condition_snps_direct_KoGES_rescue",
  "n_direct_KoGES_rescued_hits_selected",
  "EAS_n_condition_snps_direct_KoGES_rescue",
  "n_EAS_selected",
  "n_EUR_selected"
)
rename_to <- c(
  "step2_cojo_condition_status",
  "step2_EAS_primary_ready", "step2_EUR_supporting_ready",
  "step2_EAS_prior_present", "step2_EUR_prior_present", "step2_EUR_only_prior",
  "step2_EAS_clump_signal_set_complete",
  "step2_EUR_clump_signal_set_complete",
  "step2_EAS_clump_component_complete",
  "step2_all_source_clump_signal_set_complete",
  "step2_EUR_supporting_gate_failure_reason",
  "step2_prior_coverage_status",
  "step2_prior_set_coverage_status",
  "step2_prior_identity_status",
  "step2_raw_prior_reference_or_exact_KoGES_fraction",
  "step2_proxy_transferability_status",
  "step2_prior_signal_set_complete",
  "step2_prior_signal_selection_scope",
  "step2_n_raw_prior_hits_input",
  "step2_n_raw_prior_hits_invalid",
  "step2_n_raw_prior_hits_represented_1000G",
  "step2_n_raw_prior_hits_unclumpable_1000G",
  "step2_n_raw_prior_hits_direct_KoGES_rescued",
  "step2_n_raw_prior_hits_unresolved_after_direct_KoGES_rescue",
  "step2_n_external_clump_signals",
  "step2_n_clump_representatives_unrepresented",
  "step2_EAS_n_external_clump_signals",
  "step2_EAS_n_clump_representatives_unrepresented",
  "step2_EUR_n_external_clump_signals",
  "step2_EUR_n_clump_representatives_unrepresented",
  "step2_n_clump_representatives_exact_variant",
  "step2_n_clump_representatives_ld_proxy",
  "step2_n_clump_representatives_panel_wide_proxy",
  "step2_n_clump_representatives_eur_only_ld_proxy",
  "step2_strict_primary_ready",
  "step2_primary_gate_failure_reason",
  "step2_sensitivity_ready",
  "step2_prior_input_failure_reasons",
  "step2_reference_coverage_reasons",
  "step2_unclumpable_prior_reasons",
  "step2_unrepresented_clump_signal_reasons",
  "step2_n_condition_snps_available",
  "step2_n_condition_snps_ancestry_clumped",
  "step2_n_condition_snps_selected",
  "step2_n_condition_snps_from_1000G_clumps",
  "step2_n_condition_snps_direct_KoGES_rescue",
  "step2_n_direct_KoGES_rescued_hits_selected",
  "step2_EAS_n_condition_snps_direct_KoGES_rescue",
  "step2_n_EAS_selected",
  "step2_n_EUR_selected"
)
present_rename <- rename_from %in% names(step2)
setnames(step2, rename_from[present_rename], rename_to[present_rename])
if ("step2_prior_signal_set_complete" %in% names(step2)) {
  step2[, step2_prior_signal_set_complete :=
    as_logical_flag(step2_prior_signal_set_complete)]
}
if ("step2_sensitivity_ready" %in% names(step2)) {
  step2[, step2_sensitivity_ready := as_logical_flag(step2_sensitivity_ready)]
}
for (nm in c(
  "step2_EAS_primary_ready", "step2_EUR_supporting_ready",
  "step2_EAS_prior_present", "step2_EUR_prior_present", "step2_EUR_only_prior",
  "step2_EAS_clump_signal_set_complete",
  "step2_EUR_clump_signal_set_complete",
  "step2_EAS_clump_component_complete",
  "step2_all_source_clump_signal_set_complete"
)) {
  if (!nm %in% names(step2)) step2[, (nm) := FALSE]
  step2[, (nm) := as_logical_flag(get(nm))]
}

ready_ids_from_manifest <- unique(
  step2$target_id[step2$step2_EAS_primary_ready %in% TRUE]
)
eur_ready_ids_from_manifest <- unique(
  step2$target_id[step2$step2_EUR_supporting_ready %in% TRUE]
)

allowed_analysis_sets <- c("EAS_primary", "EUR_supporting", "legacy_all_sources")
unexpected_analysis_sets <- setdiff(unique(step3b$analysis_set), allowed_analysis_sets)
if (length(unexpected_analysis_sets) > 0L) {
  stop(
    "Unexpected Step3b analysis_set value(s): ",
    paste(unexpected_analysis_sets, collapse = ", ")
  )
}

if (REQUIRE_READY_TARGETS_QC || shard_mode) {
  if (file.exists(READY_TARGETS_FILE)) {
    ready_file <- fread(
      READY_TARGETS_FILE, select = "target_id", showProgress = FALSE
    )
    ready_ids_file <- unique(as.character(ready_file$target_id))
    only_manifest <- setdiff(ready_ids_from_manifest, ready_ids_file)
    only_file <- setdiff(ready_ids_file, ready_ids_from_manifest)
    if (length(only_manifest) > 0L || length(only_file) > 0L) {
      stop(
        "Step2 manifest EAS_primary_ready=TRUE disagrees with READY_TARGETS_FILE. ",
        "only_in_manifest=", length(only_manifest),
        " only_in_ready_file=", length(only_file),
        ". Examples manifest: ",
        paste(head(only_manifest, 10L), collapse = ", "),
        " | Examples file: ",
        paste(head(only_file, 10L), collapse = ", ")
      )
    }
  } else {
    warning(
      "READY_TARGETS_FILE missing (", READY_TARGETS_FILE,
      "); using Step2 manifest EAS_primary_ready=TRUE as the ready set.",
      call. = FALSE
    )
  }

  ready_ids <- ready_ids_from_manifest
  eas_result_ids <- unique(as.character(
    step3b[analysis_set == "EAS_primary", target_id]
  ))
  missing_targets <- setdiff(ready_ids, eas_result_ids)
  extra_targets <- setdiff(eas_result_ids, ready_ids)
  cat("EAS-primary ready-target QC\n")
  cat("  Expected Step2 EAS-primary targets:", length(ready_ids), "\n")
  cat("  Step3b EAS-primary targets:", length(eas_result_ids), "\n")
  if (length(missing_targets) > 0L) {
    stop(
      "Step3b EAS-primary incomplete vs Step2 ready targets: missing ",
      length(missing_targets), " target_id(s). Examples: ",
      paste(head(missing_targets, 20L), collapse = ", "),
      if (length(missing_targets) > 20L) " ..." else ""
    )
  }
  if (length(extra_targets) > 0L) {
    stop(
      "Step3b EAS-primary has ", length(extra_targets),
      " target_id(s) not in Step2 EAS-primary ready targets. Examples: ",
      paste(head(extra_targets, 20L), collapse = ", "),
      if (length(extra_targets) > 20L) " ..." else ""
    )
  }
  cat("  EAS-primary QC: OK\n")

  eur_result_ids <- unique(as.character(
    step3b[analysis_set == "EUR_supporting", target_id]
  ))
  if (REQUIRE_EUR_SUPPORTING_QC) {
    eur_missing <- setdiff(eur_ready_ids_from_manifest, eur_result_ids)
    eur_extra <- setdiff(eur_result_ids, eur_ready_ids_from_manifest)
    if (length(eur_missing) > 0L || length(eur_extra) > 0L) {
      stop(
        "EUR-supporting Step3b results disagree with Step2 readiness. missing=",
        length(eur_missing), " extra=", length(eur_extra),
        ". Examples missing: ", paste(head(eur_missing, 10L), collapse = ", "),
        " | extra: ", paste(head(eur_extra, 10L), collapse = ", ")
      )
    }
    cat("  EUR-supporting QC: OK (", length(eur_result_ids), " targets)\n", sep = "")
  } else if (length(eur_result_ids) > 0L) {
    cat("  EUR-supporting rows available (not completeness-gated): ",
        length(eur_result_ids), "\n", sep = "")
  }
}

manifest <- fread(STEP1_MANIFEST, showProgress = FALSE)
need_manifest <- c(
  "target_id", "lead_row_id", "koges_phenotype", "koges_pheno_code",
  "koges_lead_variant", "koges_chr", "koges_bp", "cojo_condition_status"
)
miss_manifest <- setdiff(need_manifest, names(manifest))
if (length(miss_manifest) > 0L) {
  stop(
    "Step1 manifest missing columns: ",
    paste(miss_manifest, collapse = ", ")
  )
}
manifest[, target_id := as.character(target_id)]
if (anyDuplicated(manifest$target_id)) {
  stop("Step1 manifest has duplicate target_id.")
}

unknown_status <- setdiff(
  unique(manifest$cojo_condition_status),
  c(STATUS_EXACT, STATUS_PRIMARY)
)
if (length(unknown_status) > 0L) {
  stop(
    "Unexpected cojo_condition_status in Step1 manifest: ",
    paste(unknown_status, collapse = ", ")
  )
}

integrated <- merge(
  manifest[, ..need_manifest],
  step2,
  by = "target_id",
  all.x = TRUE,
  sort = FALSE
)
if (anyNA(integrated$cojo_ready)) {
  stop(
    "Step1 targets missing from Step2 manifest: ",
    sum(is.na(integrated$cojo_ready)),
    ". Examples: ",
    paste(
      head(integrated$target_id[is.na(integrated$cojo_ready)], 20L),
      collapse = ", "
    )
  )
}
if (!"step2_n_clump_representatives_unrepresented" %in% names(integrated)) {
  integrated[, step2_n_clump_representatives_unrepresented := NA_integer_]
}
if (!"step2_sensitivity_ready" %in% names(integrated)) {
  integrated[, step2_sensitivity_ready := FALSE]
}
if (!"step2_prior_set_coverage_status" %in% names(integrated)) {
  integrated[, step2_prior_set_coverage_status := NA_character_]
}
if (!"step2_proxy_transferability_status" %in% names(integrated)) {
  integrated[, step2_proxy_transferability_status := NA_character_]
}

eas_step3b <- copy(step3b[analysis_set == "EAS_primary"])
eur_step3b <- copy(step3b[analysis_set == "EUR_supporting"])
if (anyDuplicated(eas_step3b$target_id)) {
  stop("More than one EAS-primary Step3b row was found for a target.")
}
if (anyDuplicated(eur_step3b$target_id)) {
  stop("More than one EUR-supporting Step3b row was found for a target.")
}

# EAS-primary results retain their original Step3b column names because they
# define final_cojo_status.  EUR results are explicitly prefixed and never
# overwrite the EAS conclusion.
integrated <- merge(
  integrated, eas_step3b,
  by = "target_id", all.x = TRUE, sort = FALSE
)
if (nrow(eur_step3b) > 0L) {
  eur_rename <- setdiff(names(eur_step3b), "target_id")
  setnames(eur_step3b, eur_rename, paste0("EUR_supporting_", eur_rename))
  integrated <- merge(
    integrated, eur_step3b,
    by = "target_id", all.x = TRUE, sort = FALSE
  )
}

# Only EAS-primary-ready targets must have an EAS-primary Step3b result.
ready_missing_step3b <- integrated[
  step2_EAS_primary_ready %in% TRUE & is.na(conditional_status), target_id
]
if ((REQUIRE_READY_TARGETS_QC || shard_mode) && length(ready_missing_step3b) > 0L) {
  stop(
    "EAS-primary-ready targets missing Step3b results: ",
    length(ready_missing_step3b), ". Examples: ",
    paste(head(ready_missing_step3b, 20L), collapse = ", "),
    if (length(ready_missing_step3b) > 20L) " ..." else ""
  )
}

n_primary_not_ready <- integrated[
  cojo_condition_status == STATUS_PRIMARY & !(step2_EAS_primary_ready %in% TRUE),
  .N
]
cat(
  "Step1 primary targets without an EAS-primary COJO test: ",
  n_primary_not_ready, "\n", sep = ""
)

integrated[, final_cojo_status := fcase(
  cojo_condition_status == STATUS_EXACT, STATUS_EXACT,
  cojo_condition_status == STATUS_PRIMARY & step2_EUR_only_prior %in% TRUE,
  STATUS_EUR_ONLY_NO_EAS_PRIMARY,
  cojo_condition_status == STATUS_PRIMARY & step2_EAS_primary_ready %in% TRUE &
    conditional_status %chin% c(
      STATUS_INDEPENDENT, STATUS_NOT_GWS, STATUS_HIGH_LD_KNOWN,
      STATUS_UNASSESSABLE, STATUS_UNASSESSABLE_COLLINEARITY
    ),
  conditional_status,
  cojo_condition_status == STATUS_PRIMARY & step2_EAS_primary_ready %in% TRUE,
  STATUS_UNASSESSABLE,
  cojo_condition_status == STATUS_PRIMARY,
  STATUS_UNASSESSABLE,
  default = STATUS_UNASSESSABLE
)]

if (!"step2_EAS_primary_ready" %in% names(integrated)) {
  integrated[, step2_EAS_primary_ready := FALSE]
}
if (!"step2_primary_gate_failure_reason" %in% names(integrated)) {
  integrated[, step2_primary_gate_failure_reason := NA_character_]
}
integrated[, step2_EAS_primary_ready := as_logical_flag(step2_EAS_primary_ready)]
integrated[, strict_primary_gate_passed := step2_EAS_primary_ready %in% TRUE]
integrated[
  final_cojo_status == STATUS_NOT_GWS,
  final_cojo_status := STATUS_NOT_GWS_PRIMARY
]
integrated[
  final_cojo_status == STATUS_INDEPENDENT,
  final_cojo_status := STATUS_INDEPENDENT_PRIMARY
]

# Keep the positional-known label separate from conditional evidence. Every
# Step1 target entered because of an external GWS hit within +/-1 Mb; the EAS
# primary test answers a narrower, ancestry-matched conditional question.
integrated[, positional_known_by_external_GWS_within_1Mb := TRUE]
integrated[, positional_known_status := "known_by_external_GWS_within_1Mb"]
integrated[, conditional_assessment_status := fcase(
  cojo_condition_status == STATUS_EXACT,
  "not_tested_exact_known_variant",
  final_cojo_status == STATUS_EUR_ONLY_NO_EAS_PRIMARY,
  "not_assessed_EUR_only_prior",
  final_cojo_status == STATUS_INDEPENDENT_PRIMARY,
  STATUS_INDEPENDENT,
  final_cojo_status == STATUS_NOT_GWS_PRIMARY,
  STATUS_NOT_GWS,
  final_cojo_status == STATUS_HIGH_LD_KNOWN,
  STATUS_HIGH_LD_KNOWN,
  final_cojo_status == STATUS_UNASSESSABLE_COLLINEARITY,
  STATUS_UNASSESSABLE_COLLINEARITY,
  default = STATUS_UNASSESSABLE
)]
integrated[, conditional_result_code := fcase(
  conditional_assessment_status == "not_tested_exact_known_variant",
  "exact_known",
  conditional_assessment_status == "not_assessed_EUR_only_prior",
  "EUR_only_no_EAS_primary",
  conditional_assessment_status == STATUS_INDEPENDENT,
  "independent_given_EAS_priors",
  conditional_assessment_status == STATUS_NOT_GWS,
  "attenuated_given_EAS_priors",
  conditional_assessment_status == STATUS_HIGH_LD_KNOWN,
  "high_LD_with_EAS_prior",
  default = "unassessable"
)]
integrated[, conditional_model_complete_for_primary :=
  strict_primary_gate_passed &
    condition_set_effectively_complete_GCTA %in% TRUE]
integrated[, conditional_evidence_grade := fcase(
  conditional_assessment_status == "not_tested_exact_known_variant",
  "not_applicable_exact_variant",
  conditional_assessment_status == "not_assessed_EUR_only_prior",
  "not_primary_EUR_only_prior",
  conditional_assessment_status %chin% c(
    STATUS_UNASSESSABLE, STATUS_UNASSESSABLE_COLLINEARITY
  ), "unassessable",
  conditional_assessment_status == STATUS_HIGH_LD_KNOWN &
    conditional_model_complete_for_primary,
  "high_LD_EAS_primary_ready",
  conditional_assessment_status == STATUS_HIGH_LD_KNOWN,
  "high_LD_EAS_model_incomplete",
  strict_primary_gate_passed, "primary_EAS_ready",
  default = "not_primary"
)]
integrated[, reviewer_response_interpretation := fcase(
  conditional_assessment_status == "not_tested_exact_known_variant",
  "external_exact_variant_match; conditional_test_not_needed",
  conditional_assessment_status == "not_assessed_EUR_only_prior",
  "EUR_only_external_prior; no_EAS_matched_conditional_primary_test",
  conditional_assessment_status == STATUS_INDEPENDENT & strict_primary_gate_passed,
  "positional_known_but_additional_KoGES_independent_association_given_EAS_prior_signals",
  conditional_assessment_status == STATUS_NOT_GWS & strict_primary_gate_passed,
  "positional_known_and_not_GWS_given_EAS_prior_signals_not_causal_proof",
  conditional_assessment_status == STATUS_HIGH_LD_KNOWN,
  "positional_known_and_high_LD_with_EAS_prior_signal_not_causal_proof",
  default = "positional_known_EAS_primary_conditional_test_unassessable"
)]

# Prefer phenotype/lead from Step1; keep Step3b lead_SNP when present.
integrated[
  is.na(phenotype) | !nzchar(as.character(phenotype)),
  phenotype := koges_phenotype
]
integrated[
  is.na(lead_SNP) | !nzchar(as.character(lead_SNP)),
  lead_SNP := koges_lead_variant
]

for (nm in c(
  "step2_prior_set_coverage_status", "step2_proxy_transferability_status",
  "step2_primary_gate_failure_reason",
  "step2_raw_prior_reference_or_exact_KoGES_fraction",
  "step2_EAS_clump_signal_set_complete",
  "step2_EAS_clump_component_complete",
  "step2_EUR_clump_signal_set_complete",
  "step2_all_source_clump_signal_set_complete"
)) {
  if (!nm %in% names(integrated)) integrated[, (nm) := NA_character_]
}
for (nm in c(
  "EUR_supporting_conditional_status", "EUR_supporting_conditional_p",
  "EUR_supporting_conditional_reason_code", "EUR_supporting_GCTA_run_status",
  "EUR_supporting_condition_set_effectively_complete_GCTA",
  "EUR_supporting_n_condition_SNP_requested",
  "EUR_supporting_n_condition_SNP_retained_by_GCTA"
)) {
  if (!nm %in% names(integrated)) integrated[, (nm) := NA]
}
for (nm in c(
  "step2_n_condition_snps_direct_KoGES_rescue",
  "step2_EAS_n_condition_snps_direct_KoGES_rescue",
  "step2_n_direct_KoGES_rescued_hits_selected",
  "step2_n_condition_snps_from_1000G_clumps",
  "step2_n_condition_snps_selected"
)) {
  if (!nm %in% names(integrated)) integrated[, (nm) := NA_integer_]
}

setorder(integrated, lead_row_id, target_id)
fwrite(integrated, OUT_FILE, sep = "\t", quote = FALSE, na = "NA")

# Compact paper-facing table: retain the full integrated audit above, but do not
# force readers to decode every execution/bookkeeping column.  EUR proxy
# transferability and EUR_supporting GCTA fields stay in the full table only.
paper_results <- integrated[, .(
  target_id,
  phenotype,
  lead_SNP,
  positional_known_status,
  conditional_result_code,
  conditional_evidence_grade,
  conditional_reason_code,
  final_cojo_status,
  conditional_p,
  primary_condition_set = "EAS_primary",
  EAS_clump_signal_set_complete = step2_EAS_clump_signal_set_complete,
  EAS_clump_component_complete = step2_EAS_clump_component_complete,
  n_condition_snps_direct_KoGES_rescue =
    step2_n_condition_snps_direct_KoGES_rescue,
  n_condition_snps_from_1000G_clumps =
    step2_n_condition_snps_from_1000G_clumps,
  n_condition_snps_selected = step2_n_condition_snps_selected,
  EAS_n_condition_snps_direct_KoGES_rescue =
    step2_EAS_n_condition_snps_direct_KoGES_rescue,
  EUR_clump_signal_set_complete = step2_EUR_clump_signal_set_complete,
  all_source_clump_signal_set_complete =
    step2_all_source_clump_signal_set_complete,
  raw_prior_reference_or_exact_KoGES_fraction =
    step2_raw_prior_reference_or_exact_KoGES_fraction,
  prior_set_coverage_status = step2_prior_set_coverage_status,
  condition_set_status_GCTA,
  strict_primary_gate_passed,
  conditional_model_complete_for_primary,
  primary_gate_failure_reason = step2_primary_gate_failure_reason,
  n_condition_SNP_requested,
  n_condition_SNP_retained_by_GCTA,
  max_lead_condition_r2_KoGES,
  lead_highLD_condition_snps,
  gwas_N_median,
  ld_reference_N,
  ld_reference_to_GWAS_N_ratio,
  ld_reference_N_warning
)]
fwrite(paper_results, OUT_PAPER, sep = "\t", quote = FALSE, na = "NA")

summary_dt <- integrated[, .(n_targets = .N), by = .(final_cojo_status)]
setorder(summary_dt, -n_targets)
fwrite(summary_dt, OUT_SUMMARY, sep = "\t", quote = FALSE, na = "NA")
ancestry_summary <- integrated[
  cojo_condition_status == STATUS_PRIMARY,
  .(n_targets = .N),
  by = .(
    step2_EAS_prior_present,
    step2_EAS_clump_component_complete,
    step2_EAS_clump_signal_set_complete,
    step2_EUR_only_prior,
    final_cojo_status
  )
]
setorder(ancestry_summary, step2_EAS_prior_present, -n_targets)
ancestry_summary_file <- sub(
  "\\.tsv$", "_by_ancestry_clump_completeness.tsv", OUT_SUMMARY
)
if (identical(ancestry_summary_file, OUT_SUMMARY)) {
  ancestry_summary_file <- paste0(
    OUT_SUMMARY, "_by_ancestry_clump_completeness.tsv"
  )
}
fwrite(ancestry_summary, ancestry_summary_file,
       sep = "\t", quote = FALSE, na = "NA")

if (!"step2_prior_identity_status" %in% names(integrated)) {
  integrated[, step2_prior_identity_status := NA_character_]
}
identity_summary <- integrated[
  cojo_condition_status == STATUS_PRIMARY,
  .(n_targets = .N),
  by = .(step2_prior_identity_status, cojo_ready, final_cojo_status)
]
setorder(identity_summary, step2_prior_identity_status, -n_targets)
identity_summary_file <- sub("\\.tsv$", "_by_1000G_identity.tsv", OUT_SUMMARY)
if (identical(identity_summary_file, OUT_SUMMARY)) {
  identity_summary_file <- paste0(OUT_SUMMARY, "_by_1000G_identity.tsv")
}
fwrite(identity_summary, identity_summary_file, sep = "\t", quote = FALSE, na = "NA")

grade_summary <- integrated[
  cojo_condition_status == STATUS_PRIMARY,
  .(n_targets = .N),
  by = .(
    conditional_evidence_grade, step2_primary_gate_failure_reason,
    conditional_assessment_status
  )
]
setorder(grade_summary, conditional_evidence_grade, -n_targets)
grade_summary_file <- sub("\\.tsv$", "_by_primary_gate.tsv", OUT_SUMMARY)
if (identical(grade_summary_file, OUT_SUMMARY)) {
  grade_summary_file <- paste0(OUT_SUMMARY, "_by_primary_gate.tsv")
}
fwrite(grade_summary, grade_summary_file, sep = "\t", quote = FALSE, na = "NA")

coverage_summary <- integrated[
  cojo_condition_status == STATUS_PRIMARY,
  .(n_targets = .N),
  by = .(
    step2_prior_set_coverage_status,
    step2_EAS_primary_ready,
    final_cojo_status
  )
]
setorder(coverage_summary, step2_prior_set_coverage_status, -n_targets)
coverage_summary_file <- sub("\\.tsv$", "_by_prior_set_coverage.tsv", OUT_SUMMARY)
if (identical(coverage_summary_file, OUT_SUMMARY)) {
  coverage_summary_file <- paste0(OUT_SUMMARY, "_by_prior_set_coverage.tsv")
}
fwrite(coverage_summary, coverage_summary_file, sep = "\t", quote = FALSE, na = "NA")

cat("Integrated full-denominator results:", OUT_FILE, "\n")
cat("Compact paper-facing results:", OUT_PAPER, "\n")
cat("Status summary:", OUT_SUMMARY, "\n")
cat("Ancestry-clump-completeness summary:", ancestry_summary_file, "\n")
cat("Identity-stratified summary:", identity_summary_file, "\n")
cat("Primary-gate-stratified summary:", grade_summary_file, "\n")
cat("Prior-set-coverage QC summary:", coverage_summary_file, "\n")
cat(
  "EAS-primary-ready with partial_unresolved_prior (QC, not a gate failure):",
  sum(
    integrated$step2_EAS_primary_ready %in% TRUE &
      integrated$step2_prior_set_coverage_status == "partial_unresolved_prior",
    na.rm = TRUE
  ),
  "\n"
)
cat("Total Step1 targets:", nrow(integrated), "\n")
for (i in seq_len(nrow(summary_dt))) {
  cat(
    "  ", summary_dt$final_cojo_status[i], ": ",
    summary_dt$n_targets[i], "\n",
    sep = ""
  )
}
cat(
  "GCTA-complete among Step2 EAS-primary-ready targets:",
  sum(
    integrated$cojo_condition_status == STATUS_PRIMARY &
      integrated$cojo_ready %in% TRUE &
      integrated$condition_set_complete_GCTA %in% TRUE,
    na.rm = TRUE
  ),
  "\n"
)
