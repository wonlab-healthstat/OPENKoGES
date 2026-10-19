#!/usr/bin/env Rscript

# ==============================================================================
# Step3b, EAS dual-mode GCTA-COJO
# ------------------------------------------------------------------------------
# Runs every row in the EAS-only Step2 ready table.  Condition SNPs are
# 1000G clump representatives in primary mode, plus direct KoGES-exact
# rescues in sensitivity mode.  Step3b does not re-apply the Step2 gate.
# EUR_supporting is not run.
#
# Do not point OUT_DIR at the joint_EAS_EUR GCTA_conditional_V4 directory.
#
# Implementation: sets EAS-only default paths, then runs the shared Step3b
# engine ../shared/step3b.run_joint_cojo_conditional_gcta_FINAL_V4.R
#
# Usage:
#   EAS_RUN_MODE=primary RUN_MODE=prepare RUN_GCTA=false \
#     Rscript step3b.run_joint_cojo_conditional_gcta_EAS_dualmode_V5.R
#   EAS_RUN_MODE=primary KOGES_BFILE=$REF KOGES_BFILE_FULL=$FULL N_PARALLEL=32 \
#     Rscript step3b.run_joint_cojo_conditional_gcta_EAS_dualmode_V5.R
# ==============================================================================

this <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", this[grep("^--file=", this)])
if (length(file_arg) != 1L || !nzchar(file_arg[1L])) {
  stop("Run this file with Rscript so --file= is available.")
}
root <- dirname(normalizePath(file_arg[1L], mustWork = TRUE))
engine <- normalizePath(
  file.path(root, "..", "shared", "step3b.run_joint_cojo_conditional_gcta_FINAL_V4.R"),
  mustWork = FALSE
)
if (!file.exists(engine)) {
  stop("Shared Step3b engine not found: ", engine)
}

base <- Sys.getenv(
  "BASE_DIR",
  "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
)

EAS_RUN_MODE <- tolower(trimws(Sys.getenv("EAS_RUN_MODE", "primary")))
if (!EAS_RUN_MODE %in% c("primary", "sensitivity")) {
  stop("EAS_RUN_MODE must be 'primary' or 'sensitivity'.")
}
EAS_RUN_SUFFIX <- if (EAS_RUN_MODE == "primary") {
  "EAS_primary_clumponly"
} else {
  "EAS_sensitivity_rescue"
}
if (!nzchar(Sys.getenv("STEP2_DIR"))) {
  Sys.setenv(STEP2_DIR = file.path(base, "COJO", "EAS-only", paste0("ancestry_clumped_1000G_V5_", EAS_RUN_SUFFIX)))
}
if (!nzchar(Sys.getenv("READY_TARGETS_FILE"))) {
  Sys.setenv(READY_TARGETS_FILE = file.path(
    Sys.getenv("STEP2_DIR"), "joint_COJO_ready_targets.tsv"
  ))
}
if (!nzchar(Sys.getenv("MA_DIR"))) {
  Sys.setenv(MA_DIR = file.path(base, "COJO", "EAS-only", paste0("GCTA_ma_V5_", EAS_RUN_SUFFIX)))
}
if (!nzchar(Sys.getenv("OUT_DIR"))) {
  Sys.setenv(OUT_DIR = file.path(base, "COJO", "EAS-only", paste0("GCTA_conditional_V5_", EAS_RUN_SUFFIX)))
}
# Do not set RESULTS_TAG here.  The shared engine already chooses:
#   non-sharded + one analysis_set → that set name (EAS_primary)
#   N_SHARDS / SHARD_INDEX        → shard0, shard1, ...
# Forcing EAS_primary in the wrapper would make concurrent shards
# overwrite the same aggregate TSV.

source(engine, local = FALSE)
