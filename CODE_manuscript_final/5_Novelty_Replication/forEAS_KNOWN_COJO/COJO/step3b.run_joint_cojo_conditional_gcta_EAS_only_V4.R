#!/usr/bin/env Rscript

# ==============================================================================
# Step3b, EAS-primary GCTA-COJO only
# ------------------------------------------------------------------------------
# Runs every row in the EAS-only Step2 EAS-primary job table.  Condition SNPs
# may be 1000G clump representatives and/or direct KoGES-exact rescues; Step3b
# does not re-apply the Step2 gate.  EUR_supporting is not run.
#
# Do not point OUT_DIR at the joint_EAS_EUR GCTA_conditional_V4 directory.
#
# Implementation: sets EAS-only default paths, then runs the shared Step3b
# engine ../shared/step3b.run_joint_cojo_conditional_gcta_FINAL_V4.R
#
# Usage:
#   RUN_MODE=prepare RUN_GCTA=false Rscript EAS-only/step3b.run_joint_cojo_conditional_gcta_EAS_only_V4.R
#   KOGES_BFILE=$REF KOGES_BFILE_FULL=$FULL N_PARALLEL=32 \
#     Rscript EAS-only/step3b.run_joint_cojo_conditional_gcta_EAS_only_V4.R
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
if (!nzchar(Sys.getenv("STEP2_DIR"))) {
  Sys.setenv(STEP2_DIR = file.path(base, "COJO", "EAS-only", "ancestry_clumped_1000G_V4"))
}
if (!nzchar(Sys.getenv("READY_TARGETS_FILE"))) {
  Sys.setenv(READY_TARGETS_FILE = file.path(
    Sys.getenv("STEP2_DIR"), "joint_COJO_EAS_primary_targets.tsv"
  ))
}
if (!nzchar(Sys.getenv("MA_DIR"))) {
  Sys.setenv(MA_DIR = file.path(base, "COJO", "EAS-only", "GCTA_ma_V4"))
}
if (!nzchar(Sys.getenv("OUT_DIR"))) {
  Sys.setenv(OUT_DIR = file.path(base, "COJO", "EAS-only", "GCTA_conditional_V4"))
}
# Do not set RESULTS_TAG here.  The shared engine already chooses:
#   non-sharded + one analysis_set → that set name (EAS_primary)
#   N_SHARDS / SHARD_INDEX        → shard0, shard1, ...
# Forcing EAS_primary in the wrapper would make concurrent shards
# overwrite the same aggregate TSV.

source(engine, local = FALSE)
