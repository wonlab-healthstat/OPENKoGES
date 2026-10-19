#!/usr/bin/env Rscript

# ==============================================================================
# Step3a, EAS-primary ready targets only
# ------------------------------------------------------------------------------
# Converts BOLT sumstats to GCTA .ma for phenotypes in the EAS-only Step2
# ready table.  Do not point MA_OUT_DIR at the joint_EAS_EUR GCTA_ma_V4 directory
# while that run is live.
#
# Implementation: sets EAS-only default paths, then runs the shared Step3a
# engine ../shared/step3a.convert_bolt_sumstats_to_gcta_ma_FINAL_4.R
#
# Usage:
#   Rscript EAS-only/step3a.convert_bolt_sumstats_to_gcta_ma_EAS_only_V4.R
# ==============================================================================

this <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", this[grep("^--file=", this)])
if (length(file_arg) != 1L || !nzchar(file_arg[1L])) {
  stop("Run this file with Rscript so --file= is available.")
}
root <- dirname(normalizePath(file_arg[1L], mustWork = TRUE))
engine <- normalizePath(
  file.path(root, "..", "shared", "step3a.convert_bolt_sumstats_to_gcta_ma_FINAL_4.R"),
  mustWork = FALSE
)
if (!file.exists(engine)) {
  stop("Shared Step3a engine not found: ", engine)
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
if (!nzchar(Sys.getenv("MA_OUT_DIR"))) {
  Sys.setenv(MA_OUT_DIR = file.path(base, "COJO", "EAS-only", "GCTA_ma_V4"))
}

source(engine, local = FALSE)
