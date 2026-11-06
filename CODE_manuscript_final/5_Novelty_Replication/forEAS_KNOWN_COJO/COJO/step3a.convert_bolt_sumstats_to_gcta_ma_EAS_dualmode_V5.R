#!/usr/bin/env Rscript

# ==============================================================================
# Step3a, EAS dual-mode (primary clump-only / rescue sensitivity)
# ------------------------------------------------------------------------------
# Converts BOLT sumstats to GCTA .ma for phenotypes in the EAS-only Step2
# ready table.  Do not point MA_OUT_DIR at the joint_EAS_EUR GCTA_ma_V4 directory
# while that run is live.
#
# Implementation: sets EAS-only default paths, then runs the shared Step3a
# engine ../shared/step3a.convert_bolt_sumstats_to_gcta_ma_FINAL_4.R
#
# Usage:
#   EAS_RUN_MODE=primary Rscript step3a.convert_bolt_sumstats_to_gcta_ma_EAS_dualmode_V5.R
# ==============================================================================


# Resolve openKoGES repo root without hard-coding host absolute paths.
.openkoges_repo_root <- function() {
  env <- Sys.getenv("OPENKOGES_ROOT", unset = "")
  if (nzchar(env)) return(normalizePath(env, winslash = "/", mustWork = FALSE))
  args <- commandArgs(trailingOnly = FALSE)
  fa <- sub("^--file=", "", args[startsWith(args, "--file=")])
  start <- if (length(fa)) dirname(normalizePath(fa[[1]], winslash = "/", mustWork = FALSE)) else getwd()
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  for (i in seq_len(12L)) {
    if (dir.exists(file.path(d, "CODE_manuscript_final")) &&
        dir.exists(file.path(d, "DATA_manuscript_final"))) {
      return(d)
    }
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate openKoGES root; set OPENKOGES_ROOT.", call. = FALSE)
}
REPO_ROOT <- .openkoges_repo_root()
REVISION_ROOT <- Sys.getenv(
  "OPENKOGES_REVISION_ROOT",
  unset = file.path(dirname(REPO_ROOT), "openKoGES_revision")
)

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
  file.path(REVISION_ROOT, "DATA_cumulative_inci", "forBOLT-LMM", "BOLT_lead_variant")
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
if (!nzchar(Sys.getenv("MA_OUT_DIR"))) {
  Sys.setenv(MA_OUT_DIR = file.path(base, "COJO", "EAS-only", paste0("GCTA_ma_V5_", EAS_RUN_SUFFIX)))
}

source(engine, local = FALSE)
