#!/usr/bin/env Rscript

# ==============================================================================
# Step1, EAS sources only
# ------------------------------------------------------------------------------
# Searches BBJ, TWB, and TPMI.  UKB / FinnGen / MVP are not read.
# known_exact_variant and the Step2 target list are therefore defined from
# EAS evidence alone.
#
# Implementation: forces SOURCE_LAYER=EAS and a distinct default OUT_DIR, then
# runs the shared Step1 engine.  Do not point this OUT_DIR at joint_EAS_EUR/step1.
#
# Usage:
#   Rscript EAS-only/step1.build_joint_cojo_condition_candidates_EAS_only_V2.R
#   N_PARALLEL=8 Rscript EAS-only/step1.build_joint_cojo_condition_candidates_EAS_only_V2.R
# ==============================================================================

this <- commandArgs(trailingOnly = FALSE)
file_arg <- sub("^--file=", "", this[grep("^--file=", this)])
if (length(file_arg) != 1L || !nzchar(file_arg[1L])) {
  stop("Run this file with Rscript so --file= is available.")
}
root <- dirname(normalizePath(file_arg[1L], mustWork = TRUE))
engine <- normalizePath(
  file.path(root, "..", "shared", "step1.build_joint_cojo_condition_candidates_FINAL_V2.R"),
  mustWork = FALSE
)
if (!file.exists(engine)) {
  stop("Shared Step1 engine not found: ", engine)
}

Sys.setenv(SOURCE_LAYER = "EAS")
if (!nzchar(Sys.getenv("OUT_DIR"))) {
  base <- Sys.getenv(
    "BASE_DIR",
    "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
  )
  Sys.setenv(OUT_DIR = file.path(base, "COJO", "EAS-only", "step1"))
}

source(engine, local = FALSE)
