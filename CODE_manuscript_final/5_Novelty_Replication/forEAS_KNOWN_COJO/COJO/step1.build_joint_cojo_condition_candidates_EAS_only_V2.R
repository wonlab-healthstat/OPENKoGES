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
    file.path(REVISION_ROOT, "DATA_cumulative_inci", "forBOLT-LMM", "BOLT_lead_variant")
  )
  Sys.setenv(OUT_DIR = file.path(base, "COJO", "EAS-only", "step1"))
}

source(engine, local = FALSE)
