#!/usr/bin/env Rscript
# BOLTLMM.R
#
# Purpose:
#   Genome-wide association analysis of continuous traits in KoGES using BOLT-LMM (--lmm). 
#
# Inputs (defaults; override with CLI flags):
#   --phenoFile  phenotype table (must contain FID, IID, and pheno columns)
#   --covarFile  covariate table (FID, IID, SEX, AGE terms, PCs)
#   --bfile      PLINK prefix
#   --modelSnps  SNP list for the mixed-model kinship
#
# Outputs:
#   --outDir/<PHENO>_BOLT_result


library(data.table)
library(parallel)

# -----------------------------------------------------------------------------
# Default paths (KoGES)
# -----------------------------------------------------------------------------
pheno_file <- "/data6/OpenGWAS/KBA160K/Phenotype/KBA_pheno_clin_RIN_152719_LDLclean_Ca.txt"
covar_file <- "/data6/OpenGWAS/KBA160K/Phenotype/KBA_pheno_covar_152719.txt"
bfile <- "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated"
model_snps <- "/data6/FROM_MR/eksong/OpenKoGES/BOLT-LMM/pruning/KBA160K_pruning_0.01_50_5_0.1.prune.in"
genetic_map <- "/home2/ymji/4.Tools/BOLT-LMM_v2.4/tables/genetic_map_hg19_withX.txt.gz"
ldscores <- "/home2/ymji/4.Tools/BOLT-LMM_v2.4/tables/LDSCORE.1000G_EAS.tab.gz"
out.dir <- "/data6/OpenGWAS/KBA160K/Result/BOLT-LMM"

# Default covariates (categorical + quantitative)
covar_cols  <- c("B_SEX_BINARY")
qcovar_cols <- c("B_AGE", "B_AGE_SQ", "B_AGE_SEX", "B_AGE_SQ_SEX", paste0("PC", 1:20))

# Default phenotype list when none are given on the CLI
pheno_list_default <- c("RIN_B_HOMOCYST")

# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------
# --njob / --threads / --phenoFile / --covarFile / --bfile / --modelSnps
# --outDir / --pheno-list=A,B,C
# remaining positional args = phenotype column names

args <- commandArgs(trailingOnly = TRUE)
njob <- 5L
bolt_threads <- 5L
pheno_list <- character(0)

parse_flag_eq <- function(args, name) {
  idx <- grep(paste0("^--", name, "="), args)
  if (length(idx) == 0) return(list(args = args, value = NULL))
  list(args = args[-idx[1]], value = sub(paste0("^--", name, "="), "", args[idx[1]]))
}

parse_flag_space <- function(args, name) {
  idx <- which(args == paste0("--", name))
  if (length(idx) == 0) return(list(args = args, value = NULL))
  if (idx[1] == length(args)) stop("--", name, " requires a value")
  list(args = args[-c(idx[1], idx[1] + 1)], value = args[idx[1] + 1])
}

parse_opt <- function(args, name) {
  p <- parse_flag_eq(args, name)
  if (!is.null(p$value)) return(p)
  parse_flag_space(args, name)
}

if (length(args) > 0) {
  p <- parse_opt(args, "njob"); args <- p$args
  if (!is.null(p$value)) {
    v <- suppressWarnings(as.integer(p$value))
    if (is.na(v) || v < 1) stop("Invalid --njob: ", p$value)
    njob <- v
  }

  p <- parse_opt(args, "threads"); args <- p$args
  if (!is.null(p$value)) {
    v <- suppressWarnings(as.integer(p$value))
    if (is.na(v) || v < 1) stop("Invalid --threads: ", p$value)
    bolt_threads <- v
  }

  p <- parse_opt(args, "phenoFile"); args <- p$args
  if (!is.null(p$value)) pheno_file <- p$value

  p <- parse_opt(args, "covarFile"); args <- p$args
  if (!is.null(p$value)) covar_file <- p$value

  p <- parse_opt(args, "bfile"); args <- p$args
  if (!is.null(p$value)) bfile <- p$value

  p <- parse_opt(args, "modelSnps"); args <- p$args
  if (!is.null(p$value)) model_snps <- p$value

  p <- parse_opt(args, "outDir"); args <- p$args
  if (!is.null(p$value)) out.dir <- p$value

  p <- parse_opt(args, "pheno-list"); args <- p$args
  if (!is.null(p$value)) {
    pheno_list <- unique(trimws(unlist(strsplit(p$value, ",", fixed = TRUE))))
    pheno_list <- pheno_list[nzchar(pheno_list)]
  }

  if (length(args) > 0) {
    extra <- unique(trimws(args))
    extra <- extra[nzchar(extra)]
    pheno_list <- unique(c(pheno_list, extra))
  }
}

if (length(pheno_list) == 0) pheno_list <- pheno_list_default

dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Validate inputs
# -----------------------------------------------------------------------------
if (!file.exists(pheno_file)) stop("Phenotype file not found: ", pheno_file)
if (!file.exists(covar_file)) stop("Covariate file not found: ", covar_file)
if (!file.exists(paste0(bfile, ".bed"))) stop("PLINK bed not found: ", bfile, ".bed")
if (!file.exists(model_snps)) stop("modelSnps file not found: ", model_snps)

pheno_hdr <- names(fread(pheno_file, nrows = 0))
need_id <- c("FID", "IID")
missing_id <- setdiff(need_id, pheno_hdr)
if (length(missing_id) > 0) {
  stop("Phenotype file missing columns: ", paste(missing_id, collapse = ", "))
}
missing_pheno <- setdiff(pheno_list, pheno_hdr)
if (length(missing_pheno) > 0) {
  stop(
    "Requested phenotype column(s) not in phenoFile:\n  ",
    paste(missing_pheno, collapse = ", ")
  )
}

cat("phenoFile :", pheno_file, "\n")
cat("covarFile :", covar_file, "\n")
cat("outDir    :", out.dir, "\n")
cat("njob      :", njob, " threads/BOLT:", bolt_threads, "\n")
cat("phenotypes:", paste(pheno_list, collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# BOLT-LMM (one trait)
# -----------------------------------------------------------------------------
run_bolt <- function(pheno) {
  out_file <- file.path(out.dir, paste0(pheno, "_BOLT_result"))
  cat("\n[BOLT]", pheno, "->", out_file, "\n")

  cmd <- paste(
    "bolt",
    "--bfile", bfile,
    "--modelSnps", model_snps,
    "--phenoFile", pheno_file,
    paste0("--phenoCol=", pheno),
    "--covarFile", covar_file,
    paste("--covarCol", covar_cols, collapse = " "),
    paste("--qCovarCol", qcovar_cols, collapse = " "),
    "--geneticMapFile", genetic_map,
    "--LDscoresFile", ldscores,
    "--lmm",
    "--verboseStats",
    paste0("--statsFile ", out_file),
    "--numThreads", bolt_threads
  )

  status <- system(cmd)
  if (!identical(as.integer(status), 0L)) {
    cat("[ERROR] BOLT failed for", pheno, "status=", status, "\n")
    return(NULL)
  }
  cat("[DONE]", pheno, "\n")
  out_file
}

# -----------------------------------------------------------------------------
# Run
# -----------------------------------------------------------------------------
if (njob <= 1L) {
  lapply(pheno_list, run_bolt)
} else {
  mclapply(pheno_list, run_bolt, mc.cores = njob)
}

cat("[DONE]\n")
