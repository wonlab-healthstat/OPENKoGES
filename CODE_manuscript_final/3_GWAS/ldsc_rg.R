#!/usr/bin/env Rscript
# ============================================================
# ldsc_rg.R
#
# Estimate pairwise SNP-based genetic correlations (rg) within
# KoGES using LD Score regression (LDSC). All phenotype pairs
# are run separately for continuous traits and for binary traits.
#
# Continuous pairs use BOLT-LMM munged sumstats; binary pairs
# use SAIGE munged sumstats. East Asian LD scores are used as
# the reference panel.
#
# Inputs:
#   Munged continuous sumstats (BOLT_LMM)
#   Munged binary sumstats (SAIGE_binary)
#   LDSC executable and East Asian LD-score reference
#
# Outputs (existing pairwise LDSC logs; cont and bin share this dir):
#   /data4/joh/openkoges/figure_new0923/intermediate/ldsc_rg_data5/*.log
#
# Usage:
#   Rscript ldsc_rg.R
# ============================================================

suppressPackageStartupMessages({
  library(parallel)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
LDSC  <- "/home2/eksong/OpenKoGES/LDSC/tools/ldsc.py"
REF   <- "/home2/eksong/OpenKoGES/LDSC/tools/eas_ldscores/"
BOLT  <- "/data5/joh/openkoges/Heritability/munge/KoGES/BOLT_LMM"
SAIGE <- "/data5/joh/openkoges/Heritability/munge/KoGES/SAIGE_binary"

OUT_RG <- "/data4/joh/openkoges/figure_new0923/intermediate/ldsc_rg_data5"
dir.create(OUT_RG, recursive = TRUE, showWarnings = FALSE)
OUT_RG_CONT <- OUT_RG
OUT_RG_BIN  <- OUT_RG

N_WORKERS <- 28L
LUNG <- c("FEV1best", "FEV1pp", "FFRbest", "FVCbest", "FVCpp")

# ------------------------------------------------------------
# phenotype lists
# ------------------------------------------------------------
CONT_PHENOS <- c(
  "ALBUMIN","ALP","ALT","APOA","APOB","AST","BMI","BUN",
  "CPEPTIDE","CREATINE","CYSTATINC","DBP","eGFR",
  "FEV1best","FEV1pp","FFRbest","FOLATE","FVCbest","FVCpp",
  "GLU0","HB","HBA1C","HCT","HDL","HEIGHT","HIP","HOMOCYST",
  "INSUL","IRON","LDL","MCH","MCHC","MCV","PHOSPHATE","PLAT",
  "PULSE","R_GTP","RBC_B","SBP","T_BIL","TCHL","TG","URICACID",
  "VITD","WAIST","WBC_B","WEIGHT","WHR"
)

BIN_PHENOS <- c(
  "HTN2","DM2","LIP","THY","ULCER","LIV","FLIV","STOMUL",
  "THYCA","NOI","MI","UTCA","CLD","GCA","COLCA","LCA","GT",
  "ANG","BRCA","CAD","OSTE","TB","COPD2","FRAC","GLAU","CVA",
  "BRON","STR","ASTH","CKD","KD","GB","UB","UT","ALLER","CATA",
  "GASTRO","ARTH","PER","POL","BPH"
)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1"
)

cont_file <- function(ph) {
  if (ph %in% LUNG) {
    return(file.path(BOLT, sprintf("B_%s_noNC.sumstats.gz", ph)))
  }
  if (ph == "LDL") {
    return(file.path(BOLT, "RIN_B_LDLclean.sumstats.gz"))
  }
  file.path(BOLT, sprintf("RIN_B_%s.sumstats.gz", ph))
}

bin_file <- function(ph) {
  file.path(SAIGE, sprintf("B_%s.sumstats.gz", ph))
}

log_ok <- function(path) {
  if (!file.exists(path) || file.info(path)$size < 200) return(FALSE)
  txt <- tryCatch(readLines(path, warn = FALSE), error = function(e) character())
  hit <- grep("Summary of Genetic Correlation Results", txt, fixed = TRUE)
  if (length(hit) == 0) return(FALSE)
  any(grepl("sumstats\\.gz", txt[hit[1]:length(txt)]))
}

pair_grid <- function(phenos, file_fn, kind, out_dir) {
  miss <- phenos[!file.exists(vapply(phenos, file_fn, character(1)))]
  if (length(miss)) {
    stop(kind, " missing munge: ", paste(miss, collapse = ", "))
  }
  comb <- t(combn(phenos, 2))
  data.frame(
    kind = kind,
    p1 = comb[, 1],
    p2 = comb[, 2],
    f1 = vapply(comb[, 1], file_fn, character(1)),
    f2 = vapply(comb[, 2], file_fn, character(1)),
    stem = paste(comb[, 1], comb[, 2], sep = "__"),
    out_dir = out_dir,
    stringsAsFactors = FALSE
  )
}

jobs <- rbind(
  pair_grid(CONT_PHENOS, cont_file, "cont", OUT_RG_CONT),
  pair_grid(BIN_PHENOS, bin_file, "bin", OUT_RG_BIN)
)

message("BOLT  ", BOLT)
message("SAIGE ", SAIGE)
message("OUT cont ", OUT_RG_CONT)
message("OUT bin  ", OUT_RG_BIN)
message("pairs total: ", nrow(jobs),
        "  cont ", sum(jobs$kind == "cont"),
        "  bin ", sum(jobs$kind == "bin"))
message("workers: ", N_WORKERS)

run_one <- function(i) {
  stem <- jobs$stem[i]
  outstem <- file.path(jobs$out_dir[i], stem)
  logfile <- paste0(outstem, ".log")
  if (log_ok(logfile)) return(c(status = "skip", stem = stem))

  cmd <- sprintf(
    "python %s --rg %s,%s --ref-ld-chr %s --w-ld-chr %s --out %s > /dev/null 2>&1",
    shQuote(LDSC),
    shQuote(jobs$f1[i]),
    shQuote(jobs$f2[i]),
    shQuote(REF),
    shQuote(REF),
    shQuote(outstem)
  )
  rc <- system(cmd)
  status <- if (log_ok(logfile)) "ok" else "fail"
  c(status = status, stem = stem, rc = as.character(rc))
}

res <- mclapply(seq_len(nrow(jobs)), run_one,
                mc.cores = N_WORKERS, mc.preschedule = FALSE)

status <- vapply(res, function(x) x[["status"]], character(1))
n_ok   <- sum(status == "ok")
n_skip <- sum(status == "skip")
n_fail <- sum(!status %in% c("ok", "skip"))

if (n_fail) {
  fail_stems <- vapply(res[status == "fail"], function(x) x[["stem"]], character(1))
  message("FAILED: ", paste(fail_stems, collapse = ", "))
}

message("DONE run_ok=", n_ok, " skip=", n_skip, " fail=", n_fail)
if (n_fail) quit(status = 2)
