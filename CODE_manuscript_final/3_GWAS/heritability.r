#!/usr/bin/env Rscript
# ============================================================
# heritability.r
#
# Estimate SNP-based heritability (h2) with LD Score regression
# for KoGES continuous and binary traits.
#
# Continuous traits: BOLT-LMM sumstats (munge + --h2).
# Binary traits: SAIGE cumulative-incidence SEfixed sumstats
# (munge + --h2 on the observed scale).
#
# Inputs:
#   Continuous trait list / sumstat paths (Excel matching sheet)
#   SAIGE binary GWAS sumstats directory
#   LDSC tools and East Asian LD-score reference
#
# Outputs (intermediate LDSC logs / munged sumstats):
#   /data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_260923/
#   /data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_binary_260924/
#
# Summary tables: Rscript SuppleTable3.r
#
# Usage:
#   Rscript heritability.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(parallel)
  library(readxl)
})

###############################################################################
# Paths (self-contained; no external config)
###############################################################################

CONT <- "/data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_260923"
BIN  <- "/data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_binary_260924"
dir.create(file.path(CONT, "munge"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONT, "h2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONT, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(CONT, "input"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BIN, "munge"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BIN, "h2"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(BIN, "logs"), recursive = TRUE, showWarnings = FALSE)

LDSC <- "/home2/eksong/OpenKoGES/LDSC/tools"
HM3  <- file.path(LDSC, "w_hm3.snplist")
EAS  <- file.path(LDSC, "eas_ldscores")
Sys.setenv(PYTHONPATH = LDSC)

XLSX <- "/data4/sjim/openKoGES_revision/sumstats/EUR/KoGES_EUR_EAS_matching_sjim_260913.xlsx"
SAIGE_DIR <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"
BOLT_LDL_SRC <- "/data6/OpenGWAS/KBA160K/Result/BOLT-LMM/RIN_B_LDLclean_BOLT_result"

bin_label <- c(
  HTN2 = "Hypertension", DM2 = "Diabetes mellitus", LIP = "Hyperlipidemia",
  THY = "Thyroid disease", GT = "Gout", OSTE = "Osteoporosis",
  BPH = "Benign prostatic hyperplasia", MI = "Myocardial infarction",
  ANG = "Angina pectoris", CAD = "Coronary artery disease",
  CVA = "Cerebrovascular disease", STR = "Stroke",
  CLD = "Chronic lung disease", TB = "Pulmonary tuberculosis",
  COPD2 = "Chronic obstructive pulmonary disease (COPD)",
  BRON = "Chronic bronchitis", ASTH = "Asthma", ULCER = "Peptic ulcer",
  LIV = "Acute liver disease", FLIV = "Fatty liver",
  STOMUL = "Chronic gastritis and peptic ulcer", GASTRO = "Chronic gastritis",
  POL = "Intestinal polyp", GB = "Gallbladder disease",
  CKD = "Chronic kidney disease", KD = "Kidney disease",
  UB = "Cystitis", UT = "Urinary tract infection",
  THYCA = "Thyroid cancer", UTCA = "Uterine cancer",
  GCA = "Gastric cancer", COLCA = "Colorectal cancer", LCA = "Lung cancer",
  BRCA = "Breast cancer", NOI = "Neurosis", FRAC = "Fracture",
  GLAU = "Glaucoma", ALLER = "Allergic disease", CATA = "Cataract",
  ARTH = "Arthritis", PER = "Periodontal disease"
)

bolt_variable <- function(path) {
  var <- basename(path)
  var <- sub("^RIN_", "", var)
  sub("_BOLT_result$", "", var)
}

run_cmd <- function(cmd, log, append = FALSE) {
  redir <- if (append) ">>" else ">"
  system(paste(cmd, redir, shQuote(log), "2>&1"))
}

munge_ldsc <- function(sumstats, out_prefix, h2_prefix, log, munge_args, append_log = FALSE) {
  munge_cmd <- paste(
    "python", shQuote(file.path(LDSC, "munge_sumstats.py")),
    "--sumstats", shQuote(sumstats),
    "--out", shQuote(out_prefix),
    munge_args,
    "--merge-alleles", shQuote(HM3)
  )
  st <- run_cmd(munge_cmd, log, append = append_log)
  cat("\nMUNGE_EXIT", st, "\n", file = log, append = TRUE)
  if (st != 0) return("munge_fail")

  ldsc_cmd <- paste(
    "python", shQuote(file.path(LDSC, "ldsc.py")),
    "--h2", shQuote(paste0(out_prefix, ".sumstats.gz")),
    "--ref-ld-chr", shQuote(paste0(EAS, "/")),
    "--w-ld-chr", shQuote(paste0(EAS, "/")),
    "--out", shQuote(h2_prefix)
  )
  st2 <- run_cmd(ldsc_cmd, log, append = TRUE)
  cat("\nLDSC_EXIT", st2, "\n", file = log, append = TRUE)
  if (st2 != 0) return("ldsc_fail")
  "ok"
}

# LDL P column can be mixed-type in munge; rewrite P in scientific notation.
clean_ldl <- function() {
  if (!file.exists(BOLT_LDL_SRC)) {
    stop("LDL BOLT sumstat not found: ", BOLT_LDL_SRC, call. = FALSE)
  }
  out <- file.path(CONT, "input", "B_LDLclean.txt")
  dt <- fread(
    BOLT_LDL_SRC,
    select = c("SNP", "ALLELE1", "ALLELE0", "BETA", "P_BOLT_LMM_INF")
  )
  dt[, P := sprintf("%.8e", as.numeric(P_BOLT_LMM_INF))]
  dt[, P_BOLT_LMM_INF := NULL]
  fwrite(dt, out, sep = "\t")
  out
}

run_continuous <- function(cores = 8) {
  if (!file.exists(XLSX)) stop("Continuous matching Excel not found: ", XLSX, call. = FALSE)
  raw <- as.data.table(read_excel(XLSX, sheet = "Continuous"))
  jobs <- raw[!is.na(KoGES...5), .(
    phenotype = KoGES...1,
    N = as.integer(KoGES_N),
    path = KoGES...5
  )]
  jobs[, variable := bolt_variable(path)]
  jobs <- unique(jobs, by = "variable")
  fwrite(jobs, file.path(CONT, "manifest.tsv"), sep = "\t")

  ldl_file <- NULL
  if ("B_LDLclean" %in% jobs$variable) {
    message("Cleaning LDL P column")
    ldl_file <- clean_ldl()
  }

  run_one <- function(i) {
    var <- jobs$variable[i]
    log <- file.path(CONT, "logs", paste0(var, ".run.log"))
    cat("START", var, "N=", jobs$N[i], "\n", jobs$path[i], "\n", file = log)
    if (var == "B_LDLclean") {
      args <- paste(
        "--snp SNP",
        "--N", jobs$N[i],
        "--a1 ALLELE1 --a2 ALLELE0",
        "--signed-sumstats BETA,0",
        "--p P"
      )
      sumstats <- ldl_file
    } else {
      args <- paste(
        "--snp SNP",
        "--N", jobs$N[i],
        "--a1 ALLELE1 --a2 ALLELE0",
        "--signed-sumstats BETA,0",
        "--p P_BOLT_LMM_INF"
      )
      sumstats <- jobs$path[i]
    }
    munge_ldsc(
      sumstats,
      file.path(CONT, "munge", var),
      file.path(CONT, "h2", var),
      log,
      args,
      append_log = TRUE
    )
  }

  status <- unlist(mclapply(seq_len(nrow(jobs)), run_one, mc.cores = cores))
  data.table(variable = jobs$variable, status = status)
}

run_binary <- function(cores = 6) {
  if (!dir.exists(SAIGE_DIR)) stop("SAIGE sumstat directory not found: ", SAIGE_DIR, call. = FALSE)
  files <- list.files(
    SAIGE_DIR,
    pattern = "_SAIGE_cumulative_inci\\.SAIGE\\.SEfixed\\.txt$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    stop("No SAIGE SEfixed sumstats found under: ", SAIGE_DIR, call. = FALSE)
  }
  jobs <- rbindlist(lapply(files, function(p) {
    var <- sub("_SAIGE_cumulative_inci\\.SAIGE\\.SEfixed\\.txt$", "", basename(p))
    code <- sub("^B_", "", var)
    head_dt <- fread(p, nrows = 1, select = c("N_case", "N_ctrl"))
    data.table(
      phenotype = unname(bin_label[code]),
      variable = var,
      N_cas = head_dt$N_case[1],
      N_con = head_dt$N_ctrl[1],
      path = p
    )
  }))
  if (any(is.na(jobs$phenotype))) {
    stop(
      "Unlabeled binary trait: ",
      paste(jobs$variable[is.na(jobs$phenotype)], collapse = ", "),
      call. = FALSE
    )
  }
  fwrite(jobs, file.path(BIN, "manifest.tsv"), sep = "\t")

  run_one <- function(i) {
    var <- jobs$variable[i]
    log <- file.path(BIN, "logs", paste0(var, ".run.log"))
    cat(
      "START", var, "cases=", jobs$N_cas[i], "controls=", jobs$N_con[i], "\n",
      jobs$path[i], "\n",
      file = log
    )
    args <- paste(
      "--snp MarkerID",
      "--N-cas", jobs$N_cas[i],
      "--N-con", jobs$N_con[i],
      "--a1 Allele2 --a2 Allele1",
      "--signed-sumstats BETA,0",
      "--p p.value"
    )
    munge_ldsc(
      jobs$path[i],
      file.path(BIN, "munge", var),
      file.path(BIN, "h2", var),
      log,
      args,
      append_log = TRUE
    )
  }

  status <- unlist(mclapply(seq_len(nrow(jobs)), run_one, mc.cores = cores))
  data.table(variable = jobs$variable, status = status)
}

###############################################################################
# Main
###############################################################################

message("Continuous LDSC")
cont_status <- run_continuous(8)
print(cont_status[status != "ok"])

message("Binary LDSC")
bin_status <- run_binary(6)
print(bin_status[status != "ok"])

message("Done (tables: Rscript SuppleTable3.r)")
