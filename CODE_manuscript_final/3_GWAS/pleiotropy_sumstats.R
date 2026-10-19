#!/usr/bin/env Rscript
# ============================================================
# pleiotropy_sumstats.R
#
# Purpose:
#   Extract lead SNPs from each phenotype's GWAS sumstats
#   (P < 5e-6), so the locus plot can show cross-trait overlap.
#
# Adapted from:
#   /data4/joh/openkoges/figure_new0923/code/figure3/01_extract_pleiotropy_sumstats.R
#
# Inputs:
#   Clumping_BOLT / Clumping  (.clumped, tag 5e-8_0.1_1000)
#   GWAS_sumstat (BOLT continuous / SAIGE binary)
#
# Outputs:
#   DATA_manuscript_final/Figure3/
#     lead_snps_continuous.tsv
#     all_lead_rsids.txt
#     sumstat_lead_continuous_all.tsv
#   DATA_manuscript_final/SuppleFigure5/
#     lead_snps_binary.tsv
#     sumstat_lead_binary_all.tsv
#
# Usage:
#   Rscript pleiotropy_sumstats.R
#   FIGURE3_CORES=8 Rscript pleiotropy_sumstats.R
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
INT_DIR_CONT <- "/home2/sjim/openKoGES/DATA_manuscript_final/Figure3"
INT_DIR_BIN  <- "/home2/sjim/openKoGES/DATA_manuscript_final/SuppleFigure5"
dir.create(INT_DIR_CONT, recursive = TRUE, showWarnings = FALSE)
dir.create(INT_DIR_BIN,  recursive = TRUE, showWarnings = FALSE)

CLUMP_DIR_CONT  <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping_BOLT"
CLUMP_DIR_BIN   <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping"
SUMSTAT_DIR     <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"

P_THR_PLOT <- 5e-6
CLUMP_TAG  <- "5e-8_0.1_1000"

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

# ------------------------------------------------------------
# file resolvers
# ------------------------------------------------------------
first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

resolve_cont_clump <- function(ph) {
  first_existing(c(
    file.path(CLUMP_DIR_CONT, sprintf("RIN_B_%s_bolt_clump_%s.clumped", ph, CLUMP_TAG)),
    file.path(CLUMP_DIR_CONT, sprintf("RIN_B_%sclean_bolt_clump_%s.clumped", ph, CLUMP_TAG)),
    file.path(CLUMP_DIR_CONT, sprintf("B_%s_noNC_bolt_clump_%s.clumped", ph, CLUMP_TAG)),
    file.path(CLUMP_DIR_CONT, sprintf("B_%s_bolt_clump_%s.clumped", ph, CLUMP_TAG))
  ))
}

resolve_cont_sumstat <- function(ph) {
  first_existing(c(
    file.path(SUMSTAT_DIR, sprintf("RIN_B_%s_BOLT_result", ph)),
    file.path(SUMSTAT_DIR, sprintf("RIN_B_%sclean_BOLT_result", ph)),
    file.path(SUMSTAT_DIR, sprintf("RIN_B_%s_noNC_BOLT_result", ph)),
    file.path(SUMSTAT_DIR, sprintf("B_%s_noNC_BOLT_result", ph))
  ))
}

resolve_bin_clump <- function(ph) {
  first_existing(c(
    file.path(CLUMP_DIR_BIN, sprintf(
      "B_%s_SAIGE_cumulative_inci.SAIGE.txt_saige_clump_%s.clumped", ph, CLUMP_TAG
    ))
  ))
}

resolve_bin_sumstat <- function(ph) {
  first_existing(c(
    file.path(SUMSTAT_DIR, sprintf("B_%s_SAIGE_cumulative_inci.SAIGE.SEfixed.txt", ph)),
    file.path(SUMSTAT_DIR, sprintf("B_%s_SAIGE_cumulative_inci.SAIGE.txt", ph))
  ))
}

read_clump_snps <- function(path, lead_pheno) {
  if (is.na(path) || !file.exists(path) || file.info(path)$size == 0) {
    message("Skip (no clump file): ", lead_pheno)
    return(NULL)
  }
  dt <- tryCatch(fread(path), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(NULL)
  snp_col <- intersect(c("SNP", "snp", "MarkerID"), names(dt))[1]
  if (is.na(snp_col)) {
    dt <- fread(path, header = FALSE)
    setnames(dt, 1, "SNP")
    snp_col <- "SNP"
  }
  out <- unique(dt[!is.na(get(snp_col)) & get(snp_col) != "",
                   .(SNP = as.character(get(snp_col)))])
  out[, lead_pheno := lead_pheno]
  out
}

n_cores <- as.integer(Sys.getenv("FIGURE3_CORES", "12"))
n_cores <- max(1L, min(n_cores, parallel::detectCores(), 16L))

# ------------------------------------------------------------
# 1. collect lead SNPs from clumping
# ------------------------------------------------------------
lead_cont <- rbindlist(
  lapply(CONT_PHENOS, function(ph) read_clump_snps(resolve_cont_clump(ph), ph)),
  fill = TRUE
)
lead_bin <- rbindlist(
  lapply(BIN_PHENOS, function(ph) read_clump_snps(resolve_bin_clump(ph), ph)),
  fill = TRUE
)

fwrite(lead_cont, file.path(INT_DIR_CONT, "lead_snps_continuous.tsv"), sep = "\t")
fwrite(lead_bin,  file.path(INT_DIR_BIN,  "lead_snps_binary.tsv"), sep = "\t")

leads_all <- unique(c(lead_cont$SNP, lead_bin$SNP))
leads_all <- leads_all[!is.na(leads_all) & nchar(leads_all) > 0]
leads_file <- file.path(INT_DIR_CONT, "all_lead_rsids.txt")
writeLines(leads_all, leads_file)

message("Lead SNPs: continuous=", nrow(lead_cont),
        " binary=", nrow(lead_bin),
        " unique=", length(leads_all))
message("Continuous dir: ", INT_DIR_CONT)
message("Binary dir: ", INT_DIR_BIN)

# ------------------------------------------------------------
# 2. extract those SNPs from each sumstat (P < 5e-6)
# ------------------------------------------------------------
extract_one <- function(ph, sumstat, kind) {
  if (is.na(sumstat) || !file.exists(sumstat)) {
    message("Missing sumstat: ", ph)
    return(NULL)
  }
  tmp_file <- tempfile(
    pattern = sprintf("sumstat_lead_%s_%s_", kind, ph),
    tmpdir = if (kind == "cont") INT_DIR_CONT else INT_DIR_BIN
  )

  snp_field <- if (kind == "cont") 1 else 3
  cmd <- sprintf(
    "awk -F'\\t' 'NR==FNR { s[$1]=1; next } FNR==1 { print; next } ($%d in s)' '%s' '%s' > '%s'",
    snp_field, leads_file, sumstat, tmp_file
  )
  status <- system(cmd, ignore.stderr = TRUE)
  if (status != 0 || !file.exists(tmp_file) || file.info(tmp_file)$size == 0) {
    message("Extract failed: ", ph)
    return(NULL)
  }

  dt <- tryCatch(fread(tmp_file), error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) {
    unlink(tmp_file)
    return(NULL)
  }

  if (kind == "cont") {
    snp_col <- intersect(c("SNP", "snp"), names(dt))[1]
    p_col <- if ("P_BOLT_LMM_INF" %in% names(dt)) "P_BOLT_LMM_INF" else
      if ("P_BOLT_LMM" %in% names(dt)) "P_BOLT_LMM" else
        if ("P" %in% names(dt)) "P" else NA_character_
    pos_col <- intersect(c("BP", "POS", "pos"), names(dt))[1]
    chr_col <- intersect(c("CHR", "Chr"), names(dt))[1]
    if (any(is.na(c(snp_col, p_col, pos_col, chr_col)))) {
      message("Unexpected BOLT columns: ", ph, " | ", paste(names(dt), collapse = ","))
      unlink(tmp_file)
      return(NULL)
    }
    out <- dt[, .(
      SNP = as.character(get(snp_col)),
      CHR = as.integer(get(chr_col)),
      POS = as.numeric(get(pos_col)),
      P   = as.numeric(get(p_col)),
      pheno = ph
    )]
  } else {
    snp_col <- intersect(c("MarkerID", "SNP", "snp"), names(dt))[1]
    p_col <- intersect(c("p.value", "P", "p"), names(dt))[1]
    pos_col <- intersect(c("POS", "BP", "pos"), names(dt))[1]
    chr_col <- intersect(c("CHR", "Chr"), names(dt))[1]
    if (any(is.na(c(snp_col, p_col, pos_col, chr_col)))) {
      message("Unexpected SAIGE columns: ", ph)
      unlink(tmp_file)
      return(NULL)
    }
    out <- dt[, .(
      SNP = as.character(get(snp_col)),
      CHR = as.integer(get(chr_col)),
      POS = as.numeric(get(pos_col)),
      P   = as.numeric(get(p_col)),
      pheno = ph
    )]
  }

  out <- out[!is.na(P) & P < P_THR_PLOT & !is.na(SNP) & !is.na(CHR) & !is.na(POS)]
  unlink(tmp_file)
  cat(sprintf("[%s] %s  rows=%d\n", kind, ph, nrow(out)), file = stderr())
  out
}

message("Extracting continuous sumstats on ", n_cores, " cores ...")
sum_cont_list <- parallel::mclapply(CONT_PHENOS, function(ph) {
  extract_one(ph, resolve_cont_sumstat(ph), "cont")
}, mc.cores = n_cores)
sum_cont <- rbindlist(sum_cont_list, fill = TRUE)
fwrite(sum_cont, file.path(INT_DIR_CONT, "sumstat_lead_continuous_all.tsv"), sep = "\t")

message("Extracting binary sumstats on ", n_cores, " cores ...")
sum_bin_list <- parallel::mclapply(BIN_PHENOS, function(ph) {
  extract_one(ph, resolve_bin_sumstat(ph), "bin")
}, mc.cores = n_cores)
sum_bin <- rbindlist(sum_bin_list, fill = TRUE)
fwrite(sum_bin, file.path(INT_DIR_BIN, "sumstat_lead_binary_all.tsv"), sep = "\t")

message("Done. continuous rows=", nrow(sum_cont), " binary rows=", nrow(sum_bin))
message("Wrote continuous under: ", INT_DIR_CONT)
message("Wrote binary under: ", INT_DIR_BIN)
