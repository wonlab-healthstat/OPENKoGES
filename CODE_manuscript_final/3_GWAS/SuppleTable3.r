#!/usr/bin/env Rscript
# ============================================================
# SuppleTable3.r
#
# Build Supplementary Table 3 (SNP heritability summaries) from
# existing LDSC h2 logs produced by heritability.r.
#
# Inputs:
#   /data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_260923/
#   /data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_binary_260924/
#
# Outputs:
#   RESULT_manuscript_final/SuppleTable3/SuppleTable3A.csv  (continuous)
#   RESULT_manuscript_final/SuppleTable3/SuppleTable3B.csv  (binary)
#
# Usage:
#   Rscript SuppleTable3.r
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

###############################################################################
# Paths
###############################################################################

CONT <- "/data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_260923"
BIN  <- "/data4/joh/openkoges/table_new0923/intermediate/suppletable3/ldsc_binary_260924"
OUT_TABLE_DIR <- "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable3"
dir.create(OUT_TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

OUT_CONT_CSV <- file.path(OUT_TABLE_DIR, "SuppleTable3_cont.csv")
OUT_BIN_CSV  <- file.path(OUT_TABLE_DIR, "SuppleTable3_cate.csv")

###############################################################################
# Helpers
###############################################################################

parse_est_se <- function(line) {
  if (is.na(line)) return(list(est = NA_real_, se = NA_real_, note = NA_character_))
  if (grepl("^Ratio < 0", line)) {
    return(list(est = NA_real_, se = NA_real_, note = "<0"))
  }
  if (grepl("^Ratio:\\s*NA", line)) {
    return(list(est = NA_real_, se = NA_real_, note = "NA"))
  }
  m <- regexec(
    ":\\s*([-+0-9.eE]+)(?:\\s*\\(([-+0-9.eE]+)\\))?",
    line,
    perl = TRUE
  )
  v <- regmatches(line, m)[[1]]
  if (length(v) < 2) return(list(est = NA_real_, se = NA_real_, note = line))
  list(
    est = as.numeric(v[2]),
    se = if (length(v) >= 3 && nzchar(v[3])) as.numeric(v[3]) else NA_real_,
    note = NA_character_
  )
}

parse_log <- function(f) {
  if (!file.exists(f)) {
    stop("Missing LDSC log: ", f, call. = FALSE)
  }
  lines <- readLines(f, warn = FALSE)
  grab <- function(pat) {
    hit <- grep(pat, lines, value = TRUE, perl = TRUE)
    if (!length(hit)) NA_character_ else hit[1]
  }
  h2_line <- grab("Total (Observed|Liability) scale h2:")
  h2 <- parse_est_se(h2_line)
  lam <- parse_est_se(grab("^Lambda GC:"))
  mean_chi <- parse_est_se(grab("^Mean Chi\\^2:"))
  intercept <- parse_est_se(grab("^Intercept:"))
  ratio <- parse_est_se(grab("^Ratio"))
  data.table(
    h2 = h2$est, h2_se = h2$se,
    mean_chi2 = mean_chi$est,
    intercept = intercept$est, intercept_se = intercept$se,
    lambda_GC = lam$est,
    attenuation_ratio = ratio$est, attenuation_ratio_se = ratio$se,
    attenuation_ratio_note = ifelse(is.na(ratio$note), "", ratio$note)
  )
}

fmt4 <- function(x) ifelse(is.na(x), "", sprintf("%.4f", x))

fmt_pair <- function(est, se, note = NULL) {
  n <- length(est)
  if (is.null(note)) note <- rep("", n)
  note <- rep_len(as.character(note), n)
  out <- ifelse(
    is.na(est), "",
    ifelse(is.na(se), fmt4(est), paste0(fmt4(est), " (", fmt4(se), ")"))
  )
  out[note == "<0"] <- "< 0"
  out[note == "NA"] <- "NA"
  out
}

format_pub_table <- function(dt) {
  out <- data.table(
    Phenotype = dt$phenotype,
    `Observed-scale h2 (SE)` = fmt_pair(dt$h2, dt$h2_se),
    `Mean chi2` = fmt4(dt$mean_chi2),
    `Intercept (SE)` = fmt_pair(dt$intercept, dt$intercept_se),
    lambdaGC = fmt4(dt$lambda_GC),
    `Attenuation ratio (SE)` = fmt_pair(
      dt$attenuation_ratio,
      dt$attenuation_ratio_se,
      dt$attenuation_ratio_note
    )
  )
  out
}

###############################################################################
# Continuous (SuppleTable3A)
###############################################################################

cont_man_path <- file.path(CONT, "manifest.tsv")
if (!file.exists(cont_man_path)) stop("Missing continuous manifest: ", cont_man_path, call. = FALSE)

cont_man <- fread(cont_man_path)
cont <- rbindlist(lapply(seq_len(nrow(cont_man)), function(i) {
  logf <- file.path(CONT, "h2", paste0(cont_man$variable[i], ".log"))
  cbind(cont_man[i], parse_log(logf))
}), fill = TRUE)

cont_pub <- format_pub_table(cont)
fwrite(cont_pub, OUT_CONT_CSV, sep = ",")
message("Wrote: ", normalizePath(OUT_CONT_CSV, winslash = "/"),
        "  n=", nrow(cont_pub))

###############################################################################
# Binary (SuppleTable3B)
###############################################################################

bin_man_path <- file.path(BIN, "manifest.tsv")
if (!file.exists(bin_man_path)) stop("Missing binary manifest: ", bin_man_path, call. = FALSE)

bin_man <- fread(bin_man_path)
binary <- rbindlist(lapply(seq_len(nrow(bin_man)), function(i) {
  logf <- file.path(BIN, "h2", paste0(bin_man$variable[i], ".log"))
  cbind(bin_man[i], parse_log(logf))
}), fill = TRUE)
setorder(binary, phenotype)

bin_pub <- format_pub_table(binary)
fwrite(bin_pub, OUT_BIN_CSV, sep = ",")
message("Wrote: ", normalizePath(OUT_BIN_CSV, winslash = "/"),
        "  n=", nrow(bin_pub))

message("Done")
