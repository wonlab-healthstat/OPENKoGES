#!/usr/bin/env Rscript
# ============================================================
# SuppleTable4.r
#
# Supplementary Table 4. LD Score Regression statistics using the
# Korea1K LD reference panel.
#
# Columns match the manuscript table:
#   Phenotype
#   Observed scale h2 (Standard Error)
#   mean chi2
#   Intercept (Standard Error)
#   lambdaGC
#   Attenuation ratio (Standard Error)
#
# Inputs:
#   Continuous: Supplementary_Table_3a_Korea1K.xlsx (sheet Continuous)
#   Binary: Korea1K LDSC logs, one per phenotype
#
# Output:
#   RESULT_manuscript_final/SuppleTable4/SuppleTable4_cont.csv
#   RESULT_manuscript_final/SuppleTable4/SuppleTable4_cate.csv
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

CONT_XLSX <- "/data4/joh/openkoges/table_new0923/output/Supplementary_Table_3a_Korea1K.xlsx"
BIN_DIR <- "/data4/joh/openkoges/Heritability_GC/KoGES_heritability_KOREA1K_qc98_unrelated/Bin"
OUT_DIR <- Sys.getenv(
  "SUPPLETABLE4_OUT",
  unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable4"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Same display names and order as Supplementary Table 3 (binary).
BIN_PHENOTYPE <- c(
  B_LIV = "Acute liver disease",
  B_ALLER = "Allergic disease",
  B_ANG = "Angina pectoris",
  B_ARTH = "Arthritis",
  B_ASTH = "Asthma",
  B_BPH = "Benign prostatic hyperplasia",
  B_BRCA = "Breast cancer",
  B_CATA = "Cataract",
  B_CVA = "Cerebrovascular disease",
  B_BRON = "Chronic bronchitis",
  B_GASTRO = "Chronic gastritis",
  B_STOMUL = "Chronic gastritis and peptic ulcer",
  B_CKD = "Chronic kidney disease",
  B_CLD = "Chronic lung disease",
  B_COPD2 = "Chronic obstructive pulmonary disease (COPD)",
  B_COLCA = "Colorectal cancer",
  B_CAD = "Coronary artery disease",
  B_UTCA = "Cystitis",
  B_DM2 = "Diabetes mellitus",
  B_FLIV = "Fatty liver",
  B_FRAC = "Fracture",
  B_GB = "Gallbladder disease",
  B_GCA = "Gastric cancer",
  B_GLAU = "Glaucoma",
  B_GT = "Gout",
  B_LIP = "Hyperlipidemia",
  B_HTN2 = "Hypertension",
  B_POL = "Intestinal polyp",
  B_KD = "Kidney disease",
  B_LCA = "Lung cancer",
  B_MI = "Myocardial infarction",
  B_NOI = "Neurosis",
  B_OSTE = "Osteoporosis",
  B_ULCER = "Peptic ulcer",
  B_PER = "Periodontal disease",
  B_TB = "Pulmonary tuberculosis",
  B_STR = "Stroke",
  B_THYCA = "Thyroid cancer",
  B_THY = "Thyroid disease",
  B_UT = "Urinary tract infection",
  B_UB = "Uterine cancer"
)

parse_est_se <- function(line) {
  if (is.na(line)) return(list(est = NA_real_, se = NA_real_, note = NA_character_))
  if (grepl("^Ratio < 0", line)) {
    return(list(est = NA_real_, se = NA_real_, note = "<0"))
  }
  if (grepl("^Ratio:\\s*NA", line)) {
    return(list(est = NA_real_, se = NA_real_, note = "NA"))
  }
  m <- regexec(":\\s*([-+0-9.eE]+)(?:\\s*\\(([-+0-9.eE]+)\\))?", line, perl = TRUE)
  v <- regmatches(line, m)[[1]]
  if (length(v) < 2) return(list(est = NA_real_, se = NA_real_, note = line))
  list(
    est = as.numeric(v[2]),
    se = if (length(v) >= 3 && nzchar(v[3])) as.numeric(v[3]) else NA_real_,
    note = NA_character_
  )
}

parse_log <- function(f) {
  if (!file.exists(f)) stop("Missing LDSC log: ", f, call. = FALSE)
  lines <- readLines(f, warn = FALSE)
  grab <- function(pat) {
    hit <- grep(pat, lines, value = TRUE, perl = TRUE)
    if (!length(hit)) NA_character_ else hit[1]
  }
  h2 <- parse_est_se(grab("Total Observed scale h2:"))
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
  data.table(
    Phenotype = dt$Phenotype,
    `Observed scale h2 (Standard Error)` = fmt_pair(dt$h2, dt$h2_se),
    `mean chi2` = fmt4(dt$mean_chi2),
    `Intercept (Standard Error)` = fmt_pair(dt$intercept, dt$intercept_se),
    lambdaGC = fmt4(dt$lambda_GC),
    `Attenuation ratio (Standard Error)` = fmt_pair(
      dt$attenuation_ratio, dt$attenuation_ratio_se, dt$attenuation_ratio_note
    )
  )
}

cont <- as.data.table(read.xlsx(CONT_XLSX, sheet = "Continuous"))
cont[, Phenotype := trimws(gsub("\\s*\\*$", "", Phenotype))]
need <- c("Phenotype", "h2", "h2_se", "mean_chi2", "intercept", "intercept_se",
          "lambda_GC", "attenuation_ratio", "attenuation_ratio_se")
missing <- setdiff(need, names(cont))
if (length(missing)) stop("Continuous Korea1K table is missing: ", paste(missing, collapse = ", "))
if (!"attenuation_ratio_note" %in% names(cont)) cont[, attenuation_ratio_note := ""]
cont_pub <- format_pub_table(cont)
cont_csv <- file.path(OUT_DIR, "SuppleTable4_cont.csv")
fwrite(cont_pub, cont_csv)
message("Wrote ", cont_csv, " n=", nrow(cont_pub))

bin <- rbindlist(lapply(names(BIN_PHENOTYPE), function(code) {
  cbind(Phenotype = unname(BIN_PHENOTYPE[[code]]), parse_log(file.path(BIN_DIR, paste0(code, ".log"))))
}))
bin_pub <- format_pub_table(bin)
cate_csv <- file.path(OUT_DIR, "SuppleTable4_cate.csv")
fwrite(bin_pub, cate_csv)
message("Wrote ", cate_csv, " n=", nrow(bin_pub))
