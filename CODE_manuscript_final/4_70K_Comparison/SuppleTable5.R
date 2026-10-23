#!/usr/bin/env Rscript
# ============================================================
# SuppleTable5.R
#
# Supplementary Table 5. Comparison of genome-wide significant loci
# in KoGES 160K and KoGES 70K across 34 shared traits.
#
# Bidirectional +/- 1 Mb overlap is recomputed from:
#   1) KoGES 70K clumped leads
#   2) KoGES 70K genome-wide significant variants
#   3) KoGES 160K longitudinal clumped leads
#   4) KoGES 160K genome-wide significant positions
#
# Manuscript columns:
#   Significant loci in KoGES160K
#   Significant loci in KoGES70K
#   Overlapping loci          = 160K loci with a 70K GWS variant within 1 Mb
#   Novel loci in KoGES160K   = 160K loci with no 70K GWS variant within 1 Mb
#   Loci only in KoGES70K     = 70K loci with no 160K GWS variant within 1 Mb
#
# Inputs (override with env vars):
#   SUPPLETABLE5_ROOT             70K vs 160K locus-compare root
#     clump/   {PHENO}_clump_5e-8_0.1_1000.clumped
#     gws/     {PHENO}.tsv
#     sig160/  {PHENO160}.sig.tsv
#   SUPPLETABLE5_CLUMP_160K_CONT  160K BOLT clump directory
#   SUPPLETABLE5_CLUMP_160K_BIN   160K SAIGE clump directory
#
# Output:
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_cont.csv
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_cate.csv
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_only70_leads.tsv
#   RESULT_manuscript_final/SuppleTable5/SuppleTable5_only70_summary.txt
#
# Usage:
#   Rscript SuppleTable5.R
# ============================================================

suppressPackageStartupMessages(library(data.table))

WINDOW_BP <- 1e6L

###############################################################################
# Paths
###############################################################################

ROOT <- Sys.getenv(
  "SUPPLETABLE5_ROOT",
  unset = "/data6/OpenGWAS/Sumstats/KBA70K_rsid_fixed/locus_compare"
)
CLUMP_70K <- file.path(ROOT, "clump")
GWS_70K <- file.path(ROOT, "gws")
GWS_160K <- file.path(ROOT, "sig160")
CLUMP_160K_CONT <- Sys.getenv(
  "SUPPLETABLE5_CLUMP_160K_CONT",
  unset = "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping_BOLT"
)
CLUMP_160K_BIN <- Sys.getenv(
  "SUPPLETABLE5_CLUMP_160K_BIN",
  unset = "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping"
)
OUT_DIR <- Sys.getenv(
  "SUPPLETABLE5_OUT",
  unset = "/home2/sjim/openKoGES/RESULT_manuscript_final/SuppleTable5"
)
SUMSTAT_160K_DIR <- Sys.getenv(
  "SUPPLETABLE5_SUMSTAT_160K",
  unset = "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"
)
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Phenotypes whose definitions changed between the 70K and 160K analyses.
# Edit this list for the manuscript "definition updated" sentence; leave empty
# if that count should not be attributed automatically.
PHENOS_DEFINITION_UPDATED <- character(0)

###############################################################################
# Trait map (25 continuous + 9 binary shared traits)
###############################################################################

traits <- data.table(
  group = c(rep("continuous", 25L), rep("binary", 9L)),
  pheno70 = c(
    "ALB", "ALT", "AST", "BMI", "BUN", "CREATINE", "DBP", "Glucose",
    "HB", "HbA1C", "HDL", "HEIGHT", "HIP", "LDL", "PLT", "PULSE",
    "GGTP", "RBC", "SBP", "TBIL", "TC", "TG", "WAIST", "WBC", "WEIGHT",
    "BRCA", "DM", "GCA", "LIP", "HTN", "LCA", "THYCA", "THYD", "UTCA"
  ),
  pheno160 = c(
    "B_ALBUMIN", "B_ALT", "B_AST", "B_BMI", "B_BUN", "B_CREATINE",
    "B_DBP", "B_GLU0", "B_HB", "B_HBA1C", "B_HDL", "B_HEIGHT",
    "B_HIP", "B_LDLclean", "B_PLAT", "B_PULSE", "B_R_GTP", "B_RBC_B",
    "B_SBP", "B_T_BIL", "B_TCHL", "B_TG", "B_WAIST", "B_WBC_B",
    "B_WEIGHT", "B_BRCA", "B_DM2", "B_GCA", "B_LIP", "B_HTN2",
    "B_LCA", "B_THYCA", "B_THY", "B_UTCA"
  ),
  phenotype = c(
    "Albumin (Alb)", "Alanine aminotransferase (ALT)",
    "Aspartate transaminase (AST)", "Body mass index (BMI)",
    "Blood urea nitrogen (BUN)", "Serum creatinine (sCr)",
    "Diastolic blood pressure (DBP)", "Glucose", "Hemoglobin (Hb)",
    "HbA1C", "HDL-cholesterol (HDL-C)", "Body height", "Hip",
    "LDL-cholesterol (LDL-C)", "Platelet count (Plt)", "Pulse",
    "Gamma glutamyl transferase (r_GTP)", "Red blood cell count (RBC)",
    "Systolic blood pressure (SBP)", "Total bilirubin (TBil)",
    "Total cholesterol (TC)", "Triglycerides (TG)", "Waist",
    "White blood cell count (WBC)", "Body weight (BW)",
    "Breast cancer", "Diabetes mellitus", "Gastric cancer",
    "Hyperlipidemia", "Hypertension", "Lung cancer", "Thyroid cancer",
    "Thyroid disease", "Uterine cancer"
  )
)

###############################################################################
# Helpers
###############################################################################

read_clump <- function(path) {
  empty <- data.table(CHR = integer(), BP = integer(), SNP = character(), P = numeric())
  if (!file.exists(path) || file.info(path)$size == 0) return(empty)
  x <- suppressWarnings(fread(path, fill = TRUE, quote = "", showProgress = FALSE))
  required <- c("CHR", "BP", "SNP", "P")
  if (!all(required %in% names(x))) stop("Missing clump columns in ", path)
  x[, `:=`(
    CHR = suppressWarnings(as.integer(CHR)),
    BP = suppressWarnings(as.integer(BP)),
    SNP = as.character(SNP),
    P = suppressWarnings(as.numeric(P))
  )]
  unique(x[CHR %in% 1:22 & !is.na(BP) & SNP != "", ..required])
}

read_70k_gws <- function(path) {
  empty <- data.table(
    rsID = character(), CHR = integer(), POS = integer(),
    ALT = character(), REF = character(), BETA = numeric(), P = numeric()
  )
  if (!file.exists(path) || file.info(path)$size == 0) return(empty)
  x <- fread(path, showProgress = FALSE)
  need <- c("rsID", "CHR", "POS")
  if (!all(need %in% names(x))) stop("Missing 70K GWS columns in ", path)
  if (!"ALT" %in% names(x)) x[, ALT := NA_character_]
  if (!"REF" %in% names(x)) x[, REF := NA_character_]
  if (!"BETA" %in% names(x)) x[, BETA := NA_real_]
  pcol <- if ("P" %in% names(x)) "P" else stop("No P in ", path)
  x <- x[, .(
    rsID = as.character(rsID),
    CHR = suppressWarnings(as.integer(CHR)),
    POS = suppressWarnings(as.integer(POS)),
    ALT = toupper(as.character(ALT)),
    REF = toupper(as.character(REF)),
    BETA = suppressWarnings(as.numeric(BETA)),
    P = suppressWarnings(as.numeric(get(pcol)))
  )]
  unique(x[CHR %in% 1:22 & !is.na(POS) & rsID != ""])
}

read_160k_gws <- function(path) {
  empty <- data.table(CHR = integer(), POS = integer())
  if (!file.exists(path) || file.info(path)$size == 0) return(empty)
  x <- fread(
    path,
    header = FALSE,
    col.names = c("CHR", "POS"),
    showProgress = FALSE
  )
  x[, `:=`(
    CHR = suppressWarnings(as.integer(CHR)),
    POS = suppressWarnings(as.integer(POS))
  )]
  unique(x[CHR %in% 1:22 & !is.na(POS)])
}

within_window <- function(leads, positions) {
  if (!nrow(leads)) return(logical())
  answer <- rep(FALSE, nrow(leads))
  if (!nrow(positions)) return(answer)

  for (chr in unique(leads$CHR)) {
    index <- which(leads$CHR == chr)
    reference_pos <- sort(unique(positions[CHR == chr, POS]))
    if (!length(reference_pos)) next
    left <- findInterval(leads$BP[index] - WINDOW_BP - 1L, reference_pos) + 1L
    valid <- left <= length(reference_pos)
    answer[index[valid]] <- (
      reference_pos[left[valid]] <= leads$BP[index[valid]] + WINDOW_BP
    )
  }
  answer
}

clump_160k_path <- function(pheno160, group) {
  if (group == "continuous") {
    file.path(
      CLUMP_160K_CONT,
      paste0("RIN_", pheno160, "_bolt_clump_5e-8_0.1_1000.clumped")
    )
  } else {
    file.path(
      CLUMP_160K_BIN,
      paste0(
        pheno160,
        "_SAIGE_cumulative_inci.SAIGE.txt_saige_clump_5e-8_0.1_1000.clumped"
      )
    )
  }
}

sumstat_160k_path <- function(pheno160, group) {
  if (group == "continuous") {
    file.path(SUMSTAT_160K_DIR, paste0("RIN_", pheno160, "_BOLT_result"))
  } else {
    file.path(
      SUMSTAT_160K_DIR,
      paste0(pheno160, "_SAIGE_cumulative_inci.SAIGE.SEfixed.txt")
    )
  }
}

analyze_trait <- function(row) {
  message("Classifying ", row$phenotype)
  lead70 <- read_clump(file.path(
    CLUMP_70K,
    paste0(row$pheno70, "_clump_5e-8_0.1_1000.clumped")
  ))
  lead160 <- read_clump(clump_160k_path(row$pheno160, row$group))
  gws70 <- read_70k_gws(file.path(GWS_70K, paste0(row$pheno70, ".tsv")))
  sig70 <- unique(gws70[, .(CHR, POS)])
  sig160 <- read_160k_gws(file.path(GWS_160K, paste0(row$pheno160, ".sig.tsv")))

  kept70 <- within_window(lead70, sig160)
  with70 <- within_window(lead160, sig70)

  only70 <- lead70[!kept70]
  if (nrow(only70)) {
    meta <- gws70[, .(
      SNP = rsID, ALT_70K = ALT, REF_70K = REF,
      BETA_70K = BETA, P_70K_gws = P
    )]
    only70 <- merge(only70, meta, by = "SNP", all.x = TRUE)
    only70[, `:=`(
      group = row$group,
      pheno70 = row$pheno70,
      pheno160 = row$pheno160,
      phenotype = row$phenotype,
      P_70K = P
    )]
    only70 <- only70[, .(
      group, pheno70, pheno160, phenotype, CHR, SNP, BP,
      P_70K, ALT_70K, REF_70K, BETA_70K
    )]
  } else {
    only70 <- data.table(
      group = character(), pheno70 = character(), pheno160 = character(),
      phenotype = character(), CHR = integer(), SNP = character(),
      BP = integer(), P_70K = numeric(), ALT_70K = character(),
      REF_70K = character(), BETA_70K = numeric()
    )
  }

  list(
    counts = data.table(
      group = row$group,
      pheno70 = row$pheno70,
      pheno160 = row$pheno160,
      phenotype = row$phenotype,
      loci70 = nrow(lead70),
      kept70 = sum(kept70),
      only70 = sum(!kept70),
      loci160 = nrow(lead160),
      with70 = sum(with70),
      only160 = sum(!with70)
    ),
    only70 = only70
  )
}

# Allele-align 160K BETA onto the 70K effect allele (ALT_70K).
align_beta_160k <- function(beta70_allele, other70, beta160, ea160, oa160) {
  beta70_allele <- toupper(as.character(beta70_allele))
  other70 <- toupper(as.character(other70))
  ea160 <- toupper(as.character(ea160))
  oa160 <- toupper(as.character(oa160))
  if (anyNA(c(beta70_allele, other70, beta160, ea160, oa160))) return(NA_real_)
  if (identical(beta70_allele, ea160) && identical(other70, oa160)) return(as.numeric(beta160))
  if (identical(beta70_allele, oa160) && identical(other70, ea160)) return(-as.numeric(beta160))
  NA_real_
}

lookup_160k_for_trait <- function(leads) {
  if (!nrow(leads)) {
    return(cbind(leads, data.table(
      EA_160K = character(), OA_160K = character(),
      BETA_160K_raw = numeric(), BETA_160K_aligned = numeric(),
      P_160K = numeric(), found_160K = logical()
    )))
  }
  row0 <- leads[1]
  path <- sumstat_160k_path(row0$pheno160, row0$group)
  if (!file.exists(path)) {
    stop("Missing 160K sumstat: ", path)
  }
  message("Looking up 160K effects for ", row0$phenotype, " (n=", nrow(leads), ")")

  if (row0$group == "continuous") {
    hdr <- names(fread(path, nrows = 0L, showProgress = FALSE))
    pcol <- if ("P_BOLT_LMM" %in% hdr) {
      "P_BOLT_LMM"
    } else if ("P_BOLT_LMM_INF" %in% hdr) {
      "P_BOLT_LMM_INF"
    } else {
      stop("No BOLT P column in ", path)
    }
    ss <- fread(
      path,
      select = c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0", "BETA", pcol),
      showProgress = FALSE
    )
    setnames(ss, pcol, "P_160K")
    hit <- merge(
      leads,
      ss[, .(
        SNP, CHR = as.integer(CHR), BP = as.integer(BP),
        EA_160K = toupper(as.character(ALLELE1)),
        OA_160K = toupper(as.character(ALLELE0)),
        BETA_160K_raw = as.numeric(BETA),
        P_160K = as.numeric(P_160K)
      )],
      by = c("SNP", "CHR", "BP"),
      all.x = TRUE
    )
  } else {
    ss <- fread(
      path,
      select = c("CHR", "POS", "Allele1", "Allele2", "BETA", "p.value"),
      showProgress = FALSE
    )
    hit <- merge(
      leads,
      ss[, .(
        CHR = as.integer(CHR), BP = as.integer(POS),
        EA_160K = toupper(as.character(Allele2)),
        OA_160K = toupper(as.character(Allele1)),
        BETA_160K_raw = as.numeric(BETA),
        P_160K = as.numeric(p.value)
      )],
      by = c("CHR", "BP"),
      all.x = TRUE
    )
  }

  hit[, BETA_160K_aligned := mapply(
    align_beta_160k, ALT_70K, REF_70K, BETA_160K_raw, EA_160K, OA_160K
  )]
  hit[, found_160K := !is.na(P_160K)]
  hit
}

###############################################################################
# Run: locus counts (Supplementary Table 5)
###############################################################################

res <- lapply(seq_len(nrow(traits)), function(i) analyze_trait(traits[i]))
counts <- rbindlist(lapply(res, `[[`, "counts"))
only70 <- rbindlist(lapply(res, `[[`, "only70"), fill = TRUE)
stopifnot(
  all(counts$loci70 == counts$kept70 + counts$only70),
  all(counts$loci160 == counts$with70 + counts$only160),
  nrow(only70) == sum(counts$only70)
)

counts[, trait_order := match(pheno70, traits$pheno70)]
setorder(counts, trait_order)
counts[, trait_order := NULL]

format_table <- function(group_name) {
  counts[group == group_name, .(
    Phenotype = phenotype,
    `Significant loci in KoGES160K` = loci160,
    `Significant loci in KoGES70K` = loci70,
    `Overlapping loci` = with70,
    `Novel loci in KoGES160K` = only160,
    `Loci only in KoGES70K` = only70
  )]
}

cont_out <- format_table("continuous")
cate_out <- format_table("binary")
if (nrow(cont_out) != 25L) stop("Continuous traits should be 25, found ", nrow(cont_out))
if (nrow(cate_out) != 9L) stop("Binary traits should be 9, found ", nrow(cate_out))

cont_csv <- file.path(OUT_DIR, "SuppleTable5_cont.csv")
cate_csv <- file.path(OUT_DIR, "SuppleTable5_cate.csv")
fwrite(cont_out, cont_csv)
fwrite(cate_out, cate_csv)
message("Wrote ", cont_csv, " n=", nrow(cont_out))
message("Wrote ", cate_csv, " n=", nrow(cate_out))

###############################################################################
# Only-in-70K loci: 70K P thresholds and 160K effect lookup
###############################################################################

if (nrow(only70)) {
  only70 <- rbindlist(
    lapply(split(only70, by = c("group", "pheno160")), lookup_160k_for_trait),
    fill = TRUE
  )

  only70[, marginal_70K := !is.na(P_70K) & P_70K < 5e-8 & P_70K >= 1e-9]
  only70[, same_direction_160K := !is.na(BETA_70K) & !is.na(BETA_160K_aligned) &
           sign(BETA_70K) == sign(BETA_160K_aligned) & BETA_70K != 0 & BETA_160K_aligned != 0]
  only70[, nominal_160K := !is.na(P_160K) & P_160K < 0.05]
  only70[, p_lt_1e5_160K := !is.na(P_160K) & P_160K < 1e-5]
  only70[, definition_updated := pheno70 %in% PHENOS_DEFINITION_UPDATED]

  only70[, trait_order := match(pheno70, traits$pheno70)]
  setorder(only70, trait_order, CHR, BP)
  only70[, trait_order := NULL]

  n_only <- nrow(only70)
  n_marginal <- sum(only70$marginal_70K, na.rm = TRUE)
  n_found <- sum(only70$found_160K, na.rm = TRUE)
  n_same_dir <- sum(only70$same_direction_160K, na.rm = TRUE)
  n_nominal <- sum(only70$nominal_160K, na.rm = TRUE)
  n_lt_1e5 <- sum(only70$p_lt_1e5_160K, na.rm = TRUE)
  n_def <- sum(only70$definition_updated, na.rm = TRUE)
  def_phenos <- only70[definition_updated == TRUE, unique(phenotype)]

  leads_tsv <- file.path(OUT_DIR, "SuppleTable5_only70_leads.tsv")
  summary_txt <- file.path(OUT_DIR, "SuppleTable5_only70_summary.txt")
  fwrite(only70, leads_tsv, sep = "\t")

  summary_lines <- c(
    "KoGES 70K loci not recovered in KoGES 160K (Supplementary Table 5)",
    paste0("Total only-in-70K loci: ", n_only),
    paste0(
      "Marginal 70K GWS (5e-8 > P >= 1e-9): ", n_marginal, " / ", n_only
    ),
    paste0("Found in 160K sumstats: ", n_found, " / ", n_only),
    paste0(
      "Same direction of effect in 160K (allele-aligned): ",
      n_same_dir, " / ", n_found
    ),
    paste0("Nominally significant in 160K (P < 0.05): ", n_nominal, " / ", n_found),
    paste0("P < 1e-5 in 160K: ", n_lt_1e5, " / ", n_found),
    paste0(
      "Loci in phenotypes with updated definitions: ", n_def, " / ", n_only,
      if (length(def_phenos)) {
        paste0(" (", paste(def_phenos, collapse = ", "), ")")
      } else if (!length(PHENOS_DEFINITION_UPDATED)) {
        " [PHENOS_DEFINITION_UPDATED is empty; set in SuppleTable5.R if needed]"
      } else {
        ""
      }
    ),
    "",
    "Per-phenotype only-in-70K counts:",
    paste(
      only70[, .N, by = .(phenotype)][order(-N), sprintf("  %s: %d", phenotype, N)],
      collapse = "\n"
    )
  )
  writeLines(summary_lines, summary_txt)

  message("Wrote ", leads_tsv, " n=", n_only)
  message("Wrote ", summary_txt)
  message(paste(summary_lines[1:8], collapse = "\n"))
} else {
  message("No only-in-70K loci; skipping 160K lookup.")
}
