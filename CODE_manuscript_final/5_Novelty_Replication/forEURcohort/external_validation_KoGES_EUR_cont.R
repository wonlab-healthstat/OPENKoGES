#!/usr/bin/env Rscript
# external_validation_KoGES_EUR_cont.R
#
# Purpose:
#   Assess locus-level novelty (+/-1 Mb GWS) and exact-variant replication of
#   KoGES BOLT-LMM continuous (Cont) lead variants against one European
#   external cohort (UKB, FinnGen, or MVP), selected via --cohort.
#   Phenotype map from META_EXCEL (BOLT stem -> external file stem).
#
# Inputs:
#   KOGES_LEAD_FILE (default):
#     DATA_manuscript_final/Extended_Data/Cont_lead_variants.tsv
#   EXT_DIR (or UKB_DIR / FINNGEN_DIR / MVP_DIR):
#     external continuous sumstats directory for the chosen cohort
#     (default under DATA_baseline/forEURcohort/{COHORT}/Cont when unset)
#
# Outputs (OUTPUT_DIR; default Extended_Data):
#   KoGES_BOLT_leads_{TAG}_novelty_replication_cont.tsv
#
# Usage:
#   Rscript external_validation_KoGES_EUR_cont.R --cohort=UKB
#   Rscript external_validation_KoGES_EUR_cont.R --cohort=FinnGen
#   Rscript external_validation_KoGES_EUR_cont.R --cohort=MVP
#   EXT_DIR=... OUTPUT_DIR=... Rscript external_validation_KoGES_EUR_cont.R --cohort=UKB
#

suppressPackageStartupMessages({
  library(data.table)
})

source("/home2/sjim/openKoGES/CODE_GWASsumstatQC/eur_koges_trait_map.R")

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
grab_eq <- function(flag) {
  idx <- grep(paste0("^", flag, "="), args)
  if (length(idx) == 0) return(NULL)
  val <- sub(paste0("^", flag, "="), "", args[idx[1]])
  args <<- args[-idx]
  val
}
grab_sp <- function(flag) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(NULL)
  if (idx[1] == length(args)) stop(flag, " requires a value")
  val <- args[idx[1] + 1]
  args <<- args[-c(idx[1], idx[1] + 1)]
  val
}
if (any(args %in% c("--help", "-h"))) {
  cat("Usage: Rscript external_validation_KoGES_EUR_cont.R --cohort=UKB|FinnGen|MVP\n")
  quit(save = "no", status = 0)
}
v <- grab_eq("--cohort"); if (is.null(v)) v <- grab_sp("--cohort")
if (is.null(v)) v <- Sys.getenv("COHORT", "")
COHORT <- trimws(as.character(v))
if (!COHORT %chin% c("UKB", "FinnGen", "MVP")) {
  stop("--cohort=UKB|FinnGen|MVP is required.")
}
if (length(args)) stop("Unknown args: ", paste(args, collapse = " "))

TAG <- COHORT
pfx <- tolower(COHORT)
RSID_COL <- switch(COHORT, UKB = "rsID_ukb", FinnGen = "rsID_finngen", MVP = "rsID_mvp")

KOGES_LEAD_FILE <- Sys.getenv(
  "KOGES_LEAD_FILE",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/Cont_lead_variants.tsv"
)

EXT_DIR <- Sys.getenv("EXT_DIR", "")
if (!nzchar(EXT_DIR)) {
  alt <- switch(
    COHORT,
    UKB = Sys.getenv("UKB_DIR", ""),
    FinnGen = Sys.getenv("FINNGEN_DIR", ""),
    MVP = Sys.getenv("MVP_DIR", "")
  )
  EXT_DIR <- if (nzchar(alt)) {
    alt
  } else {
    file.path(
      "/data4/sjim/openKoGES_revision/DATA_baseline/forEURcohort",
      COHORT, "Cont"
    )
  }
}

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
)

NOVELTY_WINDOW_BP <- as.integer(Sys.getenv("NOVELTY_WINDOW_BP", "1000000"))
GWS_P <- as.numeric(Sys.getenv("GWS_P", "5e-8"))
NOMINAL_REPLICATION_P <- as.numeric(
  Sys.getenv("NOMINAL_REPLICATION_P", "0.05")
)

# Available primary corrections:
#   global_mapped   : 0.05 / all KoGES lead variants with a mapped external phenotype
#   global_testable : 0.05 / exact lead variants testable externally
#   trait_mapped    : 0.05 / mapped KoGES lead variants within each phenotype
#   trait_testable  : 0.05 / exact testable lead variants within each phenotype
# Primary: global_testable = 0.05 / N exact leads successfully tested.
# global_mapped is a more conservative sensitivity option via PRIMARY_CORRECTION.
PRIMARY_CORRECTION <- Sys.getenv("PRIMARY_CORRECTION", "global_testable")
VALID_CORRECTIONS <- c(
  "global_mapped", "global_testable", "trait_mapped", "trait_testable"
)
if (!PRIMARY_CORRECTION %chin% VALID_CORRECTIONS) {
  stop(
    "Invalid PRIMARY_CORRECTION: ", PRIMARY_CORRECTION,
    ". Allowed: ", paste(VALID_CORRECTIONS, collapse = ", ")
  )
}

if (!is.finite(NOVELTY_WINDOW_BP) || NOVELTY_WINDOW_BP < 0L) {
  stop("NOVELTY_WINDOW_BP must be a non-negative integer.")
}
if (!is.finite(GWS_P) || GWS_P <= 0 || GWS_P > 1) {
  stop("GWS_P must be in (0,1].")
}
if (
  !is.finite(NOMINAL_REPLICATION_P) ||
  NOMINAL_REPLICATION_P <= 0 ||
  NOMINAL_REPLICATION_P > 1
) {
  stop("NOMINAL_REPLICATION_P must be in (0,1].")
}

RESULT_FILE <- file.path(
  OUTPUT_DIR,
  paste0("KoGES_BOLT_leads_", TAG, "_novelty_replication_cont.tsv")
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(KOGES_LEAD_FILE)) {
  stop("KoGES lead file not found: ", KOGES_LEAD_FILE)
}
if (!dir.exists(EXT_DIR)) {
  stop(COHORT, " directory not found: ", EXT_DIR)
}

# -----------------------------------------------------------------------------
# Phenotype map (BOLT stem -> external stem)
# -----------------------------------------------------------------------------
META_EXCEL <- META_EXCEL_DEFAULT
SHEET_NAME <- META_SHEET_DEFAULT
bolt_to_ukb_stem <- build_bolt_to_stem_from_excel(
  cohort = COHORT,
  meta_excel = META_EXCEL,
  sheet = SHEET_NAME
)
cat("[", COHORT, "] bolt_to_stem from META_EXCEL: ", length(bolt_to_ukb_stem),
    " phenotypes\n", sep = "")

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------
pick_column <- function(nms, candidates, label, required = TRUE) {
  hit <- candidates[candidates %chin% nms]
  if (length(hit) == 0L) {
    if (required) {
      stop(
        "Could not find ", label, " column. Tried: ",
        paste(candidates, collapse = ", ")
      )
    }
    return(NA_character_)
  }
  hit[[1L]]
}

norm_chr <- function(x) {
  x <- trimws(as.character(x))
  x <- sub("^chr", "", x, ignore.case = TRUE)
  suppressWarnings(as.integer(x))
}

norm_pos <- function(x) {
  suppressWarnings(as.integer(as.numeric(x)))
}

norm_allele <- function(x) {
  toupper(trimws(as.character(x)))
}

as_flag <- function(x) {
  toupper(trimws(as.character(x))) %chin% c("1", "TRUE", "T", "YES", "Y")
}

complement_allele <- function(x) {
  chartr("ACGT", "TGCA", x)
}

is_snv_pair <- function(a1, a2) {
  !is.na(a1) && !is.na(a2) &&
    nchar(a1) == 1L && nchar(a2) == 1L &&
    grepl("^[ACGT]$", a1) && grepl("^[ACGT]$", a2)
}

safe_sign_concordance <- function(beta1, beta2) {
  if (
    !is.finite(beta1) || !is.finite(beta2) ||
    beta1 == 0 || beta2 == 0
  ) {
    return(NA)
  }
  sign(beta1) == sign(beta2)
}

# Convert KoGES BOLT effect (EA vs NEA) to KoGES BIM A1 effect.
# Complement-based matching is allowed only for single-nucleotide variants.
align_koges_beta_to_a1 <- function(ea, nea, beta, a1, a2) {
  ea <- norm_allele(ea)
  nea <- norm_allele(nea)
  a1 <- norm_allele(a1)
  a2 <- norm_allele(a2)
  beta <- suppressWarnings(as.numeric(beta))
  
  fail <- list(
    ok = FALSE,
    beta_a1 = NA_real_,
    match_type = "unmatched"
  )
  
  if (
    any(is.na(c(ea, nea, a1, a2))) ||
    any(!nzchar(c(ea, nea, a1, a2))) ||
    !is.finite(beta)
  ) {
    fail$match_type <- "missing_allele_or_beta"
    return(fail)
  }
  
  if (ea == a1 && nea == a2) {
    return(list(ok = TRUE, beta_a1 = beta, match_type = "perfect"))
  }
  if (ea == a2 && nea == a1) {
    return(list(ok = TRUE, beta_a1 = -beta, match_type = "reverse"))
  }
  
  if (is_snv_pair(ea, nea) && is_snv_pair(a1, a2)) {
    ea_c <- complement_allele(ea)
    nea_c <- complement_allele(nea)
    
    if (ea_c == a1 && nea_c == a2) {
      return(list(ok = TRUE, beta_a1 = beta, match_type = "flip"))
    }
    if (ea_c == a2 && nea_c == a1) {
      return(list(ok = TRUE, beta_a1 = -beta, match_type = "reverse_flip"))
    }
  }
  
  fail
}

empty_gws_hits <- function() {
  data.table(
    lead_row_id = integer(),
    koges_phenotype = character(),
    koges_pheno_code = character(),
    koges_lead_variant = character(),
    koges_chr = integer(),
    koges_bp = integer(),
    ukb_output_stem = character(),
    ukb_chr = integer(),
    ukb_pos = integer(),
    distance_from_koges_lead_bp = integer(),
    ukb_snp = character(),
    ukb_rsid_koges = character(),
    ukb_effect_allele_orig = character(),
    ukb_other_allele_orig = character(),
    ukb_beta_orig = numeric(),
    ukb_se = numeric(),
    ukb_p = numeric(),
    ukb_eaf_orig = numeric(),
    ukb_info = numeric(),
    ukb_variant_type = character(),
    ukb_match_type = character(),
    ukb_usable_for_replication = logical()
  )
}

read_ext_cont <- function(path) {
  header <- names(fread(path, nrows = 0L))
  required <- c(
    "CHR", "POS", RSID_COL, "rsID_koges", "A1_koges", "A2_koges",
    "effect_allele_orig", "other_allele_orig", "BETA_orig",
    "BETA_aligned", "SE", "P", "effect_allele_frequency_orig",
    "effect_allele_frequency_aligned", "effect_allele_aligned",
    "other_allele_aligned", "variant_type", "match_type",
    "usable_for_replication", "usable_for_novelty"
  )
  missing <- setdiff(required, header)
  if (length(missing) > 0L) {
    stop(
      COHORT, " file is missing columns (", basename(path), "): ",
      paste(missing, collapse = ", ")
    )
  }

  optional <- c(
    "source_panel", "trait", "output_stem", "BP_GRCh38",
    "variantID_koges", "ref_koges", "alt_koges",
    "EAF_A1_koges", "MAF", "NCHROBS_koges",
    "is_indel", "is_palindromic"
  )
  dt <- fread(path, select = unique(c(required, intersect(optional, header))))

  # Keep SNP-like column for downstream reporting compatibility.
  dt[, SNP := as.character(get(RSID_COL))]

  dt[, `:=`(
    CHR = norm_chr(CHR),
    POS = norm_pos(POS),
    P = suppressWarnings(as.numeric(P)),
    BETA_orig = suppressWarnings(as.numeric(BETA_orig)),
    BETA_aligned = suppressWarnings(as.numeric(BETA_aligned)),
    SE = suppressWarnings(as.numeric(SE)),
    usable_for_replication = as_flag(usable_for_replication),
    usable_for_novelty = as_flag(usable_for_novelty),
    A1_koges = norm_allele(A1_koges),
    A2_koges = norm_allele(A2_koges),
    effect_allele_orig = norm_allele(effect_allele_orig),
    other_allele_orig = norm_allele(other_allele_orig),
    effect_allele_aligned = norm_allele(effect_allele_aligned),
    other_allele_aligned = norm_allele(other_allele_aligned)
  )]

  if (anyNA(dt[, .(CHR, POS)])) {
    stop(COHORT, " file contains missing/invalid CHR or POS: ", path)
  }
  if (uniqueN(dt, by = c("CHR", "POS")) != nrow(dt)) {
    stop(COHORT, " file is not unique by CHR:POS: ", path)
  }

  internal_bad <- dt[
    usable_for_replication == TRUE &
      (
        is.na(effect_allele_aligned) |
          is.na(other_allele_aligned) |
          effect_allele_aligned != A1_koges |
          other_allele_aligned != A2_koges |
          !is.finite(BETA_aligned)
      )
    ]
  if (nrow(internal_bad) > 0L) {
    stop(
      COHORT, " aligned-column integrity check failed in ", basename(path),
      ". First position: ", internal_bad$CHR[1L], ":", internal_bad$POS[1L]
    )
  }

  setkey(dt, CHR, POS)
  dt
}

# -----------------------------------------------------------------------------
# 3. Read and normalize KoGES lead table
# -----------------------------------------------------------------------------
lead_raw <- fread(KOGES_LEAD_FILE)
if (nrow(lead_raw) == 0L) {
  stop("KoGES lead table is empty: ", KOGES_LEAD_FILE)
}

nms <- names(lead_raw)
col_phen <- pick_column(
  nms, c("Phenotype", "Phenotype_label"), "phenotype label", required = FALSE
)
col_code <- pick_column(
  nms, c("VarNM", "Phenotype_file", "Pheno_code"), "BOLT phenotype code"
)
col_lead <- pick_column(
  nms, c("Lead_variant", "Lead variant"), "lead variant"
)
col_variant_id <- pick_column(
  nms, c("VariantID", "Variant ID"), "variant ID", required = FALSE
)
col_chr <- pick_column(nms, c("CHR", "Chromosome"), "chromosome")
col_bp <- pick_column(nms, c("BP", "Base pair position"), "base-pair position")
col_variant_type <- pick_column(
  nms, c("Variant_type", "Variant type"), "variant type", required = FALSE
)
col_gene <- pick_column(
  nms, c("Nearest_gene", "Nearest gene"), "nearest gene", required = FALSE
)
col_nea <- pick_column(
  nms, c("NEA", "Other allele", "Other Allele"), "other allele"
)
col_ea <- pick_column(
  nms, c("EA", "Effect allele", "Effect Allele"), "effect allele"
)
col_beta <- pick_column(nms, c("Beta", "BETA"), "BETA")
col_se <- pick_column(nms, c("SE", "Standard error"), "standard error")
col_p <- pick_column(nms, c("P", "P-value"), "P-value")
col_p_field <- pick_column(
  nms, c("P_field", "Association P field"), "P field", required = FALSE
)
col_eaf <- pick_column(
  nms, c("EAF", "Effect allele frequency"), "effect allele frequency",
  required = FALSE
)
col_n <- pick_column(
  nms, c("N", "Sample size"), "sample size", required = FALSE
)

lead <- data.table(
  lead_row_id = seq_len(nrow(lead_raw)),
  koges_phenotype = if (!is.na(col_phen)) {
    as.character(lead_raw[[col_phen]])
  } else {
    NA_character_
  },
  koges_pheno_code = as.character(lead_raw[[col_code]]),
  koges_lead_variant = as.character(lead_raw[[col_lead]]),
  koges_variant_id = if (!is.na(col_variant_id)) {
    as.character(lead_raw[[col_variant_id]])
  } else {
    NA_character_
  },
  koges_chr = norm_chr(lead_raw[[col_chr]]),
  koges_bp = norm_pos(lead_raw[[col_bp]]),
  koges_variant_type = if (!is.na(col_variant_type)) {
    as.character(lead_raw[[col_variant_type]])
  } else {
    NA_character_
  },
  koges_nearest_gene = if (!is.na(col_gene)) {
    as.character(lead_raw[[col_gene]])
  } else {
    NA_character_
  },
  koges_nea = norm_allele(lead_raw[[col_nea]]),
  koges_ea = norm_allele(lead_raw[[col_ea]]),
  koges_beta = suppressWarnings(as.numeric(lead_raw[[col_beta]])),
  koges_se = suppressWarnings(as.numeric(lead_raw[[col_se]])),
  koges_p = suppressWarnings(as.numeric(lead_raw[[col_p]])),
  koges_p_field = if (!is.na(col_p_field)) {
    as.character(lead_raw[[col_p_field]])
  } else {
    NA_character_
  },
  koges_eaf = if (!is.na(col_eaf)) {
    suppressWarnings(as.numeric(lead_raw[[col_eaf]]))
  } else {
    NA_real_
  },
  koges_n = if (!is.na(col_n)) {
    suppressWarnings(as.numeric(lead_raw[[col_n]]))
  } else {
    NA_real_
  }
)

lead[, koges_pheno_code := sub("\\.txt$", "", basename(koges_pheno_code))]
lead[, ukb_output_stem := unname(bolt_to_ukb_stem[koges_pheno_code])]
lead[, ukb_file := fifelse(
  !is.na(ukb_output_stem),
  file.path(EXT_DIR, paste0(ukb_output_stem, ".txt")),
  NA_character_
)]
lead[, ukb_mapping_status := fcase(
  is.na(ukb_output_stem), "no_matching_ukb_phenotype",
  !file.exists(ukb_file), "mapped_ukb_file_missing",
  default = "mapped_file_found"
)]

invalid_lead <- lead[
  is.na(koges_chr) | !(koges_chr %between% c(1L, 22L)) |
    is.na(koges_bp) | koges_bp <= 0L |
    is.na(koges_ea) | !nzchar(koges_ea) |
    is.na(koges_nea) | !nzchar(koges_nea) |
    !is.finite(koges_beta) |
    !is.finite(koges_p) | koges_p < 0 | koges_p > 1
  ]
if (nrow(invalid_lead) > 0L) {
  stop(
    "Invalid KoGES lead rows detected. First lead_row_id: ",
    invalid_lead$lead_row_id[1L]
  )
}

if (anyDuplicated(lead, by = c("koges_pheno_code", "koges_chr", "koges_bp"))) {
  stop("Duplicate KoGES lead CHR:BP within a phenotype was detected.")
}

n_non_gws_leads <- lead[koges_p > GWS_P, .N]
if (n_non_gws_leads > 0L) {
  warning(
    n_non_gws_leads,
    " lead-table rows have KoGES P > GWS_P. They are retained but should be checked."
  )
}

cat("KoGES lead variants:", nrow(lead), "\n")
cat(
  "Mapped to an existing ", COHORT, " phenotype file: ",
  lead[ukb_mapping_status == "mapped_file_found", .N], "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# 4. Process each external phenotype once
# -----------------------------------------------------------------------------
result_rows <- vector("list", nrow(lead))
gws_hit_list <- list()
gws_counter <- 0L

# First fill rows without a usable external mapping.
unmapped_idx <- lead[ukb_mapping_status != "mapped_file_found", lead_row_id]
for (i in unmapped_idx) {
  x <- lead[lead_row_id == i]
  result_rows[[i]] <- cbind(
    x,
    data.table(
      ukb_source_panel = NA_character_,
      ukb_source_trait = NA_character_,
      ukb_n_variants_in_window = NA_integer_,
      ukb_n_gws_hits_in_window = NA_integer_,
      ukb_top_window_chr = NA_integer_,
      ukb_top_window_pos = NA_integer_,
      ukb_top_window_snp = NA_character_,
      ukb_top_window_rsid_koges = NA_character_,
      ukb_top_window_p = NA_real_,
      ukb_top_window_beta_orig = NA_real_,
      ukb_top_window_distance_bp = NA_integer_,
      ukb_top_gws_chr = NA_integer_,
      ukb_top_gws_pos = NA_integer_,
      ukb_top_gws_snp = NA_character_,
      ukb_top_gws_rsid_koges = NA_character_,
      ukb_top_gws_p = NA_real_,
      ukb_top_gws_beta_orig = NA_real_,
      ukb_top_gws_distance_bp = NA_integer_,
      ukb_exact_chrpos_present = NA,
      ukb_exact_rsid_matches_koges = NA,
      ukb_exact_snp = NA_character_,
      ukb_exact_rsid_koges = NA_character_,
      ukb_exact_a1_koges = NA_character_,
      ukb_exact_a2_koges = NA_character_,
      ukb_exact_effect_allele_orig = NA_character_,
      ukb_exact_other_allele_orig = NA_character_,
      ukb_exact_beta_orig = NA_real_,
      ukb_exact_beta_aligned_a1 = NA_real_,
      ukb_exact_se = NA_real_,
      ukb_exact_p = NA_real_,
      ukb_exact_eaf_orig = NA_real_,
      ukb_exact_eaf_aligned_a1 = NA_real_,
      ukb_exact_info = NA_real_,
      ukb_exact_variant_type = NA_character_,
      ukb_exact_match_type = NA_character_,
      ukb_exact_usable_for_replication = NA,
      koges_to_a1_match_type = NA_character_,
      koges_beta_aligned_a1 = NA_real_,
      replication_testable = FALSE,
      replication_not_testable_reason = x$ukb_mapping_status,
      direction_concordant = NA
    )
  )
}

mapped_stems <- unique(
  lead[ukb_mapping_status == "mapped_file_found", ukb_output_stem]
)

for (stem in mapped_stems) {
  idx <- lead[
    ukb_mapping_status == "mapped_file_found" & ukb_output_stem == stem,
    lead_row_id
    ]
  ukb_path <- unique(lead[lead_row_id %in% idx, ukb_file])
  if (length(ukb_path) != 1L) {
    stop("Expected one ", COHORT, " file for output stem: ", stem)
  }
  
  cat("Reading ", COHORT, ": ", basename(ukb_path),
      " for ", length(idx), " KoGES lead variant(s)\n", sep = "")
  ukb <- read_ext_cont(ukb_path)
  
  for (i in idx) {
    x <- lead[lead_row_id == i]
    chr_i <- x$koges_chr
    bp_i <- x$koges_bp
    lower <- max(1L, as.integer(bp_i - NOVELTY_WINDOW_BP))
    upper <- as.integer(bp_i + NOVELTY_WINDOW_BP)
    
    window <- ukb[
      CHR == chr_i & POS >= lower & POS <= upper & usable_for_novelty == TRUE
      ]
    window <- window[is.finite(P)]
    gws <- window[P <= GWS_P]
    
    top_window <- if (nrow(window) > 0L) {
      window[which.min(P)][1L]
    } else {
      NULL
    }
    top_gws <- if (nrow(gws) > 0L) {
      gws[which.min(P)][1L]
    } else {
      NULL
    }
    
    if (nrow(gws) > 0L) {
      gws_counter <- gws_counter + 1L
      gws_hit_list[[gws_counter]] <- gws[, .(
        lead_row_id = i,
        koges_phenotype = x$koges_phenotype,
        koges_pheno_code = x$koges_pheno_code,
        koges_lead_variant = x$koges_lead_variant,
        koges_chr = chr_i,
        koges_bp = bp_i,
        ukb_output_stem = stem,
        ukb_chr = CHR,
        ukb_pos = POS,
        distance_from_koges_lead_bp = as.integer(POS - bp_i),
        ukb_snp = SNP,
        ukb_rsid_koges = rsID_koges,
        ukb_effect_allele_orig = effect_allele_orig,
        ukb_other_allele_orig = other_allele_orig,
        ukb_beta_orig = BETA_orig,
        ukb_se = SE,
        ukb_p = P,
        ukb_eaf_orig = effect_allele_frequency_orig,
        ukb_info = NA_real_,
        ukb_variant_type = variant_type,
        ukb_match_type = match_type,
        ukb_usable_for_replication = usable_for_replication
      )]
    }
    
    exact <- ukb[CHR == chr_i & POS == bp_i]
    if (nrow(exact) > 1L) {
      stop("Multiple ", COHORT, " rows at exact lead CHR:POS: ", chr_i, ":", bp_i)
    }
    
    source_panel <- if (nrow(ukb) > 0L && "source_panel" %in% names(ukb)) {
      unique(ukb$source_panel)[1L]
    } else {
      NA_character_
    }
    source_trait <- if (nrow(ukb) > 0L && "trait" %in% names(ukb)) {
      unique(ukb$trait)[1L]
    } else {
      NA_character_
    }
    
    exact_present <- nrow(exact) == 1L
    exact_usable <- exact_present && isTRUE(exact$usable_for_replication[1L])
    
    k_align <- if (exact_present) {
      align_koges_beta_to_a1(
        ea = x$koges_ea,
        nea = x$koges_nea,
        beta = x$koges_beta,
        a1 = exact$A1_koges[1L],
        a2 = exact$A2_koges[1L]
      )
    } else {
      list(ok = FALSE, beta_a1 = NA_real_, match_type = "exact_variant_absent")
    }
    
    replication_testable <- (
      exact_present &&
        exact_usable &&
        isTRUE(k_align$ok) &&
        is.finite(exact$BETA_aligned[1L]) &&
        is.finite(exact$P[1L])
    )
    
    not_testable_reason <- if (!exact_present) {
      "lead_variant_absent_after_ukb_basic_qc"
    } else if (!exact_usable) {
      paste0("ukb_not_usable_for_replication:", exact$match_type[1L])
    } else if (!isTRUE(k_align$ok)) {
      paste0("koges_ea_nea_not_aligned_to_bim_a1_a2:", k_align$match_type)
    } else if (!is.finite(exact$BETA_aligned[1L])) {
      "ukb_aligned_beta_missing"
    } else if (!is.finite(exact$P[1L])) {
      "ukb_p_missing"
    } else {
      NA_character_
    }
    
    direction <- if (replication_testable) {
      safe_sign_concordance(k_align$beta_a1, exact$BETA_aligned[1L])
    } else {
      NA
    }
    
    result_rows[[i]] <- cbind(
      x,
      data.table(
        ukb_source_panel = source_panel,
        ukb_source_trait = source_trait,
        ukb_n_variants_in_window = nrow(window),
        ukb_n_gws_hits_in_window = nrow(gws),
        ukb_top_window_chr = if (!is.null(top_window)) top_window$CHR else NA_integer_,
        ukb_top_window_pos = if (!is.null(top_window)) top_window$POS else NA_integer_,
        ukb_top_window_snp = if (!is.null(top_window)) top_window$SNP else NA_character_,
        ukb_top_window_rsid_koges = if (!is.null(top_window)) {
          top_window$rsID_koges
        } else {
          NA_character_
        },
        ukb_top_window_p = if (!is.null(top_window)) top_window$P else NA_real_,
        ukb_top_window_beta_orig = if (!is.null(top_window)) {
          top_window$BETA_orig
        } else {
          NA_real_
        },
        ukb_top_window_distance_bp = if (!is.null(top_window)) {
          as.integer(top_window$POS - bp_i)
        } else {
          NA_integer_
        },
        ukb_top_gws_chr = if (!is.null(top_gws)) top_gws$CHR else NA_integer_,
        ukb_top_gws_pos = if (!is.null(top_gws)) top_gws$POS else NA_integer_,
        ukb_top_gws_snp = if (!is.null(top_gws)) top_gws$SNP else NA_character_,
        ukb_top_gws_rsid_koges = if (!is.null(top_gws)) {
          top_gws$rsID_koges
        } else {
          NA_character_
        },
        ukb_top_gws_p = if (!is.null(top_gws)) top_gws$P else NA_real_,
        ukb_top_gws_beta_orig = if (!is.null(top_gws)) {
          top_gws$BETA_orig
        } else {
          NA_real_
        },
        ukb_top_gws_distance_bp = if (!is.null(top_gws)) {
          as.integer(top_gws$POS - bp_i)
        } else {
          NA_integer_
        },
        ukb_exact_chrpos_present = exact_present,
        ukb_exact_rsid_matches_koges = if (exact_present) {
          !is.na(exact$rsID_koges[1L]) &&
            toupper(exact$rsID_koges[1L]) == toupper(x$koges_lead_variant)
        } else {
          NA
        },
        ukb_exact_snp = if (exact_present) exact$SNP[1L] else NA_character_,
        ukb_exact_rsid_koges = if (exact_present) {
          exact$rsID_koges[1L]
        } else {
          NA_character_
        },
        ukb_exact_a1_koges = if (exact_present) exact$A1_koges[1L] else NA_character_,
        ukb_exact_a2_koges = if (exact_present) exact$A2_koges[1L] else NA_character_,
        ukb_exact_effect_allele_orig = if (exact_present) {
          exact$effect_allele_orig[1L]
        } else {
          NA_character_
        },
        ukb_exact_other_allele_orig = if (exact_present) {
          exact$other_allele_orig[1L]
        } else {
          NA_character_
        },
        ukb_exact_beta_orig = if (exact_present) exact$BETA_orig[1L] else NA_real_,
        ukb_exact_beta_aligned_a1 = if (exact_present) {
          exact$BETA_aligned[1L]
        } else {
          NA_real_
        },
        ukb_exact_se = if (exact_present) exact$SE[1L] else NA_real_,
        ukb_exact_p = if (exact_present) exact$P[1L] else NA_real_,
        ukb_exact_eaf_orig = if (exact_present) {
          exact$effect_allele_frequency_orig[1L]
        } else {
          NA_real_
        },
        ukb_exact_eaf_aligned_a1 = if (exact_present) {
          exact$effect_allele_frequency_aligned[1L]
        } else {
          NA_real_
        },
        ukb_exact_info = NA_real_,
        ukb_exact_variant_type = if (exact_present) {
          exact$variant_type[1L]
        } else {
          NA_character_
        },
        ukb_exact_match_type = if (exact_present) {
          exact$match_type[1L]
        } else {
          NA_character_
        },
        ukb_exact_usable_for_replication = if (exact_present) {
          exact$usable_for_replication[1L]
        } else {
          FALSE
        },
        koges_to_a1_match_type = k_align$match_type,
        koges_beta_aligned_a1 = k_align$beta_a1,
        replication_testable = replication_testable,
        replication_not_testable_reason = not_testable_reason,
        direction_concordant = direction
      )
    )
  }
  
  rm(ukb)
  invisible(gc())
}

result <- rbindlist(result_rows, fill = TRUE, use.names = TRUE)
setorder(result, lead_row_id)

if (nrow(result) != nrow(lead)) {
  stop("Internal error: result row count differs from lead row count.")
}

# -----------------------------------------------------------------------------
# 5. Novelty classification relative to UKB
# -----------------------------------------------------------------------------
result[, novelty_assessable := (
  ukb_mapping_status == "mapped_file_found" &
    !is.na(ukb_n_variants_in_window) &
    ukb_n_variants_in_window > 0L
)]
result[, novel_relative_to_ukb := fcase(
  novelty_assessable & ukb_n_gws_hits_in_window == 0L, TRUE,
  novelty_assessable & ukb_n_gws_hits_in_window > 0L, FALSE,
  default = NA
)]
result[, ukb_novelty_status := fcase(
  ukb_mapping_status == "no_matching_ukb_phenotype",
  "Not assessed: no matching UKB phenotype",
  ukb_mapping_status == "mapped_ukb_file_missing",
  "Not assessed: mapped UKB file missing",
  ukb_n_variants_in_window == 0L,
  "Not assessed: no UKB basic-QC variants within +/-1 Mb",
  ukb_n_gws_hits_in_window > 0L,
  "Not novel relative to UKB: GWS signal within +/-1 Mb",
  ukb_n_gws_hits_in_window == 0L,
  "Novel relative to UKB: no GWS signal within +/-1 Mb",
  default = "Not assessed"
)]

# -----------------------------------------------------------------------------
# 6. Multiple-testing thresholds and replication classification
# -----------------------------------------------------------------------------
result[, mapped_for_replication := ukb_mapping_status == "mapped_file_found"]

n_global_mapped <- result[mapped_for_replication == TRUE, .N]
n_global_testable <- result[replication_testable == TRUE, .N]

global_mapped_threshold <- if (n_global_mapped > 0L) {
  NOMINAL_REPLICATION_P / n_global_mapped
} else {
  NA_real_
}
global_testable_threshold <- if (n_global_testable > 0L) {
  NOMINAL_REPLICATION_P / n_global_testable
} else {
  NA_real_
}

result[, n_trait_mapped := sum(mapped_for_replication), by = koges_pheno_code]
result[, n_trait_testable := sum(replication_testable), by = koges_pheno_code]
result[, `:=`(
  replication_threshold_nominal = NOMINAL_REPLICATION_P,
  replication_threshold_global_mapped = global_mapped_threshold,
  replication_threshold_global_testable = global_testable_threshold,
  replication_threshold_trait_mapped = fifelse(
    n_trait_mapped > 0L,
    NOMINAL_REPLICATION_P / n_trait_mapped,
    NA_real_
  ),
  replication_threshold_trait_testable = fifelse(
    n_trait_testable > 0L,
    NOMINAL_REPLICATION_P / n_trait_testable,
    NA_real_
  )
)]

result[, replicated_nominal := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    ukb_exact_p <= replication_threshold_nominal
)]
result[, replicated_global_mapped := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    ukb_exact_p <= replication_threshold_global_mapped
)]
result[, replicated_global_testable := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    ukb_exact_p <= replication_threshold_global_testable
)]
result[, replicated_trait_mapped := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    ukb_exact_p <= replication_threshold_trait_mapped
)]
result[, replicated_trait_testable := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    ukb_exact_p <= replication_threshold_trait_testable
)]

primary_rep_col <- paste0("replicated_", PRIMARY_CORRECTION)
primary_thr_col <- paste0("replication_threshold_", PRIMARY_CORRECTION)
result[, replicated_primary := get(primary_rep_col)]
result[, replication_threshold_primary := get(primary_thr_col)]
result[, replication_primary_correction := PRIMARY_CORRECTION]

result[, replication_status := fcase(
  ukb_mapping_status == "no_matching_ukb_phenotype",
  "Not assessed: no matching UKB phenotype",
  ukb_mapping_status == "mapped_ukb_file_missing",
  "Not assessed: mapped UKB file missing",
  replication_testable == FALSE,
  paste0("Not testable: ", replication_not_testable_reason),
  replicated_primary == TRUE,
  paste0("Replicated: direction concordant and ", PRIMARY_CORRECTION,
         " Bonferroni significant"),
  direction_concordant %in% FALSE & ukb_exact_p <= NOMINAL_REPLICATION_P,
  "Not replicated: nominally significant but opposite direction",
  replicated_nominal == TRUE,
  paste0("Nominal replication only: does not pass ", PRIMARY_CORRECTION,
         " Bonferroni threshold"),
  default = "Not replicated"
)]

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
rename_cohort_prefix <- function(dt) {
  if (pfx == "ukb") return(dt)
  nms <- names(dt)
  new_nms <- nms
  new_nms <- sub("^ukb_", paste0(pfx, "_"), new_nms)
  new_nms <- gsub("ukb", pfx, new_nms, fixed = TRUE)
  if (any(duplicated(new_nms))) {
    stop("Column rename collision for cohort prefix: ", pfx)
  }
  setnames(dt, nms, new_nms)
  for (col in names(dt)) {
    if (is.character(dt[[col]])) {
      set(dt, j = col, value = gsub("ukb", pfx, dt[[col]], fixed = TRUE))
      set(dt, j = col, value = gsub("UKB", COHORT, dt[[col]], fixed = TRUE))
    }
  }
  dt
}
result <- rename_cohort_prefix(result)

fwrite(result, RESULT_FILE, sep = "\t", quote = FALSE, na = "NA")

cat("\nCompleted.\n")
cat("Cohort:", COHORT, "\n")
cat("Detailed result:", RESULT_FILE, "\n")
cat("Novelty (assessable):\n")
nov_col <- paste0("novel_relative_to_", pfx)
print(result[novelty_assessable %in% TRUE, .N, by = c(nov_col)])
cat("Replication (primary):\n")
print(result[replication_testable %in% TRUE, .N, by = replicated_primary])
