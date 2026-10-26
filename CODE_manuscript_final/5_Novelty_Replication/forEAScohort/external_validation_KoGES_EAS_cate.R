#!/usr/bin/env Rscript
# external_validation_KoGES_EAS_cate.R
#
# Purpose:
#   Assess locus-level novelty (+/-1 Mb GWS) and exact-variant replication of
#   KoGES SAIGE binary (Cate) lead variants against one East Asian external
#   cohort (BBJ or TPMI), selected via --cohort.
#
#   Matching Excel may map multiple external GWAS to one KoGES phenotype with
#   levels: direct match / subordinate disease / related biomarker /
#   underlying lesion. Direct-match results are OR-aggregated; other-match
#   levels can add "additionally known/replicated".
#
# Inputs:
#   KOGES_LEAD_FILE (default):
#     DATA_manuscript_final/Extended_Data/Cate_lead_variants.tsv
#   EXT_DIR (or BBJ_DIR / TPMI_DIR):
#     external binary sumstats from {COHORT}_Cate_modify.R
#     ({EXT_DIR}/{KoGES_VAR}.txt; default under
#     DATA_baseline/forEAScohort/{COHORT}/Cate when unset)
#   MATCHING_XLSX (default):
#     KoGES_EAS_matching_binary_SAIGE_checked.xlsx
#   Matching Excel may use English ('matching_level' / 'direct match') or
#   Korean legacy headers/labels; both are accepted.
#   defaults to that value for compatibility.
#
# Outputs (OUTPUT_DIR; default Extended_Data):
#   KoGES_SAIGE_leads_{TAG}_novelty_replication_cate.tsv
#   KoGES_SAIGE_leads_{TAG}_per_match_cate.tsv
#
# Usage:
#   Rscript external_validation_KoGES_EAS_cate.R --cohort=BBJ
#   Rscript external_validation_KoGES_EAS_cate.R --cohort=TPMI
#   EXT_DIR=... OUTPUT_DIR=... MATCHING_XLSX=... Rscript ... --cohort=BBJ
#

suppressPackageStartupMessages({
  library(data.table)
})

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
  cat("Usage: Rscript external_validation_KoGES_EAS_cate.R --cohort=BBJ|TPMI\n")
  quit(save = "no", status = 0)
}
v <- grab_eq("--cohort"); if (is.null(v)) v <- grab_sp("--cohort")
if (is.null(v)) v <- Sys.getenv("COHORT", "")
COHORT <- trimws(as.character(v))
if (!COHORT %chin% c("BBJ", "TPMI")) {
  stop("--cohort=BBJ|TPMI is required.")
}
if (length(args)) stop("Unknown args: ", paste(args, collapse = " "))

PATH_COL <- COHORT
TAG <- COHORT
pfx <- tolower(COHORT)

KOGES_LEAD_FILE <- Sys.getenv(
  "KOGES_LEAD_FILE",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/Cate_lead_variants.tsv"
)

EXT_DIR <- Sys.getenv("EXT_DIR", "")
if (!nzchar(EXT_DIR)) {
  alt <- switch(
    COHORT,
    BBJ = Sys.getenv("BBJ_DIR", ""),
    TPMI = Sys.getenv("TPMI_DIR", "")
  )
  EXT_DIR <- if (nzchar(alt)) {
    alt
  } else {
    file.path(
      "/data4/sjim/openKoGES_revision/DATA_baseline/forEAScohort",
      COHORT, "Cate"
    )
  }
}

MATCHING_XLSX <- Sys.getenv(
  "MATCHING_XLSX",
  paste0(
    "/data4/sjim/openKoGES_revision/Sumstats/EAS/",
    "KoGES_EAS_matching_binary_SAIGE_checked.xlsx"
  )
)

OUTPUT_DIR <- Sys.getenv(
  "OUTPUT_DIR",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
)

NOVELTY_WINDOW_BP <- as.integer(Sys.getenv("NOVELTY_WINDOW_BP", "1000000"))
GWS_P <- as.numeric(Sys.getenv("GWS_P", "5e-8"))
NOMINAL_REPLICATION_P <- as.numeric(
  Sys.getenv("NOMINAL_REPLICATION_P", "0.05")
)
# Direct-match label(s); Excel may still use Korean "직접 매칭".
DIRECT_MATCH_LABEL <- Sys.getenv("DIRECT_MATCH_LABEL", "direct match")
DIRECT_MATCH_LABELS <- unique(c(DIRECT_MATCH_LABEL, "direct match", "직접 매칭"))

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

EMPTY_TOKENS <- c("", "-", "0", "x", "na", "none", "nan")

RESULT_FILE <- file.path(
  OUTPUT_DIR, paste0("KoGES_SAIGE_leads_", TAG, "_novelty_replication_cate.tsv")
)
PER_MATCH_FILE <- file.path(
  OUTPUT_DIR, paste0("KoGES_SAIGE_leads_", TAG, "_per_match_cate.tsv")
)

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(KOGES_LEAD_FILE)) {
  stop("KoGES lead file not found: ", KOGES_LEAD_FILE)
}
if (!dir.exists(EXT_DIR)) {
  stop(COHORT, " directory not found: ", EXT_DIR)
}
if (!file.exists(MATCHING_XLSX)) {
  stop("Matching Excel not found: ", MATCHING_XLSX)
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("package 'readxl' is required")
}

# -----------------------------------------------------------------------------
# Helpers
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

is_empty <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | !nzchar(x) | x %chin% EMPTY_TOKENS
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

lead_stem_from_varnm <- function(x) {
  x <- sub("\\.txt$", "", basename(trimws(as.character(x))))
  x <- sub("\\.SAIGE(\\.SEfixed)?$", "", x, ignore.case = TRUE)
  sub("_SAIGE_cumulative_inci.*$", "", x)
}

assign_lead_stem <- function(koges_var, lead_stems) {
  koges_var <- trimws(as.character(koges_var))
  out <- rep(NA_character_, length(koges_var))
  for (i in seq_along(koges_var)) {
    v <- koges_var[i]
    if (is.na(v) || !nzchar(v)) next
    exact <- lead_stems[lead_stems == v]
    if (length(exact)) {
      out[i] <- exact[[1L]]
      next
    }
    pref <- lead_stems[startsWith(v, paste0(lead_stems, "_"))]
    if (length(pref)) {
      out[i] <- pref[[which.max(nchar(pref))]]
    }
  }
  out
}

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
    koges_lead_stem = character(),
    koges_lead_variant = character(),
    koges_chr = integer(),
    koges_bp = integer(),
    bbj_output_stem = character(),
    matching_level = character(),
    is_direct_match = logical(),
    bbj_chr = integer(),
    bbj_pos = integer(),
    distance_from_koges_lead_bp = integer(),
    bbj_snp = character(),
    bbj_rsid_koges = character(),
    bbj_effect_allele_orig = character(),
    bbj_other_allele_orig = character(),
    bbj_beta_orig = numeric(),
    bbj_se = numeric(),
    bbj_p = numeric(),
    bbj_eaf_orig = numeric(),
    bbj_info = numeric(),
    bbj_variant_type = character(),
    bbj_match_type = character(),
    bbj_usable_for_replication = logical()
  )
}

empty_match_stats <- function() {
  data.table(
    bbj_source_panel = NA_character_,
    bbj_source_trait = NA_character_,
    bbj_n_variants_in_window = NA_integer_,
    bbj_n_gws_hits_in_window = NA_integer_,
    bbj_top_window_chr = NA_integer_,
    bbj_top_window_pos = NA_integer_,
    bbj_top_window_snp = NA_character_,
    bbj_top_window_rsid_koges = NA_character_,
    bbj_top_window_p = NA_real_,
    bbj_top_window_beta_orig = NA_real_,
    bbj_top_window_distance_bp = NA_integer_,
    bbj_top_gws_chr = NA_integer_,
    bbj_top_gws_pos = NA_integer_,
    bbj_top_gws_snp = NA_character_,
    bbj_top_gws_rsid_koges = NA_character_,
    bbj_top_gws_p = NA_real_,
    bbj_top_gws_beta_orig = NA_real_,
    bbj_top_gws_distance_bp = NA_integer_,
    bbj_exact_chrpos_present = NA,
    bbj_exact_rsid_matches_koges = NA,
    bbj_exact_snp = NA_character_,
    bbj_exact_rsid_koges = NA_character_,
    bbj_exact_a1_koges = NA_character_,
    bbj_exact_a2_koges = NA_character_,
    bbj_exact_effect_allele_orig = NA_character_,
    bbj_exact_other_allele_orig = NA_character_,
    bbj_exact_beta_orig = NA_real_,
    bbj_exact_beta_aligned_a1 = NA_real_,
    bbj_exact_se = NA_real_,
    bbj_exact_p = NA_real_,
    bbj_exact_eaf_orig = NA_real_,
    bbj_exact_eaf_aligned_a1 = NA_real_,
    bbj_exact_info = NA_real_,
    bbj_exact_variant_type = NA_character_,
    bbj_exact_match_type = NA_character_,
    bbj_exact_usable_for_replication = NA,
    koges_to_a1_match_type = NA_character_,
    koges_beta_aligned_a1 = NA_real_,
    replication_testable = FALSE,
    replication_not_testable_reason = NA_character_,
    direction_concordant = NA
  )
}

read_ext_cate <- function(path) {
  header <- names(fread(path, nrows = 0L))
  if (COHORT == "BBJ") {
    required <- c(
      "CHR", "POS", "SNP", "rsID_koges", "A1_koges", "A2_koges",
      "effect_allele_orig", "other_allele_orig", "BETA_orig",
      "BETA_aligned", "SE", "P", "effect_allele_frequency_orig",
      "effect_allele_frequency_aligned", "effect_allele_aligned",
      "other_allele_aligned", "INFO", "variant_type", "match_type",
      "usable_for_replication", "usable_for_novelty"
    )
  } else {
    # TPMI: rsID_tpmi aliased to SNP; INFO optional
    required <- c(
      "CHR", "POS", "rsID_tpmi", "rsID_koges", "A1_koges", "A2_koges",
      "effect_allele_orig", "other_allele_orig", "BETA_orig",
      "BETA_aligned", "SE", "P", "effect_allele_frequency_orig",
      "effect_allele_frequency_aligned", "effect_allele_aligned",
      "other_allele_aligned", "variant_type", "match_type",
      "usable_for_replication", "usable_for_novelty"
    )
  }
  missing <- setdiff(required, header)
  if (length(missing) > 0L) {
    stop(
      COHORT, " file is missing columns (", basename(path), "): ",
      paste(missing, collapse = ", ")
    )
  }

  optional <- c(
    "bbj_format", "SNP", "INFO", "source_panel", "trait", "output_stem",
    "BP_GRCh38", "variantID_koges", "ref_koges", "alt_koges",
    "EAF_A1_koges", "MAF", "NCHROBS_koges",
    "is_indel", "is_palindromic"
  )
  # unique(): BBJ lists SNP/INFO in both required and optional
  dt <- fread(path, select = unique(c(required, intersect(optional, header))))
  if (COHORT == "TPMI") {
    dt[, SNP := as.character(rsID_tpmi)]
  }
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
  if ("INFO" %in% names(dt)) {
    dt[, INFO := suppressWarnings(as.numeric(INFO))]
  } else {
    dt[, INFO := NA_real_]
  }

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

assess_lead_in_bbj <- function(x, bbj, mapping_status) {
  stats <- empty_match_stats()
  gws_dt <- empty_gws_hits()
  if (mapping_status != "mapped_file_found") {
    stats$replication_not_testable_reason <- mapping_status
    return(list(stats = stats, gws = gws_dt))
  }

  chr_i <- x$koges_chr
  bp_i <- x$koges_bp
  lower <- max(1L, as.integer(bp_i - NOVELTY_WINDOW_BP))
  upper <- as.integer(bp_i + NOVELTY_WINDOW_BP)

  window <- bbj[
    CHR == chr_i & POS >= lower & POS <= upper & usable_for_novelty == TRUE
  ]
  window <- window[is.finite(P)]
  gws <- window[P <= GWS_P]

  top_window <- if (nrow(window) > 0L) window[which.min(P)][1L] else NULL
  top_gws <- if (nrow(gws) > 0L) gws[which.min(P)][1L] else NULL

  if (nrow(gws) > 0L) {
    gws_dt <- gws[, .(
      lead_row_id = x$lead_row_id,
      koges_phenotype = x$koges_phenotype,
      koges_pheno_code = x$koges_pheno_code,
      koges_lead_stem = x$koges_lead_stem,
      koges_lead_variant = x$koges_lead_variant,
      koges_chr = chr_i,
      koges_bp = bp_i,
      bbj_output_stem = x$bbj_output_stem,
      matching_level = x$matching_level,
      is_direct_match = x$is_direct_match,
      bbj_chr = CHR,
      bbj_pos = POS,
      distance_from_koges_lead_bp = as.integer(POS - bp_i),
      bbj_snp = SNP,
      bbj_rsid_koges = rsID_koges,
      bbj_effect_allele_orig = effect_allele_orig,
      bbj_other_allele_orig = other_allele_orig,
      bbj_beta_orig = BETA_orig,
      bbj_se = SE,
      bbj_p = P,
      bbj_eaf_orig = effect_allele_frequency_orig,
      bbj_info = if ("INFO" %in% names(bbj)) INFO else NA_real_,
      bbj_variant_type = variant_type,
      bbj_match_type = match_type,
      bbj_usable_for_replication = usable_for_replication
    )]
  }

  exact <- bbj[CHR == chr_i & POS == bp_i]
  if (nrow(exact) > 1L) {
    stop("Multiple ", COHORT, " rows at exact lead CHR:POS: ", chr_i, ":", bp_i)
  }

  source_panel <- if (nrow(bbj) > 0L && "source_panel" %in% names(bbj)) {
    unique(bbj$source_panel)[1L]
  } else {
    NA_character_
  }
  source_trait <- if (nrow(bbj) > 0L && "trait" %in% names(bbj)) {
    unique(bbj$trait)[1L]
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
    "lead_variant_absent_after_bbj_basic_qc"
  } else if (!exact_usable) {
    paste0("bbj_not_usable_for_replication:", exact$match_type[1L])
  } else if (!isTRUE(k_align$ok)) {
    paste0("koges_ea_nea_not_aligned_to_bim_a1_a2:", k_align$match_type)
  } else if (!is.finite(exact$BETA_aligned[1L])) {
    "bbj_aligned_beta_missing"
  } else if (!is.finite(exact$P[1L])) {
    "bbj_p_missing"
  } else {
    NA_character_
  }

  direction <- if (replication_testable) {
    safe_sign_concordance(k_align$beta_a1, exact$BETA_aligned[1L])
  } else {
    NA
  }

  stats <- data.table(
    bbj_source_panel = source_panel,
    bbj_source_trait = source_trait,
    bbj_n_variants_in_window = nrow(window),
    bbj_n_gws_hits_in_window = nrow(gws),
    bbj_top_window_chr = if (!is.null(top_window)) top_window$CHR else NA_integer_,
    bbj_top_window_pos = if (!is.null(top_window)) top_window$POS else NA_integer_,
    bbj_top_window_snp = if (!is.null(top_window)) top_window$SNP else NA_character_,
    bbj_top_window_rsid_koges = if (!is.null(top_window)) {
      top_window$rsID_koges
    } else {
      NA_character_
    },
    bbj_top_window_p = if (!is.null(top_window)) top_window$P else NA_real_,
    bbj_top_window_beta_orig = if (!is.null(top_window)) {
      top_window$BETA_orig
    } else {
      NA_real_
    },
    bbj_top_window_distance_bp = if (!is.null(top_window)) {
      as.integer(top_window$POS - bp_i)
    } else {
      NA_integer_
    },
    bbj_top_gws_chr = if (!is.null(top_gws)) top_gws$CHR else NA_integer_,
    bbj_top_gws_pos = if (!is.null(top_gws)) top_gws$POS else NA_integer_,
    bbj_top_gws_snp = if (!is.null(top_gws)) top_gws$SNP else NA_character_,
    bbj_top_gws_rsid_koges = if (!is.null(top_gws)) {
      top_gws$rsID_koges
    } else {
      NA_character_
    },
    bbj_top_gws_p = if (!is.null(top_gws)) top_gws$P else NA_real_,
    bbj_top_gws_beta_orig = if (!is.null(top_gws)) {
      top_gws$BETA_orig
    } else {
      NA_real_
    },
    bbj_top_gws_distance_bp = if (!is.null(top_gws)) {
      as.integer(top_gws$POS - bp_i)
    } else {
      NA_integer_
    },
    bbj_exact_chrpos_present = exact_present,
    bbj_exact_rsid_matches_koges = if (exact_present) {
      !is.na(exact$rsID_koges[1L]) &&
        toupper(exact$rsID_koges[1L]) == toupper(x$koges_lead_variant)
    } else {
      NA
    },
    bbj_exact_snp = if (exact_present) exact$SNP[1L] else NA_character_,
    bbj_exact_rsid_koges = if (exact_present) {
      exact$rsID_koges[1L]
    } else {
      NA_character_
    },
    bbj_exact_a1_koges = if (exact_present) exact$A1_koges[1L] else NA_character_,
    bbj_exact_a2_koges = if (exact_present) exact$A2_koges[1L] else NA_character_,
    bbj_exact_effect_allele_orig = if (exact_present) {
      exact$effect_allele_orig[1L]
    } else {
      NA_character_
    },
    bbj_exact_other_allele_orig = if (exact_present) {
      exact$other_allele_orig[1L]
    } else {
      NA_character_
    },
    bbj_exact_beta_orig = if (exact_present) exact$BETA_orig[1L] else NA_real_,
    bbj_exact_beta_aligned_a1 = if (exact_present) {
      exact$BETA_aligned[1L]
    } else {
      NA_real_
    },
    bbj_exact_se = if (exact_present) exact$SE[1L] else NA_real_,
    bbj_exact_p = if (exact_present) exact$P[1L] else NA_real_,
    bbj_exact_eaf_orig = if (exact_present) {
      exact$effect_allele_frequency_orig[1L]
    } else {
      NA_real_
    },
    bbj_exact_eaf_aligned_a1 = if (exact_present) {
      exact$effect_allele_frequency_aligned[1L]
    } else {
      NA_real_
    },
    bbj_exact_info = if (exact_present && "INFO" %in% names(bbj)) exact$INFO[1L] else NA_real_,
    bbj_exact_variant_type = if (exact_present) {
      exact$variant_type[1L]
    } else {
      NA_character_
    },
    bbj_exact_match_type = if (exact_present) {
      exact$match_type[1L]
    } else {
      NA_character_
    },
    bbj_exact_usable_for_replication = if (exact_present) {
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

  list(stats = stats, gws = gws_dt)
}

collapse_chr <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = ";")
}

# -----------------------------------------------------------------------------
# Lead table
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
  nms, c("VarNM", "Phenotype_file", "Pheno_code"), "SAIGE phenotype code"
)
col_lead <- pick_column(
  nms, c("Lead variant", "Lead_variant"), "lead variant"
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
col_eaf <- pick_column(
  nms, c("EAF", "Effect allele frequency"), "effect allele frequency",
  required = FALSE
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
  koges_eaf = if (!is.na(col_eaf)) {
    suppressWarnings(as.numeric(lead_raw[[col_eaf]]))
  } else {
    NA_real_
  }
)

lead[, koges_pheno_code := sub("\\.txt$", "", basename(koges_pheno_code))]
lead[, koges_lead_stem := lead_stem_from_varnm(koges_pheno_code)]

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

# -----------------------------------------------------------------------------
# Matching Excel map
# -----------------------------------------------------------------------------
meta <- as.data.table(readxl::read_excel(MATCHING_XLSX, sheet = 1))
match_level_col <- if ("matching_level" %chin% names(meta)) {
  "matching_level"
} else if ("매칭 수준" %chin% names(meta)) {
  "매칭 수준"
} else {
  stop(
    "Matching Excel is missing a matching-level column. ",
    "Expected 'matching_level' or Korean '매칭 수준'."
  )
}
need_meta <- c("KoGES_VAR", PATH_COL, match_level_col)
if (length(miss_meta <- setdiff(need_meta, names(meta)))) {
  stop("Excel missing columns: ", paste(miss_meta, collapse = ", "))
}

meta[, `:=`(
  KoGES_VAR = trimws(as.character(KoGES_VAR)),
  path_raw = trimws(as.character(get(PATH_COL))),
  matching_level = trimws(as.character(get(match_level_col))),
  excel_phenotype = if ("KoGES phenotype" %in% names(meta)) {
    trimws(as.character(`KoGES phenotype`))
  } else {
    NA_character_
  },
  excel_external = if ("External phenotype" %in% names(meta)) {
    trimws(as.character(`External phenotype`))
  } else {
    NA_character_
  },
  excel_data_type = if ("Data type" %in% names(meta)) {
    trimws(as.character(`Data type`))
  } else {
    NA_character_
  }
)]

meta[, is_direct_match := matching_level %chin% DIRECT_MATCH_LABELS]
# Keep direct-match rows even if the external path is empty (not assessable).
map <- meta[!is_empty(KoGES_VAR) & (!is_empty(path_raw) | is_direct_match == TRUE)]
if (!nrow(map)) stop("No Excel ", COHORT, " / direct-match rows.")

lead_stems <- unique(lead$koges_lead_stem)
map[, koges_lead_stem := assign_lead_stem(KoGES_VAR, lead_stems)]
map_unassigned <- map[is.na(koges_lead_stem)]
if (nrow(map_unassigned)) {
  warning(
    "Excel ", COHORT, " rows not assigned to a Cate lead stem: ",
    paste(unique(map_unassigned$KoGES_VAR), collapse = ", ")
  )
}

map <- map[!is.na(koges_lead_stem)]
map[, bbj_file := file.path(EXT_DIR, paste0(KoGES_VAR, ".txt"))]
map[, bbj_mapping_status := fcase(
  file.exists(bbj_file), "mapped_file_found",
  is_empty(path_raw), "no_bbj_source",
  default = "mapped_bbj_file_missing"
)]

cat("KoGES lead variants:", nrow(lead), "\n")
cat("Excel ", COHORT, " mappings assigned to lead stems: ", nrow(map), "\n", sep = "")
cat(
  "  direct match:", map[is_direct_match == TRUE, .N],
  " / other:", map[is_direct_match == FALSE, .N], "\n"
)
cat(
  COHORT, " files found:",
  map[bbj_mapping_status == "mapped_file_found", uniqueN(bbj_file)],
  "/", uniqueN(map$bbj_file), "\n"
)

# -----------------------------------------------------------------------------
# External sumstats (per-match)
# -----------------------------------------------------------------------------
match_jobs <- merge(
  lead,
  map[, .(
    koges_lead_stem,
    bbj_output_stem = KoGES_VAR,
    matching_level,
    is_direct_match,
    excel_external,
    excel_data_type,
    bbj_file,
    bbj_mapping_status
  )],
  by = "koges_lead_stem",
  all.x = TRUE,
  allow.cartesian = TRUE,
  sort = FALSE
)
setorder(match_jobs, lead_row_id, -is_direct_match, bbj_output_stem)

unmapped_leads <- lead[!koges_lead_stem %chin% map$koges_lead_stem]
cat(
  paste0("Leads with no Excel ", COHORT, " mapping:"),
  nrow(unmapped_leads),
  "variants /",
  uniqueN(unmapped_leads$koges_lead_stem),
  "stems\n"
)

bbj_cache <- new.env(parent = emptyenv())
get_bbj <- function(path) {
  if (!nzchar(path) || is.na(path)) return(NULL)
  if (exists(path, envir = bbj_cache, inherits = FALSE)) {
    return(get(path, envir = bbj_cache, inherits = FALSE))
  }
  cat("Reading ", COHORT, ": ", basename(path), "\n", sep = "")
  dt <- read_ext_cate(path)
  assign(path, dt, envir = bbj_cache)
  dt
}

per_match_list <- vector("list", nrow(match_jobs))
gws_hit_list <- list()
gws_counter <- 0L

for (j in seq_len(nrow(match_jobs))) {
  job <- match_jobs[j]
  if (is.na(job$bbj_output_stem) || !nzchar(job$bbj_output_stem)) {
    stats <- empty_match_stats()
    stats$replication_not_testable_reason <- "no_matching_bbj_phenotype"
    per_match_list[[j]] <- cbind(
      job[, .(
        lead_row_id, koges_phenotype, koges_pheno_code, koges_lead_stem,
        koges_lead_variant, koges_chr, koges_bp
      )],
      data.table(
        bbj_output_stem = NA_character_,
        matching_level = NA_character_,
        is_direct_match = NA,
        excel_external = NA_character_,
        excel_data_type = NA_character_,
        bbj_file = NA_character_,
        bbj_mapping_status = "no_matching_bbj_phenotype"
      ),
      stats
    )
    next
  }

  bbj_dt <- NULL
  if (job$bbj_mapping_status == "mapped_file_found") {
    bbj_dt <- get_bbj(job$bbj_file)
  }
  assessed <- assess_lead_in_bbj(job, bbj_dt, job$bbj_mapping_status)
  if (nrow(assessed$gws) > 0L) {
    gws_counter <- gws_counter + 1L
    gws_hit_list[[gws_counter]] <- assessed$gws
  }
  per_match_list[[j]] <- cbind(
    job[, .(
      lead_row_id, koges_phenotype, koges_pheno_code, koges_lead_stem,
      koges_lead_variant, koges_chr, koges_bp,
      bbj_output_stem, matching_level, is_direct_match,
      excel_external, excel_data_type, bbj_file, bbj_mapping_status
    )],
    assessed$stats
  )
}

per_match <- rbindlist(per_match_list, fill = TRUE, use.names = TRUE)

# -----------------------------------------------------------------------------
# Novelty
# -----------------------------------------------------------------------------
per_match[, novelty_assessable := (
  bbj_mapping_status == "mapped_file_found" &
    !is.na(bbj_n_variants_in_window) &
    bbj_n_variants_in_window > 0L
)]
per_match[, known_in_window := fcase(
  novelty_assessable & bbj_n_gws_hits_in_window > 0L, TRUE,
  novelty_assessable & bbj_n_gws_hits_in_window == 0L, FALSE,
  default = NA
)]

agg_scope <- function(pm, scope_direct) {
  sub <- if (isTRUE(scope_direct)) {
    pm[is_direct_match %in% TRUE]
  } else if (identical(scope_direct, FALSE)) {
    pm[is_direct_match %in% FALSE]
  } else {
    pm[!is.na(bbj_output_stem)]
  }
  sub[, .(
    n_mapped_files = uniqueN(bbj_output_stem[!is.na(bbj_output_stem)]),
    n_files_found = uniqueN(
      bbj_output_stem[bbj_mapping_status == "mapped_file_found"]
    ),
    mapped_stems = collapse_chr(bbj_output_stem),
    matching_levels = collapse_chr(matching_level),
    mapping_status = {
      st <- unique(bbj_mapping_status)
      if ("mapped_file_found" %chin% st) {
        "mapped_file_found"
      } else if ("mapped_bbj_file_missing" %chin% st) {
        "mapped_bbj_file_missing"
      } else if ("no_bbj_source" %chin% st) {
        "no_bbj_source"
      } else {
        "no_matching_bbj_phenotype"
      }
    },
    novelty_assessable = any(novelty_assessable %in% TRUE),
    n_gws_hits_in_window = {
      v <- bbj_n_gws_hits_in_window[novelty_assessable %in% TRUE]
      if (!length(v)) NA_integer_ else as.integer(sum(v, na.rm = TRUE))
    },
    # Multiple direct matches: any known -> known
    # else any novel -> novel (unassessable + novel still counts as novel)
    # all unassessable -> NA
    known_in_window = {
      n_known <- sum(known_in_window %in% TRUE)
      n_novel <- sum(known_in_window %in% FALSE)
      if (n_known > 0L) {
        TRUE
      } else if (n_novel > 0L) {
        FALSE
      } else {
        NA
      }
    },
    known_stems = collapse_chr(bbj_output_stem[known_in_window %in% TRUE]),
    n_variants_in_window = {
      v <- bbj_n_variants_in_window[novelty_assessable %in% TRUE]
      if (!length(v)) NA_integer_ else as.integer(max(v, na.rm = TRUE))
    },
    top_gws_p = {
      v <- bbj_top_gws_p[is.finite(bbj_top_gws_p)]
      if (!length(v)) NA_real_ else min(v)
    },
    top_gws_stem = {
      if (!any(is.finite(bbj_top_gws_p))) {
        NA_character_
      } else {
        bbj_output_stem[which.min(fifelse(is.finite(bbj_top_gws_p),
                                           bbj_top_gws_p, Inf))]
      }
    },
    exact_chrpos_present = any(bbj_exact_chrpos_present %in% TRUE),
    replication_testable = any(replication_testable %in% TRUE),
    testable_stems = collapse_chr(
      bbj_output_stem[replication_testable %in% TRUE]
    ),
    not_testable_reason = {
      if (any(replication_testable %in% TRUE)) {
        NA_character_
      } else {
        collapse_chr(unique(replication_not_testable_reason))
      }
    }
  ), by = lead_row_id]
}

agg_direct <- agg_scope(per_match, TRUE)
setnames(
  agg_direct,
  setdiff(names(agg_direct), "lead_row_id"),
  paste0(setdiff(names(agg_direct), "lead_row_id"), "_direct")
)
agg_other <- agg_scope(per_match, FALSE)
setnames(
  agg_other,
  setdiff(names(agg_other), "lead_row_id"),
  paste0(setdiff(names(agg_other), "lead_row_id"), "_other")
)
agg_any <- agg_scope(per_match, NA)
setnames(
  agg_any,
  setdiff(names(agg_any), "lead_row_id"),
  paste0(setdiff(names(agg_any), "lead_row_id"), "_any")
)

result <- merge(lead, agg_direct, by = "lead_row_id", all.x = TRUE)
result <- merge(result, agg_other, by = "lead_row_id", all.x = TRUE)
result <- merge(result, agg_any, by = "lead_row_id", all.x = TRUE)
setorder(result, lead_row_id)

fill_status <- function(x, default = "no_matching_bbj_phenotype") {
  x <- as.character(x)
  x[is.na(x) | !nzchar(x)] <- default
  x
}
result[, mapping_status_direct := fill_status(mapping_status_direct)]
result[, mapping_status_other := fill_status(mapping_status_other)]
result[, mapping_status_any := fill_status(mapping_status_any)]
result[is.na(n_mapped_files_direct), n_mapped_files_direct := 0L]
result[is.na(n_mapped_files_other), n_mapped_files_other := 0L]
result[is.na(n_mapped_files_any), n_mapped_files_any := 0L]
result[is.na(novelty_assessable_direct), novelty_assessable_direct := FALSE]
result[is.na(novelty_assessable_other), novelty_assessable_other := FALSE]
result[is.na(novelty_assessable_any), novelty_assessable_any := FALSE]
result[is.na(replication_testable_direct), replication_testable_direct := FALSE]
result[is.na(replication_testable_other), replication_testable_other := FALSE]
result[is.na(replication_testable_any), replication_testable_any := FALSE]

result[, novel_relative_to_bbj_direct := fcase(
  known_in_window_direct %in% TRUE, FALSE,
  known_in_window_direct %in% FALSE, TRUE,
  default = NA
)]
result[, novel_relative_to_bbj_other := fcase(
  known_in_window_other %in% TRUE, FALSE,
  known_in_window_other %in% FALSE, TRUE,
  default = NA
)]
result[, novel_relative_to_bbj_any := fcase(
  known_in_window_any %in% TRUE, FALSE,
  known_in_window_any %in% FALSE, TRUE,
  default = NA
)]

# NOTE: base::isTRUE() is NOT vectorized (always FALSE on length>1 columns).
# "Additionally known" = not already Known at direct match, but Known on other matching
# (includes direct Not assessed / Novel + other Known).
result[, additionally_known_from_other := (
  !(known_in_window_direct %in% TRUE) & (known_in_window_other %in% TRUE)
)]

result[, bbj_novelty_status_direct := fcase(
  known_in_window_direct %in% TRUE,
  paste0(
    "Known (direct match): GWS in ",
    fifelse(is.na(known_stems_direct), mapped_stems_direct, known_stems_direct)
  ),
  known_in_window_direct %in% FALSE,
  paste0(
    "Novel (direct match): no GWS in any direct match ",
    COHORT, " (unassessable files ignored)"
  ),
  mapping_status_direct == "no_matching_bbj_phenotype",
  paste0("Not assessed: no direct match ", COHORT, " phenotype"),
  mapping_status_direct %chin% c("mapped_bbj_file_missing", "no_bbj_source"),
  paste0("Not assessed: direct match ", COHORT, " source/file missing"),
  novelty_assessable_direct == FALSE,
  paste0(
    "Not assessed: no direct match ", COHORT,
    " basic-QC variants within +/-1 Mb"
  ),
  default = "Not assessed"
)]

result[, bbj_novelty_class := fcase(
  known_in_window_direct %in% TRUE,
  "known_direct",
  additionally_known_from_other == TRUE,
  "additionally_known_from_other",
  novel_relative_to_bbj_any %in% TRUE,
  "novel_after_all_matches",
  default = "not_assessed"
)]

result[, bbj_novelty_status := fcase(
  bbj_novelty_class == "known_direct",
  "Known at direct match",
  bbj_novelty_class == "additionally_known_from_other",
  paste0(
    "Additionally known from other matching: ",
    fifelse(is.na(known_stems_other), "", known_stems_other)
  ),
  bbj_novelty_class == "novel_after_all_matches",
  "Novel after direct match and other matching",
  mapping_status_any == "no_matching_bbj_phenotype",
  paste0("Not assessed: no matching ", COHORT, " phenotype"),
  mapping_status_any == "mapped_bbj_file_missing",
  paste0("Not assessed: mapped ", COHORT, " file missing"),
  default = "Not assessed"
)]

# -----------------------------------------------------------------------------
# Replication
# -----------------------------------------------------------------------------
n_global_mapped_direct <- result[mapping_status_direct == "mapped_file_found", .N]
n_global_testable_direct <- result[replication_testable_direct == TRUE, .N]
n_global_mapped_other <- result[mapping_status_other == "mapped_file_found", .N]
n_global_testable_other <- result[replication_testable_other == TRUE, .N]
n_global_mapped_any <- result[mapping_status_any == "mapped_file_found", .N]
n_global_testable_any <- result[replication_testable_any == TRUE, .N]

thr <- function(n) {
  if (is.na(n) || n <= 0L) NA_real_ else NOMINAL_REPLICATION_P / n
}

result[, n_trait_mapped_direct := sum(mapping_status_direct == "mapped_file_found"),
       by = koges_pheno_code]
result[, n_trait_testable_direct := sum(replication_testable_direct),
       by = koges_pheno_code]
result[, n_trait_mapped_other := sum(mapping_status_other == "mapped_file_found"),
       by = koges_pheno_code]
result[, n_trait_testable_other := sum(replication_testable_other),
       by = koges_pheno_code]
result[, n_trait_mapped_any := sum(mapping_status_any == "mapped_file_found"),
       by = koges_pheno_code]
result[, n_trait_testable_any := sum(replication_testable_any),
       by = koges_pheno_code]

result[, `:=`(
  replication_threshold_nominal = NOMINAL_REPLICATION_P,
  replication_threshold_global_mapped_direct = thr(n_global_mapped_direct),
  replication_threshold_global_testable_direct = thr(n_global_testable_direct),
  replication_threshold_trait_mapped_direct = thr(n_trait_mapped_direct),
  replication_threshold_trait_testable_direct = thr(n_trait_testable_direct),
  replication_threshold_global_mapped_other = thr(n_global_mapped_other),
  replication_threshold_global_testable_other = thr(n_global_testable_other),
  replication_threshold_trait_mapped_other = thr(n_trait_mapped_other),
  replication_threshold_trait_testable_other = thr(n_trait_testable_other)
)]

primary_thr_direct <- switch(
  PRIMARY_CORRECTION,
  global_mapped = "replication_threshold_global_mapped_direct",
  global_testable = "replication_threshold_global_testable_direct",
  trait_mapped = "replication_threshold_trait_mapped_direct",
  trait_testable = "replication_threshold_trait_testable_direct"
)
primary_thr_other <- switch(
  PRIMARY_CORRECTION,
  global_mapped = "replication_threshold_global_mapped_other",
  global_testable = "replication_threshold_global_testable_other",
  trait_mapped = "replication_threshold_trait_mapped_other",
  trait_testable = "replication_threshold_trait_testable_other"
)
result[, replication_threshold_primary_direct := get(primary_thr_direct)]
result[, replication_threshold_primary_other := get(primary_thr_other)]
result[, replication_primary_correction := PRIMARY_CORRECTION]

# attach thresholds onto per-match rows and call replication
per_match <- merge(
  per_match,
  result[, .(
    lead_row_id,
    koges_pheno_code,
    replication_threshold_nominal,
    replication_threshold_primary_direct,
    replication_threshold_primary_other
  )],
  by = "lead_row_id",
  all.x = TRUE
)
per_match[, primary_threshold := fifelse(
  is_direct_match %in% TRUE,
  replication_threshold_primary_direct,
  replication_threshold_primary_other
)]
per_match[, replicated_nominal := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    bbj_exact_p <= replication_threshold_nominal
)]
per_match[, replicated_primary := (
  replication_testable == TRUE &
    direction_concordant %in% TRUE &
    is.finite(primary_threshold) &
    bbj_exact_p <= primary_threshold
)]

rep_scope <- function(pm, scope_direct) {
  sub <- if (isTRUE(scope_direct)) {
    pm[is_direct_match %in% TRUE]
  } else {
    pm[is_direct_match %in% FALSE]
  }
  sub[, .(
    replicated_nominal = any(replicated_nominal %in% TRUE),
    replicated_primary = any(replicated_primary %in% TRUE),
    replicated_stems_nominal = collapse_chr(
      bbj_output_stem[replicated_nominal %in% TRUE]
    ),
    replicated_stems_primary = collapse_chr(
      bbj_output_stem[replicated_primary %in% TRUE]
    ),
    best_testable_p = {
      v <- bbj_exact_p[replication_testable %in% TRUE]
      if (!length(v)) NA_real_ else min(v, na.rm = TRUE)
    },
    direction_concordant = any(direction_concordant %in% TRUE)
  ), by = lead_row_id]
}

rep_direct <- rep_scope(per_match, TRUE)
setnames(
  rep_direct,
  setdiff(names(rep_direct), "lead_row_id"),
  paste0(setdiff(names(rep_direct), "lead_row_id"), "_direct")
)
rep_other <- rep_scope(per_match, FALSE)
setnames(
  rep_other,
  setdiff(names(rep_other), "lead_row_id"),
  paste0(setdiff(names(rep_other), "lead_row_id"), "_other")
)

result <- merge(result, rep_direct, by = "lead_row_id", all.x = TRUE)
result <- merge(result, rep_other, by = "lead_row_id", all.x = TRUE)
result[is.na(replicated_nominal_direct), replicated_nominal_direct := FALSE]
result[is.na(replicated_primary_direct), replicated_primary_direct := FALSE]
result[is.na(replicated_nominal_other), replicated_nominal_other := FALSE]
result[is.na(replicated_primary_other), replicated_primary_other := FALSE]
result[is.na(direction_concordant_direct), direction_concordant_direct := NA]
result[is.na(direction_concordant_other), direction_concordant_other := NA]

result[, additionally_replicated_nominal_from_other := (
  replicated_nominal_direct == FALSE & replicated_nominal_other == TRUE
)]
result[, additionally_replicated_primary_from_other := (
  replicated_primary_direct == FALSE & replicated_primary_other == TRUE
)]

result[, replication_class := fcase(
  replicated_primary_direct == TRUE,
  "replicated_direct",
  additionally_replicated_primary_from_other == TRUE,
  "additionally_replicated_from_other",
  replicated_nominal_direct == TRUE,
  "nominal_only_direct",
  additionally_replicated_nominal_from_other == TRUE,
  "additionally_replicated_nominal_from_other",
  replication_testable_direct == TRUE | replication_testable_other == TRUE,
  "not_replicated",
  default = "not_testable"
)]

result[, replication_status_direct := fcase(
  mapping_status_direct == "no_matching_bbj_phenotype",
  paste0("Not assessed: no direct match ", COHORT, " phenotype"),
  mapping_status_direct == "mapped_bbj_file_missing",
  paste0("Not assessed: direct match ", COHORT, " file missing"),
  replication_testable_direct == FALSE,
  paste0("Not testable (direct match): ", not_testable_reason_direct),
  replicated_primary_direct == TRUE,
  paste0(
    "Replicated (direct match): direction concordant and ",
    PRIMARY_CORRECTION, " Bonferroni significant"
  ),
  direction_concordant_direct %in% FALSE &
    is.finite(best_testable_p_direct) &
    best_testable_p_direct <= NOMINAL_REPLICATION_P,
  "Not replicated (direct match): nominally significant but opposite direction",
  replicated_nominal_direct == TRUE,
  paste0(
    "Nominal replication only (direct match): does not pass ",
    PRIMARY_CORRECTION, " Bonferroni threshold"
  ),
  default = "Not replicated (direct match)"
)]

result[, replication_status := fcase(
  replication_class == "replicated_direct",
  replication_status_direct,
  replication_class == "additionally_replicated_from_other",
  paste0(
    "Additionally replicated from other matching (",
    PRIMARY_CORRECTION, "): ",
    fifelse(is.na(replicated_stems_primary_other), "",
            replicated_stems_primary_other)
  ),
  replication_class == "nominal_only_direct",
  replication_status_direct,
  replication_class == "additionally_replicated_nominal_from_other",
  paste0(
    "Additionally replicated nominally from other matching: ",
    fifelse(is.na(replicated_stems_nominal_other), "",
            replicated_stems_nominal_other)
  ),
  replication_class == "not_replicated",
  "Not replicated after direct match and other matching",
  mapping_status_any == "no_matching_bbj_phenotype",
  paste0("Not assessed: no matching ", COHORT, " phenotype"),
  rep(TRUE, .N),
  paste0(
    "Not testable: ",
    fifelse(is.na(not_testable_reason_any), "", not_testable_reason_any)
  )
)]

# convenience aliases matching the Cont column names for the direct call
result[, novel_relative_to_bbj := novel_relative_to_bbj_direct]
result[, replicated_primary := replicated_primary_direct]
result[, replicated_nominal := replicated_nominal_direct]
result[, replication_testable := replication_testable_direct]
result[, mapped_for_replication := mapping_status_direct == "mapped_file_found"]

# Direct-match exact EAF on the lead-level table (Cont-compatible column name).
# At most one direct-match row per lead_row_id is expected.
exact_eaf_direct <- per_match[
  is_direct_match %in% TRUE,
  .(lead_row_id, bbj_exact_eaf_aligned_a1)
]
if (anyDuplicated(exact_eaf_direct$lead_row_id)) {
  stop(
    "Multiple direct-match rows per lead_row_id when attaching ",
    "bbj_exact_eaf_aligned_a1."
  )
}
if ("bbj_exact_eaf_aligned_a1" %chin% names(result)) {
  result[, bbj_exact_eaf_aligned_a1 := NULL]
}
result <- merge(
  result, exact_eaf_direct, by = "lead_row_id", all.x = TRUE, sort = FALSE
)

if (nrow(result) != nrow(lead)) {
  stop("Internal error: result row count differs from lead row count.")
}

# -----------------------------------------------------------------------------
# Write outputs
# -----------------------------------------------------------------------------
rename_cohort_prefix <- function(dt) {
  if (pfx == "bbj") return(dt)
  nms <- names(dt)
  new_nms <- nms
  new_nms <- sub("^bbj_", paste0(pfx, "_"), new_nms)
  new_nms <- gsub("bbj", pfx, new_nms, fixed = TRUE)
  if (any(duplicated(new_nms))) {
    stop("Column rename collision for cohort prefix: ", pfx)
  }
  setnames(dt, nms, new_nms)
  for (col in names(dt)) {
    if (is.character(dt[[col]])) {
      set(dt, j = col, value = gsub("bbj", pfx, dt[[col]], fixed = TRUE))
      set(dt, j = col, value = gsub("BBJ", COHORT, dt[[col]], fixed = TRUE))
    }
  }
  dt
}
result <- rename_cohort_prefix(result)
per_match <- rename_cohort_prefix(per_match)

stem_col <- paste0(pfx, "_output_stem")
nov_col <- paste0(pfx, "_novelty_class")
setorderv(per_match, c("lead_row_id", "is_direct_match", stem_col), c(1L, -1L, 1L))
fwrite(result, RESULT_FILE, sep = "\t", quote = FALSE, na = "NA")
fwrite(per_match, PER_MATCH_FILE, sep = "\t", quote = FALSE, na = "NA")

cat("\nCompleted.\n")
cat("Cohort:", COHORT, "\n")
cat("Lead-level result:", RESULT_FILE, "\n")
cat("Per-match result:", PER_MATCH_FILE, "\n")
cat("Novelty class:\n")
print(result[, .N, by = c(nov_col)])
cat("Replication class:\n")
print(result[, .N, by = replication_class])
