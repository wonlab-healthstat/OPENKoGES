#!/usr/bin/env Rscript

################################################################################
# Lead-level EUR panel evidence (UKB + FinnGen + MVP)
#
# External-only fixed-effect IVW beta meta for KoGES leads (RIN phenotypes).
# NOT a genome-wide meta-GWAS. KoGES discovery beta is NOT included in the meta.
#
# Inputs
#   novelty-replication TSVs (external_validation_KoGES_EUR_cont.R --cohort=UKB|FinnGen|MVP)
#   (generated here into Extended_Data/)
#     eur_external_phenotype_N.tsv
#     eur_panel_phenotype_eligibility.tsv
#
# Index conventions
#   Phenotype index: koges_pheno_code
#   Phenotype--lead index: lead_row_id (one KoGES phenotype--lead association)
#   Variant label: koges_lead_variant. The same variant in different phenotypes
#                  has distinct lead_row_id values and is a distinct test.
#   Cohort index: UKB / FinnGen / MVP. EUR_n_contributing_cohorts is the number
#                 of usable external cohorts for that one lead_row_id, not a
#                 count of phenotypes or unique variants.
#
# States
#   EUR_panel_available : k_eur >= 1 and phenotype is panel-meta eligible
#   EUR_panel_analyzed  : available AND allele anchors OK
#   EUR_panel_replication_* : TRUE/FALSE only if analyzed; else NA
#                             (blocked/unavailable is Not analyzed, never
#                              Not replicated)
#
# Cohort usable (defensive)
#   mapping_status == "mapped_file_found" AND
#   replication_testable AND exact_usable_for_replication AND
#   alleles present AND finite KoGES-aligned beta AND
#   finite beta/SE with SE > 0
#   (plus eligibility-authorized stem match after merge)
#
# Primary evidence (external-only IVW)
#   w_c = 1 / SE_c^2
#   beta_meta = sum(w_c * beta_c) / sum(w_c)
#   SE_meta   = 1 / sqrt(sum(w_c))
#   Z_meta    = beta_meta / SE_meta
#   P_meta    = 2 * Phi(-|Z_meta|)   # two-sided
#   Direction : sign(beta_KoGES) == sign(beta_meta)
# Heterogeneity (k>=2): Q, Q_P = Pr(Chi-square_df >= Q), df=k-1,
#   I^2 = max(0,(Q-df)/Q)*100
#   Q contribution: q_c = w_c*(beta_c - beta_meta)^2 ; frac = q_c/Q
#   Q-top cohort is reported only for k=3 (k=2 cannot identify a unique driver)
#   Leave-one-cohort-out (LOCO): IVW on C_i \ {c}; output includes the
#   standardized omitted-cohort residual and, for k=3, the minimum LOCO I^2
#   No k_i multiplicity correction (single pre-specified meta statistic)
#
# Sensitivity
#   sample-size weighted signed-Z meta when external N available
#
# Bonferroni
#   M_global = sum(EUR_panel_analyzed), i.e. the union of exact-testable,
#   allele-anchor-validated phenotype--lead pairs (lead_row_id) actually tested
#   in the EUR panel; it is not a count of unique variants or cohort rows.
#   alpha = 0.05 / M_global
#   Nominal: direction-concordant AND P_meta <= 0.05
#   Bonferroni: direction-concordant AND P_meta <= alpha
################################################################################
suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 0. Paths / options
# ------------------------------------------------------------------------------
ext_data <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
ukb_path <- Sys.getenv(
  "UKB_NOVELTY_PATH",
  file.path(ext_data, "KoGES_BOLT_leads_UKB_novelty_replication_cont.tsv")
)
finngen_path <- Sys.getenv(
  "FINNGEN_NOVELTY_PATH",
  file.path(ext_data, "KoGES_BOLT_leads_FinnGen_novelty_replication_cont.tsv")
)
mvp_path <- Sys.getenv(
  "MVP_NOVELTY_PATH",
  file.path(ext_data, "KoGES_BOLT_leads_MVP_novelty_replication_cont.tsv")
)
n_path <- Sys.getenv(
  "EUR_EXTERNAL_N_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eur_external_phenotype_N.tsv"
)
elig_path <- Sys.getenv(
  "EUR_PANEL_ELIGIBILITY_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eur_panel_phenotype_eligibility.tsv"
)
out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
)
out_path <- file.path(out_dir, "KoGES_BOLT_leads_EUR_panel_meta.tsv")
alpha_global <- 0.05
allele_mismatch_action <- tolower(Sys.getenv("EUR_ALLELE_MISMATCH_ACTION", "stop"))
if (!allele_mismatch_action %chin% c("stop", "warn")) {
  stop("EUR_ALLELE_MISMATCH_ACTION must be 'stop' or 'warn'.")
}
p_concordance_tol <- as.numeric(Sys.getenv("EUR_P_CONCORDANCE_TOL", "1e-3"))

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------
as_flag <- function(x) {
  toupper(trimws(as.character(x))) %chin% c("1", "TRUE", "T", "YES", "Y")
}
require_columns <- function(dt, required, label) {
  missing <- setdiff(required, names(dt))
  if (length(missing) > 0L) {
    stop(label, " missing column(s): ", paste(missing, collapse = ", "))
  }
}
safe_sign <- function(x) {
  ifelse(is.finite(x) & x != 0, sign(x), NA_real_)
}
z_to_p_twosided <- function(z) {
  ifelse(is.finite(z), 2 * pnorm(abs(z), lower.tail = FALSE), NA_real_)
}
neglog10p_from_z_twosided <- function(z) {
  ifelse(
    is.finite(z),
    -log10(2) - pnorm(abs(z), lower.tail = FALSE, log.p = TRUE) / log(10),
    NA_real_
  )
}
ivw_fixed <- function(betas, ses) {
  ok <- is.finite(betas) & is.finite(ses) & ses > 0
  betas <- betas[ok]
  ses <- ses[ok]
  k <- length(betas)
  if (k == 0L) {
    return(list(
      beta = NA_real_, se = NA_real_, z = NA_real_,
      q = NA_real_, q_p = NA_real_, df = NA_integer_, i2 = NA_real_, k = 0L
    ))
  }
  w <- 1 / (ses * ses)
  sw <- sum(w)
  b <- sum(w * betas) / sw
  se <- 1 / sqrt(sw)
  z <- b / se
  if (k >= 2L) {
    q <- sum(w * (betas - b)^2)
    df <- as.integer(k - 1L)
    q_p <- if (is.finite(q)) pchisq(q, df = df, lower.tail = FALSE) else NA_real_
    i2 <- if (is.finite(q) && q > 0) max(0, (q - df) / q) * 100 else 0
  } else {
    q <- NA_real_
    q_p <- NA_real_
    df <- NA_integer_
    i2 <- NA_real_
  }
  list(beta = b, se = se, z = z, q = q, q_p = q_p, df = df, i2 = i2, k = k)
}
# Per-cohort Q contributions + leave-one-cohort-out IVW (k >= 2)
ivw_q_contrib_loco <- function(betas, ses, labs) {
  stopifnot(length(betas) == length(ses), length(ses) == length(labs))
  empty_loco <- list(
    beta = NA_real_, se = NA_real_, z = NA_real_,
    p = NA_real_, q = NA_real_, i2 = NA_real_, k = NA_integer_,
    std_resid = NA_real_
  )
  out <- list(
    q_contrib = setNames(rep(NA_real_, length(labs)), labs),
    q_frac = setNames(rep(NA_real_, length(labs)), labs),
    q_top_cohort = NA_character_,
    loco = setNames(replicate(length(labs), empty_loco, simplify = FALSE), labs),
    loco_min_i2 = NA_real_,
    loco_min_i2_cohort = NA_character_,
    loco_i2_change_at_min = NA_real_
  )
  k <- length(betas)
  if (k < 2L) {
    return(out)
  }
  iv <- ivw_fixed(betas, ses)
  w <- 1 / (ses * ses)
  qc <- as.numeric(w * (betas - iv$beta)^2)
  names(qc) <- labs
  out$q_contrib <- qc
  if (is.finite(iv$q) && iv$q > 0) {
    out$q_frac <- qc / iv$q
  } else {
    out$q_frac <- setNames(rep(NA_real_, k), labs)
  }
  if (k >= 3L && is.finite(iv$q) && iv$q > 0) {
    out$q_top_cohort <- labs[[which.max(qc)]]
  }
  for (j in seq_len(k)) {
    omit <- labs[[j]]
    keep <- seq_len(k) != j
    ivo <- ivw_fixed(betas[keep], ses[keep])
    se_diff <- sqrt(ses[[j]]^2 + ivo$se^2)
    std_resid <- if (
      is.finite(betas[[j]]) && is.finite(ivo$beta) &&
        is.finite(se_diff) && se_diff > 0
    ) {
      (betas[[j]] - ivo$beta) / se_diff
    } else {
      NA_real_
    }
    out$loco[[omit]] <- list(
      beta = ivo$beta,
      se = ivo$se,
      z = ivo$z,
      p = z_to_p_twosided(ivo$z),
      q = ivo$q,
      i2 = ivo$i2,
      k = as.integer(ivo$k),
      std_resid = std_resid
    )
  }
  loco_i2 <- vapply(out$loco, function(x) x$i2, numeric(1))
  if (is.finite(iv$i2) && any(is.finite(loco_i2))) {
    jmin <- which.min(replace(loco_i2, !is.finite(loco_i2), Inf))
    out$loco_min_i2 <- as.numeric(loco_i2[[jmin]])
    out$loco_min_i2_cohort <- labs[[jmin]]
    out$loco_i2_change_at_min <- iv$i2 - out$loco_min_i2
  }
  out
}
# Legacy alias used by P-vs-Z concordance diagnostics
z_to_p <- z_to_p_twosided
p_display <- function(p) {
  fifelse(
    !is.finite(p),
    NA_character_,
    fifelse(
      p == 0,
      paste0("<", format(.Machine$double.xmin, scientific = TRUE, digits = 6)),
      format(p, scientific = TRUE, digits = 6)
    )
  )
}
assert_unique <- function(x, label) {
  if (anyDuplicated(x)) {
    stop(label, " is not unique.")
  }
}
prefix_to_cohort <- function(prefix) {
  switch(
    prefix,
    ukb = "UKB",
    finngen = "FinnGen",
    mvp = "MVP",
    stop("Unknown cohort prefix: ", prefix)
  )
}
pick_single <- function(use_ukb, use_finngen, use_mvp, v_ukb, v_finngen, v_mvp) {
  fifelse(
    use_ukb,
    v_ukb,
    fifelse(use_finngen, v_finngen, fifelse(use_mvp, v_mvp, NA_real_))
  )
}
pick_single_chr <- function(
    use_ukb, use_finngen, use_mvp, v_ukb, v_finngen, v_mvp) {
  fifelse(
    use_ukb,
    v_ukb,
    fifelse(use_finngen, v_finngen, fifelse(use_mvp, v_mvp, NA_character_))
  )
}
allele_ok_one <- function(use, a1, a2) {
  if (!any(use, na.rm = TRUE)) {
    return(NA)
  }
  a1u <- unique(a1[use])
  a2u <- unique(a2[use])
  length(a1u) == 1L &&
    length(a2u) == 1L &&
    !anyNA(a1u) &&
    !anyNA(a2u) &&
    nzchar(a1u[[1L]]) &&
    nzchar(a2u[[1L]])
}

# ------------------------------------------------------------------------------
# 2. Generate N + eligibility, then read + validate inputs
# ------------------------------------------------------------------------------
Sys.setenv(
  EUR_EXTERNAL_N_PATH = n_path,
  EUR_PANEL_ELIGIBILITY_PATH = elig_path
)
source(
  "/home2/sjim/openKoGES/CODE_manuscript_final/4_Novelty_Replication/forEURcohort/build_eur_external_sample_size.R",
  local = TRUE
)
n_table_present <- file.exists(n_path)
if (!n_table_present) {
  warning(
    "External N table not found after build: ", n_path,
    ". Primary IVW will run, but multi-cohort signed-Z sensitivity will be NA."
  )
}
if (!file.exists(elig_path)) {
  stop("Eligibility table not found after build: ", elig_path)
}
ukb <- fread(ukb_path, sep = "\t", na.strings = c("NA", ""))
finngen <- fread(finngen_path, sep = "\t", na.strings = c("NA", ""))
mvp <- fread(mvp_path, sep = "\t", na.strings = c("NA", ""))
n_tab <- if (n_table_present) {
  fread(n_path, sep = "\t", na.strings = c("NA", ""))
} else {
  data.table(
    output_stem = character(), ukb_n = numeric(), finngen_n = numeric(),
    mvp_n = numeric()
  )
}
elig <- fread(elig_path, sep = "\t", na.strings = c("NA", ""))
cat("UKB rows:", nrow(ukb), "\n")
cat("FinnGen rows:", nrow(finngen), "\n")
cat("MVP rows:", nrow(mvp), "\n")
key_cols <- c(
  "lead_row_id", "koges_pheno_code", "koges_lead_variant",
  "koges_chr", "koges_bp"
)
require_columns(ukb, key_cols, "UKB")
require_columns(finngen, key_cols, "FinnGen")
require_columns(mvp, key_cols, "MVP")
assert_unique(ukb$lead_row_id, "UKB lead_row_id")
assert_unique(finngen$lead_row_id, "FinnGen lead_row_id")
assert_unique(mvp$lead_row_id, "MVP lead_row_id")
assert_complete_identity <- function(dt, label, cols) {
  bad <- !complete.cases(dt[, ..cols])
  if (any(bad)) {
    stop(
      label, " has ", sum(bad), " row(s) with a missing lead identity field. ",
      "lead_row_id, phenotype, variant, chromosome, and position are required."
    )
  }
}
assert_complete_identity(ukb, "UKB", key_cols)
assert_complete_identity(finngen, "FinnGen", key_cols)
assert_complete_identity(mvp, "MVP", key_cols)
chk_ukb <- ukb[, .(
  lead_row_id,
  pheno = koges_pheno_code,
  rsid = koges_lead_variant,
  chr = koges_chr,
  bp = koges_bp
)]
chk_finngen <- finngen[, .(
  lead_row_id,
  pheno = koges_pheno_code,
  rsid = koges_lead_variant,
  chr = koges_chr,
  bp = koges_bp
)]
chk_mvp <- mvp[, .(
  lead_row_id,
  pheno = koges_pheno_code,
  rsid = koges_lead_variant,
  chr = koges_chr,
  bp = koges_bp
)]
chk <- merge(chk_ukb, chk_finngen, by = "lead_row_id",
             suffixes = c("_ukb", "_finngen"))
chk <- merge(chk, chk_mvp, by = "lead_row_id")
setnames(chk, c("pheno", "rsid", "chr", "bp"),
         c("pheno_mvp", "rsid_mvp", "chr_mvp", "bp_mvp"))
n_key_mismatch <- chk[
  pheno_ukb != pheno_finngen |
    pheno_ukb != pheno_mvp |
    rsid_ukb != rsid_finngen |
    rsid_ukb != rsid_mvp |
    chr_ukb != chr_finngen |
    chr_ukb != chr_mvp |
    bp_ukb != bp_finngen |
    bp_ukb != bp_mvp
]
if (nrow(n_key_mismatch) > 0L) {
  stop(
    "Lead identity mismatch across UKB/FinnGen/MVP for ",
    nrow(n_key_mismatch), " lead_row_id values. ",
    "Example lead_row_id: ", n_key_mismatch$lead_row_id[1L]
  )
}
if (nrow(chk) != nrow(ukb) || nrow(chk) != nrow(finngen) || nrow(chk) != nrow(mvp)) {
  stop("lead_row_id sets differ across UKB / FinnGen / MVP after join.")
}
# N table validation (wide on disk -> long lookup for joins)
require_columns(
  n_tab,
  c("output_stem", "ukb_n", "finngen_n", "mvp_n"),
  "N table"
)
assert_unique(n_tab$output_stem, "N table output_stem")
if ("koges_n" %chin% names(n_tab)) {
  stop("N table must not contain koges_n.")
}
n_tab[, `:=`(
  ukb_n = as.numeric(ukb_n),
  finngen_n = as.numeric(finngen_n),
  mvp_n = as.numeric(mvp_n)
)]
n_long <- rbindlist(list(
  n_tab[is.finite(ukb_n) & ukb_n > 0,
        .(cohort = "UKB", output_stem, external_n = ukb_n)],
  n_tab[is.finite(finngen_n) & finngen_n > 0,
        .(cohort = "FinnGen", output_stem, external_n = finngen_n)],
  n_tab[is.finite(mvp_n) & mvp_n > 0,
        .(cohort = "MVP", output_stem, external_n = mvp_n)]
))
assert_unique(
  paste(n_long$cohort, n_long$output_stem, sep = "\t"),
  "N long key (cohort, output_stem)"
)
require_columns(
  elig,
  c(
    "koges_pheno_code", "ukb_output_stem", "finngen_output_stem",
    "mvp_output_stem", "eur_panel_meta_eligible", "meta_exclusion_reason"
  ),
  "eligibility table"
)
assert_unique(elig$koges_pheno_code, "eligibility koges_pheno_code")
elig[, eur_panel_meta_eligible := as_flag(eur_panel_meta_eligible)]
elig[, meta_exclusion_reason := as.character(meta_exclusion_reason)]
default_excl_reason <-
  "No scientifically harmonized quantitative phenotype in the EUR panel"
warn_eligible_stem_missing_n <- function(stem_col, cohort_label) {
  stems <- unique(na.omit(as.character(
    elig[[stem_col]][elig$eur_panel_meta_eligible %in% TRUE]
  )))
  stems <- stems[nzchar(stems)]
  if (length(stems) == 0L) {
    return(invisible(NULL))
  }
  hit <- n_long[cohort == cohort_label & output_stem %chin% stems]
  miss <- setdiff(stems, hit$output_stem)
  if (length(miss) > 0L) {
    warning(
      "Eligible ", stem_col, " stem(s) missing finite ", cohort_label,
      " N in eur_external_phenotype_N.tsv: ", paste(miss, collapse = ", "),
      ". Primary IVW is unaffected; signed-Z sensitivity will be unavailable ",
      "for the affected multi-cohort lead(s)."
    )
  }
}
if (n_table_present) {
  warn_eligible_stem_missing_n("ukb_output_stem", "UKB")
  warn_eligible_stem_missing_n("finngen_output_stem", "FinnGen")
  warn_eligible_stem_missing_n("mvp_output_stem", "MVP")
}

# ------------------------------------------------------------------------------
# 3. Build wide table by explicit lead_row_id joins
# ------------------------------------------------------------------------------
base_cols <- c(
  "lead_row_id", "koges_phenotype", "koges_pheno_code",
  "koges_lead_variant", "koges_chr", "koges_bp"
)
require_columns(ukb, base_cols, "UKB base")
# Identity from UKB join key only. Do NOT take koges_beta_aligned_a1 from UKB
# alone: when UKB phenotype is unmapped that cell is often NA even if
# FinnGen/MVP have usable evidence for the same lead.
res <- ukb[, .(
  lead_row_id,
  koges_phenotype,
  koges_pheno_code,
  koges_lead_variant,
  koges_chr,
  koges_bp
)]
attach_cohort <- function(res, dt, prefix) {
  need <- c(
    "lead_row_id",
    "koges_beta_aligned_a1",
    paste0(prefix, "_output_stem"),
    paste0(prefix, "_mapping_status"),
    paste0(prefix, "_exact_a1_koges"),
    paste0(prefix, "_exact_a2_koges"),
    paste0(prefix, "_exact_beta_aligned_a1"),
    paste0(prefix, "_exact_se"),
    paste0(prefix, "_exact_p"),
    paste0(prefix, "_exact_eaf_aligned_a1"),
    paste0(prefix, "_exact_info"),
    paste0(prefix, "_exact_usable_for_replication"),
    "replication_testable"
  )
  require_columns(dt, need, paste0(prefix_to_cohort(prefix), " cohort block"))
  add <- dt[, ..need]
  setnames(
    add,
    old = need[-1L],
    new = c(
      "koges_beta_aligned_a1",
      "output_stem", "mapping_status", "exact_a1_koges", "exact_a2_koges",
      "beta", "se", "p", "eaf", "info", "exact_usable_for_replication",
      "replication_testable"
    )
  )
  add[, `:=`(
    mapping_status = as.character(mapping_status),
    output_stem = as.character(output_stem),
    replication_testable = as_flag(replication_testable),
    exact_usable_for_replication = as_flag(exact_usable_for_replication),
    koges_beta_aligned_a1 = as.numeric(koges_beta_aligned_a1),
    beta = as.numeric(beta),
    se = as.numeric(se),
    p = as.numeric(p),
    exact_a1_koges = toupper(trimws(as.character(exact_a1_koges))),
    exact_a2_koges = toupper(trimws(as.character(exact_a2_koges)))
  )]
  add <- merge(
    add,
    n_long[cohort == prefix_to_cohort(prefix),
           .(output_stem, ext_n = external_n)],
    by = "output_stem",
    all.x = TRUE,
    allow.cartesian = FALSE
  )
  add[, cohort_usable :=
    (mapping_status == "mapped_file_found") &
    (replication_testable %in% TRUE) &
    (exact_usable_for_replication %in% TRUE) &
    !is.na(exact_a1_koges) & nzchar(exact_a1_koges) &
    !is.na(exact_a2_koges) & nzchar(exact_a2_koges) &
    is.finite(koges_beta_aligned_a1) &
    is.finite(beta) &
    is.finite(se) &
    se > 0
  ]
  add[, z_from_beta_se := beta / se]
  add[, p_from_z := z_to_p(z_from_beta_se)]
  add[, p_vs_z_absdiff := fifelse(
    cohort_usable & is.finite(p) & is.finite(p_from_z),
    abs(p - p_from_z),
    NA_real_
  )]
  out <- add[, .(
    lead_row_id,
    koges_beta_aligned_a1,
    mapping_status,
    output_stem,
    replication_testable,
    exact_usable_for_replication,
    exact_a1_koges,
    exact_a2_koges,
    exact_beta_aligned_a1 = beta,
    exact_se = se,
    exact_p = p,
    exact_eaf_aligned_a1 = eaf,
    exact_info = info,
    external_n = ext_n,
    cohort_usable,
    z_from_beta_se,
    p_from_z,
    p_vs_z_absdiff
  )]
  setnames(
    out,
    old = setdiff(names(out), "lead_row_id"),
    new = paste0(prefix, "_", setdiff(names(out), "lead_row_id"))
  )
  merge(res, out, by = "lead_row_id", all.x = TRUE, allow.cartesian = FALSE)
}
res <- attach_cohort(res, ukb, "ukb")
res <- attach_cohort(res, finngen, "finngen")
res <- attach_cohort(res, mvp, "mvp")
assert_unique(res$lead_row_id, "EUR panel lead_row_id after cohort joins")
res <- merge(
  res,
  elig[, .(
    koges_pheno_code,
    elig_ukb_output_stem = ukb_output_stem,
    elig_finngen_output_stem = finngen_output_stem,
    elig_mvp_output_stem = mvp_output_stem,
    eur_panel_meta_eligible,
    meta_exclusion_reason
  )],
  by = "koges_pheno_code",
  all.x = TRUE,
  allow.cartesian = FALSE
)
assert_unique(res$lead_row_id, "EUR panel lead_row_id after eligibility join")
res[
  is.na(eur_panel_meta_eligible),
  `:=`(
    eur_panel_meta_eligible = FALSE,
    meta_exclusion_reason = default_excl_reason
  )
]
res[, meta_exclusion_reason := as.character(meta_exclusion_reason)]
res[
  eur_panel_meta_eligible == FALSE &
    (is.na(meta_exclusion_reason) |
       !nzchar(meta_exclusion_reason) |
       meta_exclusion_reason %chin% c("NA", "TRUE", "FALSE")),
  meta_exclusion_reason := default_excl_reason
]
res[
  eur_panel_meta_eligible == FALSE,
  `:=`(
    ukb_cohort_usable = FALSE,
    finngen_cohort_usable = FALSE,
    mvp_cohort_usable = FALSE
  )
]
stem_authorized <- function(mapped_stem, elig_stem) {
  !is.na(elig_stem) & nzchar(as.character(elig_stem)) &
    !is.na(mapped_stem) & nzchar(as.character(mapped_stem)) &
    as.character(mapped_stem) == as.character(elig_stem)
}
res[
  !stem_authorized(ukb_output_stem, elig_ukb_output_stem),
  ukb_cohort_usable := FALSE
]
res[
  !stem_authorized(finngen_output_stem, elig_finngen_output_stem),
  finngen_cohort_usable := FALSE
]
res[
  !stem_authorized(mvp_output_stem, elig_mvp_output_stem),
  mvp_cohort_usable := FALSE
]
n_usable_missing_n <- sum(
  (res$ukb_cohort_usable %in% TRUE &
     !(is.finite(res$ukb_external_n) & res$ukb_external_n > 0)) +
    (res$finngen_cohort_usable %in% TRUE &
       !(is.finite(res$finngen_external_n) & res$finngen_external_n > 0)) +
    (res$mvp_cohort_usable %in% TRUE &
       !(is.finite(res$mvp_external_n) & res$mvp_external_n > 0)),
  na.rm = TRUE
)
if (n_usable_missing_n > 0L) {
  cat(
    "Note: usable cohort rows with missing/nonpositive external N: ",
    n_usable_missing_n,
    " (IVW unaffected; blocks sensitivity signed-Z when k>=2)\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# 4. Panel aggregation
# ------------------------------------------------------------------------------
res[, EUR_n_contributing_cohorts := as.integer(
  (ukb_cohort_usable %in% TRUE) +
    (finngen_cohort_usable %in% TRUE) +
    (mvp_cohort_usable %in% TRUE)
)]
res[, EUR_contributing_cohorts := {
  labs <- c(
    if (isTRUE(ukb_cohort_usable)) "UKB",
    if (isTRUE(finngen_cohort_usable)) "FinnGen",
    if (isTRUE(mvp_cohort_usable)) "MVP"
  )
  if (length(labs) == 0L) NA_character_ else paste(labs, collapse = "+")
}, by = lead_row_id]
res[, EUR_allele_anchor_ok := {
  use <- c(
    isTRUE(ukb_cohort_usable),
    isTRUE(finngen_cohort_usable),
    isTRUE(mvp_cohort_usable)
  )
  a1 <- c(ukb_exact_a1_koges, finngen_exact_a1_koges, mvp_exact_a1_koges)
  a2 <- c(ukb_exact_a2_koges, finngen_exact_a2_koges, mvp_exact_a2_koges)
  allele_ok_one(use, a1, a2)
}, by = lead_row_id]
res[, EUR_missing_external_n := {
  miss <- FALSE
  if (isTRUE(ukb_cohort_usable)) {
    miss <- miss || !(is.finite(ukb_external_n) && ukb_external_n > 0)
  }
  if (isTRUE(finngen_cohort_usable)) {
    miss <- miss || !(is.finite(finngen_external_n) && finngen_external_n > 0)
  }
  if (isTRUE(mvp_cohort_usable)) {
    miss <- miss || !(is.finite(mvp_external_n) && mvp_external_n > 0)
  }
  miss
}, by = lead_row_id]
res[, `:=`(
  EUR_panel_available = EUR_n_contributing_cohorts >= 1L,
  EUR_panel_a1_koges = NA_character_,
  EUR_panel_a2_koges = NA_character_,
  EUR_evidence_type = NA_character_,
  EUR_panel_beta = NA_real_,
  EUR_panel_SE = NA_real_,
  EUR_panel_Z = NA_real_,
  EUR_panel_Z_rep = NA_real_,
  EUR_panel_P = NA_real_,
  EUR_panel_Q = NA_real_,
  EUR_panel_Q_P = NA_real_,
  EUR_panel_Q_df = NA_integer_,
  EUR_panel_I2 = NA_real_,
  EUR_Q_contrib_UKB = NA_real_,
  EUR_Q_contrib_FinnGen = NA_real_,
  EUR_Q_contrib_MVP = NA_real_,
  EUR_Q_frac_UKB = NA_real_,
  EUR_Q_frac_FinnGen = NA_real_,
  EUR_Q_frac_MVP = NA_real_,
  EUR_Q_top_cohort = NA_character_,
  EUR_LOCO_omit_UKB_beta = NA_real_,
  EUR_LOCO_omit_UKB_SE = NA_real_,
  EUR_LOCO_omit_UKB_Z = NA_real_,
  EUR_LOCO_omit_UKB_P = NA_real_,
  EUR_LOCO_omit_UKB_Q = NA_real_,
  EUR_LOCO_omit_UKB_I2 = NA_real_,
  EUR_LOCO_omit_UKB_std_resid = NA_real_,
  EUR_LOCO_omit_FinnGen_beta = NA_real_,
  EUR_LOCO_omit_FinnGen_SE = NA_real_,
  EUR_LOCO_omit_FinnGen_Z = NA_real_,
  EUR_LOCO_omit_FinnGen_P = NA_real_,
  EUR_LOCO_omit_FinnGen_Q = NA_real_,
  EUR_LOCO_omit_FinnGen_I2 = NA_real_,
  EUR_LOCO_omit_FinnGen_std_resid = NA_real_,
  EUR_LOCO_omit_MVP_beta = NA_real_,
  EUR_LOCO_omit_MVP_SE = NA_real_,
  EUR_LOCO_omit_MVP_Z = NA_real_,
  EUR_LOCO_omit_MVP_P = NA_real_,
  EUR_LOCO_omit_MVP_Q = NA_real_,
  EUR_LOCO_omit_MVP_I2 = NA_real_,
  EUR_LOCO_omit_MVP_std_resid = NA_real_,
  EUR_LOCO_min_I2 = NA_real_,
  EUR_LOCO_min_I2_cohort = NA_character_,
  EUR_LOCO_I2_change_at_min = NA_real_,
  EUR_panel_method_note = NA_character_,
  EUR_sensitivity_signedZ_Z = NA_real_,
  EUR_sensitivity_signedZ_P = NA_real_,
  EUR_sensitivity_signedZ_note = NA_character_
)]
# IVW needs beta/SE only; external N is not required for primary analysis
res[, EUR_panel_analyzed :=
  EUR_panel_available & (EUR_allele_anchor_ok %in% TRUE)
]
res[
  EUR_allele_anchor_ok == TRUE & ukb_cohort_usable == TRUE,
  `:=`(EUR_panel_a1_koges = ukb_exact_a1_koges,
       EUR_panel_a2_koges = ukb_exact_a2_koges)
]
res[
  EUR_allele_anchor_ok == TRUE & is.na(EUR_panel_a1_koges) &
    finngen_cohort_usable == TRUE,
  `:=`(EUR_panel_a1_koges = finngen_exact_a1_koges,
       EUR_panel_a2_koges = finngen_exact_a2_koges)
]
res[
  EUR_allele_anchor_ok == TRUE & is.na(EUR_panel_a1_koges) &
    mvp_cohort_usable == TRUE,
  `:=`(EUR_panel_a1_koges = mvp_exact_a1_koges,
       EUR_panel_a2_koges = mvp_exact_a2_koges)
]
res[
  eur_panel_meta_eligible == FALSE,
  EUR_evidence_type := paste0("Not analyzed (phenotype excluded: ",
                              meta_exclusion_reason, ")")
]
res[
  is.na(EUR_evidence_type) & EUR_n_contributing_cohorts == 0L,
  EUR_evidence_type := "Not available"
]
res[
  is.na(EUR_evidence_type) &
    EUR_panel_available == TRUE &
    EUR_allele_anchor_ok == FALSE,
  EUR_evidence_type := "Not analyzed (allele anchor mismatch)"
]
n_allele_mismatch <- sum(
  res$EUR_panel_available == TRUE & res$EUR_allele_anchor_ok == FALSE,
  na.rm = TRUE
)
if (n_allele_mismatch > 0L) {
  ex <- res[
    EUR_panel_available == TRUE & EUR_allele_anchor_ok == FALSE,
    .(lead_row_id, koges_pheno_code, koges_lead_variant,
      ukb_exact_a1_koges, ukb_exact_a2_koges,
      finngen_exact_a1_koges, finngen_exact_a2_koges,
      mvp_exact_a1_koges, mvp_exact_a2_koges)
  ][1L]
  msg <- paste0(
    "Allele anchor mismatch in ", n_allele_mismatch,
    " available lead(s). Example lead_row_id=", ex$lead_row_id,
    " pheno=", ex$koges_pheno_code,
    " variant=", ex$koges_lead_variant,
    ". Inspect allele columns before interpreting panel results."
  )
  if (identical(allele_mismatch_action, "stop")) {
    stop(msg)
  }
  warning(msg)
}
# Primary: external-only fixed-effect IVW (k==1 reduces to single-cohort IVW)
idx_ivw <- which(res$EUR_panel_analyzed == TRUE)
if (length(idx_ivw) > 0L) {
  res[
    idx_ivw,
    c(
      "EUR_evidence_type", "EUR_panel_beta", "EUR_panel_SE", "EUR_panel_Z",
      "EUR_panel_P", "EUR_panel_Q", "EUR_panel_Q_P", "EUR_panel_Q_df",
      "EUR_panel_I2",
      "EUR_panel_method_note"
    ) := {
      betas <- c(
        if (isTRUE(ukb_cohort_usable)) ukb_exact_beta_aligned_a1,
        if (isTRUE(finngen_cohort_usable)) finngen_exact_beta_aligned_a1,
        if (isTRUE(mvp_cohort_usable)) mvp_exact_beta_aligned_a1
      )
      ses <- c(
        if (isTRUE(ukb_cohort_usable)) ukb_exact_se,
        if (isTRUE(finngen_cohort_usable)) finngen_exact_se,
        if (isTRUE(mvp_cohort_usable)) mvp_exact_se
      )
      labs <- c(
        if (isTRUE(ukb_cohort_usable)) "UKB",
        if (isTRUE(finngen_cohort_usable)) "FinnGen",
        if (isTRUE(mvp_cohort_usable)) "MVP"
      )
      iv <- ivw_fixed(betas, ses)
      etype <- if (iv$k >= 2L) "Meta-analysis (IVW)" else "Single cohort (IVW)"
      note <- paste0(
        "external_IVW_fixed; cohorts=", paste(labs, collapse = "+"),
        "; k=", iv$k
      )
      list(
        etype, iv$beta, iv$se, iv$z,
        z_to_p_twosided(iv$z), iv$q, iv$q_p, iv$df, iv$i2, note
      )
    },
    by = lead_row_id
  ]
}
# Cohort Q contributions + leave-one-cohort-out (k >= 2)
idx_het <- which(
  res$EUR_panel_analyzed == TRUE & res$EUR_n_contributing_cohorts >= 2L
)
if (length(idx_het) > 0L) {
  res[
    idx_het,
    c(
      "EUR_Q_contrib_UKB", "EUR_Q_contrib_FinnGen", "EUR_Q_contrib_MVP",
      "EUR_Q_frac_UKB", "EUR_Q_frac_FinnGen", "EUR_Q_frac_MVP",
      "EUR_Q_top_cohort",
      "EUR_LOCO_omit_UKB_beta", "EUR_LOCO_omit_UKB_SE", "EUR_LOCO_omit_UKB_Z",
      "EUR_LOCO_omit_UKB_P", "EUR_LOCO_omit_UKB_Q", "EUR_LOCO_omit_UKB_I2",
      "EUR_LOCO_omit_UKB_std_resid",
      "EUR_LOCO_omit_FinnGen_beta", "EUR_LOCO_omit_FinnGen_SE",
      "EUR_LOCO_omit_FinnGen_Z", "EUR_LOCO_omit_FinnGen_P",
      "EUR_LOCO_omit_FinnGen_Q", "EUR_LOCO_omit_FinnGen_I2",
      "EUR_LOCO_omit_FinnGen_std_resid",
      "EUR_LOCO_omit_MVP_beta", "EUR_LOCO_omit_MVP_SE", "EUR_LOCO_omit_MVP_Z",
      "EUR_LOCO_omit_MVP_P", "EUR_LOCO_omit_MVP_Q", "EUR_LOCO_omit_MVP_I2",
      "EUR_LOCO_omit_MVP_std_resid",
      "EUR_LOCO_min_I2", "EUR_LOCO_min_I2_cohort",
      "EUR_LOCO_I2_change_at_min"
    ) := {
      betas <- c(
        if (isTRUE(ukb_cohort_usable)) ukb_exact_beta_aligned_a1,
        if (isTRUE(finngen_cohort_usable)) finngen_exact_beta_aligned_a1,
        if (isTRUE(mvp_cohort_usable)) mvp_exact_beta_aligned_a1
      )
      ses <- c(
        if (isTRUE(ukb_cohort_usable)) ukb_exact_se,
        if (isTRUE(finngen_cohort_usable)) finngen_exact_se,
        if (isTRUE(mvp_cohort_usable)) mvp_exact_se
      )
      labs <- c(
        if (isTRUE(ukb_cohort_usable)) "UKB",
        if (isTRUE(finngen_cohort_usable)) "FinnGen",
        if (isTRUE(mvp_cohort_usable)) "MVP"
      )
      ht <- ivw_q_contrib_loco(betas, ses, labs)
      loco_pack <- function(lab) {
        if (lab %chin% names(ht$loco)) {
          x <- ht$loco[[lab]]
          list(x$beta, x$se, x$z, x$p, x$q, x$i2, x$std_resid)
        } else {
          list(
            NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
            NA_real_
          )
        }
      }
      qv <- function(lab) {
        if (lab %chin% names(ht$q_contrib)) ht$q_contrib[[lab]] else NA_real_
      }
      qf <- function(lab) {
        if (lab %chin% names(ht$q_frac)) ht$q_frac[[lab]] else NA_real_
      }
      lu <- loco_pack("UKB")
      lf <- loco_pack("FinnGen")
      lm <- loco_pack("MVP")
      list(
        qv("UKB"), qv("FinnGen"), qv("MVP"),
        qf("UKB"), qf("FinnGen"), qf("MVP"),
        ht$q_top_cohort,
        lu[[1]], lu[[2]], lu[[3]], lu[[4]], lu[[5]], lu[[6]], lu[[7]],
        lf[[1]], lf[[2]], lf[[3]], lf[[4]], lf[[5]], lf[[6]], lf[[7]],
        lm[[1]], lm[[2]], lm[[3]], lm[[4]], lm[[5]], lm[[6]], lm[[7]],
        ht$loco_min_i2, ht$loco_min_i2_cohort, ht$loco_i2_change_at_min
      )
    },
    by = lead_row_id
  ]
}
res[, EUR_panel_koges_beta_aligned_a1 := pick_single(
  ukb_cohort_usable %in% TRUE,
  finngen_cohort_usable %in% TRUE,
  mvp_cohort_usable %in% TRUE,
  ukb_koges_beta_aligned_a1,
  finngen_koges_beta_aligned_a1,
  mvp_koges_beta_aligned_a1
)]
res[, EUR_koges_beta_consistent := {
  vals <- c(
    if (isTRUE(ukb_cohort_usable)) ukb_koges_beta_aligned_a1 else NA_real_,
    if (isTRUE(finngen_cohort_usable)) finngen_koges_beta_aligned_a1 else NA_real_,
    if (isTRUE(mvp_cohort_usable)) mvp_koges_beta_aligned_a1 else NA_real_
  )
  vals <- vals[is.finite(vals)]
  if (length(vals) <= 1L) {
    TRUE
  } else {
    diff(range(vals)) <= 1e-12
  }
}, by = lead_row_id]
n_koges_beta_mismatch <- sum(
  res$EUR_panel_analyzed == TRUE &
    res$EUR_n_contributing_cohorts >= 2L &
    res$EUR_koges_beta_consistent == FALSE,
  na.rm = TRUE
)
if (n_koges_beta_mismatch > 0L) {
  ex <- res[
    EUR_panel_analyzed == TRUE &
      EUR_n_contributing_cohorts >= 2L &
      EUR_koges_beta_consistent == FALSE,
    .(lead_row_id, koges_pheno_code, koges_lead_variant,
      ukb_koges_beta_aligned_a1, finngen_koges_beta_aligned_a1,
      mvp_koges_beta_aligned_a1)
  ][1L]
  stop(
    "Usable-cohort KoGES-aligned beta mismatch in ",
    n_koges_beta_mismatch, " analyzed multi-cohort lead(s). ",
    "Example lead_row_id=", ex$lead_row_id,
    " pheno=", ex$koges_pheno_code,
    " variant=", ex$koges_lead_variant,
    " ukb=", ex$ukb_koges_beta_aligned_a1,
    " finngen=", ex$finngen_koges_beta_aligned_a1,
    " mvp=", ex$mvp_koges_beta_aligned_a1
  )
}
n_missing_panel_koges_beta <- sum(
  res$EUR_panel_analyzed == TRUE &
    !is.finite(res$EUR_panel_koges_beta_aligned_a1),
  na.rm = TRUE
)
if (n_missing_panel_koges_beta > 0L) {
  ex <- res[
    EUR_panel_analyzed == TRUE & !is.finite(EUR_panel_koges_beta_aligned_a1),
    .(lead_row_id, koges_pheno_code, koges_lead_variant,
      EUR_contributing_cohorts,
      ukb_koges_beta_aligned_a1, finngen_koges_beta_aligned_a1,
      mvp_koges_beta_aligned_a1)
  ][1L]
  stop(
    "EUR_panel_koges_beta_aligned_a1 missing for ",
    n_missing_panel_koges_beta, " analyzed lead(s). ",
    "Example lead_row_id=", ex$lead_row_id,
    " pheno=", ex$koges_pheno_code,
    " cohorts=", ex$EUR_contributing_cohorts
  )
}
# Optional diagnostic: signed replication Z (not used for primary P)
res[
  EUR_panel_analyzed == TRUE,
  EUR_panel_Z_rep := safe_sign(EUR_panel_koges_beta_aligned_a1) * EUR_panel_Z
]
# Direction concordance required for replication
res[, EUR_panel_direction_concordant := NA]
res[
  EUR_panel_analyzed == TRUE,
  EUR_panel_direction_concordant := {
    ks <- safe_sign(EUR_panel_koges_beta_aligned_a1)
    ps <- safe_sign(EUR_panel_beta)
    is.finite(ks) & is.finite(ps) & (ks == ps)
  }
]
# LOCO direction flags are sensitivity diagnostics, not additional primary tests.
loco_direction_concordant <- function(koges_beta, loco_beta) {
  ks <- safe_sign(koges_beta)
  ls <- safe_sign(loco_beta)
  ifelse(is.finite(ks) & is.finite(ls), ks == ls, NA)
}
res[, c(
  "EUR_LOCO_omit_UKB_direction_concordant",
  "EUR_LOCO_omit_FinnGen_direction_concordant",
  "EUR_LOCO_omit_MVP_direction_concordant"
) := NA]
res[
  EUR_panel_analyzed == TRUE,
  c(
    "EUR_LOCO_omit_UKB_direction_concordant",
    "EUR_LOCO_omit_FinnGen_direction_concordant",
    "EUR_LOCO_omit_MVP_direction_concordant"
  ) := list(
    loco_direction_concordant(
      EUR_panel_koges_beta_aligned_a1, EUR_LOCO_omit_UKB_beta
    ),
    loco_direction_concordant(
      EUR_panel_koges_beta_aligned_a1, EUR_LOCO_omit_FinnGen_beta
    ),
    loco_direction_concordant(
      EUR_panel_koges_beta_aligned_a1, EUR_LOCO_omit_MVP_beta
    )
  )
]
res[, EUR_source_direction_pattern := NA_character_]
res[
  EUR_panel_analyzed == TRUE,
  EUR_source_direction_pattern := {
    ks <- safe_sign(EUR_panel_koges_beta_aligned_a1)
    conc <- c(
      if (isTRUE(ukb_cohort_usable)) {
        zs <- safe_sign(ukb_z_from_beta_se)
        is.finite(ks) & is.finite(zs) & (ks == zs)
      },
      if (isTRUE(finngen_cohort_usable)) {
        zs <- safe_sign(finngen_z_from_beta_se)
        is.finite(ks) & is.finite(zs) & (ks == zs)
      },
      if (isTRUE(mvp_cohort_usable)) {
        zs <- safe_sign(mvp_z_from_beta_se)
        is.finite(ks) & is.finite(zs) & (ks == zs)
      }
    )
    if (length(conc) == 0L || anyNA(conc)) {
      NA_character_
    } else if (all(conc)) {
      "All concordant"
    } else if (all(!conc)) {
      "All discordant"
    } else {
      "Mixed"
    }
  },
  by = lead_row_id
]
# Sensitivity: sample-size weighted signed-Z (requires external N when k>=2)
idx_sens <- which(
  res$EUR_panel_analyzed == TRUE & (
    res$EUR_n_contributing_cohorts == 1L |
      (res$EUR_n_contributing_cohorts >= 2L & !(res$EUR_missing_external_n %in% TRUE))
  )
)
if (length(idx_sens) > 0L) {
  res[
    idx_sens,
    c(
      "EUR_sensitivity_signedZ_Z",
      "EUR_sensitivity_signedZ_P",
      "EUR_sensitivity_signedZ_note"
    ) := {
      if (EUR_n_contributing_cohorts == 1L) {
        z <- pick_single(
          ukb_cohort_usable, finngen_cohort_usable, mvp_cohort_usable,
          ukb_z_from_beta_se, finngen_z_from_beta_se, mvp_z_from_beta_se
        )
        coh <- pick_single_chr(
          ukb_cohort_usable, finngen_cohort_usable, mvp_cohort_usable,
          "UKB", "FinnGen", "MVP"
        )
        note <- paste0("sensitivity_signed_Z_single; cohort=", coh)
      } else {
        num <- 0
        den <- 0
        if (isTRUE(ukb_cohort_usable)) {
          num <- num + sqrt(ukb_external_n) * ukb_z_from_beta_se
          den <- den + ukb_external_n
        }
        if (isTRUE(finngen_cohort_usable)) {
          num <- num + sqrt(finngen_external_n) * finngen_z_from_beta_se
          den <- den + finngen_external_n
        }
        if (isTRUE(mvp_cohort_usable)) {
          num <- num + sqrt(mvp_external_n) * mvp_z_from_beta_se
          den <- den + mvp_external_n
        }
        z <- num / sqrt(den)
        labs <- c(
          if (isTRUE(ukb_cohort_usable)) "UKB",
          if (isTRUE(finngen_cohort_usable)) "FinnGen",
          if (isTRUE(mvp_cohort_usable)) "MVP"
        )
        note <- paste0(
          "sensitivity_signed_Z_meta; cohorts=", paste(labs, collapse = "+")
        )
      }
      list(z, z_to_p_twosided(z), note)
    },
    by = lead_row_id
  ]
}
res[, EUR_panel_neglog10P := fifelse(
  EUR_panel_analyzed == TRUE,
  neglog10p_from_z_twosided(EUR_panel_Z),
  NA_real_
)]
res[, EUR_panel_P_display := fifelse(
  EUR_panel_analyzed == TRUE,
  p_display(EUR_panel_P),
  NA_character_
)]

# ------------------------------------------------------------------------------
# 5. Bonferroni on M_global = sum(EUR_panel_analyzed);
#    two-sided P + direction concordance
# ------------------------------------------------------------------------------
n_panel_available <- sum(res$EUR_panel_available, na.rm = TRUE)
n_panel_analyzed <- sum(res$EUR_panel_analyzed, na.rm = TRUE)
m_global <- n_panel_analyzed
if (m_global == 0L) {
  stop("No EUR_panel_analyzed leads; cannot define Bonferroni denominator M_global.")
}
alpha_bonf <- alpha_global / m_global
res[, `:=`(
  EUR_panel_N_available = as.integer(n_panel_available),
  EUR_panel_N_analyzed = as.integer(n_panel_analyzed),
  EUR_panel_M_global = as.integer(m_global),
  EUR_panel_alpha_Bonferroni = alpha_bonf,
  EUR_panel_replication_nominal = NA,
  EUR_panel_replication_Bonferroni = NA
)]
res[
  EUR_panel_analyzed == TRUE,
  `:=`(
    EUR_panel_replication_nominal =
      (EUR_panel_direction_concordant %in% TRUE) &
      is.finite(EUR_panel_P) &
      EUR_panel_P <= alpha_global,
    EUR_panel_replication_Bonferroni =
      (EUR_panel_direction_concordant %in% TRUE) &
      is.finite(EUR_panel_P) &
      EUR_panel_P <= alpha_bonf
  )
]
pdiff <- c(
  res$ukb_p_vs_z_absdiff,
  res$finngen_p_vs_z_absdiff,
  res$mvp_p_vs_z_absdiff
)
pdiff <- pdiff[is.finite(pdiff)]
if (length(pdiff) > 0L) {
  n_bad <- sum(pdiff > p_concordance_tol)
  cat(
    "P vs 2*Phi(-|beta/SE|) check: n_usable_comparisons=", length(pdiff),
    " max_absdiff=", format(max(pdiff), digits = 4),
    " n_absdiff>", p_concordance_tol, "=", n_bad, "\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# 6. Write outputs
# ------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setorder(res, lead_row_id)
fwrite(res, out_path, sep = "\t", quote = FALSE, na = "NA")
summary_dt <- data.table(
  metric = c(
    "n_leads",
    "n_EUR_panel_available",
    "n_EUR_panel_analyzed",
    "n_available_minus_analyzed",
    "M_global",
    "alpha_EUR_Bonferroni",
    "n_k0",
    "n_k1",
    "n_k2",
    "n_k3",
    "n_single_cohort_ivw",
    "n_meta_analysis_ivw",
    "n_not_available",
    "n_not_analyzed_allele_mismatch",
    "n_leads_phenotype_excluded",
    "n_phenotypes_excluded",
    "n_panel_direction_concordant",
    "n_all_source_direction_concordant",
    "n_mixed_source_direction",
    "n_nominal_replicated",
    "n_bonferroni_replicated",
    "n_analyzed_not_replicated_nominal",
    "n_sensitivity_signedZ",
    "n_signedZ_nominal_replicated",
    "n_signedZ_bonferroni_replicated",
    "n_meta_Q_P_lt_0.05",
    "n_meta_Q_P_lt_0.001",
    "n_meta_I2_ge_50",
    "n_meta_I2_ge_75",
    "n_Q_top_UKB_among_k3_I2_ge_50",
    "n_Q_top_FinnGen_among_k3_I2_ge_50",
    "n_Q_top_MVP_among_k3_I2_ge_50",
    "n_LOCO_min_I2_UKB_among_k3_I2_ge_50",
    "n_LOCO_min_I2_FinnGen_among_k3_I2_ge_50",
    "n_LOCO_min_I2_MVP_among_k3_I2_ge_50"
  ),
  value = {
    meta_rows <- res$EUR_evidence_type == "Meta-analysis (IVW)"
    i2_ge50 <- meta_rows &
      is.finite(res$EUR_panel_I2) & res$EUR_panel_I2 >= 50
    i2_ge75 <- meta_rows &
      is.finite(res$EUR_panel_I2) & res$EUR_panel_I2 >= 75
    k3_i2_ge50 <- i2_ge50 & res$EUR_n_contributing_cohorts == 3L
    phenotype_excluded <- res$eur_panel_meta_eligible %in% FALSE
    signedZ_available <- is.finite(res$EUR_sensitivity_signedZ_Z) &
      is.finite(res$EUR_sensitivity_signedZ_P)
    signedZ_dir <- {
      ks <- safe_sign(res$EUR_panel_koges_beta_aligned_a1)
      zs <- safe_sign(res$EUR_sensitivity_signedZ_Z)
      is.finite(ks) & is.finite(zs) & (ks == zs)
    }
    c(
    nrow(res),
    n_panel_available,
    n_panel_analyzed,
    n_panel_available - n_panel_analyzed,
    m_global,
    alpha_bonf,
    sum(res$EUR_n_contributing_cohorts == 0L, na.rm = TRUE),
    sum(res$EUR_n_contributing_cohorts == 1L, na.rm = TRUE),
    sum(res$EUR_n_contributing_cohorts == 2L, na.rm = TRUE),
    sum(res$EUR_n_contributing_cohorts == 3L, na.rm = TRUE),
    sum(res$EUR_evidence_type == "Single cohort (IVW)", na.rm = TRUE),
    sum(res$EUR_evidence_type == "Meta-analysis (IVW)", na.rm = TRUE),
    sum(res$EUR_evidence_type == "Not available", na.rm = TRUE),
    sum(res$EUR_evidence_type == "Not analyzed (allele anchor mismatch)", na.rm = TRUE),
    sum(phenotype_excluded, na.rm = TRUE),
    uniqueN(res$koges_pheno_code[phenotype_excluded]),
    sum(res$EUR_panel_direction_concordant %in% TRUE, na.rm = TRUE),
    sum(res$EUR_source_direction_pattern == "All concordant", na.rm = TRUE),
    sum(res$EUR_source_direction_pattern == "Mixed", na.rm = TRUE),
    sum(res$EUR_panel_replication_nominal %in% TRUE, na.rm = TRUE),
    sum(res$EUR_panel_replication_Bonferroni %in% TRUE, na.rm = TRUE),
    sum(res$EUR_panel_analyzed %in% TRUE &
          res$EUR_panel_replication_nominal %in% FALSE, na.rm = TRUE),
    sum(signedZ_available, na.rm = TRUE),
    sum(
      signedZ_available & signedZ_dir &
        res$EUR_sensitivity_signedZ_P <= alpha_global,
      na.rm = TRUE
    ),
    sum(
      signedZ_available & signedZ_dir &
        res$EUR_sensitivity_signedZ_P <= alpha_bonf,
      na.rm = TRUE
    ),
    sum(meta_rows & is.finite(res$EUR_panel_Q_P) &
          res$EUR_panel_Q_P < 0.05, na.rm = TRUE),
    sum(meta_rows & is.finite(res$EUR_panel_Q_P) &
          res$EUR_panel_Q_P < 0.001, na.rm = TRUE),
    sum(i2_ge50, na.rm = TRUE),
    sum(i2_ge75, na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_Q_top_cohort == "UKB", na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_Q_top_cohort == "FinnGen", na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_Q_top_cohort == "MVP", na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_LOCO_min_I2_cohort == "UKB", na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_LOCO_min_I2_cohort == "FinnGen", na.rm = TRUE),
    sum(k3_i2_ge50 & res$EUR_LOCO_min_I2_cohort == "MVP", na.rm = TRUE)
  )}
)
cat("Wrote:", out_path, "\n")
cat("EUR_panel_available =", n_panel_available, "\n")
cat("EUR_panel_analyzed  =", n_panel_analyzed, "\n")
cat("M_global            =", m_global, "\n")
cat("alpha_EUR_Bonferroni =", format(alpha_bonf, digits = 6), "\n")
print(summary_dt)
if (6L %in% res$lead_row_id) {
  cat("\n--- Spot check: lead_row_id = 6 ---\n")
  print(res[lead_row_id == 6L, .(
    koges_pheno_code,
    koges_lead_variant,
    EUR_n_contributing_cohorts,
    EUR_contributing_cohorts,
    EUR_panel_available,
    EUR_panel_analyzed,
    EUR_evidence_type,
    EUR_panel_koges_beta_aligned_a1,
    EUR_panel_beta,
    EUR_panel_SE,
    EUR_panel_Z,
    EUR_panel_Z_rep,
    EUR_panel_P,
    EUR_panel_Q,
    EUR_panel_Q_P,
    EUR_panel_I2,
    EUR_Q_contrib_UKB,
    EUR_Q_contrib_FinnGen,
    EUR_Q_contrib_MVP,
    EUR_Q_frac_UKB,
    EUR_Q_frac_FinnGen,
    EUR_Q_frac_MVP,
    EUR_Q_top_cohort,
    EUR_LOCO_omit_UKB_I2,
    EUR_LOCO_omit_UKB_std_resid,
    EUR_LOCO_omit_UKB_direction_concordant,
    EUR_LOCO_omit_FinnGen_I2,
    EUR_LOCO_omit_FinnGen_std_resid,
    EUR_LOCO_omit_FinnGen_direction_concordant,
    EUR_LOCO_omit_MVP_I2,
    EUR_LOCO_omit_MVP_std_resid,
    EUR_LOCO_omit_MVP_direction_concordant,
    EUR_LOCO_min_I2,
    EUR_LOCO_min_I2_cohort,
    EUR_LOCO_I2_change_at_min,
    EUR_panel_P_display,
    EUR_panel_neglog10P,
    EUR_panel_direction_concordant,
    EUR_source_direction_pattern,
    EUR_panel_replication_nominal,
    EUR_panel_replication_Bonferroni,
    EUR_panel_method_note,
    EUR_sensitivity_signedZ_Z,
    EUR_sensitivity_signedZ_P
  )])
}
