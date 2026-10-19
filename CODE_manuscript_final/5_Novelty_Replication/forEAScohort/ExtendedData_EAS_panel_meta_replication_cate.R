#!/usr/bin/env Rscript
################################################################################
# Lead-level EAS panel evidence (BBJ + TPMI), binary / Cate
#
# Cont counterpart: ExtendedData_EAS_panel_meta_replication_cont.R
# TWB has no binary matching; BBJ / TPMI only.
#
# External-only fixed-effect IVW beta meta for KoGES SAIGE leads.
# NOT a genome-wide meta-GWAS. KoGES discovery beta is NOT included.
#
# Inputs
#   BBJ / TPMI novelty TSV + per_match TSV
#     (external_validation_KoGES_EAS_cate.R --cohort=BBJ|TPMI)
#   KoGES_EAS_matching_binary_SAIGE_checked.xlsx
#     case/control (effective N) and direct-match eligibility.
#     Other matching (e.g. related diseases) is excluded from panel meta.
#   Signed-Z sensitivity (METAL sample-size method; not primary IVW):
#     Z_i = sign(beta_i) * Phi^{-1}(1 - P_i/2)  from reported exact P
#     N_eff = 4 * n_case * n_control / (n_case + n_control)
#     Z_meta = sum(sqrt(N_eff,i) * Z_i) / sqrt(sum(N_eff,i))
#   Primary evidence remains external-only fixed-effect IVW on beta/SE.
#
# Index conventions
#   Phenotype index: koges_lead_stem (B_CAD, B_BRCA, ...)
#   Phenotype--lead index: lead_row_id
#   Cohort index: BBJ / TPMI. Max k = 2.
#
# States
#   EAS_panel_available : k >= 1 and phenotype is panel-meta eligible
#   EAS_panel_analyzed  : available AND allele anchors OK
#   EAS_panel_replication_* : TRUE/FALSE only if analyzed; else NA
#
# Cohort usable
#   Direct-match rows with mapping_status == "mapped_file_found" AND
#   replication_testable AND exact_usable_for_replication AND
#   alleles present AND finite KoGES-aligned beta AND
#   finite beta/SE with SE > 0
#   (plus eligibility-authorized stem)
#
# Primary evidence (external-only IVW)
#   w_c = 1 / SE_c^2
#   beta_meta = sum(w_c * beta_c) / sum(w_c)
#   SE_meta   = 1 / sqrt(sum(w_c))
#   Z_meta    = beta_meta / SE_meta
#   P_meta    = 2 * Phi(-|Z_meta|)
#   Direction : sign(beta_KoGES) == sign(beta_meta)
# Heterogeneity (k>=2): Q, Q_P, df=k-1, I^2
#   Q-top cohort is reported only for k=3; with BBJ+TPMI max k=2 so Q-top is NA
#   LOCO: IVW on the remaining cohort when k=2
#
# Bonferroni
#   M_global = sum(EAS_panel_analyzed)
#   alpha = 0.05 / M_global
################################################################################
suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 0. Paths / options
# ------------------------------------------------------------------------------
bbj_path <- Sys.getenv(
  "BBJ_NOVELTY_PATH",
  paste0(
    "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/",
    "KoGES_SAIGE_leads_BBJ_novelty_replication_cate.tsv"
  )
)
bbj_pm_path <- Sys.getenv(
  "BBJ_PER_MATCH_PATH",
  paste0(
    "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/",
    "KoGES_SAIGE_leads_BBJ_per_match_cate.tsv"
  )
)
tpmi_path <- Sys.getenv(
  "TPMI_NOVELTY_PATH",
  paste0(
    "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/",
    "KoGES_SAIGE_leads_TPMI_novelty_replication_cate.tsv"
  )
)
tpmi_pm_path <- Sys.getenv(
  "TPMI_PER_MATCH_PATH",
  paste0(
    "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/",
    "KoGES_SAIGE_leads_TPMI_per_match_cate.tsv"
  )
)
matching_xlsx <- Sys.getenv(
  "MATCHING_XLSX",
  paste0(
    "/data4/sjim/openKoGES_revision/Sumstats/EAS/",
    "KoGES_EAS_matching_binary_SAIGE_checked.xlsx"
  )
)
n_path <- Sys.getenv(
  "EAS_EXTERNAL_N_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eas_external_phenotype_N_cate.tsv"
)
elig_path <- Sys.getenv(
  "EAS_PANEL_ELIGIBILITY_PATH",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/eas_panel_phenotype_eligibility_cate.tsv"
)
out_dir <- Sys.getenv(
  "OUTPUT_DIR",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
)
out_path <- file.path(out_dir, "KoGES_SAIGE_leads_EAS_panel_meta.tsv")
summary_path <- file.path(out_dir, "KoGES_SAIGE_leads_EAS_panel_meta_summary.tsv")
alpha_global <- 0.05
allele_mismatch_action <- tolower(Sys.getenv("EAS_ALLELE_MISMATCH_ACTION", "stop"))
if (!allele_mismatch_action %chin% c("stop", "warn")) {
  stop("EAS_ALLELE_MISMATCH_ACTION must be 'stop' or 'warn'.")
}
p_concordance_tol <- as.numeric(Sys.getenv("EAS_P_CONCORDANCE_TOL", "1e-3"))
DIRECT_MATCH_LABEL <- Sys.getenv("DIRECT_MATCH_LABEL", "direct match")
DIRECT_MATCH_LABELS <- unique(c(DIRECT_MATCH_LABEL, "direct match", "직접 매칭"))
EMPTY_TOKENS <- c("", "-", "0", "x", "na", "none", "nan")

need_files <- c(bbj_path, bbj_pm_path, tpmi_path, tpmi_pm_path, matching_xlsx)
miss_files <- need_files[!file.exists(need_files)]
if (length(miss_files)) {
  stop(
    "Required input missing:\n  ",
    paste(miss_files, collapse = "\n  ")
  )
}
if (!requireNamespace("readxl", quietly = TRUE)) {
  stop("package 'readxl' is required")
}

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------
as_flag <- function(x) {
  toupper(trimws(as.character(x))) %chin% c("1", "TRUE", "T", "YES", "Y")
}
is_empty <- function(x) {
  x <- tolower(trimws(as.character(x)))
  is.na(x) | !nzchar(x) | x %chin% EMPTY_TOKENS
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
# METAL-style signed Z from reported two-sided P and effect direction.
signed_z_from_p <- function(beta, p) {
  beta <- as.numeric(beta)
  p <- as.numeric(p)
  ok <- is.finite(beta) & beta != 0 & is.finite(p) & p >= 0 & p <= 1
  out <- rep(NA_real_, length(beta))
  if (!any(ok)) {
    return(out)
  }
  p_safe <- pmax(p[ok], .Machine$double.xmin)
  out[ok] <- sign(beta[ok]) * qnorm(p_safe / 2, lower.tail = FALSE)
  out
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
  switch(prefix, bbj = "BBJ", tpmi = "TPMI", stop("Unknown cohort prefix: ", prefix))
}
pick_single <- function(use_bbj, use_tpmi, v_bbj, v_tpmi) {
  fifelse(use_bbj, v_bbj, fifelse(use_tpmi, v_tpmi, NA_real_))
}
pick_single_chr <- function(use_bbj, use_tpmi, v_bbj, v_tpmi) {
  fifelse(use_bbj, v_bbj, fifelse(use_tpmi, v_tpmi, NA_character_))
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
collapse_chr <- function(x) {
  x <- unique(x[!is.na(x) & nzchar(as.character(x))])
  if (!length(x)) return(NA_character_)
  paste(x, collapse = ";")
}
numeric_or_na <- function(x) {
  suppressWarnings(as.numeric(x))
}
# Binary effective sample size for signed-Z weights (METAL / Willer et al.).
neff_from_case_control <- function(ncase, nctrl) {
  ncase <- numeric_or_na(ncase)
  nctrl <- numeric_or_na(nctrl)
  fifelse(
    is.finite(ncase) & is.finite(nctrl) & ncase > 0 & nctrl > 0,
    4 * ncase * nctrl / (ncase + nctrl),
    NA_real_
  )
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
stem_in_list <- function(mapped_stem, elig_stems) {
  if (is.na(mapped_stem) || !nzchar(as.character(mapped_stem))) return(FALSE)
  if (is.na(elig_stems) || !nzchar(as.character(elig_stems))) return(FALSE)
  mapped_stem %chin% strsplit(as.character(elig_stems), ";", fixed = TRUE)[[1L]]
}

# ------------------------------------------------------------------------------
# 2. Read novelty leads (identity) + per-match (exact stats)
# ------------------------------------------------------------------------------
bbj <- fread(bbj_path, sep = "\t", na.strings = c("NA", ""))
tpmi <- fread(tpmi_path, sep = "\t", na.strings = c("NA", ""))
bbj_pm <- fread(bbj_pm_path, sep = "\t", na.strings = c("NA", ""))
tpmi_pm <- fread(tpmi_pm_path, sep = "\t", na.strings = c("NA", ""))
cat("BBJ lead rows:", nrow(bbj), "\n")
cat("TPMI lead rows:", nrow(tpmi), "\n")

key_cols <- c(
  "lead_row_id", "koges_pheno_code", "koges_lead_variant",
  "koges_chr", "koges_bp", "koges_lead_stem"
)
require_columns(bbj, c(key_cols, "koges_phenotype"), "BBJ lead")
require_columns(tpmi, key_cols, "TPMI lead")
assert_unique(bbj$lead_row_id, "BBJ lead_row_id")
assert_unique(tpmi$lead_row_id, "TPMI lead_row_id")
if (anyNA(bbj[, ..key_cols]) || anyNA(tpmi[, ..key_cols])) {
  stop("Lead identity fields contain missing values.")
}
chk <- merge(
  bbj[, .(lead_row_id, pheno = koges_pheno_code, rsid = koges_lead_variant,
          chr = koges_chr, bp = koges_bp, stem = koges_lead_stem)],
  tpmi[, .(lead_row_id, pheno = koges_pheno_code, rsid = koges_lead_variant,
           chr = koges_chr, bp = koges_bp, stem = koges_lead_stem)],
  by = "lead_row_id",
  suffixes = c("_bbj", "_tpmi")
)
n_key_mismatch <- chk[
  pheno_bbj != pheno_tpmi |
    rsid_bbj != rsid_tpmi |
    chr_bbj != chr_tpmi |
    bp_bbj != bp_tpmi |
    stem_bbj != stem_tpmi
]
if (nrow(n_key_mismatch) > 0L) {
  stop(
    "Lead identity mismatch across BBJ/TPMI for ",
    nrow(n_key_mismatch), " lead_row_id values. Example: ",
    n_key_mismatch$lead_row_id[1L]
  )
}
if (nrow(chk) != nrow(bbj) || nrow(chk) != nrow(tpmi)) {
  stop("lead_row_id sets differ across BBJ / TPMI after join.")
}

# Per-match has at most one direct-match row per lead_row_id (verified).
filter_direct <- function(pm, prefix, label) {
  stem_col <- paste0(prefix, "_output_stem")
  map_col <- paste0(prefix, "_mapping_status")
  p_col <- paste0(prefix, "_exact_p")
  need <- c(
    "lead_row_id", "is_direct_match", "replication_testable",
    "koges_beta_aligned_a1",
    stem_col, map_col,
    paste0(prefix, "_exact_a1_koges"),
    paste0(prefix, "_exact_a2_koges"),
    paste0(prefix, "_exact_beta_aligned_a1"),
    paste0(prefix, "_exact_se"),
    p_col,
    paste0(prefix, "_exact_eaf_aligned_a1"),
    paste0(prefix, "_exact_info"),
    paste0(prefix, "_exact_usable_for_replication")
  )
  require_columns(pm, need, paste0(label, " per-match"))
  x <- copy(pm)
  if ("koges_pheno_code.x" %chin% names(x) && !"koges_pheno_code" %chin% names(x)) {
    setnames(x, "koges_pheno_code.x", "koges_pheno_code")
  }
  x <- x[is_direct_match %in% TRUE]
  if (!nrow(x)) {
    return(x[0L])
  }
  assert_unique(x$lead_row_id, paste0(label, " direct-match lead_row_id"))
  x
}

bbj_one <- filter_direct(bbj_pm, "bbj", "BBJ")
tpmi_one <- filter_direct(tpmi_pm, "tpmi", "TPMI")
cat("BBJ direct-match leads:", nrow(bbj_one), "/", nrow(bbj), "\n")
cat("TPMI direct-match leads:", nrow(tpmi_one), "/", nrow(tpmi), "\n")

# ------------------------------------------------------------------------------
# 3. Excel -> effective-N table + eligibility (direct match, BBJ/TPMI only)
# ------------------------------------------------------------------------------
meta <- as.data.table(readxl::read_excel(matching_xlsx, sheet = 1))
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
need_meta <- c("KoGES_VAR", match_level_col, "BBJ", "TPMI")
if (length(miss_meta <- setdiff(need_meta, names(meta)))) {
  stop("Excel missing columns: ", paste(miss_meta, collapse = ", "))
}
need_cc <- c("BBJ_case", "BBJ_control", "TPMI_case", "TPMI_control")
if (length(miss_cc <- setdiff(need_cc, names(meta)))) {
  stop("Excel missing case/control columns: ", paste(miss_cc, collapse = ", "))
}
meta[, `:=`(
  KoGES_VAR = trimws(as.character(KoGES_VAR)),
  matching_level = trimws(as.character(get(match_level_col))),
  BBJ = trimws(as.character(BBJ)),
  TPMI = trimws(as.character(TPMI)),
  bbj_case = numeric_or_na(BBJ_case),
  bbj_control = numeric_or_na(BBJ_control),
  tpmi_case = numeric_or_na(TPMI_case),
  tpmi_control = numeric_or_na(TPMI_control),
  bbj_n = neff_from_case_control(BBJ_case, BBJ_control),
  tpmi_n = neff_from_case_control(TPMI_case, TPMI_control)
)]
meta[, is_direct_match := matching_level %chin% DIRECT_MATCH_LABELS]
n_tab <- unique(meta[!is_empty(KoGES_VAR), .(
  output_stem = KoGES_VAR,
  bbj_case,
  bbj_control,
  bbj_n,
  tpmi_case,
  tpmi_control,
  tpmi_n
)])
n_tab <- n_tab[, .(
  bbj_case = {
    v <- bbj_case[is.finite(bbj_case) & bbj_case > 0]
    if (length(v)) v[[1L]] else NA_real_
  },
  bbj_control = {
    v <- bbj_control[is.finite(bbj_control) & bbj_control > 0]
    if (length(v)) v[[1L]] else NA_real_
  },
  bbj_n = {
    v <- bbj_n[is.finite(bbj_n) & bbj_n > 0]
    if (length(v)) v[[1L]] else NA_real_
  },
  tpmi_case = {
    v <- tpmi_case[is.finite(tpmi_case) & tpmi_case > 0]
    if (length(v)) v[[1L]] else NA_real_
  },
  tpmi_control = {
    v <- tpmi_control[is.finite(tpmi_control) & tpmi_control > 0]
    if (length(v)) v[[1L]] else NA_real_
  },
  tpmi_n = {
    v <- tpmi_n[is.finite(tpmi_n) & tpmi_n > 0]
    if (length(v)) v[[1L]] else NA_real_
  }
), by = output_stem]
assert_unique(n_tab$output_stem, "N table output_stem")
# external_n holds N_eff = 4*case*control/(case+control) for signed-Z weights
n_long <- rbindlist(list(
  n_tab[is.finite(bbj_n) & bbj_n > 0,
        .(cohort = "BBJ", output_stem, external_n = bbj_n)],
  n_tab[is.finite(tpmi_n) & tpmi_n > 0,
        .(cohort = "TPMI", output_stem, external_n = tpmi_n)]
))

lead_stems <- unique(bbj$koges_lead_stem)
direct_map <- meta[!is_empty(KoGES_VAR) & is_direct_match == TRUE]
direct_map[, koges_lead_stem := assign_lead_stem(KoGES_VAR, lead_stems)]
direct_map <- direct_map[!is.na(koges_lead_stem)]
elig <- direct_map[, .(
  bbj_output_stem = collapse_chr(KoGES_VAR[!is_empty(BBJ)]),
  tpmi_output_stem = collapse_chr(KoGES_VAR[!is_empty(TPMI)]),
  n_direct_rows = .N
), by = koges_lead_stem]
elig[, eas_panel_meta_eligible := (
  !is.na(bbj_output_stem) | !is.na(tpmi_output_stem)
)]
elig[, meta_exclusion_reason := fifelse(
  eas_panel_meta_eligible == TRUE,
  NA_character_,
  "Direct match but no BBJ/TPMI path"
)]
unmapped_stems <- setdiff(lead_stems, elig$koges_lead_stem)
if (length(unmapped_stems)) {
  elig <- rbindlist(list(
    elig,
    data.table(
      koges_lead_stem = unmapped_stems,
      bbj_output_stem = NA_character_,
      tpmi_output_stem = NA_character_,
      n_direct_rows = 0L,
      eas_panel_meta_eligible = FALSE,
      meta_exclusion_reason = "No scientifically harmonized binary phenotype in the EAS panel (no BBJ/TPMI direct match)"
    )
  ), use.names = TRUE, fill = TRUE)
}
assert_unique(elig$koges_lead_stem, "eligibility koges_lead_stem")
default_excl_reason <-
  "No scientifically harmonized binary phenotype in the EAS panel (no BBJ/TPMI direct match)"

dir.create(dirname(n_path), recursive = TRUE, showWarnings = FALSE)
fwrite(n_tab, n_path, sep = "\t", quote = FALSE, na = "NA")
fwrite(elig, elig_path, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote effective-N table (case/control -> N_eff):", n_path, "\n")
cat("Wrote eligibility:", elig_path, "\n")
cat(
  "Eligible lead stems:", sum(elig$eas_panel_meta_eligible),
  "/", nrow(elig), "\n"
)

# ------------------------------------------------------------------------------
# 4. Wide table by lead_row_id
# ------------------------------------------------------------------------------
res <- bbj[, .(
  lead_row_id,
  koges_phenotype,
  koges_pheno_code,
  koges_lead_stem,
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
  if (!nrow(dt)) {
    empty_names <- paste0(prefix, "_", c(
      "koges_beta_aligned_a1", "mapping_status", "output_stem",
      "replication_testable", "exact_usable_for_replication",
      "exact_a1_koges", "exact_a2_koges", "exact_beta_aligned_a1",
      "exact_se", "exact_p", "exact_eaf_aligned_a1", "exact_info",
      "external_n", "cohort_usable", "z_from_beta_se", "z_from_exact_p",
      "p_from_z", "p_vs_z_absdiff"
    ))
    empty <- data.table(lead_row_id = res$lead_row_id)
    for (nm in empty_names) empty[, (nm) := NA]
    empty[, (paste0(prefix, "_cohort_usable")) := FALSE]
    empty[, (paste0(prefix, "_mapping_status")) := paste0(
      "no_matching_", prefix, "_phenotype"
    )]
    return(merge(res, empty, by = "lead_row_id", all.x = TRUE))
  }
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
  add[, z_from_exact_p := signed_z_from_p(beta, p)]
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
    z_from_exact_p,
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
res <- attach_cohort(res, bbj_one, "bbj")
res <- attach_cohort(res, tpmi_one, "tpmi")
assert_unique(res$lead_row_id, "EAS panel lead_row_id after cohort joins")
res[is.na(bbj_mapping_status), bbj_mapping_status := "no_matching_bbj_phenotype"]
res[is.na(tpmi_mapping_status), tpmi_mapping_status := "no_matching_tpmi_phenotype"]
res[is.na(bbj_cohort_usable), bbj_cohort_usable := FALSE]
res[is.na(tpmi_cohort_usable), tpmi_cohort_usable := FALSE]

res <- merge(
  res,
  elig[, .(
    koges_lead_stem,
    elig_bbj_output_stem = bbj_output_stem,
    elig_tpmi_output_stem = tpmi_output_stem,
    eas_panel_meta_eligible,
    meta_exclusion_reason
  )],
  by = "koges_lead_stem",
  all.x = TRUE,
  allow.cartesian = FALSE
)
assert_unique(res$lead_row_id, "EAS panel lead_row_id after eligibility join")
res[
  is.na(eas_panel_meta_eligible),
  `:=`(
    eas_panel_meta_eligible = FALSE,
    meta_exclusion_reason = default_excl_reason
  )
]
res[, meta_exclusion_reason := as.character(meta_exclusion_reason)]
res[
  eas_panel_meta_eligible == FALSE,
  `:=`(
    bbj_cohort_usable = FALSE,
    tpmi_cohort_usable = FALSE
  )
]
res[
  !mapply(stem_in_list, bbj_output_stem, elig_bbj_output_stem),
  bbj_cohort_usable := FALSE
]
res[
  !mapply(stem_in_list, tpmi_output_stem, elig_tpmi_output_stem),
  tpmi_cohort_usable := FALSE
]
n_usable_missing_n <- sum(
  (res$bbj_cohort_usable %in% TRUE &
     !(is.finite(res$bbj_external_n) & res$bbj_external_n > 0)) +
    (res$tpmi_cohort_usable %in% TRUE &
       !(is.finite(res$tpmi_external_n) & res$tpmi_external_n > 0)),
  na.rm = TRUE
)
if (n_usable_missing_n > 0L) {
  cat(
    "Note: usable cohort rows with missing/nonpositive external effective N: ",
    n_usable_missing_n,
    " (IVW unaffected; blocks sensitivity signed-Z when k>=2)\n",
    sep = ""
  )
}

# ------------------------------------------------------------------------------
# 5. Panel aggregation
# ------------------------------------------------------------------------------
res[, EAS_n_contributing_cohorts := as.integer(
  (bbj_cohort_usable %in% TRUE) + (tpmi_cohort_usable %in% TRUE)
)]
res[, EAS_contributing_cohorts := {
  labs <- c(
    if (isTRUE(bbj_cohort_usable)) "BBJ",
    if (isTRUE(tpmi_cohort_usable)) "TPMI"
  )
  if (length(labs) == 0L) NA_character_ else paste(labs, collapse = "+")
}, by = lead_row_id]
res[, EAS_allele_anchor_ok := {
  use <- c(isTRUE(bbj_cohort_usable), isTRUE(tpmi_cohort_usable))
  a1 <- c(bbj_exact_a1_koges, tpmi_exact_a1_koges)
  a2 <- c(bbj_exact_a2_koges, tpmi_exact_a2_koges)
  allele_ok_one(use, a1, a2)
}, by = lead_row_id]
res[, EAS_missing_external_n := {
  miss <- FALSE
  if (isTRUE(bbj_cohort_usable)) {
    miss <- miss || !(is.finite(bbj_external_n) && bbj_external_n > 0)
  }
  if (isTRUE(tpmi_cohort_usable)) {
    miss <- miss || !(is.finite(tpmi_external_n) && tpmi_external_n > 0)
  }
  miss
}, by = lead_row_id]
res[, EAS_missing_signedZ_input := {
  miss <- FALSE
  if (isTRUE(bbj_cohort_usable)) {
    miss <- miss || !is.finite(bbj_z_from_exact_p)
  }
  if (isTRUE(tpmi_cohort_usable)) {
    miss <- miss || !is.finite(tpmi_z_from_exact_p)
  }
  miss
}, by = lead_row_id]
res[, `:=`(
  EAS_panel_available = EAS_n_contributing_cohorts >= 1L,
  EAS_panel_a1_koges = NA_character_,
  EAS_panel_a2_koges = NA_character_,
  EAS_evidence_type = NA_character_,
  EAS_panel_beta = NA_real_,
  EAS_panel_SE = NA_real_,
  EAS_panel_Z = NA_real_,
  EAS_panel_Z_rep = NA_real_,
  EAS_panel_P = NA_real_,
  EAS_panel_Q = NA_real_,
  EAS_panel_Q_P = NA_real_,
  EAS_panel_Q_df = NA_integer_,
  EAS_panel_I2 = NA_real_,
  EAS_Q_contrib_BBJ = NA_real_,
  EAS_Q_contrib_TPMI = NA_real_,
  EAS_Q_frac_BBJ = NA_real_,
  EAS_Q_frac_TPMI = NA_real_,
  EAS_Q_top_cohort = NA_character_,
  EAS_LOCO_omit_BBJ_beta = NA_real_,
  EAS_LOCO_omit_BBJ_SE = NA_real_,
  EAS_LOCO_omit_BBJ_Z = NA_real_,
  EAS_LOCO_omit_BBJ_P = NA_real_,
  EAS_LOCO_omit_BBJ_Q = NA_real_,
  EAS_LOCO_omit_BBJ_I2 = NA_real_,
  EAS_LOCO_omit_BBJ_std_resid = NA_real_,
  EAS_LOCO_omit_TPMI_beta = NA_real_,
  EAS_LOCO_omit_TPMI_SE = NA_real_,
  EAS_LOCO_omit_TPMI_Z = NA_real_,
  EAS_LOCO_omit_TPMI_P = NA_real_,
  EAS_LOCO_omit_TPMI_Q = NA_real_,
  EAS_LOCO_omit_TPMI_I2 = NA_real_,
  EAS_LOCO_omit_TPMI_std_resid = NA_real_,
  EAS_LOCO_min_I2 = NA_real_,
  EAS_LOCO_min_I2_cohort = NA_character_,
  EAS_LOCO_I2_change_at_min = NA_real_,
  EAS_panel_method_note = NA_character_,
  EAS_sensitivity_signedZ_Z = NA_real_,
  EAS_sensitivity_signedZ_P = NA_real_,
  EAS_sensitivity_signedZ_note = NA_character_
)]
res[, EAS_panel_analyzed :=
  EAS_panel_available & (EAS_allele_anchor_ok %in% TRUE)
]
res[
  EAS_allele_anchor_ok == TRUE & bbj_cohort_usable == TRUE,
  `:=`(EAS_panel_a1_koges = bbj_exact_a1_koges,
       EAS_panel_a2_koges = bbj_exact_a2_koges)
]
res[
  EAS_allele_anchor_ok == TRUE & is.na(EAS_panel_a1_koges) &
    tpmi_cohort_usable == TRUE,
  `:=`(EAS_panel_a1_koges = tpmi_exact_a1_koges,
       EAS_panel_a2_koges = tpmi_exact_a2_koges)
]
res[
  eas_panel_meta_eligible == FALSE,
  EAS_evidence_type := paste0(
    "Not analyzed (phenotype excluded: ", meta_exclusion_reason, ")"
  )
]
res[
  is.na(EAS_evidence_type) & EAS_n_contributing_cohorts == 0L,
  EAS_evidence_type := "Not available"
]
res[
  is.na(EAS_evidence_type) &
    EAS_panel_available == TRUE &
    EAS_allele_anchor_ok == FALSE,
  EAS_evidence_type := "Not analyzed (allele anchor mismatch)"
]
n_allele_mismatch <- sum(
  res$EAS_panel_available == TRUE & res$EAS_allele_anchor_ok == FALSE,
  na.rm = TRUE
)
if (n_allele_mismatch > 0L) {
  ex <- res[
    EAS_panel_available == TRUE & EAS_allele_anchor_ok == FALSE,
    .(lead_row_id, koges_pheno_code, koges_lead_variant,
      bbj_exact_a1_koges, bbj_exact_a2_koges,
      tpmi_exact_a1_koges, tpmi_exact_a2_koges)
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

idx_ivw <- which(res$EAS_panel_analyzed == TRUE)
if (length(idx_ivw) > 0L) {
  res[
    idx_ivw,
    c(
      "EAS_evidence_type", "EAS_panel_beta", "EAS_panel_SE", "EAS_panel_Z",
      "EAS_panel_P", "EAS_panel_Q", "EAS_panel_Q_P", "EAS_panel_Q_df",
      "EAS_panel_I2", "EAS_panel_method_note"
    ) := {
      betas <- c(
        if (isTRUE(bbj_cohort_usable)) bbj_exact_beta_aligned_a1,
        if (isTRUE(tpmi_cohort_usable)) tpmi_exact_beta_aligned_a1
      )
      ses <- c(
        if (isTRUE(bbj_cohort_usable)) bbj_exact_se,
        if (isTRUE(tpmi_cohort_usable)) tpmi_exact_se
      )
      labs <- c(
        if (isTRUE(bbj_cohort_usable)) "BBJ",
        if (isTRUE(tpmi_cohort_usable)) "TPMI"
      )
      iv <- ivw_fixed(betas, ses)
      etype <- if (iv$k >= 2L) "Meta-analysis (IVW)" else "Single cohort (IVW)"
      note <- paste0(
        "external_IVW_fixed; cohorts=", paste(labs, collapse = "+"),
        "; k=", iv$k, "; match=direct"
      )
      list(
        etype, iv$beta, iv$se, iv$z,
        z_to_p_twosided(iv$z), iv$q, iv$q_p, iv$df, iv$i2, note
      )
    },
    by = lead_row_id
  ]
}

idx_het <- which(
  res$EAS_panel_analyzed == TRUE & res$EAS_n_contributing_cohorts >= 2L
)
if (length(idx_het) > 0L) {
  res[
    idx_het,
    c(
      "EAS_Q_contrib_BBJ", "EAS_Q_contrib_TPMI",
      "EAS_Q_frac_BBJ", "EAS_Q_frac_TPMI",
      "EAS_Q_top_cohort",
      "EAS_LOCO_omit_BBJ_beta", "EAS_LOCO_omit_BBJ_SE", "EAS_LOCO_omit_BBJ_Z",
      "EAS_LOCO_omit_BBJ_P", "EAS_LOCO_omit_BBJ_Q", "EAS_LOCO_omit_BBJ_I2",
      "EAS_LOCO_omit_BBJ_std_resid",
      "EAS_LOCO_omit_TPMI_beta", "EAS_LOCO_omit_TPMI_SE", "EAS_LOCO_omit_TPMI_Z",
      "EAS_LOCO_omit_TPMI_P", "EAS_LOCO_omit_TPMI_Q", "EAS_LOCO_omit_TPMI_I2",
      "EAS_LOCO_omit_TPMI_std_resid",
      "EAS_LOCO_min_I2", "EAS_LOCO_min_I2_cohort",
      "EAS_LOCO_I2_change_at_min"
    ) := {
      betas <- c(
        if (isTRUE(bbj_cohort_usable)) bbj_exact_beta_aligned_a1,
        if (isTRUE(tpmi_cohort_usable)) tpmi_exact_beta_aligned_a1
      )
      ses <- c(
        if (isTRUE(bbj_cohort_usable)) bbj_exact_se,
        if (isTRUE(tpmi_cohort_usable)) tpmi_exact_se
      )
      labs <- c(
        if (isTRUE(bbj_cohort_usable)) "BBJ",
        if (isTRUE(tpmi_cohort_usable)) "TPMI"
      )
      ht <- ivw_q_contrib_loco(betas, ses, labs)
      loco_pack <- function(lab) {
        if (lab %chin% names(ht$loco)) {
          x <- ht$loco[[lab]]
          list(x$beta, x$se, x$z, x$p, x$q, x$i2, x$std_resid)
        } else {
          list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
        }
      }
      qv <- function(lab) {
        if (lab %chin% names(ht$q_contrib)) ht$q_contrib[[lab]] else NA_real_
      }
      qf <- function(lab) {
        if (lab %chin% names(ht$q_frac)) ht$q_frac[[lab]] else NA_real_
      }
      lu <- loco_pack("BBJ")
      lm <- loco_pack("TPMI")
      list(
        qv("BBJ"), qv("TPMI"),
        qf("BBJ"), qf("TPMI"),
        ht$q_top_cohort,
        lu[[1]], lu[[2]], lu[[3]], lu[[4]], lu[[5]], lu[[6]], lu[[7]],
        lm[[1]], lm[[2]], lm[[3]], lm[[4]], lm[[5]], lm[[6]], lm[[7]],
        ht$loco_min_i2, ht$loco_min_i2_cohort, ht$loco_i2_change_at_min
      )
    },
    by = lead_row_id
  ]
}

res[, EAS_panel_koges_beta_aligned_a1 := pick_single(
  bbj_cohort_usable %in% TRUE,
  tpmi_cohort_usable %in% TRUE,
  bbj_koges_beta_aligned_a1,
  tpmi_koges_beta_aligned_a1
)]
res[, EAS_koges_beta_consistent := {
  vals <- c(
    if (isTRUE(bbj_cohort_usable)) bbj_koges_beta_aligned_a1 else NA_real_,
    if (isTRUE(tpmi_cohort_usable)) tpmi_koges_beta_aligned_a1 else NA_real_
  )
  vals <- vals[is.finite(vals)]
  if (length(vals) <= 1L) {
    TRUE
  } else {
    diff(range(vals)) <= 1e-12
  }
}, by = lead_row_id]
n_koges_beta_mismatch <- sum(
  res$EAS_panel_analyzed == TRUE &
    res$EAS_n_contributing_cohorts >= 2L &
    res$EAS_koges_beta_consistent == FALSE,
  na.rm = TRUE
)
if (n_koges_beta_mismatch > 0L) {
  ex <- res[
    EAS_panel_analyzed == TRUE &
      EAS_n_contributing_cohorts >= 2L &
      EAS_koges_beta_consistent == FALSE
  ][1L]
  stop(
    "Usable-cohort KoGES-aligned beta mismatch in ",
    n_koges_beta_mismatch, " analyzed multi-cohort lead(s). ",
    "Example lead_row_id=", ex$lead_row_id
  )
}
n_missing_panel_koges_beta <- sum(
  res$EAS_panel_analyzed == TRUE &
    !is.finite(res$EAS_panel_koges_beta_aligned_a1),
  na.rm = TRUE
)
if (n_missing_panel_koges_beta > 0L) {
  stop(
    "EAS_panel_koges_beta_aligned_a1 missing for ",
    n_missing_panel_koges_beta, " analyzed lead(s)."
  )
}

res[
  EAS_panel_analyzed == TRUE,
  EAS_panel_Z_rep := safe_sign(EAS_panel_koges_beta_aligned_a1) * EAS_panel_Z
]
res[, EAS_panel_direction_concordant := NA]
res[
  EAS_panel_analyzed == TRUE,
  EAS_panel_direction_concordant := {
    ks <- safe_sign(EAS_panel_koges_beta_aligned_a1)
    ps <- safe_sign(EAS_panel_beta)
    is.finite(ks) & is.finite(ps) & (ks == ps)
  }
]
loco_direction_concordant <- function(koges_beta, loco_beta) {
  ks <- safe_sign(koges_beta)
  ls <- safe_sign(loco_beta)
  ifelse(is.finite(ks) & is.finite(ls), ks == ls, NA)
}
res[, c(
  "EAS_LOCO_omit_BBJ_direction_concordant",
  "EAS_LOCO_omit_TPMI_direction_concordant"
) := NA]
res[
  EAS_panel_analyzed == TRUE,
  c(
    "EAS_LOCO_omit_BBJ_direction_concordant",
    "EAS_LOCO_omit_TPMI_direction_concordant"
  ) := list(
    loco_direction_concordant(
      EAS_panel_koges_beta_aligned_a1, EAS_LOCO_omit_BBJ_beta
    ),
    loco_direction_concordant(
      EAS_panel_koges_beta_aligned_a1, EAS_LOCO_omit_TPMI_beta
    )
  )
]
res[, EAS_source_direction_pattern := NA_character_]
res[
  EAS_panel_analyzed == TRUE,
  EAS_source_direction_pattern := {
    ks <- safe_sign(EAS_panel_koges_beta_aligned_a1)
    conc <- c(
      if (isTRUE(bbj_cohort_usable)) {
        zs <- safe_sign(bbj_z_from_beta_se)
        is.finite(ks) & is.finite(zs) & (ks == zs)
      },
      if (isTRUE(tpmi_cohort_usable)) {
        zs <- safe_sign(tpmi_z_from_beta_se)
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

idx_sens <- which(
  res$EAS_panel_analyzed == TRUE &
    !(res$EAS_missing_signedZ_input %in% TRUE) & (
      res$EAS_n_contributing_cohorts == 1L |
        (res$EAS_n_contributing_cohorts >= 2L &
           !(res$EAS_missing_external_n %in% TRUE))
    )
)
if (length(idx_sens) > 0L) {
  res[
    idx_sens,
    c(
      "EAS_sensitivity_signedZ_Z",
      "EAS_sensitivity_signedZ_P",
      "EAS_sensitivity_signedZ_note"
    ) := {
      if (EAS_n_contributing_cohorts == 1L) {
        z <- pick_single(
          bbj_cohort_usable, tpmi_cohort_usable,
          bbj_z_from_exact_p, tpmi_z_from_exact_p
        )
        coh <- pick_single_chr(
          bbj_cohort_usable, tpmi_cohort_usable, "BBJ", "TPMI"
        )
        note <- paste0(
          "sensitivity_signed_Z_single; z_source=exact_P; cohort=", coh
        )
      } else {
        num <- 0
        den <- 0
        if (isTRUE(bbj_cohort_usable)) {
          num <- num + sqrt(bbj_external_n) * bbj_z_from_exact_p
          den <- den + bbj_external_n
        }
        if (isTRUE(tpmi_cohort_usable)) {
          num <- num + sqrt(tpmi_external_n) * tpmi_z_from_exact_p
          den <- den + tpmi_external_n
        }
        z <- num / sqrt(den)
        labs <- c(
          if (isTRUE(bbj_cohort_usable)) "BBJ",
          if (isTRUE(tpmi_cohort_usable)) "TPMI"
        )
        note <- paste0(
          "sensitivity_signed_Z_meta; z_source=exact_P;",
          " weight=effective_N;",
          " N_eff=4*case*control/(case+control);",
          " cohorts=", paste(labs, collapse = "+")
        )
      }
      list(z, z_to_p_twosided(z), note)
    },
    by = lead_row_id
  ]
}
res[, EAS_panel_neglog10P := fifelse(
  EAS_panel_analyzed == TRUE,
  neglog10p_from_z_twosided(EAS_panel_Z),
  NA_real_
)]
res[, EAS_panel_P_display := fifelse(
  EAS_panel_analyzed == TRUE,
  p_display(EAS_panel_P),
  NA_character_
)]

# ------------------------------------------------------------------------------
# 6. Bonferroni
# ------------------------------------------------------------------------------
n_panel_available <- sum(res$EAS_panel_available, na.rm = TRUE)
n_panel_analyzed <- sum(res$EAS_panel_analyzed, na.rm = TRUE)
m_global <- n_panel_analyzed
if (m_global == 0L) {
  stop("No EAS_panel_analyzed leads; cannot define Bonferroni denominator M_global.")
}
alpha_bonf <- alpha_global / m_global
res[, `:=`(
  EAS_panel_N_available = as.integer(n_panel_available),
  EAS_panel_N_analyzed = as.integer(n_panel_analyzed),
  EAS_panel_M_global = as.integer(m_global),
  EAS_panel_alpha_Bonferroni = alpha_bonf,
  EAS_panel_replication_nominal = NA,
  EAS_panel_replication_Bonferroni = NA
)]
res[
  EAS_panel_analyzed == TRUE,
  `:=`(
    EAS_panel_replication_nominal =
      (EAS_panel_direction_concordant %in% TRUE) &
      is.finite(EAS_panel_P) &
      EAS_panel_P <= alpha_global,
    EAS_panel_replication_Bonferroni =
      (EAS_panel_direction_concordant %in% TRUE) &
      is.finite(EAS_panel_P) &
      EAS_panel_P <= alpha_bonf
  )
]
pdiff <- c(res$bbj_p_vs_z_absdiff, res$tpmi_p_vs_z_absdiff)
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
# 7. Write
# ------------------------------------------------------------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
setorder(res, lead_row_id)
fwrite(res, out_path, sep = "\t", quote = FALSE, na = "NA")
phenotype_excluded <- res$eas_panel_meta_eligible %in% FALSE
signedZ_available <- is.finite(res$EAS_sensitivity_signedZ_Z) &
  is.finite(res$EAS_sensitivity_signedZ_P)
signedZ_dir <- {
  ks <- safe_sign(res$EAS_panel_koges_beta_aligned_a1)
  zs <- safe_sign(res$EAS_sensitivity_signedZ_Z)
  is.finite(ks) & is.finite(zs) & (ks == zs)
}
meta_rows <- res$EAS_evidence_type == "Meta-analysis (IVW)"
i2_ge50 <- meta_rows & is.finite(res$EAS_panel_I2) & res$EAS_panel_I2 >= 50
i2_ge75 <- meta_rows & is.finite(res$EAS_panel_I2) & res$EAS_panel_I2 >= 75
summary_dt <- data.table(
  metric = c(
    "n_leads",
    "n_EAS_panel_available",
    "n_EAS_panel_analyzed",
    "n_available_minus_analyzed",
    "M_global",
    "alpha_EAS_Bonferroni",
    "n_k0",
    "n_k1",
    "n_k2",
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
    "n_meta_I2_ge_75"
  ),
  value = c(
    nrow(res),
    n_panel_available,
    n_panel_analyzed,
    n_panel_available - n_panel_analyzed,
    m_global,
    alpha_bonf,
    sum(res$EAS_n_contributing_cohorts == 0L, na.rm = TRUE),
    sum(res$EAS_n_contributing_cohorts == 1L, na.rm = TRUE),
    sum(res$EAS_n_contributing_cohorts == 2L, na.rm = TRUE),
    sum(res$EAS_evidence_type == "Single cohort (IVW)", na.rm = TRUE),
    sum(res$EAS_evidence_type == "Meta-analysis (IVW)", na.rm = TRUE),
    sum(res$EAS_evidence_type == "Not available", na.rm = TRUE),
    sum(res$EAS_evidence_type == "Not analyzed (allele anchor mismatch)", na.rm = TRUE),
    sum(phenotype_excluded, na.rm = TRUE),
    uniqueN(res$koges_lead_stem[phenotype_excluded]),
    sum(res$EAS_panel_direction_concordant %in% TRUE, na.rm = TRUE),
    sum(res$EAS_source_direction_pattern == "All concordant", na.rm = TRUE),
    sum(res$EAS_source_direction_pattern == "Mixed", na.rm = TRUE),
    sum(res$EAS_panel_replication_nominal %in% TRUE, na.rm = TRUE),
    sum(res$EAS_panel_replication_Bonferroni %in% TRUE, na.rm = TRUE),
    sum(res$EAS_panel_analyzed %in% TRUE &
          res$EAS_panel_replication_nominal %in% FALSE, na.rm = TRUE),
    sum(signedZ_available, na.rm = TRUE),
    sum(
      signedZ_available & signedZ_dir &
        res$EAS_sensitivity_signedZ_P <= alpha_global,
      na.rm = TRUE
    ),
    sum(
      signedZ_available & signedZ_dir &
        res$EAS_sensitivity_signedZ_P <= alpha_bonf,
      na.rm = TRUE
    ),
    sum(meta_rows & is.finite(res$EAS_panel_Q_P) &
          res$EAS_panel_Q_P < 0.05, na.rm = TRUE),
    sum(meta_rows & is.finite(res$EAS_panel_Q_P) &
          res$EAS_panel_Q_P < 0.001, na.rm = TRUE),
    sum(i2_ge50, na.rm = TRUE),
    sum(i2_ge75, na.rm = TRUE)
  )
)
fwrite(summary_dt, summary_path, sep = "\t", quote = FALSE, na = "NA")
cat("Wrote:", out_path, "\n")
cat("Wrote:", summary_path, "\n")
cat("EAS_panel_available =", n_panel_available, "\n")
cat("EAS_panel_analyzed  =", n_panel_analyzed, "\n")
cat("M_global            =", m_global, "\n")
cat("alpha_EAS_Bonferroni =", format(alpha_bonf, digits = 6), "\n")
print(summary_dt)
