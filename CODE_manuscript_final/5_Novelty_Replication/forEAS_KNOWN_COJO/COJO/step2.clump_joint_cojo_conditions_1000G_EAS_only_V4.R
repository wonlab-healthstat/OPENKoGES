#!/usr/bin/env Rscript

# ==============================================================================
# EAS-only Step2: define EAS prior signals for targeted KoGES GCTA-COJO
# ------------------------------------------------------------------------------
# Reads EAS Step1 only (BBJ / TWB / TPMI).  EUR sources are not used.
# Default PRELD_DIR is COJO/EAS-only/step1.
#
# Candidate variants are the union of two paths:
#   1. Exact 1000G EAS match (raw identity or unique allele-level alias)
#      -> source-stratified PLINK clump -> clump representatives enter the
#      candidate set (KoGES exact variant or in-clump LD proxy, as in V4).
#   2. No exact 1000G EAS match, but the raw external alleles resolve to
#      exactly one KoGES BIM variant -> that KoGES SNP is added as-is
#      (DIRECT_KOGES_RESCUE_IN_CORE=TRUE).  These hits are not clumped.
#
# EUR_supporting is not produced.  Do not point OUT_DIR at the joint V4 run.
#
# Joint EAS+EUR engine (do not edit while a joint run is live):
#   ../joint_EAS_EUR/step2.clump_joint_cojo_conditions_1000G_FINAL_V4.R
# ------------------------------------------------------------------------------
# This step answers:
#   "Given EAS external GWS signals within +/- 1 Mb of this KoGES lead, does
#    the KoGES lead remain associated after joint conditioning?"
#
# It deliberately keeps two operations separate.
#
#  1. External-prior signal definition:
#     - read the source-level Step1 raw audit file (not the source-collapsed
#       candidate table);
#     - clump each external source separately with its ancestry-matched 1000G
#       reference (EAS for BBJ/TWB/TPMI; EUR for UKB/FinnGen/MVP);
#     - union the selected source representatives, then collapse pairwise
#       KoGES r^2 >= 0.9.  Two 1000G clump representatives from overlapping
#       sources never fold (those clumps already defined within-source
#       independence).  Direct-KoGES rescue SNPs were never 1000G-clumped, so
#       rescue-rescue and rescue-clump pairs may fold even when sources overlap.
#
#  2. KoGES conditional model:
#     - only after selection, check that every selected representative is in
#       the KoGES bfile;
#     - Step3b uses the KoGES bfile for GCTA --cojo-cond.
#
# This is a SOURCE-STRATIFIED, TARGET-LOCAL clumping procedure.  It is an
# appropriate targeted conditional analysis for the pre-specified +/-1 Mb
# prior-signal set.  It is not equivalent to clumping every full external GWAS
# genome-wide before locus filtering; do not describe it as such in the paper.
#
# Coverage vs COJO-readiness:
#   Input errors (invalid P / missing variant ID) remain failures.
#   A raw prior that is valid but absent from the ancestry-matched 1000G panel
#   is recorded as reference-coverage, not as an input error.  Only 1000G-
#   representable hits are clumped.  1000G-clumpability and KoGES direct
#   conditioning are separate.  Same-coordinate SNP-indel or shared-allele
#   encodings are representation_discordant, not validated identity; BIM CHR:BP
#   alone does not prove the same biological allele.  HLA is not a special-case
#   exception.  Discordant and other identity-mismatched 1000G SNPs still enter
#   PLINK clumping so independent external signals are defined.  After clumping,
#   a representative is exact_variant, an LD proxy at the pre-specified r^2 that
#   is itself an exact KoGES variant, or unrepresented.  The proxy search runs
#   in two stages: first the clump's own members (limited to variants that
#   external GWAS actually reported), then, if that fails, every variant in the
#   ancestry panel that is the same variant in KoGES (PANEL_WIDE_PROXY_RESCUE).
#   The r^2 threshold is identical in both stages; only the candidate pool
#   widens.  An unrepresented clump is not a mapping bug when an EUR variant is
#   absent from the KoGES common-imputed panel; it is an ancestry/panel-coverage
#   limitation and is reported as such.
#
# Analysis sets and gates:
#   EAS_primary       -- primary KoGES inference.  It requires EAS external GWS
#                        prior evidence (valid raw hits), not a 1000G clump.
#                        Every EAS 1000G clump that exists must be represented
#                        in KoGES.  1000G-unmatched hits that resolve to one
#                        exact KoGES variant enter the core when
#                        DIRECT_KOGES_RESCUE_IN_CORE=TRUE, so a rescue-only
#                        locus can be EAS_primary_ready.
#   EUR_supporting    -- separate supporting analysis.  It requires complete
#                        representation of EUR-defined clumps, but never alters
#                        the EAS-primary result.  EUR-only loci are explicitly
#                        reported as having no EAS-primary conditional test.
#   all_source_*      -- audit only.  It quantifies complete EAS+EUR coverage
#                        but is deliberately not the primary gate.
#
# Raw-hit 1000G availability remains a QC metric, not the binary primary gate.
# Exact KoGES matches among raw 1000G-missing hits are inserted into the EAS
# primary condition set when DIRECT_KOGES_RESCUE_IN_CORE=TRUE.  KoGES r^2
# grouping then collapses collinear rescue SNPs, including same-source
# rescue-rescue and rescue-clump pairs.  A 1000G-defined clump representative
# outranks a collinear rescue SNP as the retained covariate.  Same-source
# 1000G clump representatives are not folded.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

# ------------------------------------------------------------------------------
# 0. Parameters and paths
# ------------------------------------------------------------------------------

BASE_DIR <- Sys.getenv(
  "BASE_DIR",
  "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
)

PRELD_DIR <- Sys.getenv(
  "PRELD_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "step1")
)

# Step1's row-level source audit is required.  The source-collapsed
# joint_COJO_condition_candidates_preLD.tsv is intentionally not used for
# clumping because it loses source-specific P values and source identity.
RAW_HITS_FILE <- Sys.getenv(
  "RAW_HITS_FILE",
  file.path(PRELD_DIR, "joint_COJO_external_GWS_hits_raw.tsv")
)

PRELD_MANIFEST_FILE <- Sys.getenv(
  "PRELD_MANIFEST_FILE",
  file.path(PRELD_DIR, "joint_COJO_condition_manifest_preLD.tsv")
)

KOGES_BFILE <- Sys.getenv(
  "KOGES_BFILE",
  "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated"
)

REF_1000G_DIR <- Sys.getenv(
  "REF_1000G_DIR",
  file.path(BASE_DIR, "COJO", "1000G_LD_panel_indel")
)
EAS_BFILE <- Sys.getenv("EAS_BFILE", file.path(REF_1000G_DIR, "EAS"))
EUR_BFILE <- Sys.getenv("EUR_BFILE", file.path(REF_1000G_DIR, "EUR"))
PLINK_BIN <- Sys.getenv("PLINK_BIN", "/bin2/plink")

OUT_DIR <- Sys.getenv(
  "OUT_DIR",
  file.path(BASE_DIR, "COJO", "EAS-only", "ancestry_clumped_1000G_V4")
)

NOVELTY_WINDOW_BP <- as.integer(Sys.getenv("NOVELTY_WINDOW_BP", "1000000"))
CLUMP_R2 <- as.numeric(Sys.getenv("CLUMP_R2", "0.1"))
CLUMP_KB <- as.integer(Sys.getenv("CLUMP_KB", "1000"))
CLUMP_P1 <- as.numeric(Sys.getenv("CLUMP_P1", "5e-8"))
# CLUMP_P2 remains the external-GWS eligibility threshold used below.  Use a
# separate permissive PLINK SP2 threshold so the member audit retains every
# variant supplied to PLINK (including a rare P exactly equal to CLUMP_P2).
CLUMP_P2 <- as.numeric(Sys.getenv("CLUMP_P2", "5e-8"))
CLUMP_MEMBER_P2 <- as.numeric(Sys.getenv("CLUMP_MEMBER_P2", "1"))
# A PLINK clump member only passes CLUMP_R2 (default 0.1), which is far too
# weak to call it a proxy.  This separate, stricter threshold is used only
# when the external clump index itself is absent from KoGES.
# A proxy used as a *conditioning* covariate should be tighter than the usual
# r^2=0.8 lookup convention.  Use r^2 >= 0.9 for the primary run and repeat
# the pipeline with PROXY_MIN_R2=0.8 in a distinct OUT_DIR as sensitivity.
PROXY_MIN_R2 <- as.numeric(Sys.getenv("PROXY_MIN_R2", "0.9"))
# An external GWAS reports only its own variants, so an in-clump member proxy is
# limited to SNPs that source happened to submit.  When that fails, search the
# whole ancestry panel for the best 1000G-to-KoGES exact-identity variant in LD
# with the clump index.  The r^2 threshold is unchanged; only the candidate pool
# widens.  Set PANEL_WIDE_PROXY_RESCUE=FALSE to keep the in-clump-only behaviour.
PANEL_WIDE_PROXY_RESCUE <- tolower(
  Sys.getenv("PANEL_WIDE_PROXY_RESCUE", "true")
) %in% c("1", "true", "t", "yes", "y")
# The clump index and a candidate proxy can sit at opposite edges of the target
# window, so the index-to-proxy distance can reach 2 Mb even though both are
# within 1 Mb of the lead.  A 1000 kb LD window would silently drop those pairs.
PANEL_WIDE_PROXY_WINDOW_KB <- as.integer(Sys.getenv(
  "PANEL_WIDE_PROXY_WINDOW_KB",
  as.character(2L * (NOVELTY_WINDOW_BP %/% 1000L) + 1L)
))
# A raw external hit that the ancestry panel does not contain can still be
# resolved to one unique exact variant in the KoGES bfile.  Record this for
# raw-hit coverage QC.  When DIRECT_KOGES_RESCUE_IN_CORE=TRUE those SNPs are
# inserted into the EAS primary condition set (not audit-only).  They are not
# 1000G-clumped, so KoGES r^2 grouping may fold same-source rescue-rescue and
# rescue-clump pairs.  Same-source 1000G clump representatives still never fold.
DIRECT_KOGES_RESCUE <- tolower(
  Sys.getenv("DIRECT_KOGES_RESCUE", "true")
) %in% c("1", "true", "t", "yes", "y")
# EAS-only production: 1000G-unmatched hits that are an exact unique KoGES
# variant are inserted into the core candidate set (not audit-only).
DIRECT_KOGES_RESCUE_IN_CORE <- tolower(
  Sys.getenv("DIRECT_KOGES_RESCUE_IN_CORE", "true")
) %in% c("1", "true", "t", "yes", "y")
if (isTRUE(DIRECT_KOGES_RESCUE_IN_CORE) && !isTRUE(DIRECT_KOGES_RESCUE)) {
  stop("DIRECT_KOGES_RESCUE_IN_CORE=TRUE requires DIRECT_KOGES_RESCUE=TRUE.")
}
# Only exact 1000G EAS identity (raw rsID+alleles or unique CHR:BP alias) is
# clumped.  Representation-discordant 1000G SNPs are not clumped here.
EXACT_1000G_CLUMP_ONLY <- tolower(
  Sys.getenv("EXACT_1000G_CLUMP_ONLY", "true")
) %in% c("1", "true", "t", "yes", "y")
ANALYSIS_LAYER <- "EAS"
# In V4, EUR results are explicitly supporting and are not used to determine
# the EAS-primary conclusion.  Re-estimating EUR proxy LD in 1000G EAS is
# therefore optional provenance QC rather than a gate; it is disabled by
# default because most EUR-specific indices are absent from EAS by design.
VERIFY_EUR_PROXY_IN_EAS <- tolower(
  Sys.getenv("VERIFY_EUR_PROXY_IN_EAS", "false")
) %in% c("1", "true", "t", "yes", "y")
# After source-specific 1000G clumping, use KoGES LD to drop redundant
# covariates.  Two 1000G clump representatives from overlapping sources never
# fold.  Direct-KoGES rescue SNPs have no 1000G-defined independence, so they
# may fold into another rescue or a clump representative at the same r^2 even
# when the source sets overlap.  The cutoff is deliberately much stricter
# than the source-level clumping threshold.
CROSS_SOURCE_REP_R2 <- as.numeric(Sys.getenv("CROSS_SOURCE_REP_R2", "0.9"))
CROSS_SOURCE_LD_WINDOW_KB <- as.integer(Sys.getenv(
  "CROSS_SOURCE_LD_WINDOW_KB",
  as.character(2L * ceiling(NOVELTY_WINDOW_BP / 1000) + 1L)
))
N_PARALLEL <- as.integer(Sys.getenv("N_PARALLEL", Sys.getenv("N_THREADS", "16")))
# CPU-heavy, small working set: ancestry 1000G PLINK clump and snplist writes.
# Representative batches default to the same count; override with N_REP_PARALLEL
# if a high-core fork wave starts copying large membership tables.
N_REP_PARALLEL <- as.integer(Sys.getenv(
  "N_REP_PARALLEL",
  as.character(N_PARALLEL)
))
# Memory-sensitive forks: cross-source PLINK --r2 on the KoGES 160K bfile.
# Keep this well below N_PARALLEL when using 80 clump workers.
N_MEM_PARALLEL <- as.integer(Sys.getenv(
  "N_MEM_PARALLEL",
  as.character(max(1L, min(N_PARALLEL, 16L)))
))
# data.table OpenMP threads for the long single-process prep (fread/joins/
# aggregations). Forked PLINK/LD workers force 1 thread to avoid oversubscribe.
N_DT_THREADS <- as.integer(Sys.getenv(
  "N_DT_THREADS",
  as.character(max(1L, min(N_PARALLEL, 16L)))
))
OVERWRITE <- tolower(Sys.getenv("OVERWRITE", "false")) %in%
  c("1", "true", "t", "yes", "y")
# Skip prep + ancestry PLINK clumping; continue from written clump membership.
# If representative/cross-source tables are already in OUT_DIR, those stages
# are skipped too so a crash after the last representative batch can resume.
RESUME_AFTER_CLUMP <- tolower(Sys.getenv("RESUME_AFTER_CLUMP", "false")) %in%
  c("1", "true", "t", "yes", "y")
# Representative selection batches clumps so mclapply does not return ~1e5-1e6
# full proxy-audit tables in one fork wave (that OOM-kills every worker).
REP_BATCH_SIZE <- as.integer(Sys.getenv("REP_BATCH_SIZE", "1000"))

# Optional target subset for a smoke / pilot run.  Every downstream table then
# describes only these targets, so a pilot OUT_DIR must never be merged with a
# full run.  Same variable names as Step3b.
TARGET_IDS_ENV <- Sys.getenv("TARGET_IDS", "")
TARGET_ID_FILE <- Sys.getenv("TARGET_ID_FILE", "")
# Panel-wide stage-2 LD is normally computed once per (target, ancestry panel).
# Setting this to FALSE forces the original one-PLINK-job-per-clump path; it
# exists so the two paths can be compared on a subset and is never faster.
PANEL_WIDE_PROXY_GROUP_LD <- tolower(
  Sys.getenv("PANEL_WIDE_PROXY_GROUP_LD", "true")
) %in% c("1", "true", "t", "yes", "y")


CLUMP_SCOPE <- Sys.getenv("CLUMP_SCOPE", "target_local_source_stratified")
if (!identical(CLUMP_SCOPE, "target_local_source_stratified")) {
  stop(
    "Only CLUMP_SCOPE=target_local_source_stratified is implemented in this ",
    "script.  A genome-wide external-GWAS index analysis requires a separately ",
    "generated source-wide clump-membership map."
  )
}

if (!is.finite(NOVELTY_WINDOW_BP) || NOVELTY_WINDOW_BP <= 0L) {
  stop("NOVELTY_WINDOW_BP must be a positive integer.")
}
if (!is.finite(CLUMP_R2) || CLUMP_R2 <= 0 || CLUMP_R2 > 1) {
  stop("CLUMP_R2 must be in (0, 1].")
}
if (!is.finite(CLUMP_KB) || CLUMP_KB < 1L) {
  stop("CLUMP_KB must be a positive integer.")
}
if (!is.finite(CLUMP_P1) || CLUMP_P1 <= 0 || CLUMP_P1 > 1 ||
    !is.finite(CLUMP_P2) || CLUMP_P2 <= 0 || CLUMP_P2 > 1 ||
    !is.finite(CLUMP_MEMBER_P2) || CLUMP_MEMBER_P2 <= 0 || CLUMP_MEMBER_P2 > 1) {
  stop("CLUMP_P1, CLUMP_P2, and CLUMP_MEMBER_P2 must be in (0, 1].")
}
if (!is.finite(PROXY_MIN_R2) || PROXY_MIN_R2 < 0.8 || PROXY_MIN_R2 > 1) {
  stop("PROXY_MIN_R2 must be in [0.8, 1].")
}
minimum_panel_wide_proxy_window_kb <-
  2L * ceiling(NOVELTY_WINDOW_BP / 1000) + 1L
if (!is.finite(PANEL_WIDE_PROXY_WINDOW_KB) ||
    PANEL_WIDE_PROXY_WINDOW_KB < minimum_panel_wide_proxy_window_kb) {
  stop(
    "PANEL_WIDE_PROXY_WINDOW_KB must be at least ",
    minimum_panel_wide_proxy_window_kb,
    " so a clump index and a proxy at opposite target-window edges are evaluated."
  )
}
if (!is.finite(CROSS_SOURCE_REP_R2) || CROSS_SOURCE_REP_R2 < 0.9 ||
    CROSS_SOURCE_REP_R2 > 1) {
  stop("CROSS_SOURCE_REP_R2 must be in [0.9, 1].")
}
minimum_cross_source_window_kb <- 2L * ceiling(NOVELTY_WINDOW_BP / 1000) + 1L
if (!is.finite(CROSS_SOURCE_LD_WINDOW_KB) ||
    CROSS_SOURCE_LD_WINDOW_KB < minimum_cross_source_window_kb) {
  stop(
    "CROSS_SOURCE_LD_WINDOW_KB must be at least ",
    minimum_cross_source_window_kb,
    " so every pair within the target +/- window is evaluated."
  )
}
if (!is.finite(N_PARALLEL) || N_PARALLEL < 1L) {
  stop("N_PARALLEL must be at least 1.")
}
if (!is.finite(N_REP_PARALLEL) || N_REP_PARALLEL < 1L) {
  stop("N_REP_PARALLEL must be at least 1.")
}
if (!is.finite(N_MEM_PARALLEL) || N_MEM_PARALLEL < 1L) {
  stop("N_MEM_PARALLEL must be at least 1.")
}
if (!is.finite(N_DT_THREADS) || N_DT_THREADS < 1L) {
  stop("N_DT_THREADS must be at least 1.")
}
if (!is.finite(REP_BATCH_SIZE) || REP_BATCH_SIZE < 1L) {
  stop("REP_BATCH_SIZE must be at least 1.")
}

# Use DT threads during prep. parallel_lapply() drops to 1 thread only while
# forking PLINK/LD workers (OpenMP + fork oversubscribe is unsafe).
try(data.table::setDTthreads(N_DT_THREADS), silent = TRUE)
message(
  "Parallel settings: N_PARALLEL=", N_PARALLEL,
  " (1000G clump + snplist), N_REP_PARALLEL=", N_REP_PARALLEL,
  " (representative batches), N_MEM_PARALLEL=", N_MEM_PARALLEL,
  " (KoGES cross-source LD), N_DT_THREADS=", N_DT_THREADS,
  " (data.table OpenMP during prep), REP_BATCH_SIZE=", REP_BATCH_SIZE,
  ", RESUME_AFTER_CLUMP=", RESUME_AFTER_CLUMP,
  ", ANALYSIS_LAYER=", ANALYSIS_LAYER,
  ", DIRECT_KOGES_RESCUE_IN_CORE=", DIRECT_KOGES_RESCUE_IN_CORE,
  ", EXACT_1000G_CLUMP_ONLY=", EXACT_1000G_CLUMP_ONLY
)

check_bfile <- function(prefix, label) {
  if (is.na(prefix) || !nzchar(prefix)) {
    stop(label, " is empty. Set it to a PLINK bfile prefix.")
  }
  parts <- paste0(prefix, c(".bed", ".bim", ".fam"))
  missing <- parts[!file.exists(parts)]
  if (length(missing) > 0L) {
    stop("Missing ", label, " component(s):\n", paste(missing, collapse = "\n"))
  }
}

for (path in c(RAW_HITS_FILE, PRELD_MANIFEST_FILE)) {
  if (!file.exists(path)) stop("Missing Step1 input file: ", path)
}
check_bfile(KOGES_BFILE, "KOGES_BFILE")
check_bfile(EAS_BFILE, "EAS_BFILE (1000G EAS)")
# EUR panel is not used in the EAS-only script.

plink_found <- if (grepl("/", PLINK_BIN, fixed = TRUE)) {
  file.exists(PLINK_BIN)
} else {
  nzchar(Sys.which(PLINK_BIN))
}
if (!plink_found) stop("PLINK was not found. Set PLINK_BIN.")

if (dir.exists(OUT_DIR)) {
  old <- list.files(OUT_DIR, all.files = FALSE, no.. = TRUE)
  if (length(old) > 0L && !OVERWRITE) {
    stop(
      "OUT_DIR already contains files: ", OUT_DIR,
      "\nUse a fresh OUT_DIR or set OVERWRITE=TRUE after checking the old run."
    )
  }
} else {
  dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
}
CLUMP_WORK_DIR <- file.path(OUT_DIR, "ancestry_clump_work")
dir.create(CLUMP_WORK_DIR, recursive = TRUE, showWarnings = FALSE)
CROSS_SOURCE_LD_WORK_DIR <- file.path(OUT_DIR, "cross_source_KoGES_LD_work")
dir.create(CROSS_SOURCE_LD_WORK_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

target_columns <- c(
  "target_id", "lead_row_id", "koges_phenotype", "koges_pheno_code",
  "koges_lead_variant", "koges_chr", "koges_bp"
)

collapse_unique <- function(x, separator = ";") {
  x <- sort(unique(as.character(x[!is.na(x) & nzchar(as.character(x))])))
  if (length(x) == 0L) NA_character_ else paste(x, collapse = separator)
}

# Source-level collapse stores unique labels in clump_snp_resolution_statuses.
# EXACT_1000G_CLUMP_ONLY requires every remaining token to be an exact identity
# (raw rsID match, unique allele-level alias, or verified KoGES condition SNP).
all_exact_1000g_clump_statuses <- function(x) {
  pat <- paste0(
    "^(external_snp_raw_identity_verified_in_1000G_",
    "|exact_allele_level_alias_in_1000G_",
    "|verified_KoGES_condition_snp_alias_in_1000G_)"
  )
  x <- as.character(x)
  out <- rep(FALSE, length(x))
  ok <- !is.na(x) & nzchar(x)
  if (!any(ok)) return(out)
  toks <- strsplit(x[ok], ";", fixed = TRUE)
  out[ok] <- vapply(
    toks,
    function(v) {
      v <- v[nzchar(v)]
      length(v) > 0L && all(grepl(pat, v))
    },
    logical(1L)
  )
  out
}

collapse_source_sets <- function(x) {
  x <- x[!is.na(x) & nzchar(as.character(x))]
  if (length(x) == 0L) return(NA_character_)
  collapse_unique(unlist(strsplit(as.character(x), ";", fixed = TRUE)))
}

# Separate true input errors from ancestry-1000G representability.  Older QC
# files only have prior_input_status; newer files have both raw_input_status
# and reference_coverage_status.  Populate whichever side is missing.
harmonize_raw_prior_status_columns <- function(dt) {
  has_new <- all(c("raw_input_status", "reference_coverage_status") %in% names(dt))
  has_old <- "prior_input_status" %in% names(dt)
  if (!has_new && !has_old) {
    stop(
      "Raw prior QC lacks raw_input_status/reference_coverage_status ",
      "and prior_input_status."
    )
  }
  if (!has_new) {
    dt[, raw_input_status := fcase(
      prior_input_status %chin% c(
        "invalid_or_non_GWS_external_P",
        "no_valid_variant_identifier"
      ),
      prior_input_status,
      default = "valid"
    )]
    dt[, reference_coverage_status := fcase(
      raw_input_status != "valid",
      "not_evaluated_invalid_input",
      prior_input_status %chin% c(
        "eligible", "represented_in_ancestry_1000G"
      ),
      "represented_in_ancestry_1000G",
      default = "not_represented_in_ancestry_1000G"
    )]
  }
  if (!has_old) {
    dt[, prior_input_status := fcase(
      raw_input_status != "valid", as.character(raw_input_status),
      reference_coverage_status == "represented_in_ancestry_1000G",
      "eligible",
      default = "external_prior_not_represented_in_ancestry_1000G"
    )]
  }
  dt
}

summarize_raw_prior_input <- function(dt) {
  if ("direct_KoGES_rescue_eligible" %in% names(dt)) {
    dt[, direct_KoGES_rescue_eligible := as_bool(direct_KoGES_rescue_eligible)]
  } else {
    dt[, direct_KoGES_rescue_eligible := FALSE]
  }
  dt[, .(
    n_raw_prior_hits_input = .N,
    n_raw_prior_hits_invalid = sum(raw_input_status != "valid"),
    n_raw_prior_hits_valid = sum(raw_input_status == "valid"),
    n_raw_prior_hits_represented_1000G = sum(
      reference_coverage_status == "represented_in_ancestry_1000G"
    ),
    n_raw_prior_hits_unclumpable_1000G = sum(
      reference_coverage_status == "not_represented_in_ancestry_1000G"
    ),
    n_raw_prior_hits_direct_KoGES_rescued = sum(direct_KoGES_rescue_eligible),
    n_raw_prior_hits_unresolved_after_direct_KoGES_rescue = sum(
      raw_input_status == "valid" &
        reference_coverage_status == "not_represented_in_ancestry_1000G" &
        !direct_KoGES_rescue_eligible
    ),
    prior_input_failure_reasons = collapse_unique(
      raw_input_status[raw_input_status != "valid"]
    ),
    reference_coverage_reasons = collapse_unique(
      reference_coverage_status[
        reference_coverage_status == "not_represented_in_ancestry_1000G"
      ]
    ),
    direct_KoGES_rescue_statuses = if ("direct_KoGES_exact_status" %in% names(dt)) {
      collapse_unique(direct_KoGES_exact_status[direct_KoGES_rescue_eligible])
    } else {
      NA_character_
    }
  ), by = target_columns]
}

# A direct KoGES rescue is not a replacement 1000G clump index.  It is the
# narrow case in which a GWS external hit is absent from its matched 1000G
# panel but its *raw* CHR:BP+alleles identify one exact KoGES variant.  Keep it
# as an additional condition SNP and retain its source provenance so the final
# manifest can distinguish complete 1000G coverage from complete conditioning
# that includes this explicitly labelled rescue.
build_direct_koges_rescue_tables <- function(dt) {
  required <- c(
    target_columns, "source", "source_output_stem", "source_key",
    "evidence_layer", "external_p", "distance_from_koges_lead_bp",
    "external_snp", "direct_KoGES_exact_snp",
    "direct_KoGES_exact_status", "direct_KoGES_rescue_eligible"
  )
  if (!all(required %in% names(dt))) return(list(union = NULL, selection = NULL))
  rescue <- dt[as_bool(direct_KoGES_rescue_eligible) &
    is_valid_snp_id(direct_KoGES_exact_snp)]
  if (nrow(rescue) == 0L) return(list(union = NULL, selection = NULL))

  union <- rescue[, .(
    reported_in_sources = collapse_unique(source),
    reported_in_layers = collapse_unique(evidence_layer),
    source_output_stems = collapse_unique(source_output_stem),
    source_p_min = min_finite(external_p),
    source_p_values = format_p_values(external_p),
    source_abs_distance_min_bp = min_finite(abs(distance_from_koges_lead_bp)),
    source_abs_distance_max_bp = max_finite(abs(distance_from_koges_lead_bp)),
    source_variant_ids = collapse_unique(external_snp),
    n_raw_source_hits = .N,
    direct_KoGES_rescue = TRUE,
    n_direct_KoGES_rescued_hits = .N,
    direct_KoGES_rescue_statuses = collapse_unique(direct_KoGES_exact_status)
  ), by = c(target_columns, "direct_KoGES_exact_snp")]
  setnames(union, "direct_KoGES_exact_snp", "condition_snp")

  selection <- rescue[, .(
    selected_in_EAS = any(evidence_layer == "EAS"),
    selected_in_EUR = any(evidence_layer == "EUR"),
    selected_panels = collapse_unique(evidence_layer),
    selected_sources = collapse_unique(source),
    selected_source_output_stems = collapse_unique(source_output_stem),
    selected_source_keys = collapse_unique(source_key),
    n_distinct_selected_source_keys = as.integer(uniqueN(source_key)),
    panel_clump_p_min = NA_real_,
    external_clump_index_snps = NA_character_,
    selected_clump_member_snps = collapse_unique(direct_KoGES_exact_snp),
    representative_types = "direct_KoGES_exact_external_hit_no_1000G",
    representative_statuses = "direct_KoGES_exact_rescue_no_1000G_clump",
    proxy_r2_1000G_min = NA_real_,
    proxy_r2_1000G_values = NA_character_,
    n_external_clumps_represented = 0L,
    direct_KoGES_rescue = TRUE,
    n_direct_KoGES_rescued_hits = .N,
    direct_KoGES_rescue_statuses = collapse_unique(direct_KoGES_exact_status)
  ), by = c("target_id", "direct_KoGES_exact_snp")]
  setnames(selection, "direct_KoGES_exact_snp", "condition_snp")
  list(union = union, selection = selection)
}

format_p_values <- function(x) {
  x <- sort(unique(as.numeric(x[is.finite(x)])))
  if (length(x) == 0L) NA_character_ else {
    paste(formatC(x, format = "e", digits = 3), collapse = ";")
  }
}

format_numeric_values <- function(x) {
  x <- sort(unique(as.numeric(x[is.finite(x)])))
  if (length(x) == 0L) NA_character_ else {
    paste(formatC(x, format = "fg", digits = 6), collapse = ";")
  }
}

min_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else min(x)
}

max_finite <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0L) NA_real_ else max(x)
}

is_valid_snp_id <- function(x) {
  x <- trimws(as.character(x))
  !is.na(x) & nzchar(x) & !(x %in% c("NA", ".", "-", "0"))
}

as_bool <- function(x) {
  if (is.logical(x)) {
    y <- x
    y[is.na(y)] <- FALSE
    return(y)
  }
  tolower(trimws(as.character(x))) %in% c("1", "true", "t", "yes", "y")
}

normalize_chr <- function(x) {
  y <- gsub("^CHR", "", toupper(trimws(as.character(x))))
  y[y == "23"] <- "X"
  y[y == "24"] <- "Y"
  y[y == "25"] <- "XY"
  y[y %in% c("26", "M")] <- "MT"
  y
}

normalize_allele <- function(x) {
  y <- toupper(trimws(as.character(x)))
  y[y %chin% c("", "NA", ".", "-")] <- NA_character_
  y
}

complement_allele <- function(x) {
  x <- normalize_allele(x)
  is_snv <- !is.na(x) & grepl("^[ACGT]+$", x)
  out <- rep(NA_character_, length(x))
  # Single-base alleles dominate a whole-panel BIM comparison; keep them on a
  # vectorised path so genome-wide identity maps stay tractable.
  single <- is_snv & nchar(x) == 1L
  if (any(single)) out[single] <- chartr("ACGT", "TGCA", x[single])
  multi <- is_snv & !single
  if (any(multi)) {
    split_bases <- strsplit(x[multi], "", fixed = TRUE)
    out[multi] <- vapply(
      split_bases,
      function(bases) paste(rev(chartr("ACGT", "TGCA", bases)), collapse = ""),
      character(1L)
    )
  }
  out
}

# Exact biallelic identity: allele-order reversal and strand reversal only.
# Same CHR:BP SNP-indel or shared-allele encodings are NOT exact identity.
allele_pair_matches <- function(a1, a2, b1, b2) {
  a1 <- normalize_allele(a1)
  a2 <- normalize_allele(a2)
  b1 <- normalize_allele(b1)
  b2 <- normalize_allele(b2)
  complete <- !is.na(a1) & !is.na(a2) & !is.na(b1) & !is.na(b2)
  direct <- (a1 == b1 & a2 == b2) | (a1 == b2 & a2 == b1)
  ca1 <- complement_allele(a1)
  ca2 <- complement_allele(a2)
  strand <- (ca1 == b1 & ca2 == b2) | (ca1 == b2 & ca2 == b1)
  complete & (direct | strand)
}

allele_nchar <- function(x) {
  x <- normalize_allele(x)
  out <- rep(NA_integer_, length(x))
  ok <- !is.na(x)
  out[ok] <- nchar(x[ok], type = "chars")
  out
}

is_snp_allele_pair <- function(a1, a2) {
  n1 <- allele_nchar(a1)
  n2 <- allele_nchar(a2)
  !is.na(n1) & !is.na(n2) & n1 == 1L & n2 == 1L
}

is_extended_allele_pair <- function(a1, a2) {
  n1 <- allele_nchar(a1)
  n2 <- allele_nchar(a2)
  !is.na(n1) & !is.na(n2) & (n1 != 1L | n2 != 1L)
}

# Audit-only: 1000G SNP vs KoGES indel/MNP at the same CHR:BP.  This is
# representation discordance, not proof of the same biological allele.
snp_indel_chrbp_alias <- function(a1, a2, b1, b2) {
  (is_snp_allele_pair(a1, a2) & is_extended_allele_pair(b1, b2)) |
    (is_extended_allele_pair(a1, a2) & is_snp_allele_pair(b1, b2))
}

# Audit-only: same-position triallelic / ALT-choice (example: C/G vs T/G).
# Sharing one allele is not exact identity and is not used as a KoGES condition.
allele_pair_shares_one <- function(a1, a2, b1, b2) {
  a1 <- normalize_allele(a1)
  a2 <- normalize_allele(a2)
  b1 <- normalize_allele(b1)
  b2 <- normalize_allele(b2)
  complete <- !is.na(a1) & !is.na(a2) & !is.na(b1) & !is.na(b2)
  ca1 <- complement_allele(a1)
  ca2 <- complement_allele(a2)
  share_direct <- (a1 == b1 | a1 == b2 | a2 == b1 | a2 == b2)
  share_strand <- (!is.na(ca1) & (ca1 == b1 | ca1 == b2)) |
    (!is.na(ca2) & (ca2 == b1 | ca2 == b2))
  complete & (share_direct | share_strand)
}

EXACT_IDENTITY_STATUS <- "validated_1000G_to_KoGES_identity"
DISCORDANT_SNP_INDEL_STATUS <-
  "representation_discordant_same_coordinate_snp_indel"
DISCORDANT_SHARED_ALLELE_STATUS <-
  "representation_discordant_same_coordinate_shared_allele_multiallelic"
HARD_IDENTITY_MISMATCH_STATUSES <- c(
  "1000G_to_KoGES_chr_mismatch",
  "1000G_to_KoGES_bp_mismatch",
  "1000G_to_KoGES_allele_mismatch"
)
is_exact_identity_status <- function(x) {
  as.character(x) == EXACT_IDENTITY_STATUS
}

safe_target_label <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
}

# Comma / whitespace separated IDs, or one ID per line in a file.
resolve_target_subset <- function(ids_env, id_file) {
  ids <- character()
  if (nzchar(ids_env)) {
    ids <- trimws(unlist(strsplit(ids_env, "[,[:space:]]+"), use.names = FALSE))
  }
  if (nzchar(id_file)) {
    if (!file.exists(id_file)) {
      stop("TARGET_ID_FILE does not exist: ", id_file)
    }
    ids <- c(ids, trimws(readLines(id_file, warn = FALSE)))
  }
  unique(ids[nzchar(ids)])
}

# Keep source_key safe to write to TSV audit files.  It is length-prefixed so
# source/source_output_stem combinations remain unambiguous without embedding
# a literal tab in a data field.
make_source_key <- function(source, source_output_stem) {
  source <- as.character(source)
  source_output_stem <- as.character(source_output_stem)
  paste0(
    nchar(source), ":", source, "|",
    nchar(source_output_stem), ":", source_output_stem
  )
}

parallel_lapply <- function(X, FUN, ..., n_workers = N_PARALLEL) {
  n_workers <- as.integer(n_workers)
  if (!is.finite(n_workers) || n_workers < 1L) n_workers <- 1L
  n_workers <- min(n_workers, max(1L, length(X)))
  if (n_workers == 1L || length(X) <= 1L) return(lapply(X, FUN, ...))
  if (!requireNamespace("parallel", quietly = TRUE)) {
    warning("parallel package unavailable; running sequentially.")
    return(lapply(X, FUN, ...))
  }
  # Avoid OpenMP+fork oversubscribe in children; restore prep threads after.
  prev_dt_threads <- tryCatch(data.table::getDTthreads(), error = function(e) 1L)
  try(data.table::setDTthreads(1L), silent = TRUE)
  on.exit(try(data.table::setDTthreads(prev_dt_threads), silent = TRUE), add = TRUE)
  res <- parallel::mclapply(
    X, FUN, ...,
    mc.cores = n_workers,
    mc.preschedule = TRUE
  )
  bad <- vapply(res, inherits, logical(1L), what = "try-error")
  if (any(bad)) {
    first <- res[[which(bad)[[1L]]]]
    stop(
      "parallel::mclapply failed for ", sum(bad), "/", length(res),
      " tasks (mc.cores=", n_workers, "). First error: ",
      paste(as.character(first), collapse = " ")
    )
  }
  res
}

# Select one KoGES-usable representative per PLINK clump, in memory-safe batches.
# Returning ~3.7e5 full proxy-audit tables from one mclapply wave OOMs workers.
select_clump_representatives_batched <- function(clump_membership,
                                                batch_size = REP_BATCH_SIZE,
                                                n_workers = N_REP_PARALLEL) {
  if (nrow(clump_membership) == 0L) {
    return(list(
      representatives = empty_clump_representative_table(),
      proxy_audit = empty_member_proxy_audit(),
      panel_wide_audit = empty_panel_wide_proxy_audit()
    ))
  }
  clump_ids <- unique(as.character(clump_membership$clump_id))
  n_clumps <- length(clump_ids)
  message(
    "Selecting representatives for ", n_clumps, " clumps in batches of ",
    batch_size, " (N_REP_PARALLEL=", n_workers, ")"
  )
  rep_chunks <- vector("list", length = 0L)
  audit_chunks <- vector("list", length = 0L)
  panel_wide_chunks <- vector("list", length = 0L)
  batch_number <- 0L
  for (start in seq(1L, n_clumps, by = batch_size)) {
    end <- min(start + batch_size - 1L, n_clumps)
    ids <- clump_ids[start:end]
    sub <- clump_membership[clump_id %chin% ids]
    groups <- split(seq_len(nrow(sub)), sub$clump_id)
    batch_res <- parallel_lapply(
      groups,
      function(i) select_one_clump_representative(sub[i]),
      n_workers = n_workers
    )
    batch_number <- batch_number + 1L
    rep_chunks[[batch_number]] <- rbindlist(
      lapply(batch_res, function(x) x$representative),
      use.names = TRUE, fill = TRUE
    )
    audit_chunks[[batch_number]] <- rbindlist(
      lapply(batch_res, function(x) x$proxy_audit),
      use.names = TRUE, fill = TRUE
    )
    panel_wide_chunks[[batch_number]] <- rbindlist(
      lapply(batch_res, function(x) x$panel_wide_audit),
      use.names = TRUE, fill = TRUE
    )
    rm(sub, groups, batch_res)
    if (batch_number %% 10L == 0L) invisible(gc())
    message(
      "  representative batch ", batch_number, ": clumps ",
      start, "-", end, " / ", n_clumps
    )
  }
  list(
    representatives = rbindlist(rep_chunks, use.names = TRUE, fill = TRUE),
    proxy_audit = rbindlist(audit_chunks, use.names = TRUE, fill = TRUE),
    panel_wide_audit = rbindlist(
      panel_wide_chunks, use.names = TRUE, fill = TRUE
    )
  )
}

read_bim_index <- function(prefix) {
  x <- fread(
    paste0(prefix, ".bim"),
    header = FALSE,
    select = c(1, 2, 4, 5, 6),
    showProgress = FALSE
  )
  setnames(x, c("ref_chr", "ref_snp", "ref_bp", "ref_a1", "ref_a2"))
  x[, ref_chr := as.character(ref_chr)]
  x[, ref_snp := as.character(ref_snp)]
  x[, ref_bp := as.integer(ref_bp)]
  x[, ref_a1 := normalize_allele(ref_a1)]
  x[, ref_a2 := normalize_allele(ref_a2)]
  if (anyDuplicated(x$ref_snp)) {
    duplicate_ids <- unique(x$ref_snp[duplicated(x$ref_snp)])
    stop(
      "Duplicated SNP IDs in bfile ", prefix, ". Example: ",
      paste(head(duplicate_ids, 10L), collapse = ", ")
    )
  }
  x
}

write_snplist <- function(snps, path) {
  snps <- unique(as.character(snps[!is.na(snps) & nzchar(snps)]))
  fwrite(data.table(SNP = snps), path, sep = "\t", col.names = FALSE, quote = FALSE)
}

# PLINK's SP2 field is a comma-separated list of clump members.  With a
# single --clump file it may append a numeric file marker such as "(1)".
# Keep the raw field for audit and remove only that trailing marker for ID
# matching.  The index SNP is added separately because SP2 need not include it.
parse_clump_sp2_members <- function(x) {
  x <- trimws(as.character(x))
  if (length(x) == 0L || is.na(x) || !nzchar(x) || x %chin% c("NONE", "NA", ".")) {
    return(character())
  }
  members <- trimws(unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE))
  members <- sub("\\([0-9]+\\)$", "", members, perl = TRUE)
  members <- trimws(members)
  members <- members[is_valid_snp_id(members)]
  unique(members)
}

empty_clump_membership <- function() {
  data.table(
    target_id = character(), source_key = character(), source = character(),
    source_output_stem = character(), clump_panel = character(),
    clump_id = character(), clump_index_snp = character(),
    clump_member_snp = character(), member_is_index = logical(),
    member_reported_in_SP2 = logical(), clump_SP2_raw = character(),
    member_present_in_source_input = logical(), member_clump_p = numeric(),
    member_condition_snp = character(), member_KoGES_status = character(),
    member_identity_status = character(),
    clump_out_prefix = character()
  )
}

make_clump_membership <- function(common, clumped, sub, out_prefix = NA_character_) {
  if (is.null(clumped) || nrow(clumped) == 0L || !"SNP" %in% names(clumped)) {
    return(empty_clump_membership())
  }

  source_lookup <- unique(sub[, .(
    clump_snp = as.character(clump_snp),
    member_clump_p = as.numeric(clump_p),
    member_condition_snp = as.character(condition_snp),
    member_KoGES_status = as.character(condition_KoGES_status),
    member_identity_status = as.character(panel_to_KoGES_identity_status)
  )])

  rows <- lapply(seq_len(nrow(clumped)), function(i) {
    index_snp <- trimws(as.character(clumped$SNP[[i]]))
    if (!is_valid_snp_id(index_snp)) return(NULL)
    sp2_raw <- if ("SP2" %in% names(clumped)) as.character(clumped$SP2[[i]]) else NA_character_
    sp2_members <- parse_clump_sp2_members(sp2_raw)
    members <- unique(c(index_snp, sp2_members))
    member_match <- match(members, source_lookup$clump_snp)
    data.table(
      target_id = common$target_id,
      source_key = common$source_key,
      source = common$source,
      source_output_stem = common$source_output_stem,
      clump_panel = common$clump_panel,
      clump_id = paste(
        common$target_id, common$source, common$source_output_stem,
        common$clump_panel, index_snp, sep = "|"
      ),
      clump_index_snp = index_snp,
      clump_member_snp = members,
      member_is_index = members == index_snp,
      member_reported_in_SP2 = members %chin% sp2_members,
      clump_SP2_raw = sp2_raw,
      member_present_in_source_input = !is.na(member_match),
      member_clump_p = source_lookup$member_clump_p[member_match],
      member_condition_snp = source_lookup$member_condition_snp[member_match],
      member_KoGES_status = source_lookup$member_KoGES_status[member_match],
      member_identity_status = source_lookup$member_identity_status[member_match],
      clump_out_prefix = out_prefix
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (length(rows) == 0L) return(empty_clump_membership())
  rbindlist(rows, use.names = TRUE, fill = TRUE)
}

empty_proxy_ld_table <- function() {
  data.table(proxy_candidate_snp = character(), proxy_r2_1000G = numeric())
}

# Calculate explicit pairwise genotype r^2 in the same ancestry-matched 1000G
# panel that created the clump.  SP2 membership alone is not accepted as proxy
# evidence, since the clump threshold (CLUMP_R2) can be much lower than 0.8.
calculate_member_proxy_ld <- function(index_snp, proxy_candidates, ref_bfile,
                                      out_prefix, region_chr = NA_character_,
                                      region_from_bp = NA_real_,
                                      region_to_bp = NA_real_) {
  proxy_candidates <- unique(as.character(proxy_candidates))
  proxy_candidates <- proxy_candidates[is_valid_snp_id(proxy_candidates)]
  if (length(proxy_candidates) == 0L) {
    return(list(
      status = "no_eligible_member_candidates",
      table = empty_proxy_ld_table(),
      candidate_file = NA_character_, ld_file = NA_character_
    ))
  }

  label <- paste0(out_prefix, "__", safe_target_label(index_snp), "__member_proxy_ld")
  candidate_file <- paste0(label, ".snplist")
  ld_file <- paste0(label, ".ld")
  write_snplist(proxy_candidates, candidate_file)

  # Step1 restricted every source input to the target +/-1 Mb, so the index and
  # all of its clump members already lie inside that window.  Loading only the
  # window instead of the whole ancestry panel changes no reported pair.
  region_args <- if (!is.na(region_chr) && nzchar(region_chr) &&
                     is.finite(region_from_bp) && is.finite(region_to_bp)) {
    c(
      "--chr", as.character(region_chr),
      "--from-bp", format(region_from_bp, scientific = FALSE, trim = TRUE),
      "--to-bp", format(region_to_bp, scientific = FALSE, trim = TRUE)
    )
  } else {
    character()
  }
  args <- c(
    "--bfile", ref_bfile,
    region_args,
    "--r2",
    "--ld-snp", index_snp,
    "--ld-snp-list", candidate_file,
    # PLINK otherwise limits table output to a 10-variant / 1000-kb window.
    "--ld-window", "999999",
    "--ld-window-kb", as.character(CLUMP_KB),
    "--ld-window-r2", "0",
    "--out", label
  )
  exit_status <- suppressWarnings(
    system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)
  )
  if (is.null(exit_status) || as.integer(exit_status) != 0L) {
    return(list(
      status = "plink_member_proxy_ld_failed",
      table = empty_proxy_ld_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0L) {
    return(list(
      status = "plink_member_proxy_ld_no_output",
      table = empty_proxy_ld_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }

  ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) == 0L) {
    return(list(
      status = "plink_member_proxy_ld_unparseable_output",
      table = empty_proxy_ld_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }
  names_upper <- toupper(names(ld))
  a_col <- names(ld)[match("SNP_A", names_upper)]
  b_col <- names(ld)[match("SNP_B", names_upper)]
  r2_col <- names(ld)[match("R2", names_upper)]
  if (anyNA(c(a_col, b_col, r2_col))) {
    return(list(
      status = "plink_member_proxy_ld_missing_required_columns",
      table = empty_proxy_ld_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }

  snp_a <- as.character(ld[[a_col]])
  snp_b <- as.character(ld[[b_col]])
  member_snp <- ifelse(snp_a == index_snp, snp_b,
                       ifelse(snp_b == index_snp, snp_a, NA_character_))
  parsed <- data.table(
    proxy_candidate_snp = as.character(member_snp),
    proxy_r2_1000G = suppressWarnings(as.numeric(ld[[r2_col]]))
  )[
    proxy_candidate_snp %chin% proxy_candidates & is.finite(proxy_r2_1000G)
  ]
  if (nrow(parsed) == 0L) {
    parsed <- empty_proxy_ld_table()
  } else {
    parsed <- parsed[, .(proxy_r2_1000G = max(proxy_r2_1000G)),
                     by = proxy_candidate_snp]
  }
  list(
    status = "completed",
    table = parsed,
    candidate_file = candidate_file, ld_file = ld_file
  )
}

empty_panel_wide_proxy_audit <- function() {
  data.table(
    target_id = character(), source_key = character(), source = character(),
    source_output_stem = character(), clump_panel = character(),
    clump_id = character(),     clump_index_snp = character(),
    proxy_candidate_snp = character(), proxy_condition_snp = character(),
    proxy_r2_1000G = numeric(),
    proxy_candidate_bp = integer(), proxy_abs_distance_to_lead_bp = numeric(),
    proxy_selected = logical(), proxy_min_r2_required = numeric(),
    proxy_ld_reference_panel = character(), proxy_ld_status = character(),
    proxy_ld_file = character()
  )
}

# One row per ancestry-panel variant that is the same variant in KoGES: same
# chromosome, same position, and the same unordered allele pair (strand
# complement allowed).  Matching is done on coordinates and alleles, not on
# rsID, so a panel variant whose rsID differs from the KoGES rsID is still
# usable; only a genuinely different variant is rejected.  Pairs that are not
# 1:1 at a position are dropped rather than guessed.  This is the only pool a
# panel-wide proxy may be drawn from, so a rescued condition SNP is always an
# exact KoGES variant, never a re-encoded one.
# ref_koges_snp_alias is NA when the KoGES rsID equals the panel rsID; this
# keeps a genome-wide map small enough to share across fork workers.
# Position-and-allele index of an ancestry panel, used to find an allele-level
# alias for an external rsID the panel does not carry.
build_panel_position_index <- function(panel_bim) {
  dt <- data.table(
    alias_snp = panel_bim$ref_snp,
    chr_normalized = normalize_chr(panel_bim$ref_chr),
    bp = panel_bim$ref_bp,
    panel_a1 = panel_bim$ref_a1,
    panel_a2 = panel_bim$ref_a2
  )
  dt <- dt[is_valid_snp_id(alias_snp) & !is.na(chr_normalized) & !is.na(bp)]
  setkey(dt, chr_normalized, bp)
  dt
}

# Resolve a raw external CHR:BP+allele record to KoGES without relying on its
# rsID.  This is deliberately stricter than a coordinate-only alias: there
# must be exactly one KoGES BIM variant at the position with the same unordered
# allele pair (strand complement allowed).  It is used only to retain a raw
# external hit that the ancestry 1000G panel does not carry; it never asserts
# that the hit was part of a 1000G-defined clump.
resolve_unique_raw_to_koges_exact <- function(raw_dt, koges_position_index) {
  out <- data.table(
    row_id = seq_len(nrow(raw_dt)),
    direct_KoGES_exact_snp = NA_character_,
    direct_KoGES_exact_status = "missing_external_position_or_alleles"
  )
  if (nrow(raw_dt) == 0L) return(out)

  query <- data.table(
    row_id = seq_len(nrow(raw_dt)),
    chr_normalized = raw_dt$external_chr_normalized,
    bp = raw_dt$external_bp,
    external_effect_allele = raw_dt$external_effect_allele,
    external_other_allele = raw_dt$external_other_allele
  )
  query[, usable := !is.na(chr_normalized) & !is.na(bp) &
    !is.na(external_effect_allele) & !is.na(external_other_allele)]
  if (!any(query$usable)) return(out)

  hits <- koges_position_index[
    query[usable %in% TRUE],
    on = .(chr_normalized, bp),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  if (nrow(hits) == 0L) {
    out[query$usable %in% TRUE,
        direct_KoGES_exact_status := "no_exact_KoGES_allele_match"]
    return(out)
  }
  hits <- hits[allele_pair_matches(
    panel_a1, panel_a2, external_effect_allele, external_other_allele
  ) %in% TRUE]
  if (nrow(hits) == 0L) {
    out[query$usable %in% TRUE,
        direct_KoGES_exact_status := "no_exact_KoGES_allele_match"]
    return(out)
  }
  setorder(hits, row_id, alias_snp)
  resolved <- hits[, .(
    n_exact_matches = uniqueN(alias_snp),
    direct_KoGES_exact_snp = alias_snp[1L]
  ), by = row_id]
  out[query$usable %in% TRUE,
      direct_KoGES_exact_status := "no_exact_KoGES_allele_match"]
  out[resolved$row_id, direct_KoGES_exact_status := fifelse(
    resolved$n_exact_matches == 1L,
    "exact_unique_KoGES_allele_match",
    "ambiguous_exact_KoGES_allele_match"
  )]
  resolved_ok <- resolved[n_exact_matches == 1L]
  if (nrow(resolved_ok) > 0L) {
    out[resolved_ok$row_id,
        direct_KoGES_exact_snp := resolved_ok$direct_KoGES_exact_snp]
  }
  out
}

empty_panel_koges_exact_map <- function() {
  data.table(
    ref_snp = character(), ref_koges_snp_alias = character(),
    ref_chr_normalized = character(), ref_bp = integer()
  )
}

build_panel_koges_exact_map <- function(panel_bim, koges_bim, panel_label) {
  panel_pos <- data.table(
    ref_snp = panel_bim$ref_snp,
    ref_chr_normalized = normalize_chr(panel_bim$ref_chr),
    ref_bp = panel_bim$ref_bp,
    panel_a1 = panel_bim$ref_a1,
    panel_a2 = panel_bim$ref_a2
  )
  panel_pos <- panel_pos[
    is_valid_snp_id(ref_snp) & !is.na(ref_chr_normalized) & !is.na(ref_bp)
  ]
  koges_pos <- data.table(
    koges_snp = koges_bim$ref_snp,
    ref_chr_normalized = normalize_chr(koges_bim$ref_chr),
    ref_bp = koges_bim$ref_bp,
    koges_a1 = koges_bim$ref_a1,
    koges_a2 = koges_bim$ref_a2
  )
  koges_pos <- koges_pos[
    is_valid_snp_id(koges_snp) & !is.na(ref_chr_normalized) & !is.na(ref_bp)
  ]
  if (nrow(panel_pos) == 0L || nrow(koges_pos) == 0L) {
    message("No usable variants to map between 1000G ", panel_label, " and KoGES.")
    return(empty_panel_koges_exact_map())
  }

  joined <- panel_pos[
    koges_pos,
    on = .(ref_chr_normalized, ref_bp),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]
  rm(panel_pos, koges_pos)
  invisible(gc())
  if (nrow(joined) == 0L) {
    return(empty_panel_koges_exact_map())
  }
  joined <- joined[
    allele_pair_matches(panel_a1, panel_a2, koges_a1, koges_a2) %in% TRUE
  ]
  if (nrow(joined) == 0L) {
    return(empty_panel_koges_exact_map())
  }
  joined <- unique(
    joined[, .(ref_snp, koges_snp, ref_chr_normalized, ref_bp)]
  )
  # Ambiguity here would mean silently choosing between two candidate variants.
  n_before <- nrow(joined)
  joined[, n_koges_per_panel := .N, by = ref_snp]
  joined[, n_panel_per_koges := .N, by = koges_snp]
  joined <- joined[n_koges_per_panel == 1L & n_panel_per_koges == 1L]
  n_dropped <- n_before - nrow(joined)

  out <- joined[, .(
    ref_snp,
    ref_koges_snp_alias = fifelse(koges_snp == ref_snp, NA_character_, koges_snp),
    ref_chr_normalized,
    ref_bp
  )]
  setkey(out, ref_snp)
  message(
    "1000G ", panel_label, ": ", nrow(out),
    " variants with exact 1000G-to-KoGES identity (panel-wide proxy pool); ",
    sum(!is.na(out$ref_koges_snp_alias)), " via a different rsID; ",
    n_dropped, " non-1:1 position matches dropped."
  )
  out
}

# Ask PLINK for every panel variant in LD with the clump index above the proxy
# threshold.  The external GWAS input is irrelevant here: the rescue exists
# because a source may simply never have reported the LD partner.
calculate_panel_wide_proxy_ld <- function(index_snp, ref_bfile, work_prefix,
                                          region_chr, region_from_bp,
                                          region_to_bp,
                                          min_r2 = PROXY_MIN_R2) {
  label <- paste0(
    work_prefix, "__", safe_target_label(index_snp), "__panel_wide_proxy_ld"
  )
  ld_file <- paste0(label, ".ld")
  args <- c(
    "--bfile", ref_bfile,
    # Load only the target window.  A whole-panel r^2 scan per unrepresented
    # clump would dominate runtime, and candidates outside the window are
    # rejected anyway.
    "--chr", region_chr,
    "--from-bp", format(region_from_bp, scientific = FALSE, trim = TRUE),
    "--to-bp", format(region_to_bp, scientific = FALSE, trim = TRUE),
    "--r2",
    "--ld-snp", index_snp,
    "--ld-window", "999999",
    "--ld-window-kb", as.character(PANEL_WIDE_PROXY_WINDOW_KB),
    "--ld-window-r2", format(min_r2, scientific = FALSE, trim = TRUE),
    "--out", label
  )
  exit_status <- suppressWarnings(
    system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)
  )
  if (is.null(exit_status) || as.integer(exit_status) != 0L) {
    return(list(
      status = "plink_panel_wide_proxy_ld_failed",
      table = empty_proxy_ld_table(), ld_file = ld_file
    ))
  }
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0L) {
    return(list(
      status = "plink_panel_wide_proxy_ld_no_output",
      table = empty_proxy_ld_table(), ld_file = ld_file
    ))
  }
  ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) == 0L) {
    return(list(
      status = "plink_panel_wide_proxy_ld_unparseable_output",
      table = empty_proxy_ld_table(), ld_file = ld_file
    ))
  }
  names_upper <- toupper(names(ld))
  a_col <- names(ld)[match("SNP_A", names_upper)]
  b_col <- names(ld)[match("SNP_B", names_upper)]
  r2_col <- names(ld)[match("R2", names_upper)]
  if (anyNA(c(a_col, b_col, r2_col))) {
    return(list(
      status = "plink_panel_wide_proxy_ld_missing_required_columns",
      table = empty_proxy_ld_table(), ld_file = ld_file
    ))
  }
  snp_a <- as.character(ld[[a_col]])
  snp_b <- as.character(ld[[b_col]])
  partner <- ifelse(snp_a == index_snp, snp_b,
                    ifelse(snp_b == index_snp, snp_a, NA_character_))
  parsed <- data.table(
    proxy_candidate_snp = partner,
    proxy_r2_1000G = suppressWarnings(as.numeric(ld[[r2_col]]))
  )[
    !is.na(proxy_candidate_snp) & proxy_candidate_snp != index_snp &
      is.finite(proxy_r2_1000G)
  ]
  if (nrow(parsed) == 0L) {
    parsed <- empty_proxy_ld_table()
  } else {
    parsed <- parsed[, .(proxy_r2_1000G = max(proxy_r2_1000G)),
                     by = proxy_candidate_snp]
  }
  list(status = "completed", table = parsed, ld_file = ld_file)
}

empty_panel_wide_group_ld <- function() {
  out <- data.table(
    clump_panel = character(), target_id = character(),
    index_snp = character(), proxy_candidate_snp = character(),
    proxy_r2_1000G = numeric()
  )
  setkey(out, clump_panel, target_id, index_snp)
  out
}

empty_panel_wide_group_status <- function() {
  out <- data.table(
    clump_panel = character(), target_id = character(),
    group_status = character(), n_index_snps = integer(),
    ld_file = character()
  )
  setkey(out, clump_panel, target_id)
  out
}

# Every clump of one target and one ancestry panel needs r^2 against the same
# +/-1 Mb window, so a per-clump PLINK job recomputed the same window up to a
# hundred times.  --ld-snp-list anchors the reported pairs on the clump indices
# and --ld-window-kb / --ld-window-r2 are unchanged, so the retained pairs are
# identical to the per-clump calls; only the number of PLINK invocations drops.
run_one_panel_wide_group_ld <- function(job) {
  panel <- as.character(job$clump_panel)
  target_key <- as.character(job$target_id)
  # index_snps arrives as a one-element list column from the grouped job table.
  index_snps <- unique(as.character(unlist(job$index_snps, use.names = FALSE)))
  index_snps <- index_snps[is_valid_snp_id(index_snps)]
  label <- file.path(
    CLUMP_WORK_DIR,
    paste0(
      "panel_wide_group__", safe_target_label(panel), "__",
      safe_target_label(target_key)
    )
  )
  ld_file <- paste0(label, ".ld")
  fail <- function(status) {
    list(
      status = data.table(
        clump_panel = panel, target_id = target_key, group_status = status,
        n_index_snps = length(index_snps), ld_file = ld_file
      ),
      ld = empty_panel_wide_group_ld()
    )
  }
  if (length(index_snps) == 0L) {
    return(fail("no_index_snp_for_panel_wide_group"))
  }

  snplist_file <- paste0(label, ".snplist")
  write_snplist(index_snps, snplist_file)
  args <- c(
    "--bfile", if (identical(panel, "EAS")) EAS_BFILE else EUR_BFILE,
    # Load only the target window; candidates outside it are rejected anyway.
    "--chr", as.character(job$region_chr),
    "--from-bp", format(job$region_from_bp, scientific = FALSE, trim = TRUE),
    "--to-bp", format(job$region_to_bp, scientific = FALSE, trim = TRUE),
    "--r2",
    "--ld-snp-list", snplist_file,
    "--ld-window", "999999",
    "--ld-window-kb", as.character(PANEL_WIDE_PROXY_WINDOW_KB),
    "--ld-window-r2", format(PROXY_MIN_R2, scientific = FALSE, trim = TRUE),
    "--out", label
  )
  exit_status <- suppressWarnings(
    system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)
  )
  if (is.null(exit_status) || as.integer(exit_status) != 0L) {
    return(fail("plink_panel_wide_proxy_ld_failed"))
  }
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0L) {
    return(fail("plink_panel_wide_proxy_ld_no_output"))
  }
  ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(ld) || nrow(ld) == 0L) {
    return(fail("plink_panel_wide_proxy_ld_unparseable_output"))
  }
  names_upper <- toupper(names(ld))
  a_col <- names(ld)[match("SNP_A", names_upper)]
  b_col <- names(ld)[match("SNP_B", names_upper)]
  r2_col <- names(ld)[match("R2", names_upper)]
  if (anyNA(c(a_col, b_col, r2_col))) {
    return(fail("plink_panel_wide_proxy_ld_missing_required_columns"))
  }

  snp_a <- as.character(ld[[a_col]])
  snp_b <- as.character(ld[[b_col]])
  r2 <- suppressWarnings(as.numeric(ld[[r2_col]]))
  # Two clump indices from different sources can themselves be in high LD, so
  # a pair is emitted once per endpoint that is an index of this group.
  a_is_index <- snp_a %chin% index_snps
  b_is_index <- snp_b %chin% index_snps
  parsed <- rbindlist(list(
    data.table(
      index_snp = snp_a[a_is_index],
      proxy_candidate_snp = snp_b[a_is_index],
      proxy_r2_1000G = r2[a_is_index]
    ),
    data.table(
      index_snp = snp_b[b_is_index],
      proxy_candidate_snp = snp_a[b_is_index],
      proxy_r2_1000G = r2[b_is_index]
    )
  ), use.names = TRUE)
  parsed <- parsed[
    index_snp != proxy_candidate_snp & is.finite(proxy_r2_1000G)
  ]
  if (nrow(parsed) == 0L) {
    out_ld <- empty_panel_wide_group_ld()
  } else {
    parsed <- parsed[, .(proxy_r2_1000G = max(proxy_r2_1000G)),
                     by = .(index_snp, proxy_candidate_snp)]
    out_ld <- parsed[, .(
      clump_panel = panel, target_id = target_key,
      index_snp, proxy_candidate_snp, proxy_r2_1000G
    )]
  }
  list(
    status = data.table(
      clump_panel = panel, target_id = target_key, group_status = "completed",
      n_index_snps = length(index_snps), ld_file = ld_file
    ),
    ld = out_ld
  )
}

# Only a clump whose index is not already an exact KoGES variant can reach the
# panel-wide rescue, so LD is computed for exactly that set of indices.
build_panel_wide_group_ld <- function(clump_membership, lead_lookup,
                                      n_workers = N_PARALLEL) {
  empty <- list(
    ld = empty_panel_wide_group_ld(),
    status = empty_panel_wide_group_status()
  )
  if (!isTRUE(PANEL_WIDE_PROXY_RESCUE) || nrow(clump_membership) == 0L) {
    return(empty)
  }
  if (!isTRUE(PANEL_WIDE_PROXY_GROUP_LD)) {
    message(
      "PANEL_WIDE_PROXY_GROUP_LD=FALSE: every clump runs its own panel-wide ",
      "PLINK job (comparison mode)."
    )
    return(empty)
  }
  idx <- clump_membership[member_is_index %in% TRUE, .(
    target_id = as.character(target_id),
    clump_panel = as.character(clump_panel),
    clump_index_snp = as.character(clump_index_snp),
    index_exact = is_exact_identity_status(member_identity_status) %in% TRUE &
      member_KoGES_status %chin% "eligible_if_selected"
  )]
  idx <- unique(idx[index_exact %in% FALSE & is_valid_snp_id(clump_index_snp)])
  if (nrow(idx) == 0L) return(empty)

  win <- lead_lookup[.(idx$target_id)]
  idx[, `:=`(
    region_chr = as.character(win$target_chr_normalized),
    lead_bp = as.numeric(win$koges_bp)
  )]
  idx <- idx[!is.na(region_chr) & is.finite(lead_bp)]
  if (nrow(idx) == 0L) return(empty)

  jobs <- idx[, .(
    index_snps = list(unique(clump_index_snp)),
    region_from_bp = max(1, lead_bp[[1L]] - NOVELTY_WINDOW_BP),
    region_to_bp = lead_bp[[1L]] + NOVELTY_WINDOW_BP
  ), by = .(clump_panel, target_id, region_chr)]
  message(
    "Panel-wide proxy LD: ", nrow(jobs), " target-panel PLINK jobs for ",
    nrow(idx), " clump indices (one window per target-panel instead of one ",
    "per clump)."
  )
  res <- parallel_lapply(
    seq_len(nrow(jobs)),
    function(i) run_one_panel_wide_group_ld(as.list(jobs[i])),
    n_workers = n_workers
  )
  status <- rbindlist(lapply(res, function(x) x$status),
                      use.names = TRUE, fill = TRUE)
  ld <- rbindlist(lapply(res, function(x) x$ld), use.names = TRUE, fill = TRUE)
  if (nrow(status) == 0L) status <- empty_panel_wide_group_status()
  if (nrow(ld) == 0L) ld <- empty_panel_wide_group_ld()
  setkey(status, clump_panel, target_id)
  setkey(ld, clump_panel, target_id, index_snp)
  failed <- status[group_status != "completed"]
  if (nrow(failed) > 0L) {
    message(
      "  ", nrow(failed), " target-panel panel-wide LD job(s) did not complete; ",
      "affected clumps fall back to a per-clump PLINK call."
    )
  }
  list(ld = ld, status = status)
}

# Read the precomputed (target, panel) LD table.  A missing group record falls
# back to the original per-clump PLINK call so a grouping gap cannot silently
# drop a rescue that the per-clump path would have found.
lookup_panel_wide_proxy_ld <- function(index_snp, clump_panel, target_id,
                                       ref_bfile, work_prefix,
                                       region_chr, region_from_bp,
                                       region_to_bp) {
  per_clump <- function() {
    calculate_panel_wide_proxy_ld(
      index_snp = index_snp, ref_bfile = ref_bfile, work_prefix = work_prefix,
      region_chr = region_chr, region_from_bp = region_from_bp,
      region_to_bp = region_to_bp
    )
  }
  if (!exists("panel_wide_group_status", inherits = TRUE) ||
      nrow(panel_wide_group_status) == 0L) {
    return(per_clump())
  }
  # In x[i] the i expression is evaluated with x's columns in scope, so the
  # query values must not reuse the key column names.
  panel_key <- as.character(clump_panel)
  target_key <- as.character(target_id)
  index_key <- as.character(index_snp)
  st <- panel_wide_group_status[.(panel_key, target_key), nomatch = 0L]
  if (nrow(st) != 1L) return(per_clump())
  ld_file <- as.character(st$ld_file[[1L]])
  group_status <- as.character(st$group_status[[1L]])
  if (!identical(group_status, "completed")) {
    return(per_clump())
  }
  hit <- panel_wide_group_ld[
    .(panel_key, target_key, index_key), nomatch = 0L
  ]
  tbl <- if (nrow(hit) == 0L) {
    empty_proxy_ld_table()
  } else {
    hit[, .(proxy_candidate_snp, proxy_r2_1000G)]
  }
  list(status = "completed", table = tbl, ld_file = ld_file)
}

# A proxy chosen in the EUR panel only proves high LD in Europeans.  KoGES is
# Korean, so tagging of the same external signal is not guaranteed.  Re-estimate
# the index--proxy r^2 in 1000G EAS only when optional transferability QC is
# requested. A EUR-only proxy remains a EUR-supporting covariate; it is never
# part of the EAS-primary condition set.
empty_cross_panel_proxy_table <- function() {
  data.table(
    clump_id = character(), proxy_r2_1000G_EAS = numeric(),
    proxy_cross_panel_status = character()
  )
}

CROSS_PANEL_CONFIRMED <- "confirmed_in_EAS_r2_ge_threshold"

# Identity of a clump index between the two panels.  The proxy side is checked
# with the exact-identity maps instead, because a proxy is always an exact KoGES
# variant while a clump index need not be.
build_index_identity_across_panels <- function(index_snps, eur_bim, eas_bim) {
  index_snps <- unique(as.character(index_snps))
  index_snps <- index_snps[is_valid_snp_id(index_snps)]
  if (length(index_snps) == 0L) {
    out <- data.table(snp = character(), same_variant_in_EAS = logical())
    setkey(out, snp)
    return(out)
  }
  e <- match(index_snps, eur_bim$ref_snp)
  a <- match(index_snps, eas_bim$ref_snp)
  same <- !is.na(e) & !is.na(a) &
    normalize_chr(eur_bim$ref_chr[e]) == normalize_chr(eas_bim$ref_chr[a]) &
    eur_bim$ref_bp[e] == eas_bim$ref_bp[a] &
    allele_pair_matches(
      eur_bim$ref_a1[e], eur_bim$ref_a2[e],
      eas_bim$ref_a1[a], eas_bim$ref_a2[a]
    ) %in% TRUE
  same[is.na(same)] <- FALSE
  out <- data.table(snp = index_snps, same_variant_in_EAS = same)
  setkey(out, snp)
  out
}

unordered_pair_key <- function(a, b) {
  a <- as.character(a)
  b <- as.character(b)
  paste(pmin(a, b), pmax(a, b), sep = "\t")
}

# One PLINK job per chromosome over the union of index/proxy IDs, rather than one
# job per pair.
calculate_eas_pair_r2 <- function(pairs, work_dir) {
  if (nrow(pairs) == 0L) return(list(r2 = numeric(), status = character()))
  r2_by_pair <- numeric()
  status_by_pair <- character()
  for (ch in sort(unique(pairs$chr))) {
    sub <- pairs[chr == ch]
    ids <- unique(c(sub$index_snp, sub$proxy_snp))
    label <- file.path(work_dir, paste0("cross_panel_EAS_ld_chr", safe_target_label(ch)))
    snplist_file <- paste0(label, ".snplist")
    ld_file <- paste0(label, ".ld")
    write_snplist(ids, snplist_file)
    args <- c(
      "--bfile", EAS_BFILE,
      "--chr", ch,
      "--extract", snplist_file,
      "--r2",
      "--ld-window", "999999",
      "--ld-window-kb", as.character(PANEL_WIDE_PROXY_WINDOW_KB),
      "--ld-window-r2", "0",
      "--out", label
    )
    exit_status <- suppressWarnings(
      system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)
    )
    keys <- unordered_pair_key(sub$index_snp, sub$proxy_snp)
    if (is.null(exit_status) || as.integer(exit_status) != 0L ||
        !file.exists(ld_file) || file.info(ld_file)$size == 0L) {
      r2_by_pair[keys] <- NA_real_
      status_by_pair[keys] <- "not_confirmed_EAS_LD_calculation_failed"
      next
    }
    ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
    if (is.null(ld) || nrow(ld) == 0L) {
      r2_by_pair[keys] <- NA_real_
      status_by_pair[keys] <- "not_confirmed_EAS_LD_calculation_failed"
      next
    }
    names_upper <- toupper(names(ld))
    a_col <- names(ld)[match("SNP_A", names_upper)]
    b_col <- names(ld)[match("SNP_B", names_upper)]
    r2_col <- names(ld)[match("R2", names_upper)]
    if (anyNA(c(a_col, b_col, r2_col))) {
      r2_by_pair[keys] <- NA_real_
      status_by_pair[keys] <- "not_confirmed_EAS_LD_calculation_failed"
      next
    }
    observed <- data.table(
      pair_key = unordered_pair_key(ld[[a_col]], ld[[b_col]]),
      r2 = suppressWarnings(as.numeric(ld[[r2_col]]))
    )[is.finite(r2)]
    observed <- observed[, .(r2 = max(r2)), by = pair_key]
    hit <- observed$r2[match(keys, observed$pair_key)]
    r2_by_pair[keys] <- hit
    status_by_pair[keys] <- fcase(
      is.na(hit), "not_confirmed_pair_absent_from_EAS_LD_output",
      hit >= PROXY_MIN_R2, CROSS_PANEL_CONFIRMED,
      default = "not_confirmed_r2_below_threshold_in_EAS"
    )
  }
  list(r2 = r2_by_pair, status = status_by_pair)
}

verify_eur_ld_proxies_in_eas <- function(reps, index_identity, exact_maps,
                                         work_dir) {
  candidates <- reps[
    clump_panel == "EUR" &
      representation_evidence == "ld_proxy_r2_ge_threshold" &
      representative_usable_in_KoGES %in% TRUE &
      is_valid_snp_id(clump_index_snp) &
      is_valid_snp_id(selected_clump_member_snp),
    .(
      clump_id = as.character(clump_id),
      index_snp = as.character(clump_index_snp),
      proxy_snp = as.character(selected_clump_member_snp)
    )
  ]
  if (nrow(candidates) == 0L) return(empty_cross_panel_proxy_table())

  idx_match <- index_identity[.(candidates$index_snp)]
  candidates[, index_same_variant_in_EAS := idx_match$same_variant_in_EAS %in% TRUE]

  eur_map <- exact_maps$EUR
  eas_map <- exact_maps$EAS
  proxy_eur <- if (is.null(eur_map) || nrow(eur_map) == 0L) {
    NULL
  } else {
    eur_map[.(candidates$proxy_snp)]
  }
  proxy_eas <- if (is.null(eas_map) || nrow(eas_map) == 0L) {
    NULL
  } else {
    eas_map[.(candidates$proxy_snp)]
  }
  if (is.null(proxy_eur) || is.null(proxy_eas)) {
    candidates[, proxy_same_variant_in_EAS := FALSE]
    candidates[, chr := NA_character_]
  } else {
    candidates[, proxy_same_variant_in_EAS :=
      !is.na(proxy_eur$ref_bp) & !is.na(proxy_eas$ref_bp) &
        proxy_eur$ref_chr_normalized == proxy_eas$ref_chr_normalized &
        proxy_eur$ref_bp == proxy_eas$ref_bp]
    candidates[, chr := as.character(proxy_eas$ref_chr_normalized)]
  }

  candidates[, proxy_cross_panel_status := fcase(
    !index_same_variant_in_EAS, "not_confirmed_index_not_same_variant_in_EAS",
    !proxy_same_variant_in_EAS, "not_confirmed_proxy_not_same_variant_in_EAS",
    is.na(chr), "not_confirmed_proxy_not_same_variant_in_EAS",
    default = NA_character_
  )]
  candidates[, proxy_r2_1000G_EAS := NA_real_]

  testable <- candidates[is.na(proxy_cross_panel_status)]
  if (nrow(testable) > 0L) {
    pairs <- unique(testable[, .(chr, index_snp, proxy_snp)])
    message(
      "Confirming ", nrow(pairs), " EUR-panel LD proxy pair(s) in 1000G EAS."
    )
    res <- calculate_eas_pair_r2(pairs, work_dir)
    keys <- unordered_pair_key(candidates$index_snp, candidates$proxy_snp)
    fill <- is.na(candidates$proxy_cross_panel_status)
    candidates[fill, proxy_r2_1000G_EAS := unname(res$r2[keys[fill]])]
    candidates[fill, proxy_cross_panel_status := unname(res$status[keys[fill]])]
  }
  candidates[is.na(proxy_cross_panel_status),
    proxy_cross_panel_status := "not_confirmed_EAS_LD_calculation_failed"]

  candidates[, .(clump_id, proxy_r2_1000G_EAS, proxy_cross_panel_status)]
}

# ------------------------------------------------------------------------------
# Cross-source KoGES-LD representative selection
# ------------------------------------------------------------------------------

split_semicolon_values <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (length(x) == 0L) return(character())
  unique(unlist(strsplit(x, ";", fixed = TRUE), use.names = FALSE))
}

n_distinct_semicolon_values <- function(x) {
  length(split_semicolon_values(x))
}

# True only when both source sets are non-empty and share no source.
# Same multi-source key (e.g. BBJ;TWB vs BBJ;TWB) is NOT cross-source.
source_sets_are_cross_source <- function(x, y) {
  sx <- split_semicolon_values(x)
  sy <- split_semicolon_values(y)
  length(sx) > 0L &&
    length(sy) > 0L &&
    length(intersect(sx, sy)) == 0L
}

is_1000G_clump_representative <- function(n_clumps) {
  fcoalesce(as.integer(n_clumps), 0L) > 0L
}

koges_ld_grouping_rule <- function() {
  paste0(
    "pairwise_KoGES_r2>=",
    format(CROSS_SOURCE_REP_R2, scientific = FALSE, trim = TRUE),
    "; no_same_source_fold_between_1000G_clump_reps;",
    "direct_KoGES_rescue_may_fold_same_source"
  )
}

empty_cross_source_ld_pair_table <- function() {
  data.table(
    target_id = character(), condition_snp_a = character(),
    condition_snp_b = character(), KoGES_r2 = numeric(),
    selected_source_keys_a = character(), selected_source_keys_b = character(),
    pair_has_different_sources = logical(),
    pair_a_is_1000G_clump_rep = logical(),
    pair_b_is_1000G_clump_rep = logical(),
    pair_blocked_same_source_clump_reps = logical(),
    pair_meets_cross_source_r2_threshold = logical(),
    pair_eligible_for_ld_collapse = logical(),
    cross_source_ld_status = character(), cross_source_ld_candidate_file = character(),
    cross_source_ld_file = character()
  )
}

empty_cross_source_ld_member_audit <- function() {
  data.table(
    target_id = character(), condition_snp = character(),
    selected_source_keys = character(), selected_sources = character(),
    selected_panels = character(), representative_types = character(),
    source_p_min = numeric(), n_external_clumps_represented = integer(),
    direct_KoGES_rescue = logical(),
    cross_source_ld_group_id = character(), cross_source_ld_group_size = integer(),
    cross_source_ld_group_member_snps = character(),
    cross_source_ld_is_representative = logical(),
    cross_source_ld_representative_condition_snp = character(),
    cross_source_ld_r2_to_representative = numeric(),
    cross_source_ld_grouping_status = character(),
    cross_source_ld_grouping_rule = character()
  )
}

empty_cross_source_ld_target_summary <- function() {
  data.table(
    target_id = character(),
    n_condition_snps_pre_cross_source_KoGES_LD = integer(),
    n_condition_snps_post_cross_source_KoGES_LD = integer(),
    n_condition_snps_collapsed_cross_source_KoGES_LD = integer(),
    n_cross_source_KoGES_LD_representative_groups = integer(),
    n_cross_source_KoGES_high_LD_pairs = integer(),
    n_KoGES_LD_collapse_eligible_pairs = integer(),
    cross_source_KoGES_LD_grouping_status = character(),
    cross_source_KoGES_LD_candidate_file = character(),
    cross_source_KoGES_LD_file = character()
  )
}

# Calculate all pairwise LD values among the final source-level representatives
# for one target.  Only SNP IDs already confirmed in the KoGES BIM are passed
# here.  The LD window is deliberately wide enough to cover both edges of the
# target +/- NOVELTY_WINDOW_BP interval.
calculate_cross_source_KoGES_ld <- function(
    target_id, candidates, analysis_label = "all_sources") {
  candidate_snps <- sort(unique(as.character(candidates$condition_snp)))
  candidate_snps <- candidate_snps[is_valid_snp_id(candidate_snps)]
  label <- file.path(
    CROSS_SOURCE_LD_WORK_DIR,
    paste0(
      safe_target_label(target_id), "__",
      gsub("[^A-Za-z0-9._-]+", "_", analysis_label),
      "__cross_source_KoGES_LD"
    )
  )
  candidate_file <- paste0(label, ".snplist")
  ld_file <- paste0(label, ".ld")
  write_snplist(candidate_snps, candidate_file)

  args <- c(
    "--bfile", KOGES_BFILE,
    "--extract", candidate_file,
    "--r2",
    "--ld-window", "999999",
    "--ld-window-kb", as.character(CROSS_SOURCE_LD_WINDOW_KB),
    "--ld-window-r2", "0",
    "--out", label
  )
  exit_status <- suppressWarnings(
    system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)
  )
  if (is.null(exit_status) || as.integer(exit_status) != 0L) {
    return(list(
      status = "plink_cross_source_KoGES_LD_failed",
      table = empty_cross_source_ld_pair_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }
  if (!file.exists(ld_file) || file.info(ld_file)$size == 0L) {
    return(list(
      status = "plink_cross_source_KoGES_LD_no_output",
      table = empty_cross_source_ld_pair_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }

  ld <- tryCatch(fread(ld_file, showProgress = FALSE), error = function(e) NULL)
  if (is.null(ld)) {
    return(list(
      status = "plink_cross_source_KoGES_LD_unparseable_output",
      table = empty_cross_source_ld_pair_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }
  names_upper <- toupper(names(ld))
  a_col <- names(ld)[match("SNP_A", names_upper)]
  b_col <- names(ld)[match("SNP_B", names_upper)]
  r2_col <- names(ld)[match("R2", names_upper)]
  if (anyNA(c(a_col, b_col, r2_col))) {
    return(list(
      status = "plink_cross_source_KoGES_LD_missing_required_columns",
      table = empty_cross_source_ld_pair_table(),
      candidate_file = candidate_file, ld_file = ld_file
    ))
  }

  pairs <- data.table(
    condition_snp_a = as.character(ld[[a_col]]),
    condition_snp_b = as.character(ld[[b_col]]),
    KoGES_r2 = suppressWarnings(as.numeric(ld[[r2_col]]))
  )[
    condition_snp_a %chin% candidate_snps &
      condition_snp_b %chin% candidate_snps &
      condition_snp_a != condition_snp_b &
      is.finite(KoGES_r2)
  ]
  if (nrow(pairs) > 0L) {
    pairs[, `:=`(
      snp_low = pmin(condition_snp_a, condition_snp_b),
      snp_high = pmax(condition_snp_a, condition_snp_b)
    )]
    pairs <- pairs[, .(KoGES_r2 = max(KoGES_r2)), by = .(snp_low, snp_high)]
    setnames(pairs, c("snp_low", "snp_high"), c("condition_snp_a", "condition_snp_b"))
  } else {
    pairs <- data.table(
      condition_snp_a = character(), condition_snp_b = character(),
      KoGES_r2 = numeric()
    )
  }

  source_a_match <- match(pairs$condition_snp_a, candidates$condition_snp)
  source_b_match <- match(pairs$condition_snp_b, candidates$condition_snp)
  n_clumps <- if ("n_external_clumps_represented" %in% names(candidates)) {
    fcoalesce(as.integer(candidates$n_external_clumps_represented), 0L)
  } else {
    integer(nrow(candidates))
  }
  pairs[, `:=`(
    target_id = rep(as.character(target_id), nrow(pairs)),
    selected_source_keys_a = candidates$selected_source_keys[source_a_match],
    selected_source_keys_b = candidates$selected_source_keys[source_b_match],
    pair_a_is_1000G_clump_rep = is_1000G_clump_representative(n_clumps[source_a_match]),
    pair_b_is_1000G_clump_rep = is_1000G_clump_representative(n_clumps[source_b_match])
  )]
  if (nrow(pairs) == 0L) {
    pairs[, pair_has_different_sources := logical()]
  } else {
    pairs[, pair_has_different_sources := mapply(
      source_sets_are_cross_source,
      selected_source_keys_a, selected_source_keys_b,
      USE.NAMES = FALSE
    )]
  }
  pairs[, pair_blocked_same_source_clump_reps :=
    pair_a_is_1000G_clump_rep %in% TRUE &
      pair_b_is_1000G_clump_rep %in% TRUE &
      !(pair_has_different_sources %in% TRUE)]
  pairs[, pair_meets_cross_source_r2_threshold :=
    pair_has_different_sources & KoGES_r2 >= CROSS_SOURCE_REP_R2]
  pairs[, pair_eligible_for_ld_collapse :=
    is.finite(KoGES_r2) &
      KoGES_r2 >= CROSS_SOURCE_REP_R2 &
      !(pair_blocked_same_source_clump_reps %in% TRUE)]
  pairs[, `:=`(
    cross_source_ld_status = rep("completed", nrow(pairs)),
    cross_source_ld_candidate_file = rep(candidate_file, nrow(pairs)),
    cross_source_ld_file = rep(ld_file, nrow(pairs))
  )]
  setcolorder(pairs, names(empty_cross_source_ld_pair_table()))
  list(
    status = "completed",
    table = pairs,
    candidate_file = candidate_file, ld_file = ld_file
  )
}

# The selected representative is deterministic but does not rank P values from
# different studies as though they were one meta-analysis.  A 1000G-defined
# clump representative is always preferred over a direct-KoGES rescue SNP:
# rescue has no ancestry-panel LD definition and is only an auxiliary
# covariate.  Within 1000G-defined SNPs, a direct clump index outranks an
# LD proxy; then more clumps/sources; P is only a final tie-breaker.
# Order: direct 1000G clump index > 1000G clump proxy > direct KoGES rescue.
# A non-representative is dropped only when it has a direct pairwise KoGES
# r^2 >= cutoff to that chosen representative.  Two 1000G clump representatives
# whose source sets overlap never fold: 1000G already defined them as distinct
# within-source signals.  Direct-KoGES rescue SNPs have no such definition, so
# rescue-rescue and rescue-clump pairs may fold at the same r^2 even when
# sources overlap.  There is still no transitive graph-only collapse:
# membership requires direct r^2 to the chosen representative.
assign_cross_source_ld_groups <- function(
    candidates, pair_table, grouping_status, analysis_label = "all_sources") {
  candidates <- copy(candidates)
  if (!"n_external_clumps_represented" %in% names(candidates)) {
    candidates[, n_external_clumps_represented := 0L]
  }
  candidates[, `:=`(
    .cross_is_1000G_clump_rep = is_1000G_clump_representative(
      n_external_clumps_represented
    ),
    .cross_is_direct_clump_index = !is.na(representative_types) & grepl(
      "(^|;)direct_clump_index(;|$)",
      representative_types, perl = TRUE
    ),
    .cross_n_source_keys = vapply(
      selected_source_keys, n_distinct_semicolon_values, integer(1L)
    ),
    .cross_n_external_clumps = fcoalesce(
      as.integer(n_external_clumps_represented), 0L
    ),
    .cross_source_p_order = fifelse(
      is.finite(source_p_min), as.numeric(source_p_min), Inf
    )
  )]
  setorder(
    candidates,
    -.cross_is_1000G_clump_rep,
    -.cross_is_direct_clump_index,
    -.cross_n_external_clumps,
    -.cross_n_source_keys,
    .cross_source_p_order,
    condition_snp
  )

  candidates[, `:=`(
    cross_source_ld_group_id = NA_character_,
    cross_source_ld_group_size = NA_integer_,
    cross_source_ld_group_member_snps = NA_character_,
    cross_source_ld_is_representative = FALSE,
    cross_source_ld_representative_condition_snp = NA_character_,
    cross_source_ld_r2_to_representative = NA_real_,
    cross_source_ld_grouping_status = grouping_status,
    cross_source_ld_grouping_rule = koges_ld_grouping_rule()
  )]

  pair_lookup <- if (nrow(pair_table) == 0L) {
    numeric()
  } else {
    setNames(
      as.numeric(pair_table$KoGES_r2),
      paste(pair_table$condition_snp_a, pair_table$condition_snp_b, sep = "\t")
    )
  }
  lookup_pair_r2 <- function(anchor_snp, other_snps) {
    keys <- paste(pmin(anchor_snp, other_snps), pmax(anchor_snp, other_snps), sep = "\t")
    as.numeric(pair_lookup[keys])
  }

  unassigned <- rep(TRUE, nrow(candidates))
  group_number <- 0L
  while (any(unassigned)) {
    anchor_index <- which(unassigned)[[1L]]
    remaining_indices <- which(unassigned)
    remaining_snps <- candidates$condition_snp[remaining_indices]
    r2_to_anchor <- lookup_pair_r2(
      candidates$condition_snp[[anchor_index]], remaining_snps
    )
    # Grow the group greedily in current priority order.  Direct r^2 to the
    # anchor is required.  Same-source fold is blocked only when the candidate
    # and at least one already-included member are both 1000G clump
    # representatives with overlapping sources.  A rescue SNP (0 represented
    # 1000G clumps) may join even if it shares BBJ/TWB/TPMI with the group.
    group_indices <- anchor_index
    for (cand_pos in seq_along(remaining_indices)) {
      cand_index <- remaining_indices[[cand_pos]]
      if (cand_index == anchor_index) next
      r2_val <- r2_to_anchor[[cand_pos]]
      if (!is.finite(r2_val) || r2_val < CROSS_SOURCE_REP_R2) next
      cand_sources <- split_semicolon_values(
        candidates$selected_source_keys[[cand_index]]
      )
      if (length(cand_sources) == 0L) next
      if (candidates$.cross_is_1000G_clump_rep[[cand_index]]) {
        blocked_same_source_clump <- FALSE
        for (gi in group_indices) {
          if (!candidates$.cross_is_1000G_clump_rep[[gi]]) next
          group_member_sources <- split_semicolon_values(
            candidates$selected_source_keys[[gi]]
          )
          if (length(intersect(group_member_sources, cand_sources)) > 0L) {
            blocked_same_source_clump <- TRUE
            break
          }
        }
        if (blocked_same_source_clump) next
      }
      group_indices <- c(group_indices, cand_index)
    }
    group_r2 <- r2_to_anchor[match(group_indices, remaining_indices)]
    group_r2[group_indices == anchor_index] <- 1
    group_number <- group_number + 1L
    group_id <- paste0(
      candidates$target_id[[anchor_index]],
      "|", analysis_label, "|cross_source_KoGES_LD_group_",
      sprintf("%03d", group_number)
    )
    group_member_snps <- collapse_unique(candidates$condition_snp[group_indices])
    candidates[group_indices, `:=`(
      cross_source_ld_group_id = group_id,
      cross_source_ld_group_size = as.integer(length(group_indices)),
      cross_source_ld_group_member_snps = group_member_snps,
      cross_source_ld_representative_condition_snp =
        candidates$condition_snp[[anchor_index]],
      cross_source_ld_r2_to_representative = as.numeric(group_r2),
      cross_source_ld_is_representative =
        condition_snp == candidates$condition_snp[[anchor_index]]
    )]
    unassigned[group_indices] <- FALSE
  }

  candidates[, c(
    ".cross_is_1000G_clump_rep", ".cross_is_direct_clump_index",
    ".cross_n_source_keys", ".cross_n_external_clumps",
    ".cross_source_p_order"
  ) := NULL]
  candidates
}

collapse_one_cross_source_ld_group <- function(group_candidates) {
  group_candidates <- copy(group_candidates)
  representative <- group_candidates[
    cross_source_ld_is_representative %in% TRUE
  ]
  if (nrow(representative) != 1L) {
    stop("Each cross-source KoGES-LD group must have exactly one representative.")
  }
  representative <- copy(representative[1L])
  representative[, `:=`(
    reported_in_sources = collapse_source_sets(group_candidates$reported_in_sources),
    reported_in_layers = collapse_source_sets(group_candidates$reported_in_layers),
    source_output_stems = collapse_source_sets(group_candidates$source_output_stems),
    source_p_min = min_finite(group_candidates$source_p_min),
    source_p_values = collapse_source_sets(group_candidates$source_p_values),
    source_abs_distance_min_bp = min_finite(
      group_candidates$source_abs_distance_min_bp
    ),
    source_abs_distance_max_bp = max_finite(
      group_candidates$source_abs_distance_max_bp
    ),
    source_variant_ids = collapse_source_sets(group_candidates$source_variant_ids),
    n_raw_source_hits = as.integer(sum(group_candidates$n_raw_source_hits, na.rm = TRUE)),
    selected_in_EAS = any(group_candidates$selected_in_EAS %in% TRUE),
    selected_in_EUR = any(group_candidates$selected_in_EUR %in% TRUE),
    selected_panels = collapse_source_sets(group_candidates$selected_panels),
    selected_sources = collapse_source_sets(group_candidates$selected_sources),
    selected_source_output_stems = collapse_source_sets(
      group_candidates$selected_source_output_stems
    ),
    selected_source_keys = collapse_source_sets(group_candidates$selected_source_keys),
    n_distinct_selected_source_keys = as.integer(n_distinct_semicolon_values(
      group_candidates$selected_source_keys
    )),
    panel_clump_p_min = min_finite(group_candidates$panel_clump_p_min),
    external_clump_index_snps = collapse_source_sets(
      group_candidates$external_clump_index_snps
    ),
    selected_clump_member_snps = collapse_source_sets(
      group_candidates$selected_clump_member_snps
    ),
    representative_types = collapse_source_sets(group_candidates$representative_types),
    representative_statuses = collapse_source_sets(
      group_candidates$representative_statuses
    ),
    proxy_r2_1000G_min = min_finite(group_candidates$proxy_r2_1000G_min),
    proxy_r2_1000G_values = collapse_source_sets(
      group_candidates$proxy_r2_1000G_values
    ),
    n_external_clumps_represented = as.integer(sum(
      group_candidates$n_external_clumps_represented, na.rm = TRUE
    )),
    direct_KoGES_rescue = any(group_candidates$direct_KoGES_rescue %in% TRUE),
    n_direct_KoGES_rescued_hits = as.integer(sum(
      group_candidates$n_direct_KoGES_rescued_hits, na.rm = TRUE
    )),
    direct_KoGES_rescue_statuses = collapse_source_sets(
      group_candidates$direct_KoGES_rescue_statuses
    )
  )]
  representative
}

make_cross_source_ld_member_audit <- function(grouped_candidates) {
  if (nrow(grouped_candidates) == 0L) return(empty_cross_source_ld_member_audit())
  if (!"direct_KoGES_rescue" %in% names(grouped_candidates)) {
    grouped_candidates[, direct_KoGES_rescue := FALSE]
  }
  grouped_candidates[, .(
    target_id, condition_snp, selected_source_keys, selected_sources,
    selected_panels, representative_types, source_p_min,
    n_external_clumps_represented, direct_KoGES_rescue,
    cross_source_ld_group_id, cross_source_ld_group_size,
    cross_source_ld_group_member_snps,
    cross_source_ld_is_representative,
    cross_source_ld_representative_condition_snp,
    cross_source_ld_r2_to_representative,
    cross_source_ld_grouping_status, cross_source_ld_grouping_rule
  )]
}

run_one_cross_source_ld_grouping <- function(
    candidates, analysis_label = "all_sources") {
  candidates <- copy(candidates)
  if (nrow(candidates) == 0L) {
    return(list(
      collapsed_candidates = candidates,
      member_audit = empty_cross_source_ld_member_audit(),
      pair_audit = empty_cross_source_ld_pair_table(),
      target_summary = empty_cross_source_ld_target_summary()
    ))
  }
  target_id <- as.character(candidates$target_id[[1L]])
  if (uniqueN(candidates$target_id) != 1L) {
    stop("A cross-source KoGES-LD job must contain exactly one target_id.")
  }

  if (nrow(candidates) == 1L) {
    grouped <- assign_cross_source_ld_groups(
      candidates, empty_cross_source_ld_pair_table(), "not_needed_singleton",
      analysis_label = analysis_label
    )
    collapsed <- collapse_one_cross_source_ld_group(grouped)
    return(list(
      collapsed_candidates = collapsed,
      member_audit = make_cross_source_ld_member_audit(grouped),
      pair_audit = empty_cross_source_ld_pair_table(),
      target_summary = data.table(
        target_id = target_id,
        n_condition_snps_pre_cross_source_KoGES_LD = 1L,
        n_condition_snps_post_cross_source_KoGES_LD = 1L,
        n_condition_snps_collapsed_cross_source_KoGES_LD = 0L,
        n_cross_source_KoGES_LD_representative_groups = 1L,
        n_cross_source_KoGES_high_LD_pairs = 0L,
        n_KoGES_LD_collapse_eligible_pairs = 0L,
        cross_source_KoGES_LD_grouping_status = "not_needed_singleton",
        cross_source_KoGES_LD_candidate_file = NA_character_,
        cross_source_KoGES_LD_file = NA_character_
      )
    ))
  }

  ld_result <- calculate_cross_source_KoGES_ld(
    target_id, candidates, analysis_label = analysis_label
  )
  if (!identical(ld_result$status, "completed")) {
    candidates[, `:=`(
      cross_source_ld_group_id = NA_character_,
      cross_source_ld_group_size = NA_integer_,
      cross_source_ld_group_member_snps = NA_character_,
      cross_source_ld_is_representative = NA,
      cross_source_ld_representative_condition_snp = NA_character_,
      cross_source_ld_r2_to_representative = NA_real_,
      cross_source_ld_grouping_status = ld_result$status,
      cross_source_ld_grouping_rule = koges_ld_grouping_rule()
    )]
    return(list(
      collapsed_candidates = candidates,
      member_audit = make_cross_source_ld_member_audit(candidates),
      pair_audit = ld_result$table,
      target_summary = data.table(
        target_id = target_id,
        n_condition_snps_pre_cross_source_KoGES_LD = as.integer(nrow(candidates)),
        n_condition_snps_post_cross_source_KoGES_LD = as.integer(nrow(candidates)),
        n_condition_snps_collapsed_cross_source_KoGES_LD = 0L,
        n_cross_source_KoGES_LD_representative_groups = NA_integer_,
        n_cross_source_KoGES_high_LD_pairs = NA_integer_,
        n_KoGES_LD_collapse_eligible_pairs = NA_integer_,
        cross_source_KoGES_LD_grouping_status = ld_result$status,
        cross_source_KoGES_LD_candidate_file = ld_result$candidate_file,
        cross_source_KoGES_LD_file = ld_result$ld_file
      )
    ))
  }

  grouped <- assign_cross_source_ld_groups(
    candidates, ld_result$table, "completed", analysis_label = analysis_label
  )
  groups <- split(seq_len(nrow(grouped)), grouped$cross_source_ld_group_id)
  collapsed <- rbindlist(
    lapply(groups, function(i) collapse_one_cross_source_ld_group(grouped[i])),
    use.names = TRUE, fill = TRUE
  )
  list(
    collapsed_candidates = collapsed,
    member_audit = make_cross_source_ld_member_audit(grouped),
    pair_audit = ld_result$table,
    target_summary = data.table(
      target_id = target_id,
      n_condition_snps_pre_cross_source_KoGES_LD = as.integer(nrow(grouped)),
      n_condition_snps_post_cross_source_KoGES_LD = as.integer(nrow(collapsed)),
      n_condition_snps_collapsed_cross_source_KoGES_LD = as.integer(
        nrow(grouped) - nrow(collapsed)
      ),
      n_cross_source_KoGES_LD_representative_groups = as.integer(nrow(collapsed)),
      n_cross_source_KoGES_high_LD_pairs = as.integer(sum(
        ld_result$table$pair_meets_cross_source_r2_threshold %in% TRUE
      )),
      n_KoGES_LD_collapse_eligible_pairs = as.integer(sum(
        ld_result$table$pair_eligible_for_ld_collapse %in% TRUE
      )),
      cross_source_KoGES_LD_grouping_status = "completed",
      cross_source_KoGES_LD_candidate_file = ld_result$candidate_file,
      cross_source_KoGES_LD_file = ld_result$ld_file
    )
  )
}

# A job is one target x original external-GWAS source.  Do not pool BBJ/TWB/
# TPMI or UKB/FinnGen/MVP before clumping: their P values are not one common
# GWAS ranking.
run_one_clump_job <- function(job_row, panel_candidates_dt) {
  this_target_id <- as.character(job_row$target_id)
  this_source_key <- as.character(job_row$source_key)
  this_panel <- as.character(job_row$clump_panel)
  ref_bfile <- if (this_panel == "EAS") EAS_BFILE else EUR_BFILE

  sub <- panel_candidates_dt[
    target_id == this_target_id &
      source_key == this_source_key &
      clump_panel == this_panel &
      clump_reference_eligibility == "eligible"
  ]

  common <- list(
    target_id = this_target_id,
    source_key = this_source_key,
    source = as.character(job_row$source),
    source_output_stem = as.character(job_row$source_output_stem),
    clump_panel = this_panel
  )

  if (nrow(sub) == 0L) {
    return(c(common, list(
      status = "no_1000G_clumpable_candidates",
      selected = character(),
      membership = empty_clump_membership(),
      assoc_file = NA_character_,
      clumped_file = NA_character_,
      out_prefix = NA_character_
    )))
  }

  if (nrow(sub) == 1L) {
    singleton <- data.table(SNP = as.character(sub$clump_snp[[1L]]))
    return(c(common, list(
      status = "singleton_retained",
      selected = as.character(sub$clump_snp[[1L]]),
      membership = make_clump_membership(common, singleton, sub),
      assoc_file = NA_character_,
      clumped_file = NA_character_,
      out_prefix = NA_character_
    )))
  }

  label <- paste(
    safe_target_label(this_target_id),
    safe_target_label(this_source_key),
    this_panel,
    sep = "__"
  )
  assoc_file <- file.path(CLUMP_WORK_DIR, paste0(label, ".assoc.tsv"))
  out_prefix <- file.path(CLUMP_WORK_DIR, label)
  clumped_file <- paste0(out_prefix, ".clumped")

  assoc <- unique(sub[, .(
    SNP = as.character(clump_snp),
    P = as.numeric(clump_p)
  )])
  fwrite(assoc, assoc_file, sep = "\t", quote = FALSE, na = "NA")

  args <- c(
    "--bfile", ref_bfile,
    "--clump", assoc_file,
    "--clump-snp-field", "SNP",
    "--clump-field", "P",
    "--clump-p1", format(CLUMP_P1, scientific = TRUE, trim = TRUE),
    "--clump-p2", format(CLUMP_MEMBER_P2, scientific = TRUE, trim = TRUE),
    "--clump-r2", format(CLUMP_R2, scientific = FALSE, trim = TRUE),
    "--clump-kb", as.character(CLUMP_KB),
    "--out", out_prefix
  )
  exit_status <- system2(PLINK_BIN, args = args, stdout = FALSE, stderr = FALSE)

  if (is.null(exit_status) || as.integer(exit_status) != 0L) {
    return(c(common, list(
      status = "plink_clump_failed",
      selected = character(),
      membership = empty_clump_membership(),
      assoc_file = assoc_file,
      clumped_file = clumped_file,
      out_prefix = out_prefix
    )))
  }
  if (!file.exists(clumped_file) || file.info(clumped_file)$size == 0L) {
    return(c(common, list(
      status = "plink_clump_no_output",
      selected = character(),
      membership = empty_clump_membership(),
      assoc_file = assoc_file,
      clumped_file = clumped_file,
      out_prefix = out_prefix
    )))
  }

  clumped <- tryCatch(
    fread(clumped_file, fill = TRUE, showProgress = FALSE),
    error = function(e) NULL
  )
  if (is.null(clumped) || nrow(clumped) == 0L || !"SNP" %in% names(clumped)) {
    return(c(common, list(
      status = "plink_clump_unparseable_output",
      selected = character(),
      membership = empty_clump_membership(),
      assoc_file = assoc_file,
      clumped_file = clumped_file,
      out_prefix = out_prefix
    )))
  }

  selected <- unique(as.character(clumped$SNP))
  selected <- selected[!is.na(selected) & nzchar(selected)]
  if (length(selected) == 0L) {
    return(c(common, list(
      status = "plink_clump_zero_representatives",
      selected = character(),
      membership = empty_clump_membership(),
      assoc_file = assoc_file,
      clumped_file = clumped_file,
      out_prefix = out_prefix
    )))
  }

  c(common, list(
    status = "completed",
    selected = selected,
    membership = make_clump_membership(common, clumped, sub, out_prefix),
    assoc_file = assoc_file,
    clumped_file = clumped_file,
    out_prefix = out_prefix
  ))
}

empty_clump_representative_table <- function() {
  data.table(
    target_id = character(), source_key = character(), source = character(),
    source_output_stem = character(), clump_panel = character(),
    clump_id = character(), clump_index_snp = character(),
    selected_clump_member_snp = character(), selected_condition_snp = character(),
    representative_type = character(), representation_evidence = character(),
    representative_status = character(), representative_usable_in_KoGES = logical(),
    index_KoGES_status = character(), representative_KoGES_status = character(),
    proxy_r2_1000G = numeric(), proxy_ld_reference_panel = character(),
    proxy_min_r2_required = numeric(), n_clump_members = integer(),
    n_eligible_member_proxies = integer(), n_high_ld_member_proxies = integer(),
    n_panel_wide_proxy_candidates = integer(),
    clump_member_snps = character(), clump_SP2_raw = character(),
    proxy_ld_status = character(), proxy_ld_candidate_file = character(),
    proxy_ld_file = character(),
    panel_wide_proxy_ld_status = character(),
    panel_wide_proxy_ld_file = character()
  )
}

empty_member_proxy_audit <- function() {
  data.table(
    target_id = character(), source_key = character(), source = character(),
    source_output_stem = character(), clump_panel = character(),
    clump_id = character(), clump_index_snp = character(),
    proxy_candidate_snp = character(), candidate_condition_snp = character(),
    candidate_is_index = logical(),
    candidate_reported_in_SP2 = logical(), candidate_present_in_source_input = logical(),
    candidate_clump_p = numeric(), candidate_KoGES_status = character(),
    candidate_identity_status = character(),
    proxy_eligibility_status = character(), proxy_r2_1000G = numeric(),
    proxy_meets_min_r2 = logical(), proxy_selected = logical(),
    proxy_selection_status = character(), proxy_ld_reference_panel = character(),
    proxy_min_r2_required = numeric(), proxy_ld_status = character(),
    proxy_ld_candidate_file = character(), proxy_ld_file = character()
  )
}

# Select an exact KoGES variant when the clump index itself has validated
# 1000G-to-KoGES allele identity.  Otherwise permit a substitute only from the
# same PLINK clump, only if that member is an exact KoGES variant, and only
# after re-estimating index--member LD in that source's 1000G panel at
# PROXY_MIN_R2.  Same-coordinate SNP-indel or shared-allele encodings are not
# exact identity and cannot be used as the KoGES condition SNP.
select_one_clump_representative <- function(clump_members) {
  if (nrow(clump_members) == 0L) {
    return(list(
      representative = empty_clump_representative_table(),
      proxy_audit = empty_member_proxy_audit(),
      panel_wide_audit = empty_panel_wide_proxy_audit()
    ))
  }

  clump_members <- unique(copy(clump_members))
  if (!"member_identity_status" %in% names(clump_members)) {
    clump_members[, member_identity_status := NA_character_]
  }
  index_rows <- clump_members[member_is_index %in% TRUE]
  if (nrow(index_rows) != 1L) {
    stop("Each clump must contain exactly one index member. Clump: ",
         clump_members$clump_id[[1L]])
  }
  index_row <- index_rows[1L]
  index_snp <- as.character(index_row$clump_index_snp)
  index_KoGES_status <- as.character(index_row$member_KoGES_status)
  index_identity_status <- as.character(index_row$member_identity_status)
  index_is_exact <- is_exact_identity_status(index_identity_status) &
    identical(index_KoGES_status, "eligible_if_selected")
  ref_bfile <- if (index_row$clump_panel == "EAS") EAS_BFILE else EUR_BFILE
  # The target window is needed by both proxy stages: to restrict the member
  # PLINK load and to define the panel-wide rescue region.
  target_key <- as.character(index_row$target_id)
  target_row <- target_lead_lookup[.(target_key), nomatch = 0L]
  target_window_ok <- nrow(target_row) == 1L &&
    is.finite(target_row$koges_bp[[1L]]) &&
    !is.na(target_row$target_chr_normalized[[1L]])
  if (target_window_ok) {
    target_bp <- as.numeric(target_row$koges_bp[[1L]])
    target_chr <- as.character(target_row$target_chr_normalized[[1L]])
    lead_variant <- as.character(target_row$koges_lead_variant[[1L]])
    region_from_bp <- max(1, target_bp - NOVELTY_WINDOW_BP)
    region_to_bp <- target_bp + NOVELTY_WINDOW_BP
  } else {
    target_bp <- NA_real_
    target_chr <- NA_character_
    lead_variant <- NA_character_
    region_from_bp <- NA_real_
    region_to_bp <- NA_real_
  }
  member_ids <- unique(as.character(clump_members$clump_member_snp))
  member_string <- collapse_unique(member_ids)
  sp2_raw <- collapse_unique(clump_members$clump_SP2_raw)

  proxy_audit <- clump_members[, .(
    target_id, source_key, source, source_output_stem, clump_panel,
    clump_id, clump_index_snp,
    proxy_candidate_snp = clump_member_snp,
    candidate_condition_snp = member_condition_snp,
    candidate_is_index = member_is_index,
    candidate_reported_in_SP2 = member_reported_in_SP2,
    candidate_present_in_source_input = member_present_in_source_input,
    candidate_clump_p = member_clump_p,
    candidate_KoGES_status = member_KoGES_status,
    candidate_identity_status = member_identity_status
  )]
  proxy_audit[, proxy_eligibility_status := fcase(
    candidate_is_index, "index_not_proxy_candidate",
    !candidate_present_in_source_input, "member_not_in_source_input",
    !is_valid_snp_id(candidate_condition_snp),
    "member_has_no_resolved_KoGES_condition_id",
    candidate_KoGES_status != "eligible_if_selected",
    "member_not_usable_in_KoGES",
    !is_exact_identity_status(candidate_identity_status),
    "member_not_exact_1000G_to_KoGES_identity",
    default = "eligible_member_for_proxy_LD"
  )]
  proxy_audit[, `:=`(
    proxy_r2_1000G = NA_real_,
    proxy_meets_min_r2 = FALSE,
    proxy_selected = FALSE,
    proxy_selection_status = NA_character_,
    proxy_ld_reference_panel = as.character(index_row$clump_panel),
    proxy_min_r2_required = PROXY_MIN_R2,
    proxy_ld_status = NA_character_,
    proxy_ld_candidate_file = NA_character_,
    proxy_ld_file = NA_character_
  )]

  make_rep <- function(selected_clump_member_snp = NA_character_,
                       selected_condition_snp = NA_character_,
                       representative_type,
                       representation_evidence,
                       representative_status,
                       representative_usable_in_KoGES,
                       representative_KoGES_status = NA_character_,
                       proxy_r2_1000G = NA_real_,
                       n_eligible_member_proxies = 0L,
                       n_high_ld_member_proxies = 0L,
                       n_panel_wide_proxy_candidates = 0L,
                       proxy_ld_status = NA_character_,
                       proxy_ld_candidate_file = NA_character_,
                       proxy_ld_file = NA_character_,
                       panel_wide_proxy_ld_status = NA_character_,
                       panel_wide_proxy_ld_file = NA_character_) {
    data.table(
      target_id = as.character(index_row$target_id),
      source_key = as.character(index_row$source_key),
      source = as.character(index_row$source),
      source_output_stem = as.character(index_row$source_output_stem),
      clump_panel = as.character(index_row$clump_panel),
      clump_id = as.character(index_row$clump_id),
      clump_index_snp = index_snp,
      selected_clump_member_snp = selected_clump_member_snp,
      selected_condition_snp = selected_condition_snp,
      representative_type = representative_type,
      representation_evidence = representation_evidence,
      representative_status = representative_status,
      representative_usable_in_KoGES = representative_usable_in_KoGES,
      index_KoGES_status = index_KoGES_status,
      representative_KoGES_status = representative_KoGES_status,
      proxy_r2_1000G = proxy_r2_1000G,
      proxy_ld_reference_panel = as.character(index_row$clump_panel),
      proxy_min_r2_required = PROXY_MIN_R2,
      n_clump_members = as.integer(length(member_ids)),
      n_eligible_member_proxies = as.integer(n_eligible_member_proxies),
      n_high_ld_member_proxies = as.integer(n_high_ld_member_proxies),
      n_panel_wide_proxy_candidates = as.integer(n_panel_wide_proxy_candidates),
      clump_member_snps = member_string,
      clump_SP2_raw = sp2_raw,
      proxy_ld_status = proxy_ld_status,
      proxy_ld_candidate_file = proxy_ld_candidate_file,
      proxy_ld_file = proxy_ld_file,
      panel_wide_proxy_ld_status = panel_wide_proxy_ld_status,
      panel_wide_proxy_ld_file = panel_wide_proxy_ld_file
    )
  }

  work_prefix <- as.character(index_row$clump_out_prefix)
  if (is.na(work_prefix) || !nzchar(work_prefix)) {
    work_prefix <- file.path(
      CLUMP_WORK_DIR, safe_target_label(as.character(index_row$clump_id))
    )
  }

  make_panel_wide_audit <- function(candidates, selected_snp, status, ld_file) {
    if (is.null(candidates) || nrow(candidates) == 0L) {
      return(empty_panel_wide_proxy_audit())
    }
    data.table(
      target_id = as.character(index_row$target_id),
      source_key = as.character(index_row$source_key),
      source = as.character(index_row$source),
      source_output_stem = as.character(index_row$source_output_stem),
      clump_panel = as.character(index_row$clump_panel),
      clump_id = as.character(index_row$clump_id),
      clump_index_snp = index_snp,
      proxy_candidate_snp = as.character(candidates$proxy_candidate_snp),
      proxy_condition_snp = as.character(candidates$proxy_condition_snp),
      proxy_r2_1000G = as.numeric(candidates$proxy_r2_1000G),
      proxy_candidate_bp = as.integer(candidates$ref_bp),
      proxy_abs_distance_to_lead_bp = as.numeric(candidates$abs_distance_to_lead_bp),
      proxy_selected = as.character(candidates$proxy_candidate_snp) %chin% selected_snp,
      proxy_min_r2_required = PROXY_MIN_R2,
      proxy_ld_reference_panel = as.character(index_row$clump_panel),
      proxy_ld_status = status,
      proxy_ld_file = ld_file
    )
  }

  # The in-clump pool is limited to variants the external GWAS itself reported.
  # A panel-wide search asks the separate question of whether the ancestry panel
  # has any exact 1000G-to-KoGES variant in LD with this clump index.
  panel_wide_rescue <- function(unrepresented_status,
                                n_eligible_member_proxies = 0L,
                                n_high_ld_member_proxies = 0L,
                                member_ld_status = NA_character_,
                                member_ld_candidate_file = NA_character_,
                                member_ld_file = NA_character_) {
    give_up <- function(panel_wide_status, ld_file = NA_character_,
                        candidates = NULL, n_candidates = 0L) {
      list(
        representative = make_rep(
          representative_type = "none",
          representation_evidence = "unrepresented",
          representative_status = unrepresented_status,
          representative_usable_in_KoGES = FALSE,
          representative_KoGES_status = index_KoGES_status,
          n_eligible_member_proxies = n_eligible_member_proxies,
          n_high_ld_member_proxies = n_high_ld_member_proxies,
          n_panel_wide_proxy_candidates = n_candidates,
          proxy_ld_status = member_ld_status,
          proxy_ld_candidate_file = member_ld_candidate_file,
          proxy_ld_file = member_ld_file,
          panel_wide_proxy_ld_status = panel_wide_status,
          panel_wide_proxy_ld_file = ld_file
        ),
        proxy_audit = proxy_audit,
        panel_wide_audit = make_panel_wide_audit(
          candidates, character(), panel_wide_status, ld_file
        )
      )
    }

    if (!isTRUE(PANEL_WIDE_PROXY_RESCUE)) {
      return(give_up("not_attempted_panel_wide_rescue_disabled"))
    }
    exact_map <- panel_koges_exact_maps[[as.character(index_row$clump_panel)]]
    if (is.null(exact_map) || nrow(exact_map) == 0L) {
      return(give_up("no_exact_identity_pool_for_panel"))
    }
    if (!isTRUE(target_window_ok)) {
      return(give_up("target_window_unavailable"))
    }

    ld <- lookup_panel_wide_proxy_ld(
      index_snp = index_snp,
      clump_panel = as.character(index_row$clump_panel),
      target_id = target_key,
      ref_bfile = ref_bfile, work_prefix = work_prefix,
      region_chr = target_chr,
      region_from_bp = region_from_bp,
      region_to_bp = region_to_bp
    )
    if (!identical(ld$status, "completed")) {
      return(give_up(ld$status, ld$ld_file))
    }
    if (nrow(ld$table) == 0L) {
      return(give_up("no_panel_variant_at_proxy_R2_threshold", ld$ld_file))
    }

    in_ld <- merge(
      ld$table, exact_map,
      by.x = "proxy_candidate_snp", by.y = "ref_snp", sort = FALSE
    )
    if (nrow(in_ld) == 0L) {
      return(give_up("no_exact_identity_variant_in_LD", ld$ld_file))
    }
    in_ld[, abs_distance_to_lead_bp := abs(as.numeric(ref_bp) - target_bp)]
    in_ld[, proxy_condition_snp := fifelse(
      is.na(ref_koges_snp_alias), proxy_candidate_snp, ref_koges_snp_alias
    )]
    candidates <- in_ld[
      ref_chr_normalized == target_chr &
        abs_distance_to_lead_bp <= NOVELTY_WINDOW_BP &
        proxy_condition_snp != lead_variant &
        proxy_candidate_snp != index_snp &
        is.finite(proxy_r2_1000G) & proxy_r2_1000G >= PROXY_MIN_R2
    ]
    if (nrow(candidates) == 0L) {
      # Keep the rejected near-misses so an unrepresented clump is auditable.
      return(give_up(
        "no_exact_identity_variant_in_target_window", ld$ld_file,
        candidates = in_ld, n_candidates = 0L
      ))
    }

    setorder(candidates, -proxy_r2_1000G, abs_distance_to_lead_bp, proxy_candidate_snp)
    chosen <- candidates[1L]
    chosen_snp <- as.character(chosen$proxy_candidate_snp)
    list(
      representative = make_rep(
        selected_clump_member_snp = chosen_snp,
        selected_condition_snp = as.character(chosen$proxy_condition_snp),
        representative_type = "panel_wide_proxy_high_LD_1000G",
        representation_evidence = "ld_proxy_r2_ge_threshold",
        representative_status = "panel_wide_ld_proxy_selected_r2_ge_threshold",
        representative_usable_in_KoGES = TRUE,
        representative_KoGES_status = "eligible_if_selected",
        proxy_r2_1000G = as.numeric(chosen$proxy_r2_1000G),
        n_eligible_member_proxies = n_eligible_member_proxies,
        n_high_ld_member_proxies = n_high_ld_member_proxies,
        n_panel_wide_proxy_candidates = nrow(candidates),
        proxy_ld_status = member_ld_status,
        proxy_ld_candidate_file = member_ld_candidate_file,
        proxy_ld_file = member_ld_file,
        panel_wide_proxy_ld_status = ld$status,
        panel_wide_proxy_ld_file = ld$ld_file
      ),
      proxy_audit = proxy_audit,
      panel_wide_audit = make_panel_wide_audit(
        candidates, chosen_snp, ld$status, ld$ld_file
      )
    )
  }

  if (isTRUE(index_is_exact)) {
    proxy_audit[, proxy_selection_status := "exact_variant_available_in_KoGES"]
    proxy_audit[candidate_is_index %in% TRUE, proxy_r2_1000G := 1]
    return(list(
      representative = make_rep(
        selected_clump_member_snp = index_snp,
        selected_condition_snp = as.character(index_row$member_condition_snp),
        representative_type = "direct_clump_index",
        representation_evidence = "exact_variant",
        representative_status = "exact_variant_available_in_KoGES",
        representative_usable_in_KoGES = TRUE,
        representative_KoGES_status = index_KoGES_status,
        proxy_r2_1000G = 1,
        proxy_ld_status = "not_needed_exact_variant",
        panel_wide_proxy_ld_status = "not_needed_exact_variant"
      ),
      proxy_audit = proxy_audit,
      panel_wide_audit = empty_panel_wide_proxy_audit()
    ))
  }

  eligible_proxy_rows <- proxy_audit[
    proxy_eligibility_status == "eligible_member_for_proxy_LD"
  ]
  n_eligible <- nrow(eligible_proxy_rows)
  if (n_eligible == 0L) {
    proxy_audit[, proxy_selection_status := "no_eligible_member_proxy"]
    return(panel_wide_rescue(
      unrepresented_status = "unrepresented_no_exact_variant_or_ld_proxy",
      n_eligible_member_proxies = 0L,
      member_ld_status = "no_eligible_member_candidates"
    ))
  }

  ld_result <- calculate_member_proxy_ld(
    index_snp = index_snp,
    proxy_candidates = eligible_proxy_rows$proxy_candidate_snp,
    ref_bfile = ref_bfile,
    out_prefix = as.character(index_row$clump_out_prefix),
    region_chr = target_chr,
    region_from_bp = region_from_bp,
    region_to_bp = region_to_bp
  )
  if (nrow(ld_result$table) > 0L) {
    r2_match <- match(
      proxy_audit$proxy_candidate_snp,
      ld_result$table$proxy_candidate_snp
    )
    proxy_audit[, proxy_r2_1000G := ld_result$table$proxy_r2_1000G[r2_match]]
  }
  proxy_audit[, `:=`(
    proxy_ld_status = ld_result$status,
    proxy_ld_candidate_file = ld_result$candidate_file,
    proxy_ld_file = ld_result$ld_file,
    proxy_meets_min_r2 = proxy_eligibility_status == "eligible_member_for_proxy_LD" &
      is.finite(proxy_r2_1000G) & proxy_r2_1000G >= PROXY_MIN_R2
  )]

  if (!identical(ld_result$status, "completed")) {
    proxy_audit[, proxy_selection_status := "member_proxy_LD_calculation_failed"]
    return(panel_wide_rescue(
      unrepresented_status = "unrepresented_member_proxy_LD_failed",
      n_eligible_member_proxies = n_eligible,
      member_ld_status = ld_result$status,
      member_ld_candidate_file = ld_result$candidate_file,
      member_ld_file = ld_result$ld_file
    ))
  }

  high_ld <- proxy_audit[proxy_meets_min_r2 %in% TRUE]
  n_high_ld <- nrow(high_ld)
  if (n_high_ld == 0L) {
    proxy_audit[, proxy_selection_status := "no_member_meets_proxy_r2_threshold"]
    return(panel_wide_rescue(
      unrepresented_status = "unrepresented_no_member_at_proxy_R2_threshold",
      n_eligible_member_proxies = n_eligible,
      n_high_ld_member_proxies = 0L,
      member_ld_status = ld_result$status,
      member_ld_candidate_file = ld_result$candidate_file,
      member_ld_file = ld_result$ld_file
    ))
  }

  setorder(high_ld, -proxy_r2_1000G, candidate_clump_p, proxy_candidate_snp)
  selected_proxy <- high_ld[1L]
  proxy_audit[, proxy_selection_status := fifelse(
    proxy_candidate_snp == selected_proxy$proxy_candidate_snp,
    "selected_high_LD_member_proxy",
    "not_selected_after_proxy_r2_ranking"
  )]
  proxy_audit[
    proxy_candidate_snp == selected_proxy$proxy_candidate_snp,
    proxy_selected := TRUE
  ]
  list(
    representative = make_rep(
      selected_clump_member_snp = as.character(selected_proxy$proxy_candidate_snp),
      selected_condition_snp = as.character(selected_proxy$candidate_condition_snp),
      representative_type = "member_proxy_high_LD_1000G",
      representation_evidence = "ld_proxy_r2_ge_threshold",
      representative_status = "ld_proxy_selected_r2_ge_threshold",
      representative_usable_in_KoGES = TRUE,
      representative_KoGES_status = as.character(selected_proxy$candidate_KoGES_status),
      proxy_r2_1000G = as.numeric(selected_proxy$proxy_r2_1000G),
      n_eligible_member_proxies = n_eligible,
      n_high_ld_member_proxies = n_high_ld,
      proxy_ld_status = ld_result$status,
      proxy_ld_candidate_file = ld_result$candidate_file,
      proxy_ld_file = ld_result$ld_file,
      panel_wide_proxy_ld_status = "not_needed_member_proxy_available"
    ),
    proxy_audit = proxy_audit,
    panel_wide_audit = empty_panel_wide_proxy_audit()
  )
}

merge_into_manifest <- function(base, add) {
  if (nrow(add) == 0L) return(base)
  by_cols <- intersect(names(add), target_columns)
  if (length(by_cols) == 0L) by_cols <- "target_id"
  merge(base, add, by = by_cols, all.x = TRUE, sort = FALSE)
}

# ------------------------------------------------------------------------------
# 2. Read Step1 source-level prior evidence
# ------------------------------------------------------------------------------

raw_input_qc_file <- file.path(OUT_DIR, "joint_COJO_raw_prior_input_QC.tsv")
panel_reference_qc_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_level_candidates_1000G_reference_QC.tsv"
)
source_candidates_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_level_condition_candidates_preclump.tsv"
)
clump_summary_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_ancestry_clump_job_summary.tsv"
)
# The index/member map is the largest audit table in the run (tens of GB
# uncompressed).  It is written gzipped; fread reads either form, so a resume
# from an earlier uncompressed run still works.
clump_member_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_clump_members.tsv.gz"
)
clump_member_audit_file_legacy <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_clump_members.tsv"
)
resolve_clump_member_audit_file <- function() {
  if (file.exists(clump_member_audit_file)) return(clump_member_audit_file)
  if (file.exists(clump_member_audit_file_legacy)) {
    return(clump_member_audit_file_legacy)
  }
  clump_member_audit_file
}
clump_representative_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_clump_representatives.tsv"
)
member_proxy_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_member_proxy_LD_audit.tsv"
)
panel_wide_proxy_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_panel_wide_proxy_LD_audit.tsv"
)
representative_collision_file <- file.path(
  OUT_DIR,
  "joint_COJO_within_source_representative_collisions.tsv"
)
final_candidates_pre_cross_source_ld_file <- file.path(
  OUT_DIR,
  "joint_COJO_condition_candidates_pre_cross_source_KoGES_LD.tsv"
)
cross_source_ld_member_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_group_members.tsv"
)
cross_source_ld_pair_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_pairs.tsv"
)
cross_source_ld_summary_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_grouping_summary.tsv"
)
final_candidates_file <- file.path(
  OUT_DIR,
  "joint_COJO_condition_candidates_ancestry_clumped.tsv"
)
post_clump_selection_files <- c(
  clump_representative_file,
  representative_collision_file,
  cross_source_ld_summary_file,
  final_candidates_file
)

skip_post_clump_selection <- FALSE
if (isTRUE(RESUME_AFTER_CLUMP)) {
  resume_required <- c(
    PRELD_MANIFEST_FILE, raw_input_qc_file, panel_reference_qc_file,
    clump_summary_file, resolve_clump_member_audit_file()
  )
  missing_resume <- resume_required[!file.exists(resume_required)]
  if (length(missing_resume) > 0L) {
    stop(
      "RESUME_AFTER_CLUMP=TRUE but missing required file(s):\n",
      paste(missing_resume, collapse = "\n")
    )
  }
  message(
    "RESUME_AFTER_CLUMP=TRUE: skipping prep/PLINK clumping; ",
    "loading panel candidates + clump membership from OUT_DIR."
  )
  manifest <- fread(PRELD_MANIFEST_FILE, sep = "\t", showProgress = FALSE)
  required_manifest_columns <- c(target_columns, "cojo_condition_status")
  missing_manifest <- setdiff(required_manifest_columns, names(manifest))
  if (length(missing_manifest) > 0L) {
    stop("Step1 manifest missing: ", paste(missing_manifest, collapse = ", "))
  }
  if (anyDuplicated(manifest$target_id)) stop("Step1 manifest has duplicate target_id.")
  manifest[, target_id := as.character(target_id)]

  # A resumed pilot must describe the same target set as the run that produced
  # the clump membership, otherwise the manifest reports thousands of targets
  # that were never clumped.
  target_subset <- resolve_target_subset(TARGET_IDS_ENV, TARGET_ID_FILE)
  if (length(target_subset) > 0L) {
    unknown_targets <- setdiff(target_subset, manifest$target_id)
    if (length(unknown_targets) > 0L) {
      stop(
        "TARGET_IDS / TARGET_ID_FILE contains target_id absent from the Step1 ",
        "manifest: ", paste(head(unknown_targets, 20L), collapse = ", ")
      )
    }
    manifest <- manifest[target_id %chin% target_subset]
    message(
      "TARGET subset active on resume: ", nrow(manifest), " target(s). ",
      "Every output in this OUT_DIR describes only this subset."
    )
  }

  koges_bim <- read_bim_index(KOGES_BFILE)
  lead_match <- match(manifest$koges_lead_variant, koges_bim$ref_snp)
  manifest[, lead_reference_chr := koges_bim$ref_chr[lead_match]]
  manifest[, lead_reference_bp := koges_bim$ref_bp[lead_match]]
  manifest[, target_chr_normalized := normalize_chr(koges_chr)]
  manifest[, lead_reference_chr_normalized := normalize_chr(lead_reference_chr)]
  manifest[, lead_reference_status := fcase(
    is.na(lead_reference_bp), "missing_target_lead_SNP_in_KoGES_bim",
    lead_reference_chr_normalized != target_chr_normalized,
    "target_lead_chr_mismatch_in_KoGES_bim",
    lead_reference_bp != koges_bp, "target_lead_bp_mismatch_in_KoGES_bim",
    default = "eligible"
  )]

  message("Loading panel candidates: ", panel_reference_qc_file)
  panel_candidates <- fread(panel_reference_qc_file, sep = "\t", showProgress = FALSE)
  if (
    "clump_reference_eligibility" %in% names(panel_candidates) &&
      any(panel_candidates$clump_reference_eligibility ==
        "1000G_to_KoGES_variant_identity_mismatch")
  ) {
    message(
      "RESUME_AFTER_CLUMP loaded a 1000G QC table that excluded identity-",
      "discordant SNPs from clumping.  Re-run Step2 without RESUME_AFTER_CLUMP ",
      "so every 1000G-representable prior is clumped, then represented by ",
      "exact_variant / LD proxy / unrepresented."
    )
  } else if (
    "panel_to_KoGES_identity_status" %in% names(panel_candidates) &&
      any(panel_candidates$panel_to_KoGES_identity_status %chin%
        c(
          "validated_1000G_to_KoGES_chrbp_snp_indel_alias",
          "validated_1000G_to_KoGES_chrbp_shared_allele_multiallelic"
        ))
  ) {
    message(
      "RESUME_AFTER_CLUMP loaded QC that treated same-coordinate SNP-indel / ",
      "shared-allele encodings as validated identity.  Re-run representative ",
      "selection (do not skip post-clump outputs) so those clumps require ",
      "exact KoGES identity or an r^2 proxy."
    )
  }
  message("Rebuilding raw-input summary from: ", raw_input_qc_file)
  raw_primary <- fread(raw_input_qc_file, sep = "\t", showProgress = FALSE)
  raw_primary <- harmonize_raw_prior_status_columns(raw_primary)
  if (DIRECT_KOGES_RESCUE && !all(c(
    "direct_KoGES_exact_snp", "direct_KoGES_exact_status",
    "direct_KoGES_rescue_eligible"
  ) %in% names(raw_primary))) {
    stop(
      "RESUME_AFTER_CLUMP cannot enable DIRECT_KOGES_RESCUE from an older raw ",
      "QC file. Re-run Step2 without RESUME_AFTER_CLUMP once so raw external ",
      "CHR:BP+alleles can be resolved against the current KoGES BIM."
    )
  }
  direct_koges_rescue_bundle <- build_direct_koges_rescue_tables(raw_primary)
  direct_koges_rescue_union <- direct_koges_rescue_bundle$union
  direct_koges_rescue_selection <- direct_koges_rescue_bundle$selection
  rm(direct_koges_rescue_bundle)
  raw_input_summary <- summarize_raw_prior_input(raw_primary)
  rm(raw_primary)
  invisible(gc())

  if (file.exists(source_candidates_file)) {
    message("Loading source-level candidates: ", source_candidates_file)
    source_candidates <- fread(source_candidates_file, sep = "\t", showProgress = FALSE)
  } else {
    message(
      "Source-level preclump candidate file missing; ",
      "using panel candidates for source-level counts."
    )
    source_candidates <- copy(panel_candidates)
  }

  clump_result_summary <- fread(clump_summary_file, sep = "\t", showProgress = FALSE)
  bad_jobs <- clump_result_summary[
    !(clump_status %chin% c(
      "completed", "singleton_retained", "no_1000G_clumpable_candidates"
    ))
  ]
  if (nrow(bad_jobs) > 0L) {
    stop(
      "RESUME_AFTER_CLUMP cannot continue: clump job summary has failures. See: ",
      clump_summary_file
    )
  }

  # V4 writes distinct EAS-primary and EUR-supporting representative sets.
  # Rebuild post-clump selection even when legacy all-source audit files exist;
  # otherwise a resume could silently reuse a non-stratified condition set.
  skip_post_clump_selection <- FALSE
  if (all(file.exists(post_clump_selection_files))) {
    message(
      "RESUME_AFTER_CLUMP: rebuilding post-clump selection to create the ",
      "EAS-primary and EUR-supporting condition sets."
    )
  }
  if (isTRUE(skip_post_clump_selection)) {
    message(
      "RESUME_AFTER_CLUMP: representative/cross-source outputs already exist; ",
      "skipping clump-membership reload, representative selection, and KoGES-LD grouping."
    )
    clump_membership <- empty_clump_membership()
  } else {
    clump_member_audit_resume_file <- resolve_clump_member_audit_file()
    message("Loading clump membership (large file): ",
            clump_member_audit_resume_file)
    clump_membership <- fread(
      clump_member_audit_resume_file, sep = "\t", showProgress = FALSE
    )
    for (col in c(
      "member_is_index", "member_reported_in_SP2", "member_present_in_source_input"
    )) {
      if (col %in% names(clump_membership)) {
        clump_membership[, (col) := as.logical(get(col))]
      }
    }
    message(
      "Resume loaded ", nrow(panel_candidates), " panel candidates and ",
      nrow(clump_membership), " clump-member rows (",
      uniqueN(clump_membership$clump_id), " clumps)."
    )
  }
} else {
  skip_post_clump_selection <- FALSE

raw <- fread(RAW_HITS_FILE, sep = "\t", showProgress = FALSE)
manifest <- fread(PRELD_MANIFEST_FILE, sep = "\t", showProgress = FALSE)

required_manifest_columns <- c(target_columns, "cojo_condition_status")
required_raw_columns <- c(
  target_columns,
  "source", "evidence_layer", "source_ancestry", "source_output_stem",
  "condition_snp", "condition_snp_mappable", "is_known_exact_variant",
  "external_snp", "external_p", "distance_from_koges_lead_bp",
  # Needed to resolve a clumping ID at allele level.  Without the raw external
  # position and alleles, a KoGES rsID alias cannot be verified as the same
  # variant, and using it as the clump index would silently redefine the signal.
  "external_chr", "external_bp",
  "external_effect_allele", "external_other_allele"
)
missing_manifest <- setdiff(required_manifest_columns, names(manifest))
missing_raw <- setdiff(required_raw_columns, names(raw))
if (length(missing_manifest) > 0L) {
  stop("Step1 manifest missing: ", paste(missing_manifest, collapse = ", "))
}
if (length(missing_raw) > 0L) {
  stop("Step1 raw-hit file missing: ", paste(missing_raw, collapse = ", "))
}
if (anyDuplicated(manifest$target_id)) stop("Step1 manifest has duplicate target_id.")

manifest[, target_id := as.character(target_id)]
raw[, target_id := as.character(target_id)]
raw[, source := as.character(source)]
raw[, evidence_layer := as.character(evidence_layer)]
raw[, source_ancestry := as.character(source_ancestry)]
raw[, source_output_stem := as.character(source_output_stem)]
raw[, condition_snp := trimws(as.character(condition_snp))]
raw[, external_snp := trimws(as.character(external_snp))]
raw[, external_chr_normalized := normalize_chr(external_chr)]
raw[, external_bp := suppressWarnings(as.integer(external_bp))]
raw[, external_effect_allele := normalize_allele(external_effect_allele)]
raw[, external_other_allele := normalize_allele(external_other_allele)]
raw[, external_p := suppressWarnings(as.numeric(external_p))]
raw[, distance_from_koges_lead_bp := suppressWarnings(as.numeric(distance_from_koges_lead_bp))]
raw[is.finite(external_p) & external_p == 0, external_p := .Machine$double.xmin]

if (any(!raw$target_id %chin% manifest$target_id)) {
  stop("Step1 raw-hit file contains target_id absent from its manifest.")
}

target_subset <- resolve_target_subset(TARGET_IDS_ENV, TARGET_ID_FILE)
if (length(target_subset) > 0L) {
  unknown_targets <- setdiff(target_subset, manifest$target_id)
  if (length(unknown_targets) > 0L) {
    stop(
      "TARGET_IDS / TARGET_ID_FILE contains target_id absent from the Step1 ",
      "manifest: ", paste(head(unknown_targets, 20L), collapse = ", ")
    )
  }
  manifest <- manifest[target_id %chin% target_subset]
  raw <- raw[target_id %chin% target_subset]
  message(
    "TARGET subset active: ", nrow(manifest), " target(s), ", nrow(raw),
    " raw prior-hit row(s). Every output in this OUT_DIR describes only ",
    "this subset; do not merge it with a full run."
  )
}

STATUS_PRIMARY <- "preLD_candidates_available"
primary_ids <- manifest[cojo_condition_status == STATUS_PRIMARY, target_id]

# Catalog remains annotation-only.  Crucially, do NOT require condition_snp to be
# mappable in KoGES here.  Source-specific external signals must be defined in
# ancestry-matched 1000G LD before KoGES representability is evaluated.
raw_primary <- raw[
  target_id %chin% primary_ids &
    evidence_layer == ANALYSIS_LAYER &
    !as_bool(is_known_exact_variant)
]
if (nrow(raw) > 0L && any(raw$evidence_layer == "EUR")) {
  message(
    "EAS-only Step2: dropping ",
    sum(raw$evidence_layer == "EUR"),
    " EUR raw-hit row(s) from the Step1 file."
  )
}

if (nrow(raw_primary) > 0L) {
  bad_source <- raw_primary[
    is.na(source) | !nzchar(source) |
      is.na(source_output_stem) | !nzchar(source_output_stem)
  ]
  if (nrow(bad_source) > 0L) {
    stop(
      "Source-level raw hits lack source or source_output_stem. Example target: ",
      bad_source$target_id[[1L]]
    )
  }
}
raw_primary[, source_key := make_source_key(source, source_output_stem)]
raw_primary[, external_snp_is_valid := is_valid_snp_id(external_snp)]
raw_primary[, condition_snp_is_valid := is_valid_snp_id(condition_snp)]

# ------------------------------------------------------------------------------
# 3. Resolve an ancestry-panel clumping ID WITHOUT using KoGES mappability as a
#    gate, and WITHOUT ever substituting a different variant for the external
#    signal.  The clump index defines what "the external signal" means, so an
#    alias is accepted only at allele level: the original external rsID if the
#    panel has it, else the unique panel variant at the raw external CHR:BP with
#    the same unordered allele pair, else the Step1 KoGES rsID but only when that
#    rsID is verified to be the same variant in both KoGES and the panel.  A
#    same-coordinate SNP-indel or a different ALT is not an alias.
# ------------------------------------------------------------------------------

koges_bim <- read_bim_index(KOGES_BFILE)
eas_bim <- read_bim_index(EAS_BFILE)
eur_bim <- eas_bim[0]
koges_position_index <- build_panel_position_index(koges_bim)

lead_match <- match(manifest$koges_lead_variant, koges_bim$ref_snp)
manifest[, lead_reference_chr := koges_bim$ref_chr[lead_match]]
manifest[, lead_reference_bp := koges_bim$ref_bp[lead_match]]
manifest[, target_chr_normalized := normalize_chr(koges_chr)]
manifest[, lead_reference_chr_normalized := normalize_chr(lead_reference_chr)]
manifest[, lead_reference_status := fcase(
  is.na(lead_reference_bp), "missing_target_lead_SNP_in_KoGES_bim",
  lead_reference_chr_normalized != target_chr_normalized,
  "target_lead_chr_mismatch_in_KoGES_bim",
  lead_reference_bp != koges_bp, "target_lead_bp_mismatch_in_KoGES_bim",
  default = "eligible"
)]

raw_primary[, `:=`(
  external_snp_in_1000G = FALSE,
  condition_snp_in_1000G = FALSE,
  clump_snp = NA_character_,
  clump_snp_resolution_status = NA_character_,
  n_allele_level_alias_matches_1000G = NA_integer_,
  external_to_condition_snp_identity_verified = FALSE
)]

resolve_panel_ids <- function(idx, panel_bim, panel_position_index, panel_label) {
  if (length(idx) == 0L) return(invisible(NULL))

  # A panel rsID by itself is not sufficient: rsID reuse / multiallelic ALT
  # drift would otherwise silently turn the external hit into a different
  # variant.  Require the original external CHR:BP+alleles to match this exact
  # panel record before using it as the clump index.
  ext_panel_match <- match(raw_primary$external_snp[idx], panel_bim$ref_snp)
  ext_raw_identity <- raw_primary$external_snp_is_valid[idx] &
    !is.na(ext_panel_match) &
    !is.na(raw_primary$external_chr_normalized[idx]) &
    !is.na(raw_primary$external_bp[idx]) &
    !is.na(raw_primary$external_effect_allele[idx]) &
    !is.na(raw_primary$external_other_allele[idx]) &
    normalize_chr(panel_bim$ref_chr[ext_panel_match]) ==
      raw_primary$external_chr_normalized[idx] &
    panel_bim$ref_bp[ext_panel_match] == raw_primary$external_bp[idx] &
    allele_pair_matches(
      panel_bim$ref_a1[ext_panel_match], panel_bim$ref_a2[ext_panel_match],
      raw_primary$external_effect_allele[idx],
      raw_primary$external_other_allele[idx]
    ) %in% TRUE
  ext_ok <- ext_raw_identity
  cond_in_panel <- raw_primary$condition_snp_is_valid[idx] &
    raw_primary$condition_snp[idx] %chin% panel_bim$ref_snp
  raw_primary[idx, external_snp_in_1000G := ext_ok]
  raw_primary[idx, condition_snp_in_1000G := cond_in_panel]

  chosen <- rep(NA_character_, length(idx))
  status <- rep(NA_character_, length(idx))
  n_alias <- rep(NA_integer_, length(idx))
  verified_condition <- rep(FALSE, length(idx))

  chosen[ext_ok] <- raw_primary$external_snp[idx][ext_ok]
  status[ext_ok] <- paste0("external_snp_raw_identity_verified_in_1000G_", panel_label)

  query <- data.table(
    slot = seq_along(idx),
    chr_normalized = raw_primary$external_chr_normalized[idx],
    bp = raw_primary$external_bp[idx],
    external_effect_allele = raw_primary$external_effect_allele[idx],
    external_other_allele = raw_primary$external_other_allele[idx]
  )
  query[, alias_query_usable := !ext_ok & !is.na(chr_normalized) & !is.na(bp) &
    !is.na(external_effect_allele) & !is.na(external_other_allele)]
  status[!ext_ok & !query$alias_query_usable] <- paste0(
    "missing_external_position_or_alleles_in_1000G_", panel_label
  )
  status[!ext_ok & !is.na(ext_panel_match) & query$alias_query_usable] <- paste0(
    "external_rsid_present_but_raw_identity_mismatch_in_1000G_", panel_label
  )

  if (any(query$alias_query_usable)) {
    hits <- panel_position_index[
      query[alias_query_usable %in% TRUE],
      on = .(chr_normalized, bp),
      nomatch = 0L,
      allow.cartesian = TRUE
    ]
    if (nrow(hits) > 0L) {
      hits <- hits[allele_pair_matches(
        panel_a1, panel_a2, external_effect_allele, external_other_allele
      ) %in% TRUE]
    }
    if (nrow(hits) > 0L) {
      setorder(hits, slot, alias_snp)
      alias <- hits[, .(
        n_matches = uniqueN(alias_snp),
        alias_snp = alias_snp[1L]
      ), by = slot]
      n_alias[alias$slot] <- as.integer(alias$n_matches)
      unique_alias <- alias[n_matches == 1L]
      if (nrow(unique_alias) > 0L) {
        chosen[unique_alias$slot] <- unique_alias$alias_snp
        status[unique_alias$slot] <- paste0(
          "exact_allele_level_alias_in_1000G_", panel_label
        )
      }
    }
    n_alias[query$alias_query_usable %in% TRUE & is.na(n_alias)] <- 0L
  }

  # Last resort: the Step1 KoGES rsID, accepted only when that rsID is the same
  # variant as the raw external hit in BOTH KoGES and the panel.
  fallback_slots <- which(
    !ext_ok & is.na(chosen) & cond_in_panel & query$alias_query_usable %in% TRUE
  )
  if (length(fallback_slots) > 0L) {
    rows <- idx[fallback_slots]
    cond_snp <- raw_primary$condition_snp[rows]
    ea <- raw_primary$external_effect_allele[rows]
    oa <- raw_primary$external_other_allele[rows]
    ext_chr <- raw_primary$external_chr_normalized[rows]
    ext_bp <- raw_primary$external_bp[rows]

    k <- match(cond_snp, koges_bim$ref_snp)
    koges_same_variant <- !is.na(k) &
      normalize_chr(koges_bim$ref_chr[k]) == ext_chr &
      koges_bim$ref_bp[k] == ext_bp &
      allele_pair_matches(koges_bim$ref_a1[k], koges_bim$ref_a2[k], ea, oa) %in% TRUE

    p <- match(cond_snp, panel_bim$ref_snp)
    panel_same_variant <- !is.na(p) &
      normalize_chr(panel_bim$ref_chr[p]) == ext_chr &
      panel_bim$ref_bp[p] == ext_bp &
      allele_pair_matches(panel_bim$ref_a1[p], panel_bim$ref_a2[p], ea, oa) %in% TRUE

    verified <- koges_same_variant & panel_same_variant
    verified[is.na(verified)] <- FALSE
    verified_condition[fallback_slots] <- verified
    chosen[fallback_slots[verified]] <- cond_snp[verified]
    status[fallback_slots[verified]] <- paste0(
      "verified_KoGES_condition_snp_alias_in_1000G_", panel_label
    )
  }

  still_missing <- !ext_ok & is.na(chosen) & query$alias_query_usable %in% TRUE
  status[still_missing] <- fifelse(
    !is.na(n_alias[still_missing]) & n_alias[still_missing] > 1L,
    paste0("ambiguous_allele_level_alias_in_1000G_", panel_label),
    paste0("no_exact_allele_level_alias_in_1000G_", panel_label)
  )

  raw_primary[idx, clump_snp := chosen]
  raw_primary[idx, clump_snp_resolution_status := status]
  raw_primary[idx, n_allele_level_alias_matches_1000G := n_alias]
  raw_primary[idx, external_to_condition_snp_identity_verified := verified_condition]
  invisible(NULL)
}

eas_position_index <- build_panel_position_index(eas_bim)
eur_position_index <- build_panel_position_index(eur_bim)
resolve_panel_ids(
  which(raw_primary$evidence_layer == "EAS"), eas_bim, eas_position_index, "EAS"
)
resolve_panel_ids(
  which(raw_primary$evidence_layer == "EUR"), eur_bim, eur_position_index, "EUR"
)
rm(eas_position_index, eur_position_index)
invisible(gc())

message(
  "1000G clump-index resolution:\n",
  paste(
    capture.output(print(
      raw_primary[, .N, by = clump_snp_resolution_status][order(-N)]
    )),
    collapse = "\n"
  )
)

raw_primary[, external_position_alleles_available :=
  !is.na(external_chr_normalized) & !is.na(external_bp) &
    !is.na(external_effect_allele) & !is.na(external_other_allele)]
raw_primary[, raw_input_status := fcase(
  !is.finite(external_p) | external_p <= 0 | external_p > CLUMP_P2,
  "invalid_or_non_GWS_external_P",
  # An rsID-free hit is still identifiable at allele level, which is what the
  # clump index actually needs.
  !external_snp_is_valid & !condition_snp_is_valid &
    !external_position_alleles_available,
  "no_valid_variant_identifier",
  default = "valid"
)]
raw_primary[, reference_coverage_status := fcase(
  raw_input_status != "valid",
  "not_evaluated_invalid_input",
  is_valid_snp_id(clump_snp),
  "represented_in_ancestry_1000G",
  default = "not_represented_in_ancestry_1000G"
)]

# A 1000G-absent hit cannot participate in source-specific 1000G clumping, but
# it need not be discarded if its raw external alleles resolve to exactly one
# KoGES BIM variant.  This is an *additional* condition SNP, not a proxy and
# not evidence that the missing hit belongs to a particular 1000G clump.
direct_koges_resolution <- resolve_unique_raw_to_koges_exact(
  raw_primary, koges_position_index
)
raw_primary[, direct_KoGES_exact_snp := direct_koges_resolution$direct_KoGES_exact_snp]
raw_primary[, direct_KoGES_exact_status := direct_koges_resolution$direct_KoGES_exact_status]
raw_primary[, direct_KoGES_rescue_eligible := (
  DIRECT_KOGES_RESCUE &
    raw_input_status == "valid" &
    reference_coverage_status == "not_represented_in_ancestry_1000G" &
    !is.na(clump_snp_resolution_status) &
    grepl(
      "^no_exact_allele_level_alias_in_1000G_",
      clump_snp_resolution_status
    ) &
    is_valid_snp_id(direct_KoGES_exact_snp)
)]
rm(direct_koges_resolution)
invisible(gc())
# Compatibility column for older QC readers / resume from a previous V3 file.
raw_primary[, prior_input_status := fcase(
  raw_input_status != "valid", as.character(raw_input_status),
  reference_coverage_status == "represented_in_ancestry_1000G",
  "eligible",
  default = "external_prior_not_represented_in_ancestry_1000G"
)]

# For each raw prior, resolve the SNP ID that would be used later in KoGES if
# this prior becomes a selected clump member.  condition_snp is preferred because
# Step1 already harmonized source rsIDs to the KoGES SNP namespace; if unavailable,
# the 1000G clump ID itself may be used when present in KoGES.
raw_primary[, KoGES_condition_snp_row := fifelse(
  condition_snp_is_valid & condition_snp %chin% koges_bim$ref_snp,
  condition_snp,
  fifelse(is_valid_snp_id(clump_snp) & clump_snp %chin% koges_bim$ref_snp,
          clump_snp, NA_character_)
)]

direct_koges_rescue_bundle <- build_direct_koges_rescue_tables(raw_primary)
direct_koges_rescue_union <- direct_koges_rescue_bundle$union
direct_koges_rescue_selection <- direct_koges_rescue_bundle$selection
rm(direct_koges_rescue_bundle)
rm(koges_position_index)
invisible(gc())

raw_input_qc_file <- file.path(OUT_DIR, "joint_COJO_raw_prior_input_QC.tsv")
fwrite(raw_primary, raw_input_qc_file, sep = "\t", quote = FALSE, na = "NA")

raw_input_summary <- summarize_raw_prior_input(raw_primary)

# Clump only 1000G-representable valid hits.  Unavailable hits stay in the
# raw QC audit for coverage accounting and are not dropped from the target.
clumpable_raw <- raw_primary[
  raw_input_status == "valid" &
    reference_coverage_status == "represented_in_ancestry_1000G"
]

# Deduplicate only within one original external GWAS source and one 1000G clump
# SNP.  This keeps source-specific P ordering while allowing old/new rsID aliases
# to resolve to one ancestry-panel variant.
source_candidates <- clumpable_raw[, {
  mapped_ids <- unique(KoGES_condition_snp_row[
    is_valid_snp_id(KoGES_condition_snp_row)
  ])
  list(
    condition_snp = if (length(mapped_ids) == 1L) mapped_ids[[1L]] else NA_character_,
    n_KoGES_condition_id_candidates = length(mapped_ids),
    KoGES_condition_id_candidates = if (length(mapped_ids)) {
      paste(sort(mapped_ids), collapse = ";")
    } else NA_character_,
    condition_mapping_collision = length(mapped_ids) > 1L,
    clump_p = min_finite(external_p),
    source_p_min = min_finite(external_p),
    source_p_values = format_p_values(external_p),
    source_abs_distance_min_bp = min_finite(abs(distance_from_koges_lead_bp)),
    source_abs_distance_max_bp = max_finite(abs(distance_from_koges_lead_bp)),
    source_variant_ids = collapse_unique(external_snp),
    clump_snp_resolution_statuses = collapse_unique(clump_snp_resolution_status),
    n_raw_source_hits = .N
  )
}, by = c(
  target_columns,
  "source", "source_output_stem", "source_key",
  "evidence_layer", "source_ancestry", "clump_snp"
)]

# Drop the multi-GB row-level tables before forking PLINK workers. They are
# already written to the QC audit and collapsed into source_candidates.
rm(raw, raw_primary, clumpable_raw)
invisible(gc())
message(
  "Prep tables collapsed to ", nrow(source_candidates),
  " source-level candidates; freed raw prior tables before parallel clumping."
)

if (nrow(source_candidates) > 0L) {
  source_panel_check <- source_candidates[, .(
    n_panel = uniqueN(evidence_layer),
    panels = collapse_unique(evidence_layer)
  ), by = source_key]
  if (any(source_panel_check$n_panel != 1L)) {
    bad <- source_panel_check[n_panel != 1L, source_key]
    stop(
      "A source_output_stem maps to more than one ancestry panel: ",
      paste(head(bad, 10L), collapse = ", ")
    )
  }

  candidate_match <- match(source_candidates$condition_snp, koges_bim$ref_snp)
  source_candidates[, condition_KoGES_chr := koges_bim$ref_chr[candidate_match]]
  source_candidates[, condition_KoGES_bp := koges_bim$ref_bp[candidate_match]]
  source_candidates[, condition_KoGES_a1 := koges_bim$ref_a1[candidate_match]]
  source_candidates[, condition_KoGES_a2 := koges_bim$ref_a2[candidate_match]]
  source_candidates[, target_chr_normalized := normalize_chr(koges_chr)]
  source_candidates[, condition_KoGES_chr_normalized := normalize_chr(condition_KoGES_chr)]
  source_candidates[, condition_KoGES_status := fcase(
    condition_mapping_collision,
    "multiple_KoGES_condition_ids_for_same_1000G_variant",
    is_valid_snp_id(condition_snp) & condition_snp == koges_lead_variant,
    "condition_is_target_lead",
    !is.na(condition_KoGES_bp) &
      condition_KoGES_chr_normalized != target_chr_normalized,
    "KoGES_condition_chr_different_from_target",
    !is.na(condition_KoGES_bp) &
      abs(condition_KoGES_bp - koges_bp) > NOVELTY_WINDOW_BP,
    "KoGES_condition_outside_target_window",
    is.na(condition_KoGES_bp), "missing_condition_in_KoGES_bim",
    default = "eligible_if_selected"
  )]
} else {
  source_candidates[, condition_KoGES_chr := character()]
  source_candidates[, condition_KoGES_bp := integer()]
  source_candidates[, condition_KoGES_a1 := character()]
  source_candidates[, condition_KoGES_a2 := character()]
  source_candidates[, target_chr_normalized := character()]
  source_candidates[, condition_KoGES_chr_normalized := character()]
  source_candidates[, condition_KoGES_status := character()]
}

fwrite(source_candidates, source_candidates_file, sep = "\t", quote = FALSE, na = "NA")

candidate_reference_qc_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_level_candidates_KoGES_reference_QC.tsv"
)
lead_reference_qc_file <- file.path(
  OUT_DIR,
  "joint_COJO_target_leads_KoGES_reference_QC.tsv"
)
fwrite(source_candidates, candidate_reference_qc_file, sep = "\t", quote = FALSE, na = "NA")
fwrite(
  manifest[, c(
    target_columns, "lead_reference_chr", "lead_reference_bp",
    "lead_reference_status"
  ), with = FALSE],
  lead_reference_qc_file, sep = "\t", quote = FALSE, na = "NA"
)

# ------------------------------------------------------------------------------
# 4. Confirm source-specific clump SNPs in the appropriate ancestry-matched
#    1000G panel.  KoGES availability is intentionally NOT a clumping gate.
# ------------------------------------------------------------------------------

panel_candidates <- copy(source_candidates)
if (nrow(panel_candidates) > 0L) {
  panel_candidates[, clump_panel := evidence_layer]
  panel_candidates[, panel_reference_chr := NA_character_]
  panel_candidates[, panel_reference_bp := NA_integer_]
  panel_candidates[, panel_reference_a1 := NA_character_]
  panel_candidates[, panel_reference_a2 := NA_character_]

  eas_i <- which(panel_candidates$clump_panel == "EAS")
  if (length(eas_i) > 0L) {
    m <- match(panel_candidates$clump_snp[eas_i], eas_bim$ref_snp)
    panel_candidates[eas_i, panel_reference_chr := eas_bim$ref_chr[m]]
    panel_candidates[eas_i, panel_reference_bp := eas_bim$ref_bp[m]]
    panel_candidates[eas_i, panel_reference_a1 := eas_bim$ref_a1[m]]
    panel_candidates[eas_i, panel_reference_a2 := eas_bim$ref_a2[m]]
  }
  eur_i <- which(panel_candidates$clump_panel == "EUR")
  if (length(eur_i) > 0L) {
    m <- match(panel_candidates$clump_snp[eur_i], eur_bim$ref_snp)
    panel_candidates[eur_i, panel_reference_chr := eur_bim$ref_chr[m]]
    panel_candidates[eur_i, panel_reference_bp := eur_bim$ref_bp[m]]
    panel_candidates[eur_i, panel_reference_a1 := eur_bim$ref_a1[m]]
    panel_candidates[eur_i, panel_reference_a2 := eur_bim$ref_a2[m]]
  }

  panel_candidates[, panel_reference_chr_normalized := normalize_chr(panel_reference_chr)]
  panel_candidates[, panel_reference_bp_matches_KoGES :=
    !is.na(condition_KoGES_bp) & !is.na(panel_reference_bp) &
      condition_KoGES_bp == panel_reference_bp]
  panel_candidates[, panel_reference_alleles_match_KoGES := allele_pair_matches(
    panel_reference_a1, panel_reference_a2,
    condition_KoGES_a1, condition_KoGES_a2
  )]
  panel_candidates[, panel_reference_snp_indel_alias := snp_indel_chrbp_alias(
    panel_reference_a1, panel_reference_a2,
    condition_KoGES_a1, condition_KoGES_a2
  )]
  panel_candidates[, panel_reference_shared_allele := allele_pair_shares_one(
    panel_reference_a1, panel_reference_a2,
    condition_KoGES_a1, condition_KoGES_a2
  )]
  panel_candidates[, panel_to_KoGES_identity_status := fcase(
    condition_mapping_collision,
    "ambiguous_multiple_KoGES_condition_ids",
    !is_valid_snp_id(condition_snp),
    "no_resolved_KoGES_condition_id_proxy_may_rescue",
    is.na(condition_KoGES_bp),
    "missing_resolved_condition_in_KoGES_bim_proxy_may_rescue",
    is.na(panel_reference_bp),
    "missing_clump_snp_in_1000G_reference",
    !is.na(condition_KoGES_bp) &
      panel_reference_chr_normalized != condition_KoGES_chr_normalized,
    "1000G_to_KoGES_chr_mismatch",
    !panel_reference_bp_matches_KoGES,
    "1000G_to_KoGES_bp_mismatch",
    panel_reference_alleles_match_KoGES,
    EXACT_IDENTITY_STATUS,
    panel_reference_bp_matches_KoGES & panel_reference_snp_indel_alias,
    DISCORDANT_SNP_INDEL_STATUS,
    panel_reference_bp_matches_KoGES & panel_reference_shared_allele,
    DISCORDANT_SHARED_ALLELE_STATUS,
    default = "1000G_to_KoGES_allele_mismatch"
  )]
  panel_candidates[, clump_reference_eligibility := fcase(
    !is_valid_snp_id(clump_snp), "missing_1000G_clump_snp",
    is.na(panel_reference_bp),
    fifelse(clump_panel == "EAS", "missing_in_1000G_EAS", "missing_in_1000G_EUR"),
    panel_reference_chr_normalized != target_chr_normalized,
    "1000G_clump_snp_chr_different_from_target",
    abs(panel_reference_bp - koges_bp) > NOVELTY_WINDOW_BP,
    "1000G_clump_snp_outside_target_window",
    isTRUE(EXACT_1000G_CLUMP_ONLY) &
      !all_exact_1000g_clump_statuses(clump_snp_resolution_statuses),
    "not_exact_1000G_identity",
    default = "eligible"
  )]
} else {
  panel_candidates[, clump_panel := character()]
  panel_candidates[, panel_reference_chr := character()]
  panel_candidates[, panel_reference_bp := integer()]
  panel_candidates[, panel_reference_a1 := character()]
  panel_candidates[, panel_reference_a2 := character()]
  panel_candidates[, panel_reference_chr_normalized := character()]
  panel_candidates[, clump_reference_eligibility := character()]
  panel_candidates[, panel_reference_bp_matches_KoGES := logical()]
  panel_candidates[, panel_reference_alleles_match_KoGES := logical()]
  panel_candidates[, panel_reference_snp_indel_alias := logical()]
  panel_candidates[, panel_reference_shared_allele := logical()]
  panel_candidates[, panel_to_KoGES_identity_status := character()]
}

panel_reference_qc_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_level_candidates_1000G_reference_QC.tsv"
)
fwrite(panel_candidates, panel_reference_qc_file, sep = "\t", quote = FALSE, na = "NA")

n_representation_discordant <- panel_candidates[
  panel_to_KoGES_identity_status %chin% c(
    DISCORDANT_SNP_INDEL_STATUS, DISCORDANT_SHARED_ALLELE_STATUS
  ),
  .N
]
n_hard_identity_mismatch <- panel_candidates[
  panel_to_KoGES_identity_status %chin% HARD_IDENTITY_MISMATCH_STATUSES,
  .N
]
if (n_representation_discordant + n_hard_identity_mismatch > 0L) {
  message(
    n_representation_discordant,
    " source-level priors are representation-discordant (same-coordinate ",
    "SNP-indel or shared-allele) and ",
    n_hard_identity_mismatch,
    " have CHR/BP/disjoint-allele mismatch. ",
    "They remain in 1000G clumping. KoGES conditioning requires an exact ",
    "identity match or an in-clump 1000G LD proxy at r^2 >= ",
    PROXY_MIN_R2, "."
  )
}

# ------------------------------------------------------------------------------
# 5. Source-stratified ancestry-matched PLINK clumping
# ------------------------------------------------------------------------------

clump_jobs <- unique(panel_candidates[, .(
  target_id, source_key, source, source_output_stem, clump_panel
)])
setorder(clump_jobs, target_id, source_key, clump_panel)

cat(
  "Running ", nrow(clump_jobs),
  " target-source ancestry-clump jobs with N_PARALLEL=", N_PARALLEL, "\n",
  "(parallel PLINK clump starts now; prep above was data.table-threaded, not forked)\n",
  sep = ""
)

clump_results <- if (nrow(clump_jobs) == 0L) {
  list()
} else {
  parallel_lapply(
    seq_len(nrow(clump_jobs)),
    function(i) run_one_clump_job(as.list(clump_jobs[i]), panel_candidates),
    n_workers = N_PARALLEL
  )
}

clump_result_summary <- if (length(clump_results) == 0L) {
  data.table(
    target_id = character(), source_key = character(), source = character(),
    source_output_stem = character(), clump_panel = character(),
    clump_status = character(), n_selected = integer(), selected_snps = character(),
    n_clump_indices = integer(), n_clump_members = integer(),
    assoc_file = character(), clumped_file = character(), out_prefix = character()
  )
} else {
  rbindlist(lapply(clump_results, function(x) {
    data.table(
      target_id = x$target_id,
      source_key = x$source_key,
      source = x$source,
      source_output_stem = x$source_output_stem,
      clump_panel = x$clump_panel,
      clump_status = x$status,
      n_selected = length(x$selected),
      selected_snps = if (length(x$selected)) {
        paste(x$selected, collapse = ";")
      } else {
        NA_character_
      },
      n_clump_indices = uniqueN(x$membership$clump_id),
      n_clump_members = nrow(x$membership),
      assoc_file = x$assoc_file,
      clumped_file = x$clumped_file,
      out_prefix = x$out_prefix
    )
  }), use.names = TRUE, fill = TRUE)
}

clump_summary_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_ancestry_clump_job_summary.tsv"
)
fwrite(clump_result_summary, clump_summary_file, sep = "\t", quote = FALSE, na = "NA")

bad_jobs <- clump_result_summary[
  !(clump_status %chin% c(
    "completed", "singleton_retained", "no_1000G_clumpable_candidates"
  ))
]
if (nrow(bad_jobs) > 0L) {
  stop(
    "One or more source-specific PLINK clump jobs failed. See: ",
    clump_summary_file
  )
}

clump_membership <- if (length(clump_results) == 0L) {
  empty_clump_membership()
} else {
  rbindlist(lapply(clump_results, function(x) x$membership),
            use.names = TRUE, fill = TRUE)
}
if (nrow(clump_membership) > 0L) {
  clump_membership <- unique(clump_membership)
}

# A successful clump must have an auditable index/member map.  Failing here is
# preferable to later treating a missing SP2 field as though it were a complete
# prior-signal definition.
successful_job_keys <- clump_result_summary[
  clump_status %chin% c("completed", "singleton_retained"),
  paste(target_id, source_key, clump_panel, sep = "\t")
]
member_job_keys <- if (nrow(clump_membership) == 0L) character() else unique(
  paste(clump_membership$target_id, clump_membership$source_key,
        clump_membership$clump_panel, sep = "\t")
)
missing_member_map <- setdiff(successful_job_keys, member_job_keys)
if (length(missing_member_map) > 0L) {
  stop(
    "A successful PLINK clump job has no auditable index/member map. Example: ",
    missing_member_map[[1L]]
  )
}

fwrite(clump_membership, clump_member_audit_file,
       sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

} # end of !RESUME_AFTER_CLUMP (prep + ancestry PLINK clumping)

if (isTRUE(skip_post_clump_selection)) {
  message(
    "Loading already-written post-clump selection tables; ",
    "not re-running representative batches or cross-source KoGES LD."
  )
  clump_representatives <- fread(
    clump_representative_file, sep = "\t", showProgress = FALSE
  )
  within_source_rep_collisions <- fread(
    representative_collision_file, sep = "\t", showProgress = FALSE
  )
  cross_source_ld_target_summary <- fread(
    cross_source_ld_summary_file, sep = "\t", showProgress = FALSE
  )
  final_candidates <- fread(
    final_candidates_file, sep = "\t", showProgress = FALSE
  )
  if ("representative_usable_in_KoGES" %in% names(clump_representatives)) {
    clump_representatives[, representative_usable_in_KoGES :=
      as.logical(representative_usable_in_KoGES)]
  }
  if (!"representation_evidence" %in% names(clump_representatives)) {
    clump_representatives[, representation_evidence := fcase(
      representative_type == "direct_clump_index", "exact_variant",
      representative_type %chin% c(
        "member_proxy_high_LD_1000G", "panel_wide_proxy_high_LD_1000G"
      ),
      "ld_proxy_r2_ge_threshold",
      default = "unrepresented"
    )]
  }
  if (!"proxy_cross_panel_status" %in% names(clump_representatives)) {
    clump_representatives[, proxy_cross_panel_status := "not_evaluated"]
  }
  if (!"proxy_r2_1000G_EAS" %in% names(clump_representatives)) {
    clump_representatives[, proxy_r2_1000G_EAS := NA_real_]
  }
  if (!"representative_primary_eligible" %in% names(clump_representatives)) {
    message(
      "Loaded clump representatives predate EAS confirmation of EUR-panel LD ",
      "proxies.  Treating every EUR-panel proxy as transferability-unconfirmed ",
      "until Step2 re-runs representative selection. This affects the evidence ",
      "grade, not COJO execution or EAS-primary eligibility."
    )
    clump_representatives[, representative_primary_eligible := fcase(
      representation_evidence == "exact_variant", TRUE,
      representation_evidence != "ld_proxy_r2_ge_threshold", FALSE,
      clump_panel == "EAS", TRUE,
      default = FALSE
    )]
  } else {
    clump_representatives[, representative_primary_eligible :=
      as.logical(representative_primary_eligible)]
  }
  for (col in c("selected_in_EAS", "selected_in_EUR", "is_retained_for_cojo_cond")) {
    if (col %in% names(final_candidates)) {
      final_candidates[, (col) := as.logical(get(col))]
    }
  }
  if (!"direct_KoGES_rescue" %in% names(final_candidates)) {
    final_candidates[, direct_KoGES_rescue := FALSE]
  }
  if (!"n_direct_KoGES_rescued_hits" %in% names(final_candidates)) {
    final_candidates[, n_direct_KoGES_rescued_hits := 0L]
  }
  if (!"direct_KoGES_rescue_statuses" %in% names(final_candidates)) {
    final_candidates[, direct_KoGES_rescue_statuses := NA_character_]
  }
} else {

# Panel-wide proxy rescue needs the target window and the set of ancestry-panel
# variants that are the same variant in KoGES.  Build both before forking, then
# release the panel BIMs so the representative workers stay small.
target_lead_lookup <- manifest[, .(
  target_id = as.character(target_id),
  koges_lead_variant = as.character(koges_lead_variant),
  target_chr_normalized = normalize_chr(koges_chr),
  koges_bp = as.numeric(koges_bp)
)]
setkey(target_lead_lookup, target_id)

panel_koges_exact_maps <- list()
eur_index_identity_in_eas <- build_index_identity_across_panels(
  character(), data.table(), data.table()
)
if (nrow(clump_membership) > 0L) {
  if (!exists("koges_bim", inherits = FALSE)) {
    koges_bim <- read_bim_index(KOGES_BFILE)
  }
  built_panel_bims <- character()
  if (!exists("eas_bim", inherits = FALSE)) {
    eas_bim <- read_bim_index(EAS_BFILE)
    built_panel_bims <- c(built_panel_bims, "eas_bim")
  }
  if (!exists("eur_bim", inherits = FALSE)) {
    eur_bim <- if (exists("eas_bim", inherits = FALSE)) eas_bim[0] else {
      empty <- read_bim_index(EAS_BFILE)
      empty[0]
    }
    built_panel_bims <- c(built_panel_bims, "eur_bim")
  }
  panel_koges_exact_maps$EAS <- build_panel_koges_exact_map(
    eas_bim, koges_bim, "EAS"
  )
  panel_koges_exact_maps$EUR <- build_panel_koges_exact_map(
    eur_bim, koges_bim, "EUR"
  )
  if (VERIFY_EUR_PROXY_IN_EAS) {
    eur_index_identity_in_eas <- build_index_identity_across_panels(
      clump_membership[clump_panel == "EUR" & member_is_index %in% TRUE,
                       clump_index_snp],
      eur_bim, eas_bim
    )
  }
  if (length(built_panel_bims) > 0L) {
    rm(list = built_panel_bims)
    invisible(gc())
  }
}

# Compute panel-wide proxy LD once per (target, ancestry panel) before the
# representative fork wave.  Both tables are read-only in the workers, so fork
# copy-on-write shares them instead of duplicating per clump.
panel_wide_group <- build_panel_wide_group_ld(
  clump_membership, target_lead_lookup
)
panel_wide_group_ld <- panel_wide_group$ld
panel_wide_group_status <- panel_wide_group$status
rm(panel_wide_group)
invisible(gc())

clump_rep_bundle <- select_clump_representatives_batched(clump_membership)
clump_representatives <- clump_rep_bundle$representatives
member_proxy_audit <- clump_rep_bundle$proxy_audit
panel_wide_proxy_audit <- clump_rep_bundle$panel_wide_audit
rm(clump_rep_bundle)
invisible(gc())
if (nrow(clump_representatives) != uniqueN(clump_membership$clump_id) &&
    !(nrow(clump_membership) == 0L && nrow(clump_representatives) == 0L)) {
  stop(
    "Failed to create exactly one representative record per PLINK clump. Got ",
    nrow(clump_representatives), " representatives for ",
    uniqueN(clump_membership$clump_id), " clumps."
  )
}

cross_panel_proxy <- if (VERIFY_EUR_PROXY_IN_EAS) {
  verify_eur_ld_proxies_in_eas(
    clump_representatives, eur_index_identity_in_eas, panel_koges_exact_maps,
    CLUMP_WORK_DIR
  )
} else {
  data.table()
}
if (nrow(cross_panel_proxy) > 0L) {
  cross_match <- match(clump_representatives$clump_id, cross_panel_proxy$clump_id)
  clump_representatives[, proxy_r2_1000G_EAS :=
    cross_panel_proxy$proxy_r2_1000G_EAS[cross_match]]
  clump_representatives[, proxy_cross_panel_status :=
    cross_panel_proxy$proxy_cross_panel_status[cross_match]]
} else {
  clump_representatives[, proxy_r2_1000G_EAS := NA_real_]
  clump_representatives[, proxy_cross_panel_status := fcase(
    clump_panel == "EUR" & representation_evidence == "ld_proxy_r2_ge_threshold" &
      !VERIFY_EUR_PROXY_IN_EAS,
    "not_run_not_required_for_EUR_supporting_analysis",
    default = NA_character_
  )]
}
clump_representatives[is.na(proxy_cross_panel_status),
  proxy_cross_panel_status := fcase(
    representation_evidence == "exact_variant", "not_applicable_exact_variant",
    representation_evidence == "unrepresented", "not_applicable_unrepresented",
    clump_panel == "EAS", "not_applicable_EAS_panel_LD_proxy",
    default = "not_evaluated"
  )]
# An exact KoGES variant needs no LD argument.  An EAS-panel proxy is already
# estimated in the ancestry closest to KoGES.  EUR proxy transferability to EAS
# is optional QC only; it is not an execution, EAS-primary, or EUR-supporting
# exclusion gate in the ancestry-stratified design.
clump_representatives[, representative_primary_eligible := fcase(
  representation_evidence == "exact_variant", TRUE,
  representation_evidence != "ld_proxy_r2_ge_threshold", FALSE,
  clump_panel == "EAS", TRUE,
  VERIFY_EUR_PROXY_IN_EAS & proxy_cross_panel_status %chin% CROSS_PANEL_CONFIRMED,
  TRUE,
  default = FALSE
)]

clump_representative_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_clump_representatives.tsv"
)
member_proxy_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_member_proxy_LD_audit.tsv"
)
panel_wide_proxy_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_stratified_panel_wide_proxy_LD_audit.tsv"
)
fwrite(clump_representatives, clump_representative_file,
       sep = "\t", quote = FALSE, na = "NA")
fwrite(member_proxy_audit, member_proxy_audit_file,
       sep = "\t", quote = FALSE, na = "NA")
fwrite(panel_wide_proxy_audit, panel_wide_proxy_audit_file,
       sep = "\t", quote = FALSE, na = "NA")

# A representative collision within one source/panel would mean two distinct
# source-level clumps have been collapsed to the same KoGES condition SNP.  Do
# not call that a complete joint condition set.  The same SNP across *different*
# external sources is intentionally unioned and is not a collision.
successful_clump_reps <- clump_representatives[
  representative_usable_in_KoGES %in% TRUE &
    is_valid_snp_id(selected_condition_snp)
]
if (nrow(successful_clump_reps) > 0L) {
  # A panel-wide proxy is not an external hit, so its source provenance and
  # clump P value come from the clump index it represents.
  successful_clump_reps[, provenance_via_clump_index :=
    representative_type == "panel_wide_proxy_high_LD_1000G"]
  successful_clump_reps[, provenance_clump_snp := fifelse(
    provenance_via_clump_index,
    as.character(clump_index_snp),
    as.character(selected_clump_member_snp)
  )]
  representative_key <- paste(
    successful_clump_reps$target_id, successful_clump_reps$source_key,
    successful_clump_reps$clump_panel,
    successful_clump_reps$provenance_clump_snp, sep = "\t"
  )
  panel_candidate_key <- paste(
    panel_candidates$target_id, panel_candidates$source_key,
    panel_candidates$clump_panel, panel_candidates$clump_snp, sep = "\t"
  )
  representative_match <- match(representative_key, panel_candidate_key)
  successful_clump_reps[, selected_clump_p :=
    panel_candidates$clump_p[representative_match]]
} else {
  successful_clump_reps[, provenance_via_clump_index := logical()]
  successful_clump_reps[, provenance_clump_snp := character()]
  successful_clump_reps[, selected_clump_p := numeric()]
}
within_source_rep_collisions <- successful_clump_reps[, .(
  n_distinct_clump_indices = uniqueN(clump_index_snp),
  clump_index_snps = collapse_unique(clump_index_snp)
), by = .(target_id, source_key, clump_panel, selected_condition_snp)][
  n_distinct_clump_indices > 1L
]
representative_collision_file <- file.path(
  OUT_DIR,
  "joint_COJO_within_source_representative_collisions.tsv"
)
fwrite(within_source_rep_collisions, representative_collision_file,
       sep = "\t", quote = FALSE, na = "NA")

selected_keys <- successful_clump_reps[, .(
  target_id, source_key, clump_panel,
  clump_snp = as.character(provenance_clump_snp),
  condition_snp = as.character(selected_condition_snp),
  provenance_via_clump_index = provenance_via_clump_index %in% TRUE
)]
selected_keys <- unique(selected_keys)

panel_candidates[, selected_as_source_representative := FALSE]
if (nrow(selected_keys) > 0L && nrow(panel_candidates) > 0L) {
  selected_key <- paste(
    selected_keys$target_id, selected_keys$source_key,
    selected_keys$clump_panel, selected_keys$clump_snp, sep = "\t"
  )
  candidate_key <- paste(
    panel_candidates$target_id, panel_candidates$source_key,
    panel_candidates$clump_panel, panel_candidates$clump_snp, sep = "\t"
  )
  panel_candidates[
    candidate_key %chin% selected_key,
    selected_as_source_representative := TRUE
  ]
}
if (any(
  panel_candidates$selected_as_source_representative %in% TRUE &
    !(panel_candidates$clump_reference_eligibility %chin% "eligible")
)) {
  stop("A non-eligible 1000G-to-KoGES candidate was selected as a clump representative.")
}

panel_clump_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_source_level_candidates_with_clump_selection.tsv"
)
fwrite(panel_candidates, panel_clump_audit_file, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------------------
# 6. Union usable direct/proxy representatives, then apply the pre-specified
#    cross-source KoGES-LD representative rule
# ------------------------------------------------------------------------------

# The condition ID can be the external clump index or a rigorously selected
# in-clump proxy.  Preserve both IDs and the exact source-panel r^2 evidence in
# the final candidate file so the reviewer-facing conclusion is auditable.
selection_by_snp <- if (nrow(successful_clump_reps) == 0L) {
  data.table(
    target_id = character(), condition_snp = character(),
    selected_in_EAS = logical(), selected_in_EUR = logical(),
    selected_panels = character(), selected_sources = character(),
    selected_source_output_stems = character(), selected_source_keys = character(),
    n_distinct_selected_source_keys = integer(), panel_clump_p_min = numeric(),
    external_clump_index_snps = character(), selected_clump_member_snps = character(),
    representative_types = character(), representative_statuses = character(),
    proxy_r2_1000G_min = numeric(),
    proxy_r2_1000G_values = character(), n_external_clumps_represented = integer()
  )
} else {
  successful_clump_reps[, .(
    selected_in_EAS = any(clump_panel == "EAS"),
    selected_in_EUR = any(clump_panel == "EUR"),
    selected_panels = collapse_unique(clump_panel),
    selected_sources = collapse_unique(source),
    selected_source_output_stems = collapse_unique(source_output_stem),
    selected_source_keys = collapse_unique(source_key),
    n_distinct_selected_source_keys = as.integer(uniqueN(source_key)),
    panel_clump_p_min = min_finite(selected_clump_p),
    external_clump_index_snps = collapse_unique(clump_index_snp),
    selected_clump_member_snps = collapse_unique(selected_clump_member_snp),
    representative_types = collapse_unique(representative_type),
    representative_statuses = collapse_unique(representative_status),
    proxy_r2_1000G_min = min_finite(proxy_r2_1000G),
    proxy_r2_1000G_values = format_numeric_values(proxy_r2_1000G),
    n_external_clumps_represented = .N
  ), by = .(target_id, condition_snp = selected_condition_snp)]
}

# Build provenance only from source candidates that were actually selected as
# direct indices / high-LD clump-member representatives.  This avoids attributing
# a source to a final condition SNP merely because that SNP was present but not
# selected in that source's clumping result.
provenance_key_columns <- c("target_id", "source_key", "clump_panel", "clump_snp")
selected_source_rows_direct <- merge(
  selected_keys[
    provenance_via_clump_index %in% FALSE,
    c(provenance_key_columns, "condition_snp"), with = FALSE
  ],
  panel_candidates,
  by = c(provenance_key_columns, "condition_snp"),
  all.x = TRUE,
  sort = FALSE
)
# A panel-wide proxy has no panel_candidates row of its own; attribute it to the
# single clump-index row it represents, then keep the proxy as the condition SNP.
selected_source_rows_panel_wide <- if (
  nrow(selected_keys[provenance_via_clump_index %in% TRUE]) == 0L
) {
  NULL
} else {
  matched <- merge(
    selected_keys[
      provenance_via_clump_index %in% TRUE,
      c(provenance_key_columns, "condition_snp"), with = FALSE
    ],
    unique(panel_candidates, by = provenance_key_columns)[
      , setdiff(names(panel_candidates), "condition_snp"), with = FALSE
    ],
    by = provenance_key_columns,
    all.x = TRUE,
    sort = FALSE
  )
  matched
}
selected_source_rows <- rbindlist(
  list(selected_source_rows_direct, selected_source_rows_panel_wide),
  use.names = TRUE, fill = TRUE
)

candidate_union <- if (nrow(selected_source_rows) == 0L) {
  data.table(
    target_id = character(), lead_row_id = integer(), koges_phenotype = character(),
    koges_pheno_code = character(), koges_lead_variant = character(),
    koges_chr = integer(), koges_bp = integer(), condition_snp = character(),
    reported_in_sources = character(), reported_in_layers = character(),
    source_output_stems = character(), source_p_min = numeric(),
    source_p_values = character(), source_abs_distance_min_bp = numeric(),
    source_abs_distance_max_bp = numeric(), source_variant_ids = character(),
    n_raw_source_hits = integer()
  )
} else {
  selected_source_rows[, .(
    reported_in_sources = collapse_unique(source),
    reported_in_layers = collapse_unique(evidence_layer),
    source_output_stems = collapse_unique(source_output_stem),
    source_p_min = min_finite(source_p_min),
    source_p_values = collapse_source_sets(source_p_values),
    source_abs_distance_min_bp = min_finite(source_abs_distance_min_bp),
    source_abs_distance_max_bp = max_finite(source_abs_distance_max_bp),
    source_variant_ids = collapse_source_sets(source_variant_ids),
    n_raw_source_hits = sum(n_raw_source_hits)
  ), by = c(target_columns, "condition_snp")]
}

# DIRECT_KOGES_RESCUE_IN_CORE=TRUE: union 1000G-unmatched KoGES-exact rescues
# into the candidate set.  They still lack 1000G-defined independence, so the
# later KoGES r^2 grouping is allowed to fold same-source rescue SNPs. The
# rescue audit table remains the row-level provenance of those additions.
if (!"direct_KoGES_rescue" %in% names(candidate_union)) {
  candidate_union[, direct_KoGES_rescue := FALSE]
}
if (!"n_direct_KoGES_rescued_hits" %in% names(candidate_union)) {
  candidate_union[, n_direct_KoGES_rescued_hits := 0L]
}
if (!"direct_KoGES_rescue_statuses" %in% names(candidate_union)) {
  candidate_union[, direct_KoGES_rescue_statuses := NA_character_]
}
if (!"direct_KoGES_rescue" %in% names(selection_by_snp)) {
  selection_by_snp[, direct_KoGES_rescue := FALSE]
}
if (!"n_direct_KoGES_rescued_hits" %in% names(selection_by_snp)) {
  selection_by_snp[, n_direct_KoGES_rescued_hits := 0L]
}
if (!"direct_KoGES_rescue_statuses" %in% names(selection_by_snp)) {
  selection_by_snp[, direct_KoGES_rescue_statuses := NA_character_]
}

if (DIRECT_KOGES_RESCUE_IN_CORE && !is.null(direct_koges_rescue_union)) {
  candidate_union <- rbindlist(
    list(candidate_union, direct_koges_rescue_union), use.names = TRUE, fill = TRUE
  )[, .(
    reported_in_sources = collapse_unique(reported_in_sources),
    reported_in_layers = collapse_unique(reported_in_layers),
    source_output_stems = collapse_unique(source_output_stems),
    source_p_min = min_finite(source_p_min),
    source_p_values = collapse_source_sets(source_p_values),
    source_abs_distance_min_bp = min_finite(source_abs_distance_min_bp),
    source_abs_distance_max_bp = max_finite(source_abs_distance_max_bp),
    source_variant_ids = collapse_source_sets(source_variant_ids),
    n_raw_source_hits = sum(n_raw_source_hits),
    direct_KoGES_rescue = any(direct_KoGES_rescue %in% TRUE),
    n_direct_KoGES_rescued_hits = sum(n_direct_KoGES_rescued_hits),
    direct_KoGES_rescue_statuses = collapse_unique(direct_KoGES_rescue_statuses)
  ), by = c(target_columns, "condition_snp")]
}
if (DIRECT_KOGES_RESCUE_IN_CORE && !is.null(direct_koges_rescue_selection)) {
  selection_by_snp <- rbindlist(
    list(selection_by_snp, direct_koges_rescue_selection),
    use.names = TRUE, fill = TRUE
  )[, .(
    selected_in_EAS = any(selected_in_EAS %in% TRUE),
    selected_in_EUR = any(selected_in_EUR %in% TRUE),
    selected_panels = collapse_unique(selected_panels),
    selected_sources = collapse_unique(selected_sources),
    selected_source_output_stems = collapse_unique(selected_source_output_stems),
    selected_source_keys = collapse_source_sets(selected_source_keys),
    n_distinct_selected_source_keys = NA_integer_,
    panel_clump_p_min = min_finite(panel_clump_p_min),
    external_clump_index_snps = collapse_unique(external_clump_index_snps),
    selected_clump_member_snps = collapse_unique(selected_clump_member_snps),
    representative_types = collapse_unique(representative_types),
    representative_statuses = collapse_unique(representative_statuses),
    proxy_r2_1000G_min = min_finite(proxy_r2_1000G_min),
    proxy_r2_1000G_values = collapse_source_sets(proxy_r2_1000G_values),
    n_external_clumps_represented = sum(n_external_clumps_represented),
    direct_KoGES_rescue = any(direct_KoGES_rescue %in% TRUE),
    n_direct_KoGES_rescued_hits = sum(n_direct_KoGES_rescued_hits),
    direct_KoGES_rescue_statuses = collapse_unique(direct_KoGES_rescue_statuses)
  ), by = .(target_id, condition_snp)]
  selection_by_snp[, n_distinct_selected_source_keys := as.integer(
    vapply(selected_source_keys, n_distinct_semicolon_values, integer(1L))
  )]
}

final_candidates <- merge(
  candidate_union, selection_by_snp,
  by = c("target_id", "condition_snp"),
  all = FALSE, sort = FALSE
)

if (nrow(final_candidates) > 0L) {
  selected_match <- match(final_candidates$condition_snp, koges_bim$ref_snp)
  final_candidates[, condition_reference_chr := koges_bim$ref_chr[selected_match]]
  final_candidates[, condition_reference_bp := koges_bim$ref_bp[selected_match]]
  final_candidates[, target_chr_normalized := normalize_chr(koges_chr)]
  final_candidates[, condition_reference_chr_normalized :=
    normalize_chr(condition_reference_chr)]
  final_candidates[, selected_condition_KoGES_status := fcase(
    is.na(condition_reference_bp), "missing_selected_condition_in_KoGES_bim",
    condition_reference_chr_normalized != target_chr_normalized,
    "selected_condition_chr_different_from_target",
    abs(condition_reference_bp - koges_bp) > NOVELTY_WINDOW_BP,
    "selected_condition_outside_target_window",
    condition_snp == koges_lead_variant, "selected_condition_is_target_lead",
    default = "eligible"
  )]
} else {
  final_candidates[, condition_reference_chr := character()]
  final_candidates[, condition_reference_bp := integer()]
  final_candidates[, target_chr_normalized := character()]
  final_candidates[, condition_reference_chr_normalized := character()]
  final_candidates[, selected_condition_KoGES_status := character()]
}

if (anyDuplicated(final_candidates[, .(target_id, condition_snp)])) {
  stop("Duplicate target-condition pair after union of source representatives.")
}

# Preserve the source-level union before any KoGES-LD representative selection;
# this is the audit trail needed to distinguish a source-specific external
# signal from a cross-source duplicate signal.
final_candidates_pre_cross_source_ld <- copy(final_candidates)
final_candidates_pre_cross_source_ld_file <- file.path(
  OUT_DIR,
  "joint_COJO_condition_candidates_pre_cross_source_KoGES_LD.tsv"
)
fwrite(
  final_candidates_pre_cross_source_ld,
  final_candidates_pre_cross_source_ld_file,
  sep = "\t", quote = FALSE, na = "NA"
)

cross_source_ld_groups <- if (nrow(final_candidates) == 0L) {
  list()
} else {
  split(seq_len(nrow(final_candidates)), final_candidates$target_id)
}
cross_source_ld_results <- parallel_lapply(
  cross_source_ld_groups,
  function(i) run_one_cross_source_ld_grouping(final_candidates[i]),
  n_workers = N_MEM_PARALLEL
)

cross_source_ld_member_audit <- if (length(cross_source_ld_results) == 0L) {
  empty_cross_source_ld_member_audit()
} else {
  rbindlist(
    lapply(cross_source_ld_results, function(x) x$member_audit),
    use.names = TRUE, fill = TRUE
  )
}
cross_source_ld_pair_audit <- if (length(cross_source_ld_results) == 0L) {
  empty_cross_source_ld_pair_table()
} else {
  rbindlist(
    lapply(cross_source_ld_results, function(x) x$pair_audit),
    use.names = TRUE, fill = TRUE
  )
}
cross_source_ld_target_summary <- if (length(cross_source_ld_results) == 0L) {
  empty_cross_source_ld_target_summary()
} else {
  rbindlist(
    lapply(cross_source_ld_results, function(x) x$target_summary),
    use.names = TRUE, fill = TRUE
  )
}

if (nrow(final_candidates_pre_cross_source_ld) > 0L &&
    uniqueN(cross_source_ld_target_summary$target_id) !=
      uniqueN(final_candidates_pre_cross_source_ld$target_id)) {
  stop("Missing cross-source KoGES-LD grouping summary for one or more targets.")
}
nonrepresentative_cross_source_rows <- cross_source_ld_member_audit[
  cross_source_ld_grouping_status == "completed" &
    cross_source_ld_is_representative %in% FALSE
]
if (nrow(nonrepresentative_cross_source_rows) > 0L && any(
  !is.finite(nonrepresentative_cross_source_rows$cross_source_ld_r2_to_representative) |
    nonrepresentative_cross_source_rows$cross_source_ld_r2_to_representative <
      CROSS_SOURCE_REP_R2
)) {
  stop("A collapsed cross-source condition SNP lacks the required direct KoGES LD evidence.")
}
cross_source_group_source_check <- cross_source_ld_member_audit[
  !is.na(cross_source_ld_group_id),
  .(
    n_members = .N,
    n_distinct_source_keys = n_distinct_semicolon_values(selected_source_keys),
    # Members of one multi-SNP group must have pairwise-disjoint source sets.
    has_overlapping_member_sources = {
      source_lists <- lapply(selected_source_keys, split_semicolon_values)
      overlap <- FALSE
      if (length(source_lists) >= 2L) {
        for (i in seq_len(length(source_lists) - 1L)) {
          for (j in seq.int(i + 1L, length(source_lists))) {
            if (length(intersect(source_lists[[i]], source_lists[[j]])) > 0L) {
              overlap <- TRUE
              break
            }
          }
          if (overlap) break
        }
      }
      overlap
    }
  ),
  by = .(target_id, cross_source_ld_group_id)
]
if (any(
  cross_source_group_source_check$n_members > 1L &
    cross_source_group_source_check$n_distinct_source_keys < 2L
)) {
  stop("A cross-source KoGES-LD group contains more than one representative from one source only.")
}
if (any(cross_source_group_source_check$has_overlapping_member_sources %in% TRUE)) {
  stop("A cross-source KoGES-LD group contains members with overlapping source sets.")
}

if (length(cross_source_ld_results) > 0L) {
  final_candidates <- rbindlist(
    lapply(cross_source_ld_results, function(x) x$collapsed_candidates),
    use.names = TRUE, fill = TRUE
  )
} else {
  final_candidates[, `:=`(
    cross_source_ld_group_id = character(),
    cross_source_ld_group_size = integer(),
    cross_source_ld_group_member_snps = character(),
    cross_source_ld_is_representative = logical(),
    cross_source_ld_representative_condition_snp = character(),
    cross_source_ld_r2_to_representative = numeric(),
    cross_source_ld_grouping_status = character(),
    cross_source_ld_grouping_rule = character()
  )]
}

cross_source_ld_member_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_group_members.tsv"
)
cross_source_ld_pair_audit_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_pairs.tsv"
)
cross_source_ld_summary_file <- file.path(
  OUT_DIR,
  "joint_COJO_cross_source_KoGES_LD_grouping_summary.tsv"
)
fwrite(cross_source_ld_member_audit, cross_source_ld_member_audit_file,
       sep = "\t", quote = FALSE, na = "NA")
fwrite(cross_source_ld_pair_audit, cross_source_ld_pair_audit_file,
       sep = "\t", quote = FALSE, na = "NA")
fwrite(cross_source_ld_target_summary, cross_source_ld_summary_file,
       sep = "\t", quote = FALSE, na = "NA")

if (nrow(final_candidates) > 0L) {
  setorder(
    final_candidates, lead_row_id, source_p_min,
    source_abs_distance_min_bp, condition_snp
  )
  final_candidates[, condition_rank := seq_len(.N), by = target_id]
  final_candidates[, is_retained_for_cojo_cond := TRUE]
  final_candidates[, conditioning_strategy := fifelse(
    direct_KoGES_rescue %in% TRUE,
    paste0(
      "source_stratified_1000G_clump_plus_exact_direct_KoGES_rescue_then_",
      "cross_source_pairwise_KoGES_r2_ge_",
      format(CROSS_SOURCE_REP_R2, scientific = FALSE, trim = TRUE),
      "_representative_then_KoGES_COJO"
    ),
    paste0(
      "source_stratified_target_local_1000G_clump_then_high_LD_member_proxy_if_index_absent_then_",
      "cross_source_pairwise_KoGES_r2_ge_",
      format(CROSS_SOURCE_REP_R2, scientific = FALSE, trim = TRUE),
      "_representative_then_KoGES_COJO"
    )
  )]
} else {
  final_candidates[, condition_rank := integer()]
  final_candidates[, is_retained_for_cojo_cond := logical()]
  final_candidates[, conditioning_strategy := character()]
}

if (anyDuplicated(final_candidates[, .(target_id, condition_snp)])) {
  stop("Duplicate target-condition pair after cross-source KoGES-LD grouping.")
}

final_candidates_file <- file.path(
  OUT_DIR,
  "joint_COJO_condition_candidates_ancestry_clumped.tsv"
)
fwrite(final_candidates, final_candidates_file, sep = "\t", quote = FALSE, na = "NA")

# ------------------------------------------------------------------------------
# 6b. Build ancestry-stratified core condition sets
# ------------------------------------------------------------------------------
# The all-source table above is retained as an audit of the union of EAS and EUR
# evidence.  It is not the primary inference set: EUR index variants are often
# absent from the KoGES common-imputed panel for biological/panel-coverage
# reasons, not because of an ID or build error.  Build EAS and EUR condition
# sets independently from their own source-level clump representatives, then
# (when DIRECT_KOGES_RESCUE_IN_CORE=TRUE) rbind 1000G-unmatched KoGES-exact
# rescues into that panel's core.  Do not rebuild EAS core from clumps only.
if (!exists("selected_source_rows", inherits = FALSE)) {
  stop(
    "EAS/EUR-stratified V4 selection requires selected_source_rows. ",
    "Re-run post-clump selection rather than loading legacy final candidates."
  )
}
if (!exists("koges_bim", inherits = FALSE)) {
  koges_bim <- read_bim_index(KOGES_BFILE)
}

empty_panel_core_candidates <- function() {
  copy(final_candidates[0L])
}

build_panel_core_candidates <- function(panel) {
  if (!panel %chin% c("EAS", "EUR")) stop("panel must be EAS or EUR.")
  panel_in_eas <- identical(panel, "EAS")
  panel_rows <- selected_source_rows[clump_panel == panel]
  panel_reps <- successful_clump_reps[
    clump_panel == panel & representative_usable_in_KoGES %in% TRUE &
      is_valid_snp_id(selected_condition_snp)
  ]

  clump_candidates <- empty_panel_core_candidates()
  if (nrow(panel_rows) > 0L && nrow(panel_reps) > 0L) {
    panel_union <- panel_rows[, .(
      reported_in_sources = collapse_unique(source),
      reported_in_layers = collapse_unique(evidence_layer),
      source_output_stems = collapse_unique(source_output_stem),
      source_p_min = min_finite(source_p_min),
      source_p_values = collapse_source_sets(source_p_values),
      source_abs_distance_min_bp = min_finite(source_abs_distance_min_bp),
      source_abs_distance_max_bp = max_finite(source_abs_distance_max_bp),
      source_variant_ids = collapse_source_sets(source_variant_ids),
      n_raw_source_hits = as.integer(sum(n_raw_source_hits, na.rm = TRUE))
    ), by = c(target_columns, "condition_snp")]

    panel_selection <- panel_reps[, .(
      selected_in_EAS = panel_in_eas,
      selected_in_EUR = !panel_in_eas,
      selected_panels = panel,
      selected_sources = collapse_unique(source),
      selected_source_output_stems = collapse_unique(source_output_stem),
      selected_source_keys = collapse_unique(source_key),
      n_distinct_selected_source_keys = as.integer(uniqueN(source_key)),
      panel_clump_p_min = min_finite(selected_clump_p),
      external_clump_index_snps = collapse_unique(clump_index_snp),
      selected_clump_member_snps = collapse_unique(selected_clump_member_snp),
      representative_types = collapse_unique(representative_type),
      representative_statuses = collapse_unique(representative_status),
      proxy_r2_1000G_min = min_finite(proxy_r2_1000G),
      proxy_r2_1000G_values = format_numeric_values(proxy_r2_1000G),
      n_external_clumps_represented = as.integer(.N)
    ), by = .(target_id, condition_snp = selected_condition_snp)]

    clump_candidates <- merge(
      panel_union, panel_selection,
      by = c("target_id", "condition_snp"), all = FALSE, sort = FALSE
    )
    if (nrow(clump_candidates) == 0L) {
      stop("Selected ", panel, " clump representatives could not be joined to provenance.")
    }
    clump_candidates[, `:=`(
      direct_KoGES_rescue = FALSE,
      n_direct_KoGES_rescued_hits = 0L,
      direct_KoGES_rescue_statuses = NA_character_
    )]
  }

  rescue_candidates <- empty_panel_core_candidates()
  if (
    isTRUE(DIRECT_KOGES_RESCUE_IN_CORE) &&
      !is.null(direct_koges_rescue_union) &&
      !is.null(direct_koges_rescue_selection)
  ) {
    rescue_sel <- copy(if (panel_in_eas) {
      direct_koges_rescue_selection[selected_in_EAS %in% TRUE]
    } else {
      direct_koges_rescue_selection[selected_in_EUR %in% TRUE]
    })
    if (nrow(rescue_sel) > 0L) {
      rescue_keys <- rescue_sel[, paste(target_id, condition_snp)]
      rescue_union <- direct_koges_rescue_union[
        paste(target_id, condition_snp) %chin% rescue_keys
      ]
      overlap <- setdiff(
        intersect(names(rescue_union), names(rescue_sel)),
        c("target_id", "condition_snp")
      )
      if (length(overlap) > 0L) rescue_sel[, (overlap) := NULL]
      rescue_candidates <- merge(
        rescue_union, rescue_sel,
        by = c("target_id", "condition_snp"), all = FALSE, sort = FALSE
      )
    }
  }

  candidates <- rbindlist(
    list(clump_candidates, rescue_candidates),
    use.names = TRUE, fill = TRUE
  )
  if (nrow(candidates) == 0L) {
    return(empty_panel_core_candidates())
  }
  if (anyDuplicated(candidates, by = c("target_id", "condition_snp"))) {
    candidates <- candidates[, .(
      reported_in_sources = collapse_unique(reported_in_sources),
      reported_in_layers = collapse_unique(reported_in_layers),
      source_output_stems = collapse_unique(source_output_stems),
      source_p_min = min_finite(source_p_min),
      source_p_values = collapse_source_sets(source_p_values),
      source_abs_distance_min_bp = min_finite(source_abs_distance_min_bp),
      source_abs_distance_max_bp = max_finite(source_abs_distance_max_bp),
      source_variant_ids = collapse_source_sets(source_variant_ids),
      n_raw_source_hits = as.integer(sum(n_raw_source_hits, na.rm = TRUE)),
      selected_in_EAS = any(selected_in_EAS %in% TRUE),
      selected_in_EUR = any(selected_in_EUR %in% TRUE),
      selected_panels = collapse_unique(selected_panels),
      selected_sources = collapse_unique(selected_sources),
      selected_source_output_stems = collapse_unique(selected_source_output_stems),
      selected_source_keys = collapse_source_sets(selected_source_keys),
      n_distinct_selected_source_keys = NA_integer_,
      panel_clump_p_min = min_finite(panel_clump_p_min),
      external_clump_index_snps = collapse_unique(external_clump_index_snps),
      selected_clump_member_snps = collapse_unique(selected_clump_member_snps),
      representative_types = collapse_unique(representative_types),
      representative_statuses = collapse_unique(representative_statuses),
      proxy_r2_1000G_min = min_finite(proxy_r2_1000G_min),
      proxy_r2_1000G_values = collapse_source_sets(proxy_r2_1000G_values),
      n_external_clumps_represented = as.integer(sum(
        n_external_clumps_represented, na.rm = TRUE
      )),
      direct_KoGES_rescue = any(direct_KoGES_rescue %in% TRUE),
      n_direct_KoGES_rescued_hits = as.integer(sum(
        n_direct_KoGES_rescued_hits, na.rm = TRUE
      )),
      direct_KoGES_rescue_statuses = collapse_unique(direct_KoGES_rescue_statuses)
    ), by = c(target_columns, "condition_snp")]
  }
  candidates[, n_distinct_selected_source_keys := as.integer(
    vapply(selected_source_keys, n_distinct_semicolon_values, integer(1L))
  )]

  candidate_match <- match(candidates$condition_snp, koges_bim$ref_snp)
  candidates[, `:=`(
    condition_reference_chr = koges_bim$ref_chr[candidate_match],
    condition_reference_bp = koges_bim$ref_bp[candidate_match],
    target_chr_normalized = normalize_chr(koges_chr),
    condition_reference_chr_normalized = normalize_chr(
      koges_bim$ref_chr[candidate_match]
    )
  )]
  candidates[, selected_condition_KoGES_status := fcase(
    is.na(condition_reference_bp), "missing_selected_condition_in_KoGES_bim",
    condition_reference_chr_normalized != target_chr_normalized,
    "selected_condition_chr_different_from_target",
    abs(condition_reference_bp - koges_bp) > NOVELTY_WINDOW_BP,
    "selected_condition_outside_target_window",
    condition_snp == koges_lead_variant, "selected_condition_is_target_lead",
    default = "eligible"
  )]
  candidates
}

collapse_panel_core_candidates <- function(candidates, analysis_label) {
  if (nrow(candidates) == 0L) {
    return(list(
      candidates = candidates,
      member_audit = empty_cross_source_ld_member_audit(),
      pair_audit = empty_cross_source_ld_pair_table(),
      target_summary = empty_cross_source_ld_target_summary()
    ))
  }
  groups <- split(seq_len(nrow(candidates)), candidates$target_id)
  results <- parallel_lapply(
    groups,
    function(i) run_one_cross_source_ld_grouping(
      candidates[i], analysis_label = analysis_label
    ),
    n_workers = N_MEM_PARALLEL
  )
  list(
    candidates = rbindlist(
      lapply(results, function(x) x$collapsed_candidates),
      use.names = TRUE, fill = TRUE
    ),
    member_audit = rbindlist(
      lapply(results, function(x) x$member_audit),
      use.names = TRUE, fill = TRUE
    ),
    pair_audit = rbindlist(
      lapply(results, function(x) x$pair_audit),
      use.names = TRUE, fill = TRUE
    ),
    target_summary = rbindlist(
      lapply(results, function(x) x$target_summary),
      use.names = TRUE, fill = TRUE
    )
  )
}

finalize_panel_core_candidates <- function(candidates, analysis_label) {
  if (nrow(candidates) == 0L) return(candidates)
  setorder(
    candidates, lead_row_id, source_p_min,
    source_abs_distance_min_bp, condition_snp
  )
  candidates[, condition_rank := seq_len(.N), by = target_id]
  candidates[, is_retained_for_cojo_cond := TRUE]
  candidates[, conditioning_strategy := paste0(
    analysis_label,
    "_source_stratified_target_local_1000G_clump_then_high_LD_proxy",
    if (isTRUE(DIRECT_KOGES_RESCUE_IN_CORE)) {
      "_plus_direct_KoGES_exact_rescue"
    } else {
      ""
    },
    "_then_cross_source_pairwise_KoGES_r2_ge_",
    format(CROSS_SOURCE_REP_R2, scientific = FALSE, trim = TRUE),
    "_representative_then_KoGES_COJO"
  )]
  candidates
}

summarize_panel_core_candidates <- function(candidates, panel) {
  if (!panel %chin% c("EAS", "EUR")) stop("panel must be EAS or EUR.")
  if (nrow(candidates) == 0L) return(data.table(target_id = character()))
  if (!"direct_KoGES_rescue" %in% names(candidates)) {
    candidates[, direct_KoGES_rescue := FALSE]
  }
  if (!"n_external_clumps_represented" %in% names(candidates)) {
    candidates[, n_external_clumps_represented := 0L]
  }
  prefix <- paste0(panel, "_")
  out <- candidates[, .(
    n_condition_snps_selected = .N,
    n_condition_snps_from_1000G_clumps = as.integer(sum(
      fcoalesce(as.integer(n_external_clumps_represented), 0L) > 0L
    )),
    n_condition_snps_direct_KoGES_rescue =
      as.integer(sum(direct_KoGES_rescue %in% TRUE)),
    n_condition_snps_selected_in_KoGES =
      sum(selected_condition_KoGES_status == "eligible"),
    n_selected_condition_snps_not_usable_in_KoGES =
      sum(selected_condition_KoGES_status != "eligible"),
    selected_condition_KoGES_failure_reasons = collapse_unique(
      selected_condition_KoGES_status[
        selected_condition_KoGES_status != "eligible"
      ]
    ),
    condition_snp_ids_selected = collapse_unique(condition_snp),
    condition_snp_ids_selected_in_KoGES = collapse_unique(
      condition_snp[selected_condition_KoGES_status == "eligible"]
    ),
    condition_sources_selected = collapse_source_sets(selected_sources),
    condition_source_output_stems_selected = collapse_source_sets(
      selected_source_output_stems
    )
  ), by = target_id]
  setnames(
    out,
    setdiff(names(out), "target_id"),
    paste0(prefix, setdiff(names(out), "target_id"))
  )
  out
}

summarize_panel_clumps <- function(panel) {
  if (!panel %chin% c("EAS", "EUR")) stop("panel must be EAS or EUR.")
  reps <- clump_representatives[clump_panel == panel]
  if (nrow(reps) == 0L) return(data.table(target_id = character()))
  prefix <- paste0(panel, "_")
  out <- reps[, .(
    n_external_clump_signals = .N,
    n_clump_representatives_exact_variant =
      sum(representation_evidence == "exact_variant"),
    n_clump_representatives_ld_proxy =
      sum(representation_evidence == "ld_proxy_r2_ge_threshold"),
    n_clump_representatives_panel_wide_proxy =
      sum(representative_type == "panel_wide_proxy_high_LD_1000G"),
    n_clump_representatives_unrepresented = sum(
      representation_evidence == "unrepresented" |
        !(representative_usable_in_KoGES %in% TRUE)
    ),
    unrepresented_clump_signal_reasons = collapse_unique(
      representative_status[
        representation_evidence == "unrepresented" |
          !(representative_usable_in_KoGES %in% TRUE)
      ]
    )
  ), by = target_id]
  setnames(
    out,
    setdiff(names(out), "target_id"),
    paste0(prefix, setdiff(names(out), "target_id"))
  )
  out
}

summarize_panel_collisions <- function(panel) {
  collisions <- within_source_rep_collisions[clump_panel == panel]
  if (nrow(collisions) == 0L) return(data.table(target_id = character()))
  prefix <- paste0(panel, "_")
  out <- collisions[, .(
    n_within_source_representative_collisions = .N,
    within_source_representative_collision_details = paste(
      paste(source_key, clump_panel, selected_condition_snp,
            clump_index_snps, sep = "|"),
      collapse = ";"
    )
  ), by = target_id]
  setnames(
    out,
    setdiff(names(out), "target_id"),
    paste0(prefix, setdiff(names(out), "target_id"))
  )
  out
}

panel_core_bundles <- list()
for (panel in c("EAS", "EUR")) {
  analysis_label <- if (panel == "EAS") "EAS_primary" else "EUR_supporting"
  bundle <- collapse_panel_core_candidates(
    build_panel_core_candidates(panel), analysis_label = analysis_label
  )
  bundle$candidates <- finalize_panel_core_candidates(
    bundle$candidates, analysis_label = analysis_label
  )
  panel_core_bundles[[panel]] <- bundle

  candidate_file <- file.path(
    OUT_DIR, paste0("joint_COJO_", analysis_label, "_condition_candidates.tsv")
  )
  member_file <- file.path(
    OUT_DIR, paste0("joint_COJO_", analysis_label, "_cross_source_KoGES_LD_group_members.tsv")
  )
  pair_file <- file.path(
    OUT_DIR, paste0("joint_COJO_", analysis_label, "_cross_source_KoGES_LD_pairs.tsv")
  )
  summary_file <- file.path(
    OUT_DIR, paste0("joint_COJO_", analysis_label, "_cross_source_KoGES_LD_grouping_summary.tsv")
  )
  fwrite(bundle$candidates, candidate_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(bundle$member_audit, member_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(bundle$pair_audit, pair_file, sep = "\t", quote = FALSE, na = "NA")
  fwrite(bundle$target_summary, summary_file, sep = "\t", quote = FALSE, na = "NA")
}

EAS_core_candidates <- panel_core_bundles$EAS$candidates
EUR_core_candidates <- panel_core_bundles$EUR$candidates
EAS_core_candidate_summary <- summarize_panel_core_candidates(
  EAS_core_candidates, "EAS"
)
EUR_core_candidate_summary <- summarize_panel_core_candidates(
  EUR_core_candidates, "EUR"
)
EAS_clump_summary <- summarize_panel_clumps("EAS")
EUR_clump_summary <- summarize_panel_clumps("EUR")
EAS_collision_summary <- summarize_panel_collisions("EAS")
EUR_collision_summary <- summarize_panel_collisions("EUR")
EAS_cross_source_summary <- panel_core_bundles$EAS$target_summary
EUR_cross_source_summary <- panel_core_bundles$EUR$target_summary
if (nrow(EAS_cross_source_summary) > 0L) {
  setnames(
    EAS_cross_source_summary,
    setdiff(names(EAS_cross_source_summary), "target_id"),
    paste0("EAS_", setdiff(names(EAS_cross_source_summary), "target_id"))
  )
}
if (nrow(EUR_cross_source_summary) > 0L) {
  setnames(
    EUR_cross_source_summary,
    setdiff(names(EUR_cross_source_summary), "target_id"),
    paste0("EUR_", setdiff(names(EUR_cross_source_summary), "target_id"))
  )
}

} # end of !skip_post_clump_selection

# ------------------------------------------------------------------------------
# 7. Completeness accounting and final Step2 manifest
# ------------------------------------------------------------------------------

source_candidate_counts <- source_candidates[, .(
  n_source_condition_candidates_preclump = .N
), by = target_id]

panel_counts <- panel_candidates[, .(
  n_panel_candidates_preclump = .N,
  n_panel_candidates_clumpable_1000G =
    sum(clump_reference_eligibility == "eligible"),
  n_panel_candidates_exact_identity = sum(
    panel_to_KoGES_identity_status == EXACT_IDENTITY_STATUS
  ),
  n_panel_candidates_identity_mismatch = sum(
    panel_to_KoGES_identity_status %chin% HARD_IDENTITY_MISMATCH_STATUSES
  ),
  n_panel_candidates_snp_indel_alias = sum(
    panel_to_KoGES_identity_status == DISCORDANT_SNP_INDEL_STATUS
  ),
  n_panel_candidates_shared_allele_multiallelic = sum(
    panel_to_KoGES_identity_status == DISCORDANT_SHARED_ALLELE_STATUS
  ),
  n_panel_candidates_representation_discordant = sum(
    panel_to_KoGES_identity_status %chin% c(
      DISCORDANT_SNP_INDEL_STATUS, DISCORDANT_SHARED_ALLELE_STATUS
    )
  ),
  n_panel_candidates_unclumpable_1000G =
    sum(clump_reference_eligibility != "eligible"),
  unclumpable_prior_reasons = collapse_unique(
    clump_reference_eligibility[clump_reference_eligibility != "eligible"]
  ),
  identity_mismatch_prior_reasons = collapse_unique(
    panel_to_KoGES_identity_status[
      panel_to_KoGES_identity_status %chin% c(
        HARD_IDENTITY_MISMATCH_STATUSES,
        DISCORDANT_SNP_INDEL_STATUS,
        DISCORDANT_SHARED_ALLELE_STATUS
      )
    ]
  )
), by = target_id]

# Count external clump signals before deduplicating their KoGES condition IDs.
# This is what protects completeness when several source-level clumps happen
# to union to the same physical condition SNP across different external sources.
clump_representative_summary <- clump_representatives[, .(
  n_external_clump_signals = .N,
  n_clump_representatives_direct_index =
    sum(representative_type == "direct_clump_index"),
  n_clump_representatives_member_proxy =
    sum(representative_type == "member_proxy_high_LD_1000G"),
  n_clump_representatives_panel_wide_proxy =
    sum(representative_type == "panel_wide_proxy_high_LD_1000G"),
  n_clump_representatives_exact_variant =
    sum(representation_evidence == "exact_variant"),
  n_clump_representatives_ld_proxy =
    sum(representation_evidence == "ld_proxy_r2_ge_threshold"),
  n_clump_representatives_unrepresented =
    sum(representation_evidence == "unrepresented" |
      !(representative_usable_in_KoGES %in% TRUE)),
  n_clump_representatives_eur_only_ld_proxy = sum(
    representation_evidence == "ld_proxy_r2_ge_threshold" &
      representative_usable_in_KoGES %in% TRUE &
      !(representative_primary_eligible %in% TRUE)
  ),
  eur_only_ld_proxy_reasons = collapse_unique(
    proxy_cross_panel_status[
      representation_evidence == "ld_proxy_r2_ge_threshold" &
        representative_usable_in_KoGES %in% TRUE &
        !(representative_primary_eligible %in% TRUE)
    ]
  ),
  unrepresented_clump_signal_reasons = collapse_unique(
    representative_status[representation_evidence == "unrepresented" |
      !(representative_usable_in_KoGES %in% TRUE)]
  ),
  selected_proxy_r2_1000G_min = min_finite(
    proxy_r2_1000G[representation_evidence == "ld_proxy_r2_ge_threshold"]
  ),
  selected_proxy_r2_1000G_values = format_numeric_values(
    proxy_r2_1000G[representation_evidence == "ld_proxy_r2_ge_threshold"]
  )
), by = target_id]

representative_collision_summary <- within_source_rep_collisions[, .(
  n_within_source_representative_collisions = .N,
  within_source_representative_collision_details = paste(
    paste(source_key, clump_panel, selected_condition_snp,
          clump_index_snps, sep = "|"),
    collapse = ";"
  )
), by = target_id]

selected_summary <- final_candidates[, .(
  n_condition_snps_selected_external = .N,
  n_condition_snps_direct_KoGES_rescue = sum(direct_KoGES_rescue %in% TRUE),
  n_direct_KoGES_rescued_hits_selected = sum(
    n_direct_KoGES_rescued_hits, na.rm = TRUE
  ),
  direct_KoGES_rescue_statuses_selected = collapse_source_sets(
    direct_KoGES_rescue_statuses
  ),
  n_condition_snps_selected_in_KoGES =
    sum(selected_condition_KoGES_status == "eligible"),
  n_selected_condition_snps_not_usable_in_KoGES =
    sum(selected_condition_KoGES_status != "eligible"),
  selected_condition_KoGES_failure_reasons = collapse_unique(
    selected_condition_KoGES_status[selected_condition_KoGES_status != "eligible"]
  ),
  condition_snp_ids_selected_external = collapse_unique(condition_snp),
  condition_snp_ids_selected_in_KoGES = collapse_unique(
    condition_snp[selected_condition_KoGES_status == "eligible"]
  ),
  condition_selected_panels = collapse_unique(selected_panels),
  condition_sources_selected = collapse_source_sets(selected_sources),
  condition_source_output_stems_selected =
    collapse_source_sets(selected_source_output_stems),
  n_EAS_selected = sum(selected_in_EAS),
  n_EUR_selected = sum(selected_in_EUR)
), by = target_columns]

final_manifest <- copy(manifest)
final_manifest[, step1_cojo_condition_status := cojo_condition_status]
final_manifest <- merge_into_manifest(final_manifest, raw_input_summary)
final_manifest <- merge_into_manifest(final_manifest, source_candidate_counts)
final_manifest <- merge_into_manifest(final_manifest, panel_counts)
final_manifest <- merge_into_manifest(final_manifest, clump_representative_summary)
final_manifest <- merge_into_manifest(final_manifest, representative_collision_summary)
final_manifest <- merge_into_manifest(final_manifest, cross_source_ld_target_summary)
final_manifest <- merge_into_manifest(final_manifest, selected_summary)
final_manifest <- merge_into_manifest(final_manifest, EAS_clump_summary)
final_manifest <- merge_into_manifest(final_manifest, EUR_clump_summary)
final_manifest <- merge_into_manifest(final_manifest, EAS_collision_summary)
final_manifest <- merge_into_manifest(final_manifest, EUR_collision_summary)
final_manifest <- merge_into_manifest(final_manifest, EAS_core_candidate_summary)
final_manifest <- merge_into_manifest(final_manifest, EUR_core_candidate_summary)
final_manifest <- merge_into_manifest(final_manifest, EAS_cross_source_summary)
final_manifest <- merge_into_manifest(final_manifest, EUR_cross_source_summary)

count_columns <- c(
  "n_raw_prior_hits_input", "n_raw_prior_hits_invalid", "n_raw_prior_hits_valid",
  "n_raw_prior_hits_represented_1000G", "n_raw_prior_hits_unclumpable_1000G",
  "n_raw_prior_hits_direct_KoGES_rescued",
  "n_raw_prior_hits_unresolved_after_direct_KoGES_rescue",
  "n_source_condition_candidates_preclump", "n_panel_candidates_preclump",
  "n_panel_candidates_clumpable_1000G", "n_panel_candidates_identity_mismatch",
  "n_panel_candidates_exact_identity",
  "n_panel_candidates_snp_indel_alias",
  "n_panel_candidates_shared_allele_multiallelic",
  "n_panel_candidates_representation_discordant",
  "n_panel_candidates_unclumpable_1000G",
  "n_external_clump_signals", "n_clump_representatives_direct_index",
  "n_clump_representatives_member_proxy",
  "n_clump_representatives_panel_wide_proxy",
  "n_clump_representatives_exact_variant", "n_clump_representatives_ld_proxy",
  "n_clump_representatives_unrepresented",
  "n_clump_representatives_eur_only_ld_proxy",
  "n_within_source_representative_collisions",
  "n_condition_snps_pre_cross_source_KoGES_LD",
  "n_condition_snps_post_cross_source_KoGES_LD",
  "n_condition_snps_collapsed_cross_source_KoGES_LD",
  "n_cross_source_KoGES_LD_representative_groups",
  "n_cross_source_KoGES_high_LD_pairs",
  "n_KoGES_LD_collapse_eligible_pairs",
  "n_condition_snps_selected_external", "n_condition_snps_selected_in_KoGES",
  "n_condition_snps_selected",
  "n_condition_snps_from_1000G_clumps",
  "n_condition_snps_direct_KoGES_rescue",
  "n_direct_KoGES_rescued_hits_selected",
  "n_selected_condition_snps_not_usable_in_KoGES", "n_EAS_selected",
  "n_EUR_selected",
  "EAS_n_external_clump_signals",
  "EAS_n_clump_representatives_exact_variant",
  "EAS_n_clump_representatives_ld_proxy",
  "EAS_n_clump_representatives_panel_wide_proxy",
  "EAS_n_clump_representatives_unrepresented",
  "EAS_n_within_source_representative_collisions",
  "EAS_n_condition_snps_selected",
  "EAS_n_condition_snps_from_1000G_clumps",
  "EAS_n_condition_snps_direct_KoGES_rescue",
  "EAS_n_condition_snps_selected_in_KoGES",
  "EAS_n_selected_condition_snps_not_usable_in_KoGES",
  "EUR_n_external_clump_signals",
  "EUR_n_clump_representatives_exact_variant",
  "EUR_n_clump_representatives_ld_proxy",
  "EUR_n_clump_representatives_panel_wide_proxy",
  "EUR_n_clump_representatives_unrepresented",
  "EUR_n_within_source_representative_collisions",
  "EUR_n_condition_snps_selected",
  "EUR_n_condition_snps_selected_in_KoGES",
  "EUR_n_selected_condition_snps_not_usable_in_KoGES"
)
for (nm in count_columns) {
  if (!nm %in% names(final_manifest)) final_manifest[, (nm) := 0L]
  final_manifest[is.na(get(nm)), (nm) := 0L]
}
character_columns <- c(
  "prior_input_failure_reasons", "reference_coverage_reasons",
  "direct_KoGES_rescue_statuses", "direct_KoGES_rescue_statuses_selected",
  "unclumpable_prior_reasons",
  "identity_mismatch_prior_reasons",
  "unrepresented_clump_signal_reasons", "eur_only_ld_proxy_reasons",
  "selected_proxy_r2_1000G_values",
  "within_source_representative_collision_details",
  "cross_source_KoGES_LD_grouping_status",
  "cross_source_KoGES_LD_candidate_file", "cross_source_KoGES_LD_file",
  "selected_condition_KoGES_failure_reasons",
  "condition_snp_ids_selected_external", "condition_snp_ids_selected_in_KoGES",
  "condition_selected_panels", "condition_sources_selected",
  "condition_source_output_stems_selected",
  "EAS_unrepresented_clump_signal_reasons",
  "EAS_within_source_representative_collision_details",
  "EAS_selected_condition_KoGES_failure_reasons",
  "EAS_condition_snp_ids_selected",
  "EAS_condition_snp_ids_selected_in_KoGES",
  "EAS_condition_sources_selected",
  "EAS_condition_source_output_stems_selected",
  "EAS_cross_source_KoGES_LD_grouping_status",
  "EAS_cross_source_KoGES_LD_candidate_file",
  "EAS_cross_source_KoGES_LD_file",
  "EUR_unrepresented_clump_signal_reasons",
  "EUR_within_source_representative_collision_details",
  "EUR_selected_condition_KoGES_failure_reasons",
  "EUR_condition_snp_ids_selected",
  "EUR_condition_snp_ids_selected_in_KoGES",
  "EUR_condition_sources_selected",
  "EUR_condition_source_output_stems_selected",
  "EUR_cross_source_KoGES_LD_grouping_status",
  "EUR_cross_source_KoGES_LD_candidate_file",
  "EUR_cross_source_KoGES_LD_file"
)
for (nm in character_columns) {
  if (!nm %in% names(final_manifest)) final_manifest[, (nm) := NA_character_]
}

final_manifest[, prior_signal_selection_scope := CLUMP_SCOPE]
final_manifest[, prior_coverage_status := fcase(
  n_raw_prior_hits_input == 0L,
  NA_character_,
  n_raw_prior_hits_represented_1000G > 0L &
    n_raw_prior_hits_unclumpable_1000G == 0L,
  "full_1000G_coverage",
  n_raw_prior_hits_represented_1000G > 0L &
    n_raw_prior_hits_unclumpable_1000G > 0L,
  "partial_1000G_coverage",
  default = "no_representable_prior"
)]
final_manifest[, prior_set_coverage_status := fcase(
  n_raw_prior_hits_input == 0L,
  NA_character_,
  n_raw_prior_hits_invalid > 0L,
  "invalid_external_prior_input",
  n_raw_prior_hits_unresolved_after_direct_KoGES_rescue == 0L &
    n_raw_prior_hits_direct_KoGES_rescued == 0L,
  "complete_1000G_coverage",
  n_raw_prior_hits_unresolved_after_direct_KoGES_rescue == 0L &
    n_raw_prior_hits_direct_KoGES_rescued > 0L,
  "complete_with_direct_KoGES_rescue",
  n_raw_prior_hits_direct_KoGES_rescued > 0L,
  "partial_with_direct_KoGES_rescue",
  default = "partial_unresolved_prior"
)]
final_manifest[, prior_identity_status := fcase(
  n_panel_candidates_preclump == 0L,
  NA_character_,
  n_panel_candidates_identity_mismatch == 0L &
    n_panel_candidates_representation_discordant == 0L,
  "exact_1000G_to_KoGES_identity",
  n_panel_candidates_identity_mismatch == 0L,
  "representation_discordant_present",
  default = "identity_mismatch_present"
)]
# Raw-hit availability is intentionally continuous QC, not a binary gate.
# It describes whether a reported row is in its source-matched 1000G panel or
# has an exact KoGES identity; it does not claim that every row is an
# independent signal or that every exact rescue was included as a covariate.
final_manifest[, raw_prior_reference_or_exact_KoGES_fraction := fifelse(
  n_raw_prior_hits_valid > 0L,
  pmin(
    1,
    (n_raw_prior_hits_represented_1000G +
       n_raw_prior_hits_direct_KoGES_rescued) / n_raw_prior_hits_valid
  ),
  NA_real_
)]

EAS_grouping_ok <- final_manifest$EAS_cross_source_KoGES_LD_grouping_status %chin% c(
  "completed", "not_needed_singleton"
)
EUR_grouping_ok <- final_manifest$EUR_cross_source_KoGES_LD_grouping_status %chin% c(
  "completed", "not_needed_singleton"
)
all_source_grouping_ok <- final_manifest$cross_source_KoGES_LD_grouping_status %chin% c(
  "completed", "not_needed_singleton"
)

# EAS_prior_present is EAS external GWS prior evidence, not "a 1000G clump
# exists".  In this EAS-only script that is n_raw_prior_hits_valid > 0.
final_manifest[, `:=`(
  EAS_prior_present = n_raw_prior_hits_valid > 0L,
  EUR_prior_present = EUR_n_external_clump_signals > 0L
)]
final_manifest[, EUR_only_prior := !EAS_prior_present & EUR_prior_present]

# Completeness of the *clump component* does not require any clump to exist.
# Rescue-only targets have 0 EAS 1000G clumps; that is allowed.  If clumps
# do exist, every one must have a KoGES representative and no within-source
# representative collision.
final_manifest[, EAS_clump_component_complete := (
  EAS_n_external_clump_signals == 0L |
    (
      EAS_n_clump_representatives_unrepresented == 0L &
        EAS_n_within_source_representative_collisions == 0L
    )
)]
# Legacy name: the EAS condition set is complete for primary COJO.  Clumps
# are optional; selected KoGES-usable SNPs (clump reps and/or direct rescue)
# and successful cross-source grouping are not.
final_manifest[, EAS_clump_signal_set_complete := (
  EAS_clump_component_complete &
    EAS_n_condition_snps_selected > 0L &
    EAS_n_selected_condition_snps_not_usable_in_KoGES == 0L &
    EAS_grouping_ok
)]
final_manifest[, EUR_clump_signal_set_complete := (
  EUR_prior_present &
    EUR_n_clump_representatives_unrepresented == 0L &
    EUR_n_within_source_representative_collisions == 0L &
    EUR_n_condition_snps_selected > 0L &
    EUR_n_selected_condition_snps_not_usable_in_KoGES == 0L &
    EUR_grouping_ok
)]
final_manifest[, all_source_clump_signal_set_complete := (
  n_external_clump_signals > 0L &
    n_clump_representatives_unrepresented == 0L &
    n_within_source_representative_collisions == 0L &
    n_condition_snps_selected_external > 0L &
    n_selected_condition_snps_not_usable_in_KoGES == 0L &
    all_source_grouping_ok
)]

# Legacy name retained for audit compatibility; it now means all *clumps* are
# represented, not that every raw external hit is individually available.
final_manifest[, prior_signal_set_complete := all_source_clump_signal_set_complete]

# EAS-primary is the only primary result for the Korean cohort.  EUR runs use
# a separate condition set and are retained as supporting evidence, including
# for EUR-only positional prior loci.  Rescue-only EAS loci are eligible:
# they have prior evidence even when EAS_n_external_clump_signals == 0.
final_manifest[, EAS_primary_ready := (
  step1_cojo_condition_status != "known_exact_variant" &
    lead_reference_status == "eligible" &
    n_raw_prior_hits_invalid == 0L &
    EAS_prior_present &
    EAS_clump_component_complete &
    EAS_n_condition_snps_selected > 0L &
    EAS_n_selected_condition_snps_not_usable_in_KoGES == 0L &
    EAS_grouping_ok
)]
final_manifest[, EUR_supporting_ready := (
  step1_cojo_condition_status != "known_exact_variant" &
    lead_reference_status == "eligible" &
    n_raw_prior_hits_invalid == 0L &
    EUR_clump_signal_set_complete
)]

# Legacy fields consumed by existing downstream scripts now point to the
# EAS-primary model.  All-source completeness remains available explicitly.
final_manifest[, cojo_ready := EAS_primary_ready]
final_manifest[, strict_primary_ready := EAS_primary_ready]
final_manifest[, sensitivity_ready := EUR_supporting_ready]
final_manifest[, proxy_transferability_status :=
  "not_applicable_EAS_only_analysis"]
final_manifest[, primary_gate_failure_reason := fcase(
  strict_primary_ready %in% TRUE, NA_character_,
  step1_cojo_condition_status == "known_exact_variant", NA_character_,
  n_raw_prior_hits_invalid > 0L, "invalid_external_prior_input",
  lead_reference_status != "eligible", "target_lead_not_usable_in_KoGES_reference",
  EUR_only_prior, "EUR_only_prior_no_EAS_primary_test",
  !EAS_prior_present, "no_EAS_prior_evidence",
  !EAS_clump_component_complete & EAS_n_clump_representatives_unrepresented > 0L,
  "unrepresented_EAS_clump_signal",
  !EAS_clump_component_complete & EAS_n_within_source_representative_collisions > 0L,
  "within_source_EAS_representative_collision",
  EAS_n_condition_snps_selected == 0L,
  "no_EAS_condition_snp_from_clump_or_KoGES_rescue",
  EAS_n_selected_condition_snps_not_usable_in_KoGES > 0L,
  "selected_EAS_condition_not_usable_in_KoGES",
  !EAS_grouping_ok, "EAS_cross_source_KoGES_LD_grouping_failed",
  default = "EAS_primary_clump_set_incomplete"
)]
final_manifest[, EUR_supporting_gate_failure_reason := fcase(
  EUR_supporting_ready %in% TRUE, NA_character_,
  step1_cojo_condition_status == "known_exact_variant", NA_character_,
  n_raw_prior_hits_invalid > 0L, "invalid_external_prior_input",
  lead_reference_status != "eligible", "target_lead_not_usable_in_KoGES_reference",
  !EUR_prior_present, "no_EUR_prior_clump",
  EUR_n_clump_representatives_unrepresented > 0L,
  "unrepresented_EUR_clump_signal",
  EUR_n_within_source_representative_collisions > 0L,
  "within_source_EUR_representative_collision",
  EUR_n_condition_snps_selected == 0L,
  "no_EUR_condition_representative_after_source_clumping",
  EUR_n_selected_condition_snps_not_usable_in_KoGES > 0L,
  "selected_EUR_condition_not_usable_in_KoGES",
  !EUR_grouping_ok, "EUR_cross_source_KoGES_LD_grouping_failed",
  default = "EUR_supporting_clump_set_incomplete"
)]

# Keep the legacy names used by Step3a/3b, bound to the EAS-primary condition
# set.  n_condition_snps_ancestry_clumped is the *total* selected count
# (1000G clump representatives plus direct rescue).  It is not clump-only.
final_manifest[, n_condition_snps_selected := EAS_n_condition_snps_selected]
final_manifest[, n_condition_snps_from_1000G_clumps :=
  EAS_n_condition_snps_from_1000G_clumps]
final_manifest[, n_condition_snps_available := EAS_n_condition_snps_selected]
final_manifest[, n_condition_snps_ancestry_clumped := EAS_n_condition_snps_selected]
final_manifest[, n_condition_snps_KoGES_mappable :=
  EAS_n_condition_snps_selected_in_KoGES]
final_manifest[, condition_snp_ids_available := fifelse(
  cojo_ready,
  EAS_condition_snp_ids_selected,
  NA_character_
)]
final_manifest[, condition_snp_ids_sensitivity := fifelse(
  sensitivity_ready,
  EUR_condition_snp_ids_selected,
  NA_character_
)]

final_manifest[, cojo_condition_status := fcase(
  step1_cojo_condition_status == "known_exact_variant",
  "direct_overlap_only_no_conditional_test_needed",
  lead_reference_status != "eligible",
  "target_lead_not_usable_in_KoGES_reference",
  n_raw_prior_hits_input == 0L,
  "no_EAS_EUR_condition_candidate",
  n_raw_prior_hits_invalid > 0L,
  "incomplete_prior_set_invalid_external_input",
  EAS_primary_ready,
  "ready_for_EAS_primary_COJO",
  EUR_only_prior,
  "EUR_only_prior_no_EAS_primary_COJO",
  EAS_prior_present & EAS_n_clump_representatives_unrepresented > 0L,
  "EAS_prior_clump_unrepresented_in_KoGES",
  EAS_prior_present & EAS_n_within_source_representative_collisions > 0L,
  "EAS_prior_within_source_representative_collision",
  EAS_prior_present & EAS_n_condition_snps_selected == 0L,
  "EAS_primary_no_condition_snp_from_clump_or_KoGES_rescue",
  EAS_prior_present & !EAS_grouping_ok,
  "EAS_prior_cross_source_KoGES_LD_grouping_failed",
  EAS_prior_present,
  "EAS_primary_not_ready",
  default = "conditionally_unassessable_step2"
)]

final_manifest[, cojo_condition_snplist_file := NA_character_]
final_manifest[, cojo_extract_snplist_file := NA_character_]
final_manifest[, sensitivity_cojo_condition_snplist_file := NA_character_]
final_manifest[, sensitivity_cojo_extract_snplist_file := NA_character_]
final_manifest[, EAS_primary_cojo_condition_snplist_file := NA_character_]
final_manifest[, EAS_primary_cojo_extract_snplist_file := NA_character_]
final_manifest[, EUR_supporting_cojo_condition_snplist_file := NA_character_]
final_manifest[, EUR_supporting_cojo_extract_snplist_file := NA_character_]
EAS_primary_snplist_dir <- file.path(
  OUT_DIR, "per_target_EAS_primary_cojo_cond_snplists"
)
dir.create(EAS_primary_snplist_dir, recursive = TRUE, showWarnings = FALSE)
EUR_supporting_snplist_dir <- file.path(
  OUT_DIR, "per_target_EUR_supporting_cojo_cond_snplists"
)
dir.create(EUR_supporting_snplist_dir, recursive = TRUE, showWarnings = FALSE)

write_target_snplists <- function(row_indices, id_col, out_dir, suffix) {
  if (length(row_indices) == 0L) return(list())
  parallel_lapply(
    row_indices,
    function(row_index) {
      snps <- strsplit(
        final_manifest[[id_col]][row_index],
        ";", fixed = TRUE
      )[[1L]]
      snps <- unique(snps[nzchar(snps)])
      cond_file <- file.path(
        out_dir,
        paste0(final_manifest$target_id[row_index], suffix, "_cojo_cond.snplist")
      )
      extract_file <- file.path(
        out_dir,
        paste0(final_manifest$target_id[row_index], suffix, "_cojo_extract.snplist")
      )
      write_snplist(snps, cond_file)
      write_snplist(
        unique(c(final_manifest$koges_lead_variant[row_index], snps)),
        extract_file
      )
      list(row_index = row_index, cond = cond_file, extract = extract_file)
    },
    n_workers = N_PARALLEL
  )
}

ready_rows <- which(final_manifest$cojo_ready)
snplist_paths <- write_target_snplists(
  ready_rows, "condition_snp_ids_available", EAS_primary_snplist_dir,
  "_EAS_primary"
)
for (res in snplist_paths) {
  final_manifest$cojo_condition_snplist_file[res$row_index] <- res$cond
  final_manifest$cojo_extract_snplist_file[res$row_index] <- res$extract
  final_manifest$EAS_primary_cojo_condition_snplist_file[res$row_index] <- res$cond
  final_manifest$EAS_primary_cojo_extract_snplist_file[res$row_index] <- res$extract
}

sensitivity_rows <- which(final_manifest$sensitivity_ready)
sensitivity_snplist_paths <- write_target_snplists(
  sensitivity_rows, "condition_snp_ids_sensitivity",
  EUR_supporting_snplist_dir, "_EUR_supporting"
)
for (res in sensitivity_snplist_paths) {
  final_manifest$sensitivity_cojo_condition_snplist_file[res$row_index] <-
    res$cond
  final_manifest$sensitivity_cojo_extract_snplist_file[res$row_index] <-
    res$extract
  final_manifest$EUR_supporting_cojo_condition_snplist_file[res$row_index] <-
    res$cond
  final_manifest$EUR_supporting_cojo_extract_snplist_file[res$row_index] <-
    res$extract
}

setorder(final_manifest, lead_row_id, target_id)
final_manifest_file <- file.path(
  OUT_DIR,
  "joint_COJO_condition_manifest_ancestry_clumped.tsv"
)
fwrite(final_manifest, final_manifest_file, sep = "\t", quote = FALSE, na = "NA")

EAS_primary_jobs <- final_manifest[
  EAS_primary_ready %in% TRUE,
  .(
    analysis_set = "EAS_primary",
    analysis_id = paste(target_id, "EAS_primary", sep = "__"),
    analysis_role = "primary_EAS_condition_set_complete",
    target_id, lead_row_id, koges_phenotype, koges_pheno_code,
    koges_lead_variant, koges_chr, koges_bp,
    cojo_ready = TRUE, EAS_primary_ready, EUR_supporting_ready,
    EAS_prior_present, EUR_prior_present, EUR_only_prior,
    prior_coverage_status, prior_set_coverage_status, prior_identity_status,
    raw_prior_reference_or_exact_KoGES_fraction,
    proxy_transferability_status, cojo_condition_status,
    prior_signal_set_complete, all_source_clump_signal_set_complete,
    EAS_clump_signal_set_complete, EUR_clump_signal_set_complete,
    strict_primary_ready = TRUE,
    primary_gate_failure_reason,
    analysis_gate_failure_reason = primary_gate_failure_reason,
    n_condition_snps_available = EAS_n_condition_snps_selected,
    n_condition_snps_selected = EAS_n_condition_snps_selected,
    n_condition_snps_from_1000G_clumps = EAS_n_condition_snps_from_1000G_clumps,
    n_condition_snps_ancestry_clumped = EAS_n_condition_snps_selected,
    n_condition_snps_selected_external = EAS_n_condition_snps_selected,
    n_EAS_selected = EAS_n_condition_snps_selected,
    n_EUR_selected = 0L,
    n_clump_representatives_exact_variant =
      EAS_n_clump_representatives_exact_variant,
    n_clump_representatives_ld_proxy = EAS_n_clump_representatives_ld_proxy,
    n_clump_representatives_panel_wide_proxy =
      EAS_n_clump_representatives_panel_wide_proxy,
    n_clump_representatives_unrepresented =
      EAS_n_clump_representatives_unrepresented,
    n_condition_snps_direct_KoGES_rescue = EAS_n_condition_snps_direct_KoGES_rescue,
    n_direct_KoGES_rescued_hits_selected = n_direct_KoGES_rescued_hits_selected,
    cojo_condition_snplist_file = EAS_primary_cojo_condition_snplist_file,
    cojo_extract_snplist_file = EAS_primary_cojo_extract_snplist_file
  )
]
EAS_primary_jobs_file <- file.path(OUT_DIR, "joint_COJO_EAS_primary_targets.tsv")
fwrite(EAS_primary_jobs, EAS_primary_jobs_file, sep = "\t", quote = FALSE, na = "NA")
# Backward-compatible default consumed by Step3b when READY_TARGETS_FILE is not
# supplied explicitly.  It is always the EAS-primary job set in V4.
cojo_jobs_file <- file.path(OUT_DIR, "joint_COJO_ready_targets.tsv")
fwrite(EAS_primary_jobs, cojo_jobs_file, sep = "\t", quote = FALSE, na = "NA")

EUR_supporting_jobs <- final_manifest[
  EUR_supporting_ready %in% TRUE,
  .(
    analysis_set = "EUR_supporting",
    analysis_id = paste(target_id, "EUR_supporting", sep = "__"),
    analysis_role = "supporting_EUR_clump_complete",
    target_id, lead_row_id, koges_phenotype, koges_pheno_code,
    koges_lead_variant, koges_chr, koges_bp,
    cojo_ready = FALSE, EAS_primary_ready, EUR_supporting_ready = TRUE,
    EAS_prior_present, EUR_prior_present, EUR_only_prior,
    prior_coverage_status, prior_set_coverage_status, prior_identity_status,
    raw_prior_reference_or_exact_KoGES_fraction,
    proxy_transferability_status, cojo_condition_status,
    prior_signal_set_complete, all_source_clump_signal_set_complete,
    EAS_clump_signal_set_complete, EUR_clump_signal_set_complete,
    strict_primary_ready = FALSE,
    primary_gate_failure_reason,
    analysis_gate_failure_reason = EUR_supporting_gate_failure_reason,
    n_condition_snps_available = EUR_n_condition_snps_selected,
    n_condition_snps_ancestry_clumped = EUR_n_condition_snps_selected,
    n_condition_snps_selected_external = EUR_n_condition_snps_selected,
    n_EAS_selected = 0L,
    n_EUR_selected = EUR_n_condition_snps_selected,
    n_clump_representatives_exact_variant =
      EUR_n_clump_representatives_exact_variant,
    n_clump_representatives_ld_proxy = EUR_n_clump_representatives_ld_proxy,
    n_clump_representatives_panel_wide_proxy =
      EUR_n_clump_representatives_panel_wide_proxy,
    n_clump_representatives_unrepresented =
      EUR_n_clump_representatives_unrepresented,
    n_condition_snps_direct_KoGES_rescue = 0L,
    n_direct_KoGES_rescued_hits_selected = 0L,
    cojo_condition_snplist_file = EUR_supporting_cojo_condition_snplist_file,
    cojo_extract_snplist_file = EUR_supporting_cojo_extract_snplist_file
  )
]
EUR_supporting_jobs_file <- file.path(
  OUT_DIR, "joint_COJO_EUR_supporting_targets.tsv"
)
fwrite(EUR_supporting_jobs, EUR_supporting_jobs_file,
       sep = "\t", quote = FALSE, na = "NA")
# Compatibility alias: this is explicitly EUR supporting, not a reduced
# all-source primary model.
sensitivity_jobs_file <- file.path(OUT_DIR, "joint_COJO_sensitivity_targets.tsv")
fwrite(EUR_supporting_jobs, sensitivity_jobs_file,
       sep = "\t", quote = FALSE, na = "NA")

all_source_complete_audit <- final_manifest[
  step1_cojo_condition_status != "known_exact_variant",
  .(
    target_id, lead_row_id, all_source_clump_signal_set_complete,
    EAS_clump_signal_set_complete, EUR_clump_signal_set_complete,
    EAS_prior_present, EUR_prior_present, EUR_only_prior,
    n_external_clump_signals, n_clump_representatives_unrepresented,
    EAS_n_external_clump_signals, EAS_n_clump_representatives_unrepresented,
    EUR_n_external_clump_signals, EUR_n_clump_representatives_unrepresented,
    raw_prior_reference_or_exact_KoGES_fraction
  )
]
all_source_complete_audit_file <- file.path(
  OUT_DIR, "joint_COJO_ancestry_stratified_clump_completeness_audit.tsv"
)
fwrite(all_source_complete_audit, all_source_complete_audit_file,
       sep = "\t", quote = FALSE, na = "NA")

direct_koges_rescue_audit_file <- file.path(
  OUT_DIR, "joint_COJO_direct_KoGES_exact_rescue_audit.tsv"
)
if (is.null(direct_koges_rescue_union)) {
  fwrite(data.table(), direct_koges_rescue_audit_file,
         sep = "\t", quote = FALSE, na = "NA")
} else {
  fwrite(direct_koges_rescue_union, direct_koges_rescue_audit_file,
         sep = "\t", quote = FALSE, na = "NA")
}

# ------------------------------------------------------------------------------
# 8. Run README
# ------------------------------------------------------------------------------

readme_file <- file.path(OUT_DIR, "README_joint_COJO_1000G_ancestry_clumping.txt")
writeLines(c(
  "Targeted external-prior signal definition for KoGES GCTA-COJO",
  "",
  "Design:",
  "- Input is Step1's source-level raw external GWS audit, not a P-value pool across cohorts.",
  "- Direct KoGES mapping is NOT an upstream eligibility filter.",
  "- clump_snp is resolved in the ancestry-matched 1000G panel at allele level,",
  "  because the clump index is what defines the external signal:",
  "    1) the original external rsID only when its raw CHR:BP+alleles match the panel;",
  "    2) otherwise the UNIQUE panel variant at the raw external CHR:BP whose",
  "       unordered allele pair matches external effect/other allele;",
  "    3) otherwise the Step1 harmonized KoGES rsID, but only when that rsID is",
  "       verified to be the same variant as the raw external hit in BOTH KoGES",
  "       and the panel;",
  "    4) otherwise not_represented_in_ancestry_1000G.",
  "  A same-coordinate SNP-indel, a different ALT, or a coordinate-only rsID",
  "  alias is never substituted for the external index.",
  "- A valid hit absent from its ancestry 1000G panel is not clumped. If its raw",
  "  CHR:BP+alleles identify exactly one KoGES BIM variant and",
  paste0("  DIRECT_KOGES_RESCUE_IN_CORE=", DIRECT_KOGES_RESCUE_IN_CORE, ","),
  "  that KoGES SNP is added to the EAS-primary condition set.  The rescue",
  "  audit table is the row-level provenance of those additions, not a",
  "  substitute for the core covariates.  Rescue SNPs have no 1000G-defined",
  "  within-source independence, so KoGES r^2 grouping may fold same-source",
  "  rescue-rescue and rescue-clump pairs at the same cutoff used for",
  "  cross-source clump representatives.  Two 1000G clump representatives",
  "  from overlapping sources still never fold.",
  "- 1000G clumpability and KoGES direct conditioning are separate.  Every",
  "  1000G-representable prior is clumped, including same-coordinate SNP-indel",
  "  and shared-allele (multiallelic) encodings.  Those are recorded as",
  "  representation_discordant, not as validated identity.  BIM CHR:BP alone",
  "  does not prove the same biological allele.  HLA is not a special-case",
  "  exception.",
  "- After clumping, each clump is exact_variant (validated allele identity in",
  "  KoGES), ld_proxy_r2_ge_threshold, or unrepresented.",
  paste0(
    "- Proxy search stage 1: in-clump members only, 1000G r^2 >= ", PROXY_MIN_R2,
    ", and the member must itself be an exact KoGES variant."
  ),
  paste0(
    "- Proxy search stage 2 (PANEL_WIDE_PROXY_RESCUE=", PANEL_WIDE_PROXY_RESCUE,
    "): an external GWAS reports only its own variants, so stage 1 can fail"
  ),
  "  merely because the LD partner was never submitted.  Stage 2 asks PLINK for",
  paste0(
    "  every ancestry-panel variant within ", PANEL_WIDE_PROXY_WINDOW_KB,
    " kb of the clump index at r^2 >= ", PROXY_MIN_R2, ", keeps only variants"
  ),
  "  with exact 1000G-to-KoGES identity inside the target window, excludes the",
  "  target lead, and takes the highest r^2 (ties by distance to lead, then ID).",
  "  The r^2 threshold is never relaxed; only the candidate pool widens.",
  "  The identity pool is matched on CHR:BP and allele pair, not on rsID, so a",
  "  panel variant carrying a different rsID than KoGES is still usable; only",
  "  non-1:1 position matches are dropped.",
  paste0(
    "  The LD window is ", PANEL_WIDE_PROXY_WINDOW_KB,
    " kb so an index and a proxy at opposite edges of the target window are",
    " still evaluated."
  ),
  "  Stage-2 LD is computed once per (target, ancestry panel) with",
  "  --ld-snp-list over all of that target's clump indices, instead of once per",
  "  clump.  The window, --ld-window-kb and --ld-window-r2 are unchanged, so the",
  "  retained pairs are identical; only the number of PLINK jobs falls (roughly",
  "  one per target-panel rather than one per unrepresented clump).  A missing",
  "  group record falls back to the original per-clump call.",
  "  Stage-1 member LD also loads only the target window; Step1 already limited",
  "  every source input to the lead +/-1 Mb, so no reported pair changes.",
  "- A proxy selected in the EUR panel proves LD in Europeans only.  The same",
  paste0(
  "  index-proxy pair can optionally be re-estimated in 1000G EAS and is recorded as a",
  "  transferability-evidence grade (not an EAS-primary execution gate). An unconfirmed"
  ),
  "  EUR-only proxy is used only in the separately labelled EUR-supporting",
  "  condition set, never to determine the EAS-primary conclusion.",
  "- An unrepresented EAS clump makes the EAS-primary analysis unassessable.",
  "  An unrepresented EUR clump makes the EUR-supporting analysis unavailable,",
  "  but does not invalidate an otherwise complete EAS-primary model.",
  "- In the observed data, unrepresented EUR indices are almost entirely true",
  "  KoGES-panel absences rather than build/rsID/indel mapping errors. Do not",
  "  present liftover or left-normalization as the expected rescue route; the",
  "  panel-wide proxy audit records the relevant available rescue path.",
  "- Each source is clumped separately: BBJ/TWB/TPMI with 1000G EAS; UKB/FinnGen/MVP with 1000G EUR.",
  "- PLINK index SNPs and all SP2 clump members are retained in an audit map.",
  "  That map is the largest table in the run and is written gzipped; a resume",
  "  reads either the .tsv.gz or an older uncompressed .tsv.",
  paste0(
    "- After union, pairwise KoGES r^2 >= ",
    CROSS_SOURCE_REP_R2,
    " can collapse two condition SNPs to one representative."
  ),
  "  Two 1000G clump representatives collapse only when selected_source_keys",
  "  are disjoint (shared sources never fold those 1000G-defined signals).",
  "  Direct-KoGES rescue SNPs may fold at the same r^2 even when sources",
  "  overlap, including rescue-rescue and rescue-clump pairs.",
  "  Within a group, a candidate still needs direct r^2 to the chosen",
  "  representative (no transitive graph-only collapse).",
  "  Representative priority: 1000G-defined clump SNP first (direct clump",
  "  index, then proxy), then direct KoGES rescue.  Rescue never outranks a",
  "  collinear 1000G clump representative.  Among rescue-only SNPs, source",
  "  support then external P is the tie-breaker.",
  "  The member-to-representative pairwise r^2, both source sets, and",
  "  pair_eligible_for_ld_collapse are retained in the audit files.",
  "- KoGES LD is not used for source-specific 1000G clumping.  It is used for",
  "  this redundancy collapse and for the final GCTA --cojo-cond model in Step3b.",
  "",
  "Scope:",
  "- CLUMP_SCOPE=target_local_source_stratified.",
  "- This is a targeted +/-1 Mb conditional analysis, not a genome-wide independent-signal",
  "  analysis of the complete external GWAS. Do not claim source-wide GWAS clumping unless",
  "  source-wide external clump-index/member files were generated separately.",
  "",
  "Coverage vs ancestry-stratified COJO readiness:",
  "- Invalid P / missing variant IDs remain input failures.",
  "- Raw-hit 1000G availability and raw_prior_reference_or_exact_KoGES_fraction",
  "  are continuous QC metrics. They are not the binary primary gate.",
  "- EAS_primary_ready=TRUE when EAS external GWS prior evidence exists",
  "  (valid raw hits), every EAS 1000G clump that exists is represented in",
  "  KoGES, and at least one usable EAS condition SNP is selected (1000G",
  "  clump representative and/or direct KoGES-exact rescue) with successful",
  "  grouping.  A 1000G clump is not required.  Unresolved raw hits",
  "  (prior_set_coverage_status=partial_unresolved_prior) do not fail the",
  "  gate: raw GWS rows are not independent signals.  Do not write that all",
  "  reported EAS variants were conditioned.  Write that all representable",
  "  EAS prior signals defined by the conditioning workflow were accounted for.",
  "  Show partial_unresolved_prior as a separate QC slice.",
  "- EUR_supporting_ready is evaluated separately with EUR-defined clumps. It",
  "  supplies supporting evidence only and never changes the EAS-primary label.",
  "- EUR_only_prior is a valid, explicit category: it has no EAS-primary test",
  "  because the positional prior is EUR-only, not because the pipeline failed.",
  "- all_source_clump_signal_set_complete is retained as an audit of the small",
  "  subset with complete EAS+EUR clump representation; it is not the primary",
  "  Korean-cohort gate.",
  "- Same-coordinate SNP-indel / shared-allele encodings stay in clumping as",
  "  representation_discordant (prior_identity_status). They are not treated as",
  "  the external known variant without exact identity or a high-LD proxy.",
  "",
  "Parameters:",
  paste0("- RAW_HITS_FILE: ", RAW_HITS_FILE),
  paste0("- PRELD_MANIFEST_FILE: ", PRELD_MANIFEST_FILE),
  paste0("- KOGES_BFILE: ", KOGES_BFILE),
  paste0("- EAS_BFILE: ", EAS_BFILE),
  paste0("- EUR_BFILE: ", EUR_BFILE),
  paste0("- PLINK_BIN: ", PLINK_BIN),
  paste0("- clump-r2: ", CLUMP_R2),
  paste0("- clump-kb: ", CLUMP_KB),
  paste0("- clump-p1: ", CLUMP_P1),
  paste0("- external-GWS input P cutoff (CLUMP_P2): ", CLUMP_P2),
  paste0("- PLINK SP2 member-output P cutoff (CLUMP_MEMBER_P2): ", CLUMP_MEMBER_P2),
  paste0("- proxy minimum 1000G r^2 (both stages): ", PROXY_MIN_R2),
  paste0("- DIRECT_KOGES_RESCUE: ", DIRECT_KOGES_RESCUE),
  paste0("- DIRECT_KOGES_RESCUE_IN_CORE: ", DIRECT_KOGES_RESCUE_IN_CORE),
  paste0("- EXACT_1000G_CLUMP_ONLY: ", EXACT_1000G_CLUMP_ONLY),
  paste0("- ANALYSIS_LAYER: ", ANALYSIS_LAYER),
  paste0("- VERIFY_EUR_PROXY_IN_EAS (optional QC only): ", VERIFY_EUR_PROXY_IN_EAS),
  paste0("- PANEL_WIDE_PROXY_RESCUE: ", PANEL_WIDE_PROXY_RESCUE),
  paste0("- PANEL_WIDE_PROXY_WINDOW_KB: ", PANEL_WIDE_PROXY_WINDOW_KB),
  paste0(
    "- PANEL_WIDE_PROXY_GROUP_LD: ", PANEL_WIDE_PROXY_GROUP_LD,
    " (FALSE = one PLINK job per clump; comparison only)"
  ),
  paste0(
    "- target subset: ",
    if (length(resolve_target_subset(TARGET_IDS_ENV, TARGET_ID_FILE)) > 0L) {
      paste0(
        nrow(final_manifest),
        " target(s) via TARGET_IDS / TARGET_ID_FILE -- PILOT OUTPUT, ",
        "not a full run"
      )
    } else {
      "none (all Step1 targets)"
    }
  ),
  paste0("- N_PARALLEL (1000G clump + snplist): ", N_PARALLEL),
  paste0("- N_REP_PARALLEL (representative batches): ", N_REP_PARALLEL),
  paste0("- N_MEM_PARALLEL (KoGES cross-source LD): ", N_MEM_PARALLEL),
  paste0("- N_DT_THREADS (prep OpenMP): ", N_DT_THREADS),
  paste0("- REP_BATCH_SIZE: ", REP_BATCH_SIZE),
  "",
  "Key outputs:",
  paste0("- source-level preclump candidates: ", source_candidates_file),
  paste0("- raw prior input QC: ", raw_input_qc_file),
  paste0("- 1000G reference QC: ", panel_reference_qc_file),
  paste0("- source-level clump job summary: ", clump_summary_file),
  paste0("- clump index/member audit: ", clump_member_audit_file),
  paste0("- clump representative audit: ", clump_representative_file),
  paste0("- member-proxy LD audit: ", member_proxy_audit_file),
  paste0("- panel-wide proxy LD audit: ", panel_wide_proxy_audit_file),
  paste0("- within-source representative collisions: ", representative_collision_file),
  paste0("- pre-cross-source-LD candidate union: ", final_candidates_pre_cross_source_ld_file),
  paste0("- cross-source KoGES-LD group members: ", cross_source_ld_member_audit_file),
  paste0("- cross-source KoGES-LD pair audit: ", cross_source_ld_pair_audit_file),
  paste0("- cross-source KoGES-LD grouping summary: ", cross_source_ld_summary_file),
  paste0("- selected condition candidates: ", final_candidates_file),
  "- ancestry-stratified candidate/grouping files: joint_COJO_EAS_primary_* and",
  "  joint_COJO_EUR_supporting_* under OUT_DIR",
  paste0("- ancestry-stratified completeness audit: ", all_source_complete_audit_file),
  paste0("- exact direct-KoGES rescue provenance (core covariates when ",
         "DIRECT_KOGES_RESCUE_IN_CORE=", DIRECT_KOGES_RESCUE_IN_CORE, "): ",
         direct_koges_rescue_audit_file),
  paste0("- final manifest: ", final_manifest_file),
  paste0("- EAS-primary targets (default Step3b input): ", EAS_primary_jobs_file),
  paste0("- EUR-supporting targets: ", EUR_supporting_jobs_file),
  "",
  "Step3b should use the full phenotype .ma as --cojo-file and the generated",
  "per-target lead+condition --extract and --cojo-cond snplists."
), readme_file)

cat("Source-level valid condition candidates:", nrow(source_candidates), "\n")
cat("All-source audit candidate pairs:", nrow(final_candidates), "\n")
cat("EAS-primary candidate pairs:", nrow(EAS_core_candidates), "\n")
cat("EUR-supporting candidate pairs:", nrow(EUR_core_candidates), "\n")
cat("EAS-primary ready targets:", nrow(EAS_primary_jobs), "\n")
cat("EUR-supporting ready targets:", nrow(EUR_supporting_jobs), "\n")
cat("EUR-only prior targets:", sum(final_manifest$EUR_only_prior %in% TRUE), "\n")
cat(
  "All-source clump-complete audit targets:",
  sum(final_manifest$all_source_clump_signal_set_complete %in% TRUE), "\n"
)
cat(
  "EAS clump representation: exact_variant",
  sum(final_manifest$EAS_n_clump_representatives_exact_variant, na.rm = TRUE),
  "/ LD proxy",
  sum(final_manifest$EAS_n_clump_representatives_ld_proxy, na.rm = TRUE),
  "(panel-wide rescue",
  sum(final_manifest$EAS_n_clump_representatives_panel_wide_proxy, na.rm = TRUE),
  ") / unrepresented",
  sum(final_manifest$EAS_n_clump_representatives_unrepresented, na.rm = TRUE),
  "\n"
)
cat(
  "EUR clump representation: exact_variant",
  sum(final_manifest$EUR_n_clump_representatives_exact_variant, na.rm = TRUE),
  "/ LD proxy",
  sum(final_manifest$EUR_n_clump_representatives_ld_proxy, na.rm = TRUE),
  "(panel-wide rescue",
  sum(final_manifest$EUR_n_clump_representatives_panel_wide_proxy, na.rm = TRUE),
  ") / unrepresented",
  sum(final_manifest$EUR_n_clump_representatives_unrepresented, na.rm = TRUE),
  "\n"
)
cat("Final manifest:", final_manifest_file, "\n")
cat("EAS-primary targets:", EAS_primary_jobs_file, "\n")
cat("EUR-supporting targets:", EUR_supporting_jobs_file, "\n")
