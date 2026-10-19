#!/usr/bin/env Rscript
# ==============================================================================
# build_ExtendedData_conditionally_independent_172.R
# ------------------------------------------------------------------------------
# Build Extended Data for the 172 KoGES lead associations that remain
# conditionally independent after EAS COJO conditioning on all candidate
# prior variants (rescue-augmented sensitivity mode).
#
# Denominator (same rule as EUR-only Step0 / README):
#   robust_additional_independent_association == TRUE
#   AND step2_n_raw_prior_hits_unresolved_after_direct_KoGES_rescue == 0
#   (n must be 172)
#
# Conditional beta / SE / P come from the EAS sensitivity COJO run
# (1000G EAS clump-defined signals + direct KoGES-exact rescue), i.e. after
# conditioning on the full candidate set used in dual-mode sensitivity.
#
# Inputs (defaults under BASE_DIR):
#   COJO/EAS-only/GCTA_conditional_V5_EAS_dualmode_comparison.tsv
#   COJO/EAS-only/GCTA_conditional_V5_EAS_sensitivity_rescue/
#     joint_COJO_conditional_results_integrated.tsv
#   DATA_manuscript_final/Extended_Data/Cont_lead_variants.tsv
# Optional (preferred if present; same 172 IDs):
#   COJO/EUR-only/strict_novel_172/strict_novel_172_EUR_prior_inventory.tsv
#
# Outputs:
#   RESULT_manuscript_final/ExtendedData/
#     ExtendedData_conditionally_independent_172.csv
#     ExtendedData_conditionally_independent_172.tsv
#     ExtendedData_conditionally_independent_172_column_dictionary.tsv
#
# Usage:
#   Rscript --vanilla build_ExtendedData_conditionally_independent_172.R
#   BASE_DIR=/path OUT_DIR=/path Rscript --vanilla build_ExtendedData_...R
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

BASE_DIR <- Sys.getenv(
  "BASE_DIR",
  "/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant"
)
COMPARISON_FILE <- Sys.getenv(
  "COMPARISON_FILE",
  file.path(
    BASE_DIR, "COJO", "EAS-only",
    "GCTA_conditional_V5_EAS_dualmode_comparison.tsv"
  )
)
SENSITIVITY_RESULTS <- Sys.getenv(
  "SENSITIVITY_RESULTS",
  file.path(
    BASE_DIR, "COJO", "EAS-only",
    "GCTA_conditional_V5_EAS_sensitivity_rescue",
    "joint_COJO_conditional_results_integrated.tsv"
  )
)
INVENTORY_FILE <- Sys.getenv(
  "INVENTORY_FILE",
  file.path(
    BASE_DIR, "COJO", "EUR-only", "strict_novel_172",
    "strict_novel_172_EUR_prior_inventory.tsv"
  )
)
LEAD_ANNOT_FILE <- Sys.getenv(
  "LEAD_ANNOT_FILE",
  "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/Cont_lead_variants.tsv"
)
OUT_DIR <- Sys.getenv(
  "OUT_DIR",
  "/home2/sjim/openKoGES/RESULT_manuscript_final/ExtendedData"
)

STATUS_INDEPENDENT <- "conditionally_independent_given_EAS_prior_signals"
EXPECTED_N <- 172L

as_bool <- function(x) {
  if (is.null(x)) return(logical())
  if (is.logical(x)) return(!is.na(x) & x)
  if (is.numeric(x)) return(is.finite(x) & x != 0)
  toupper(trimws(as.character(x))) %chin% c("TRUE", "T", "1", "YES", "Y")
}

for (path in c(COMPARISON_FILE, SENSITIVITY_RESULTS, LEAD_ANNOT_FILE)) {
  if (!file.exists(path)) stop("Missing required file: ", path)
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Resolve the 172 target IDs
# ------------------------------------------------------------------------------
cmp <- fread(COMPARISON_FILE, sep = "\t", showProgress = FALSE)
needed_cmp <- c(
  "target_id",
  "robust_additional_independent_association",
  "step2_n_raw_prior_hits_unresolved_after_direct_KoGES_rescue",
  "dualmode_robustness_status",
  "final_cojo_status",
  "sensitivity_final_cojo_status"
)
missing_cmp <- setdiff(needed_cmp, names(cmp))
if (length(missing_cmp) > 0L) {
  stop("Comparison table missing column(s): ", paste(missing_cmp, collapse = ", "))
}

cmp[, target_id := as.character(target_id)]
cmp[, robust_additional_independent_association :=
      as_bool(robust_additional_independent_association)]
cmp[, n_unresolved := as.integer(
  step2_n_raw_prior_hits_unresolved_after_direct_KoGES_rescue
)]
cmp[is.na(n_unresolved), n_unresolved := 0L]

strict_from_cmp <- cmp[
  robust_additional_independent_association == TRUE & n_unresolved == 0L
]
if (nrow(strict_from_cmp) != EXPECTED_N) {
  stop(
    "Expected ", EXPECTED_N, " strictly novel leads from comparison; found ",
    nrow(strict_from_cmp), "."
  )
}

if (file.exists(INVENTORY_FILE)) {
  inv <- fread(INVENTORY_FILE, select = "target_id", showProgress = FALSE)
  inv[, target_id := as.character(target_id)]
  if (nrow(inv) != EXPECTED_N) {
    stop("INVENTORY_FILE does not contain ", EXPECTED_N, " rows: ", nrow(inv))
  }
  if (!setequal(inv$target_id, strict_from_cmp$target_id)) {
    stop(
      "Inventory target_id set differs from comparison-derived 172 set. ",
      "missing_in_inventory=",
      length(setdiff(strict_from_cmp$target_id, inv$target_id)),
      "; extra_in_inventory=",
      length(setdiff(inv$target_id, strict_from_cmp$target_id))
    )
  }
  target_ids <- inv$target_id
  cat("Using inventory IDs (validated against comparison): ", INVENTORY_FILE, "\n", sep = "")
} else {
  target_ids <- strict_from_cmp$target_id
  cat("Inventory not found; using comparison-derived 172 IDs.\n")
}

# ------------------------------------------------------------------------------
# 2. Sensitivity COJO results (conditional on all candidates)
# ------------------------------------------------------------------------------
sens <- fread(SENSITIVITY_RESULTS, sep = "\t", showProgress = FALSE)
sens[, target_id := as.character(target_id)]
sens_sub <- sens[target_id %chin% target_ids]
if (nrow(sens_sub) != EXPECTED_N) {
  stop(
    "Sensitivity integrated results cover ", nrow(sens_sub),
    " / ", EXPECTED_N, " target IDs."
  )
}
if (anyDuplicated(sens_sub$target_id)) {
  stop("Duplicate target_id in sensitivity results for the 172 set.")
}
if (!all(sens_sub$final_cojo_status == STATUS_INDEPENDENT)) {
  bad <- sens_sub[final_cojo_status != STATUS_INDEPENDENT, .N, by = final_cojo_status]
  stop(
    "Not all 172 sensitivity rows are ", STATUS_INDEPENDENT, ":\n",
    paste(capture.output(print(bad)), collapse = "\n")
  )
}
if (!all(is.finite(sens_sub$conditional_beta) & is.finite(sens_sub$conditional_p))) {
  stop("Non-finite conditional_beta / conditional_p in sensitivity results for the 172.")
}

# ------------------------------------------------------------------------------
# 3. Lead annotation (nearest gene, alleles, variant ID)
# ------------------------------------------------------------------------------
lead <- fread(LEAD_ANNOT_FILE, sep = "\t", showProgress = FALSE)
lead_keep <- c(
  "Phenotype", "Lead_variant", "VariantID", "CHR", "BP", "Variant_type",
  "Nearest_gene", "NEA", "EA", "REF", "ALT", "MAF", "EAF", "N", "P_field"
)
missing_lead <- setdiff(lead_keep, names(lead))
if (length(missing_lead) > 0L) {
  stop("Lead annotation missing column(s): ", paste(missing_lead, collapse = ", "))
}
lead_anno <- unique(lead[, ..lead_keep])
if (anyDuplicated(lead_anno, by = c("Phenotype", "Lead_variant"))) {
  stop("Lead annotation is not unique on Phenotype + Lead_variant.")
}

# ------------------------------------------------------------------------------
# 4. Assemble Extended Data table
# ------------------------------------------------------------------------------
meta <- strict_from_cmp[
  ,
  .(
    target_id,
    dualmode_robustness_status,
    primary_final_cojo_status = final_cojo_status,
    sensitivity_final_cojo_status
  )
]

dt <- merge(sens_sub, meta, by = "target_id", all.x = FALSE, sort = FALSE)
dt <- merge(
  dt,
  lead_anno,
  by.x = c("koges_phenotype", "koges_lead_variant"),
  by.y = c("Phenotype", "Lead_variant"),
  all.x = TRUE,
  sort = FALSE
)
if (any(is.na(dt$Nearest_gene) | is.na(dt$VariantID))) {
  stop(
    "Failed to annotate nearest gene / VariantID for ",
    sum(is.na(dt$Nearest_gene) | is.na(dt$VariantID)),
    " / ", nrow(dt), " rows."
  )
}

# Coordinate QC vs Cont lead table
coord_mismatch <- dt[
  !is.na(CHR) & !is.na(BP) &
    (as.integer(CHR) != as.integer(koges_chr) | as.integer(BP) != as.integer(koges_bp))
]
if (nrow(coord_mismatch) > 0L) {
  stop(
    "CHR/BP mismatch between COJO results and Cont lead annotation for ",
    nrow(coord_mismatch), " row(s). Examples: ",
    paste(head(coord_mismatch$target_id, 10L), collapse = ", ")
  )
}

setorder(dt, koges_phenotype, koges_chr, koges_bp, koges_lead_variant)

out <- data.table(
  `Phenotype` = dt$koges_phenotype,
  `Lead variant (rsID)` = dt$koges_lead_variant,
  `Variant ID` = dt$VariantID,
  `Chromosome` = as.integer(dt$koges_chr),
  `Base pair position (GRCh37)` = as.integer(dt$koges_bp),
  `Variant type` = dt$Variant_type,
  `Nearest gene` = dt$Nearest_gene,
  `Effect allele` = dt$EA,
  `Other allele` = dt$NEA,
  `REF` = dt$REF,
  `ALT` = dt$ALT,
  `EAF` = dt$EAF,
  `MAF` = dt$MAF,
  `N` = as.integer(dt$N),
  `Marginal beta` = dt$marginal_beta,
  `Marginal SE` = dt$marginal_se,
  `Marginal P` = dt$marginal_p,
  `Conditional beta` = dt$conditional_beta,
  `Conditional SE` = dt$conditional_se,
  `Conditional P` = dt$conditional_p,
  `n condition SNPs requested` = as.integer(dt$n_condition_SNP_requested),
  `n condition SNPs retained by GCTA` = as.integer(dt$n_condition_SNP_retained_by_GCTA),
  `n condition SNPs from 1000G clumps` = as.integer(dt$step2_n_condition_snps_from_1000G_clumps),
  `n condition SNPs from direct KoGES rescue` = as.integer(
    dt$step2_n_condition_snps_direct_KoGES_rescue
  ),
  `Condition SNPs retained by GCTA` = dt$condition_SNP_retained_by_GCTA_ids,
  `Condition set effectively complete` = as_bool(
    dt$condition_set_effectively_complete_GCTA
  ),
  `Final COJO status (sensitivity)` = dt$final_cojo_status,
  `Dual-mode robustness status` = dt$dualmode_robustness_status,
  `target_id` = dt$target_id,
  `lead_row_id` = as.integer(dt$lead_row_id),
  `Phenotype file stem` = dt$koges_pheno_code
)

if (nrow(out) != EXPECTED_N) {
  stop("Output row count is ", nrow(out), ", expected ", EXPECTED_N, ".")
}
if (anyDuplicated(out$target_id) || anyDuplicated(out[, .(`Phenotype`, `Lead variant (rsID)`)])) {
  stop("Output contains duplicate target / phenotype-lead rows.")
}

dictionary <- data.table(
  column = names(out),
  description = c(
    "KoGES continuous phenotype display name",
    "KoGES lead variant rsID",
    "Variant ID as chr:pos:ref:alt (GRCh37)",
    "Chromosome (GRCh37)",
    "Base-pair position (GRCh37)",
    "Variant class (SNV / INDEL / MNV)",
    "Nearest gene from ANNOVAR annotation of the lead variant",
    "Effect allele (EA) used for beta sign",
    "Other / non-effect allele (NEA)",
    "Genomic reference allele",
    "Genomic alternate allele",
    "Effect-allele frequency in KoGES",
    "Minor-allele frequency in KoGES",
    "GWAS sample size for the phenotype",
    "Marginal (unconditional) BOLT-LMM beta",
    "Marginal (unconditional) BOLT-LMM SE",
    "Marginal (unconditional) association P",
    "GCTA-COJO conditional beta after conditioning on all EAS prior candidate variants (sensitivity: 1000G clump + direct KoGES-exact rescue)",
    "GCTA-COJO conditional SE",
    "GCTA-COJO conditional P (genome-wide significance threshold 5e-8)",
    "Number of condition SNPs requested for GCTA --cojo-cond",
    "Number of condition SNPs retained in the GCTA conditional model",
    "Number of condition SNPs from 1000G EAS clump-defined signals",
    "Number of condition SNPs from direct KoGES-exact rescue",
    "Semicolon-separated rsIDs of condition SNPs retained by GCTA",
    "Whether the GCTA condition set was effectively complete",
    "Paper-facing COJO status from the EAS sensitivity run",
    "Dual-mode robustness label from primary vs sensitivity comparison",
    "Internal COJO target identifier",
    "Internal lead row index",
    "BOLT result file stem / phenotype code"
  )
)

out_csv <- file.path(OUT_DIR, "ExtendedData_conditionally_independent_172.csv")
out_tsv <- file.path(OUT_DIR, "ExtendedData_conditionally_independent_172.tsv")
out_dict <- file.path(
  OUT_DIR, "ExtendedData_conditionally_independent_172_column_dictionary.tsv"
)

fwrite(out, out_csv, sep = ",", quote = TRUE, na = "")
fwrite(out, out_tsv, sep = "\t", quote = FALSE, na = "NA")
fwrite(dictionary, out_dict, sep = "\t", quote = FALSE, na = "NA")

cat("=== Extended Data: 172 conditionally independent associations ===\n")
cat("Comparison: ", COMPARISON_FILE, "\n", sep = "")
cat("Sensitivity COJO: ", SENSITIVITY_RESULTS, "\n", sep = "")
cat("Lead annotation: ", LEAD_ANNOT_FILE, "\n", sep = "")
cat("n rows: ", nrow(out), "\n", sep = "")
cat(
  "Phenotypes: ", uniqueN(out$Phenotype),
  " | Conditional P < 5e-8: ",
  sum(out[["Conditional P"]] < 5e-8, na.rm = TRUE),
  "\n",
  sep = ""
)
cat(
  "Conditional beta range: [",
  format(min(out[["Conditional beta"]]), digits = 4), ", ",
  format(max(out[["Conditional beta"]]), digits = 4), "]\n",
  sep = ""
)
cat("Wrote: ", out_csv, "\n", sep = "")
cat("Wrote: ", out_tsv, "\n", sep = "")
cat("Wrote: ", out_dict, "\n", sep = "")
