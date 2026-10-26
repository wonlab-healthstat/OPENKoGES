#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   Build final merged Cont Supplementary Data 3 workbook from novelty TSVs,
#   panel-meta TSVs, and GWAS Catalog novelty TSV — in memory (no intermediate
#   Supple_Data3_6_cohort_detailed.xlsx or Extended_Data_external_replication.xlsx).
#
# Run order
#   1. external_validation (EAS/EUR Cont)
#   2. panel meta scripts (EAS/EUR Cont)
#   3. Cont_lead_variants_GWAScatalog_novel.tsv available
#   4. Rscript ExtendedData2_cont.R
#
# Inputs
#   - KoGES_BOLT_leads_{EAS,EUR}_panel_meta.tsv (panel meta; Extended_Data)
#   - Cont_lead_variants_GWAScatalog_novel.tsv (catalog novelty; Extended_Data)
#   - EAS Cont novelty TSVs under DATA_manuscript_final/Extended_Data/
#       KoGES_BOLT_leads_{BBJ,TWB,TPMI}_novelty_replication_cont.tsv
#   - EUR Cont novelty TSVs under DATA_manuscript_final/Extended_Data/
#       KoGES_BOLT_leads_{UKB,FinnGen,MVP}_novelty_replication_cont.tsv
#
# Outputs
#   - DATA_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx
#     (multi-sheet analysis workbook used by Figures / SuppleTables)
#   - Publication sheets (run after this script):
#       python3 CODE_manuscript_final/ExtendedData_export/export_ExtendedData_publication.py
#     → RESULT_manuscript_final/ExtendedData/Extended Data2A.xlsx
#     → RESULT_manuscript_final/ExtendedData/Extended Data4A.xlsx / Extended Data4B.xlsx
#   - manuscript summary CSV is built by SuppleTable7.R → RESULT_manuscript_final/SuppleTable7/SuppleTable7_cont.csv
#
# IMPORTANT
#   - One row = one KoGES phenotype--lead association (lead_row_id).
#   - Excel sheet names are capped at 31 characters.
#   - Lives with Cate under forEAS_EUR_merge/ (shared helpers in
#     ExtendedData2_common.R; builders in ExtendedData2_builders_cont.R).
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})

# ------------------------------------------------------------------------------
# 0. Paths + builders
# ------------------------------------------------------------------------------
catalog_tsv <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data/Cont_lead_variants_GWAScatalog_novel.tsv"
output_xlsx <- "/home2/sjim/openKoGES/DATA_manuscript_final/ExtendedData/ExtendedData2_continuous.xlsx"

ext_data_dir <- Sys.getenv("EXT_DATA_DIR", unset = "")
if (!nzchar(ext_data_dir)) {
  ext_data_dir <- Sys.getenv(
    "EXTENDED_DATA_DIR",
    "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
  )
}
eas_panel_tsv <- Sys.getenv(
  "EAS_PANEL_META_PATH",
  file.path(ext_data_dir, "KoGES_BOLT_leads_EAS_panel_meta.tsv")
)
eur_panel_tsv <- Sys.getenv(
  "EUR_PANEL_META_PATH",
  file.path(ext_data_dir, "KoGES_BOLT_leads_EUR_panel_meta.tsv")
)

merge_dir <- "/home2/sjim/openKoGES/CODE_manuscript_final/4_Novelty_Replication/forEAS_EUR_merge"
code_root <- dirname(merge_dir)
common_r <- file.path(merge_dir, "ExtendedData2_common.R")
builders_r <- file.path(merge_dir, "ExtendedData2_builders_cont.R")
panel_repl_r <- file.path(code_root, "panel_replication_from_tsv.R")
if (!file.exists(common_r)) stop("Common helpers not found: ", common_r)
if (!file.exists(builders_r)) stop("Builder library not found: ", builders_r)
if (!file.exists(panel_repl_r)) stop("Panel replication library not found: ", panel_repl_r)
source(common_r)
source(builders_r)
source(panel_repl_r)

eas_novelty_files <- c(
  BBJ = file.path(ext_data_dir, "KoGES_BOLT_leads_BBJ_novelty_replication_cont.tsv"),
  TWB = file.path(ext_data_dir, "KoGES_BOLT_leads_TWB_novelty_replication_cont.tsv"),
  TPMI = file.path(ext_data_dir, "KoGES_BOLT_leads_TPMI_novelty_replication_cont.tsv")
)
eur_novelty_files <- c(
  UKB = file.path(ext_data_dir, "KoGES_BOLT_leads_UKB_novelty_replication_cont.tsv"),
  FinnGen = file.path(ext_data_dir, "KoGES_BOLT_leads_FinnGen_novelty_replication_cont.tsv"),
  MVP = file.path(ext_data_dir, "KoGES_BOLT_leads_MVP_novelty_replication_cont.tsv")
)

required_inputs <- c(
  eas_panel_tsv, eur_panel_tsv, catalog_tsv,
  unname(eas_novelty_files), unname(eur_novelty_files)
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs) > 0L) {
  stop("Input file not found:\n  ", paste(missing_inputs, collapse = "\n  "))
}

# ------------------------------------------------------------------------------
# 2. Build panel replication tables from panel-meta TSVs (in memory)
# ------------------------------------------------------------------------------
message("Building Cont panel replication tables from TSV ...")
panel_tabs <- build_panel_replication_tables(eas_panel_tsv, eur_panel_tsv)
eas_meta <- panel_tabs$eas
eur_meta <- panel_tabs$eur

meta_required <- c(
  "lead_row_id", "koges_pheno_code", "koges_phenotype",
  "koges_lead_variant", "koges_chr", "koges_bp",
  "meta_exclusion_reason", "panel_analyzed", "EAS_EUR_both_testable",
  "contributing_cohorts", "koges_beta_aligned_a1",
  "external_panel_beta", "external_panel_SE", "external_panel_P",
  "direction_concordant_with_KoGES",
  "replication_nominal", "replication_Bonferroni",
  "signedZ_P", "signedZ_replication_Bonferroni"
)

require_cols(eas_meta, meta_required, "EAS_Replication")
require_cols(eur_meta, meta_required, "EUR_Replication")

eas_meta[, lead_row_id := as.integer(lead_row_id)]
eur_meta[, lead_row_id := as.integer(lead_row_id)]

if (anyDuplicated(eas_meta$lead_row_id)) stop("Duplicate lead_row_id in EAS_Replication")
if (anyDuplicated(eur_meta$lead_row_id)) stop("Duplicate lead_row_id in EUR_Replication")
if (!setequal(eas_meta$lead_row_id, eur_meta$lead_row_id)) {
  stop("EAS_Replication and EUR_Replication have different lead_row_id sets.")
}

setorder(eas_meta, lead_row_id)
setorder(eur_meta, lead_row_id)

identity_cols <- c(
  "lead_row_id", "koges_pheno_code", "koges_phenotype",
  "koges_lead_variant", "koges_chr", "koges_bp"
)

for (nm in identity_cols) {
  a <- as.character(eas_meta[[nm]])
  b <- as.character(eur_meta[[nm]])
  if (!identical(a, b)) stop("EAS/EUR identity mismatch in column: ", nm)
}

# Common-testable flag must also agree between EAS and EUR sheets.
common_e <- as_flag(eas_meta$EAS_EUR_both_testable)
common_u <- as_flag(eur_meta$EAS_EUR_both_testable)
common_equal <-
  (is.na(common_e) & is.na(common_u)) |
  (!is.na(common_e) & !is.na(common_u) & common_e == common_u)
if (!all(common_equal)) {
  stop("EAS_EUR_both_testable differs between EAS and EUR sheets.")
}

id_map <- eas_meta[, .(
  lead_row_id,
  Phenotype = as.character(koges_phenotype),
  `KoGES phenotype code` = as.character(koges_pheno_code),
  Chromosome = norm_chr(koges_chr),
  `Base pair position (GRCh37)` = norm_pos(koges_bp),
  rsID = as.character(koges_lead_variant)
)]

id_map[, identity_key := make_identity_key(
  Phenotype, Chromosome, `Base pair position (GRCh37)`, rsID
)]

if (anyDuplicated(id_map$identity_key)) {
  stop("Phenotype + chr + bp + rsID identity key is not unique in meta workbook.")
}

# ------------------------------------------------------------------------------
# 3. Build EAS/EUR summary + cohort-detail tables in memory from novelty TSVs
# ------------------------------------------------------------------------------
message("Reading Cont novelty TSVs ...")
bbj_nov <- fread(eas_novelty_files[["BBJ"]], sep = "\t", na.strings = c("NA", ""))
twb_nov <- fread(eas_novelty_files[["TWB"]], sep = "\t", na.strings = c("NA", ""))
tpmi_nov <- fread(eas_novelty_files[["TPMI"]], sep = "\t", na.strings = c("NA", ""))
ukb_nov <- fread(eur_novelty_files[["UKB"]], sep = "\t", na.strings = c("NA", ""))
finngen_nov <- fread(eur_novelty_files[["FinnGen"]], sep = "\t", na.strings = c("NA", ""))
mvp_nov <- fread(eur_novelty_files[["MVP"]], sep = "\t", na.strings = c("NA", ""))

# Same KoGES lead rows across all six Cont cohorts
ref_key <- bbj_nov[, .(lead_row_id, koges_phenotype, koges_lead_variant, koges_chr, koges_bp)]
for (nm in c("twb_nov", "tpmi_nov", "ukb_nov", "finngen_nov", "mvp_nov")) {
  obs <- get(nm)[, .(lead_row_id, koges_phenotype, koges_lead_variant, koges_chr, koges_bp)]
  if (!identical(ref_key, obs)) {
    stop("KoGES lead-association rows differ vs BBJ for: ", nm)
  }
}

raw_eas <- build_summary_eas_cont(bbj_nov, twb_nov, tpmi_nov)
raw_eur <- build_summary_eur_cont(ukb_nov, finngen_nov, mvp_nov)

cohort_detail_tables <- list(
  BBJ = build_cohort_detail_cont(bbj_nov, "bbj"),
  TWB = build_cohort_detail_cont(twb_nov, "twb"),
  TPMI = build_cohort_detail_cont(tpmi_nov, "tpmi"),
  UKB = build_cohort_detail_cont(ukb_nov, "ukb"),
  FinnGen = build_cohort_detail_cont(finngen_nov, "finngen"),
  MVP = build_cohort_detail_cont(mvp_nov, "mvp")
)
summary_tables <- list(Summary_EAS = raw_eas, Summary_EUR = raw_eur)

eas_req <- c(
  "Phenotype", "Chromosome", "Base pair position", "rsID",
  "Minimum P-value within ±1 Mb in BBJ",
  "Minimum P-value within ±1 Mb in TWB",
  "Minimum P-value within ±1 Mb in TPMI",
  "Novel relative to assessed EastAsian-ancestry cohorts"
)
eur_req <- c(
  "Phenotype", "Chromosome", "Base pair position", "rsID",
  "Minimum P-value within ±1 Mb in UKB",
  "Minimum P-value within ±1 Mb in FinnGen",
  "Minimum P-value within ±1 Mb in MVP",
  "Novel relative to assessed European-ancestry cohorts"
)
require_cols(raw_eas, eas_req, "Summary_EAS (built)")
require_cols(raw_eur, eur_req, "Summary_EUR (built)")

attach_lead_id <- function(raw, label) {
  key <- make_identity_key(
    raw[["Phenotype"]],
    raw[["Chromosome"]],
    raw[["Base pair position"]],
    raw[["rsID"]]
  )
  if (anyDuplicated(key)) stop(label, " identity key is not unique.")

  map <- id_map[, .(lead_row_id, identity_key)]
  idx <- match(key, map$identity_key)
  if (anyNA(idx)) {
    stop(
      label, ": ", sum(is.na(idx)),
      " row(s) could not be matched to lead_row_id."
    )
  }

  out <- copy(raw)
  out[, lead_row_id := map$lead_row_id[idx]]
  if (anyDuplicated(out$lead_row_id)) {
    stop(label, " produced duplicate lead_row_id after identity matching.")
  }
  out
}

raw_eas <- attach_lead_id(raw_eas, "Summary_EAS")
raw_eur <- attach_lead_id(raw_eur, "Summary_EUR")

eas_summary <- raw_eas[, .(
  lead_row_id,
  `Minimum P value within ±1 Mb in BBJ` =
    as_num(`Minimum P-value within ±1 Mb in BBJ`),
  `Minimum P value within ±1 Mb in TWB` =
    as_num(`Minimum P-value within ±1 Mb in TWB`),
  `Minimum P value within ±1 Mb in TPMI` =
    as_num(`Minimum P-value within ±1 Mb in TPMI`),
  `EAS locus novelty` =
    normalize_novelty(`Novel relative to assessed EastAsian-ancestry cohorts`)
)]

eur_summary <- raw_eur[, .(
  lead_row_id,
  `Minimum P value within ±1 Mb in UKB` =
    as_num(`Minimum P-value within ±1 Mb in UKB`),
  `Minimum P value within ±1 Mb in FinnGen` =
    as_num(`Minimum P-value within ±1 Mb in FinnGen`),
  `Minimum P value within ±1 Mb in MVP` =
    as_num(`Minimum P-value within ±1 Mb in MVP`),
  `EUR locus novelty` =
    normalize_novelty(`Novel relative to assessed European-ancestry cohorts`)
)]

# ------------------------------------------------------------------------------
# 4. Read GWAS Catalog novelty file and validate row_id == lead_row_id
# ------------------------------------------------------------------------------
catalog <- fread(
  catalog_tsv,
  sep = "\t",
  na.strings = c("NA", ""),
  data.table = TRUE
)

catalog_req <- c(
  "row_id", "Phenotype", "VarNM", "Lead variant", "VariantID",
  "CHR", "BP_GRCh37", "Nearest_gene", "NEA", "EA",
  "Beta", "SE", "P", "N",
  "efo_ids", "novel"
)
require_cols(catalog, catalog_req, "Cont_lead_variants_GWAScatalog_novel.tsv")

catalog[, lead_row_id := as.integer(row_id)]
if (anyDuplicated(catalog$lead_row_id)) stop("Duplicate row_id in GWAS Catalog TSV.")
if (!setequal(catalog$lead_row_id, id_map$lead_row_id)) {
  stop("GWAS Catalog row_id set does not match meta lead_row_id set.")
}

# Validate strong identity fields by lead_row_id.
cat_check <- merge(
  catalog[, .(
    lead_row_id,
    catalog_phenotype = as.character(Phenotype),
    catalog_code = as.character(VarNM),
    catalog_rsid = as.character(`Lead variant`),
    catalog_chr = norm_chr(CHR),
    catalog_bp37 = norm_pos(BP_GRCh37)
  )],
  id_map[, .(
    lead_row_id,
    meta_phenotype = Phenotype,
    meta_code = `KoGES phenotype code`,
    meta_rsid = rsID,
    meta_chr = Chromosome,
    meta_bp37 = `Base pair position (GRCh37)`
  )],
  by = "lead_row_id",
  sort = FALSE
)

bad_identity <- cat_check[
  catalog_phenotype != meta_phenotype |
    catalog_code != meta_code |
    catalog_rsid != meta_rsid |
    catalog_chr != meta_chr
]

# BP_GRCh37 is missing for a small number of Catalog rows; only validate where present.
bad_bp <- cat_check[
  !is.na(catalog_bp37) & nzchar(catalog_bp37) &
    catalog_bp37 != meta_bp37
]

if (nrow(bad_identity) > 0L) {
  stop("GWAS Catalog identity mismatch for ", nrow(bad_identity), " lead(s).")
}
if (nrow(bad_bp) > 0L) {
  stop("GWAS Catalog GRCh37 position mismatch for ", nrow(bad_bp), " lead(s).")
}

catalog_novel_flag <- as_flag(catalog$novel)

catalog_block <- catalog[, .(
  lead_row_id,
  `Variant ID` = as.character(VariantID),
  `Nearest gene` = as.character(Nearest_gene),
  `Effect allele` = as.character(EA),
  `Other allele` = as.character(NEA),
  `KoGES beta` = as_num(Beta),
  `KoGES standard error` = as_num(SE),
  `KoGES P value` = as_num(P),
  `KoGES N` = as_num(N),
  # Kept only to build Phenotype_GWAS_Catalog_map; dropped from main summary.
  `GWAS Catalog EFO ID` = as.character(efo_ids),
  `GWAS Catalog novelty` =
    fifelse(
      catalog_novel_flag %in% TRUE, "Novel",
      fifelse(catalog_novel_flag %in% FALSE, "Known", "Not assessed")
    )
)]

# ------------------------------------------------------------------------------
# 5. Prepare panel-level meta-analysis replication blocks
# ------------------------------------------------------------------------------
eas_meta_block <- make_meta_block(eas_meta, "EAS")
eur_meta_block <- make_meta_block(eur_meta, "EUR")

# ------------------------------------------------------------------------------
# 6. Merge all evidence layers
# ------------------------------------------------------------------------------
out <- Reduce(
  function(x, y) merge(x, y, by = "lead_row_id", all.x = TRUE, sort = FALSE),
  list(
    id_map[, setdiff(names(id_map), "identity_key"), with = FALSE],
    catalog_block,
    eas_summary,
    eur_summary,
    eas_meta_block,
    eur_meta_block
  )
)

if (nrow(out) != nrow(id_map)) {
  stop("Merged output row count differs from master lead count.")
}
if (anyDuplicated(out$lead_row_id)) {
  stop("Merged output has duplicate lead_row_id.")
}

# Integrated novelty fields.
# Strict novelty across EAS/EUR/GWAS Catalog: 3-level classification.
#   Known            : any evidence layer is Known
#   Strict level     : all three layers are Novel
#   Not strict level : otherwise (e.g. incomplete / not assessed layers)
out[, `Strict novelty across EAS, EUR, and GWAS Catalog` := {
  known <-
    `EAS locus novelty` == "Known" |
    `EUR locus novelty` == "Known" |
    `GWAS Catalog novelty` == "Known"
  novel_all <-
    `EAS locus novelty` == "Novel" &
    `EUR locus novelty` == "Novel" &
    `GWAS Catalog novelty` == "Novel"
  fifelse(known, "Known", fifelse(novel_all, "Strict level", "Not strict level"))
}]

out[, `Novelty evidence pattern` := paste0(
  "EAS=", `EAS locus novelty`,
  "; EUR=", `EUR locus novelty`,
  "; GWAS Catalog=", `GWAS Catalog novelty`
)]

# Integrated novelty: Strict novel + Incomplete assessment → Novel;
# any Known evidence layer → Known in ≥1 evidence source.
out[, `Integrated novelty status` := {
  known <-
    `EAS locus novelty` == "Known" |
    `EUR locus novelty` == "Known" |
    `GWAS Catalog novelty` == "Known"
  fifelse(known, "Known in ≥1 evidence source", "Novel")
}]

# Unique Phenotype ↔ GWAS Catalog EFO ID map (separate sheet).
pheno_catalog_map <- unique(out[, .(
  Phenotype,
  `GWAS Catalog EFO ID`
)])
setorder(pheno_catalog_map, Phenotype, `GWAS Catalog EFO ID`)

# Reorder columns: KoGES -> novelty -> panel meta replication.
# GWAS Catalog EFO ID / status are omitted here (EFO lives on
# Phenotype_GWAS_Catalog_map).
final_cols <- c(
  # KoGES association
  "lead_row_id",
  "Phenotype",
  "KoGES phenotype code",
  "Chromosome",
  "Base pair position (GRCh37)",
  "rsID",
  "Variant ID",
  "Nearest gene",
  "Effect allele",
  "Other allele",
  "KoGES beta",
  "KoGES standard error",
  "KoGES P value",
  "KoGES N",

  # Novelty
  "Minimum P value within ±1 Mb in BBJ",
  "Minimum P value within ±1 Mb in TWB",
  "Minimum P value within ±1 Mb in TPMI",
  "EAS locus novelty",
  "Minimum P value within ±1 Mb in UKB",
  "Minimum P value within ±1 Mb in FinnGen",
  "Minimum P value within ±1 Mb in MVP",
  "EUR locus novelty",
  "GWAS Catalog novelty",
  "Integrated novelty status",
  "Strict novelty across EAS, EUR, and GWAS Catalog",
  "Novelty evidence pattern",

  # Panel meta-analysis replication
  "EAS panel meta analyzed",
  "EAS panel contributing cohorts",
  "EAS panel meta beta",
  "EAS panel meta standard error",
  "EAS panel meta P value",
  "EAS meta direction vs KoGES",
  "EAS meta nominal replication",
  "EAS meta Bonferroni replication",
  "EAS signed-Z Bonferroni sensitivity",
  "EUR panel meta analyzed",
  "EUR panel contributing cohorts",
  "EUR panel meta beta",
  "EUR panel meta standard error",
  "EUR panel meta P value",
  "EUR meta direction vs KoGES",
  "EUR meta nominal replication",
  "EUR meta Bonferroni replication",
  "EUR signed-Z Bonferroni sensitivity"
)

require_cols(out, final_cols, "Merged output")
out <- out[, ..final_cols]
setorder(out, lead_row_id)

# ------------------------------------------------------------------------------
# 7. Build definitions and thresholds
# ------------------------------------------------------------------------------
# Meta-analysis Bonferroni thresholds.
n_eas_meta <- sum(as_flag(eas_meta$panel_analyzed) %in% TRUE, na.rm = TRUE)
n_eur_meta <- sum(as_flag(eur_meta$panel_analyzed) %in% TRUE, na.rm = TRUE)

panel_test_info <- rbindlist(list(
  data.table(
    Analysis = "EAS panel meta-analysis replication",
    N_tests = n_eas_meta,
    Type = "Panel meta-analysis"
  ),
  data.table(
    Analysis = "EUR panel meta-analysis replication",
    N_tests = n_eur_meta,
    Type = "Panel meta-analysis"
  )
))

panel_test_info[, Bonferroni_threshold := 0.05 / N_tests]

definitions <- data.table(
  Field = c(
    "lead_row_id",
    "Minimum P value within ±1 Mb in [cohort]",
    "EAS locus novelty",
    "EUR locus novelty",
    "GWAS Catalog novelty",
    "Integrated novelty status",
    "Strict novelty across EAS, EUR, and GWAS Catalog",
    "Novelty evidence pattern",
    "EAS / EUR panel meta analyzed",
    "EAS / EUR panel contributing cohorts",
    "EAS / EUR panel meta beta / SE / P value",
    "EUR panel IVW cohort composition",
    "EAS / EUR meta direction vs KoGES",
    "EAS / EUR meta nominal replication",
    "EAS / EUR meta Bonferroni replication",
    "EAS / EUR signed-Z Bonferroni sensitivity",
    "Not assessed / Not analyzed",
    "Phenotype_GWAS_Catalog_map (separate sheet)"
  ),
  Definition = c(
    "Unique phenotype–lead association identifier. One row is one KoGES phenotype–lead test; the same SNP in different phenotypes is a different lead_row_id.",
    "Lowest P value among genome-wide significant external signals within ±1 Mb of the KoGES lead in that cohort (BBJ, TWB, TPMI for EAS; UKB, FinnGen, MVP for EUR). Used for locus-level novelty assessment.",
    "Locus-level novelty vs assessed East Asian cohorts (BBJ, TWB, TPMI). Levels: Novel (no GWS signal within ±1 Mb in assessed cohorts), Known (GWS signal present in ≥1 assessed cohort), Not assessed (no evaluable cohort).",
    "Locus-level novelty vs assessed European cohorts (UKB, FinnGen, MVP). Levels: Novel, Known, Not assessed (same rules as EAS).",
    "Trait-matched novelty from GWAS Catalog mapping (Cont_lead_variants_GWAScatalog_novel.tsv; ±1 Mb). Levels: Novel, Known, Not assessed. EFO IDs are not in the main summary sheet; see Phenotype_GWAS_Catalog_map.",
    "Two-level summary across EAS, EUR, and GWAS Catalog. Novel = none of the three layers is Known (covers both full three-layer novelty and incomplete assessment). Known in ≥1 evidence source = any of EAS, EUR, or GWAS Catalog is Known.",
    "Three-level classification across EAS, EUR, and GWAS Catalog. Known = any layer is Known. Strict level = all three layers are Novel. Not strict level = otherwise (e.g. some layers Not assessed, and none Known).",
    "Compact string of the three layer labels, e.g. EAS=Novel; EUR=Known; GWAS Catalog=Novel.",
    "Whether the phenotype–lead was included in that ancestry panel’s external-only meta-analysis (Analyzed / Not analyzed).",
    "Cohorts that contributed effect estimates to the PRIMARY panel IVW for that lead. EAS pools usable BBJ/TWB/TPMI; EUR pools usable UKB/FinnGen/MVP. Composition varies by lead depending on which cohorts are usable.",
    "Primary external-panel summary statistics from fixed-effect IVW of the contributing cohorts listed for that lead (effect size, SE, and two-sided P). All usable cohorts in that ancestry panel enter the primary IVW (EAS and EUR use the same rule).",
    "EUR primary IVW includes every usable cohort among UKB, FinnGen, and MVP (same logic as EAS BBJ/TWB/TPMI). FinnGen is not preferentially dropped when UKB and/or MVP are usable. Per-cohort FinnGen columns and locus-novelty scans also use FinnGen. Signed-Z sensitivity uses the same usable-cohort set as primary IVW.",
    "Whether the panel meta effect direction is concordant with KoGES (Concordant / Discordant / Not available / Not analyzed).",
    "Panel-level nominal (P < 0.05) replication status requiring concordant direction with KoGES (Replicated / Not replicated / Not available / Not analyzed).",
    "Primary panel-level replication: concordant direction with KoGES and panel-specific Bonferroni significance (0.05 / N_tests for that panel). Levels: Replicated / Not replicated / Not available / Not analyzed.",
    "Sensitivity using sample-size-weighted signed-Z meta-analysis and the same panel Bonferroni family. Uses the same usable cohorts as the primary panel IVW (EAS: BBJ/TWB/TPMI; EUR: UKB/FinnGen/MVP). Levels: Replicated / Not replicated / Not available / Not analyzed.",
    "Not assessed = novelty evidence unavailable or not eligible for that layer. Not analyzed = no panel-level replication statistic was produced for that lead.",
    "Unique Phenotype / GWAS Catalog EFO ID rows. Phenotype–Catalog EFO mapping is stored here instead of repeating EFO columns on every lead-level summary row."
  )
)

# ------------------------------------------------------------------------------
# 8. Verification messages
# ------------------------------------------------------------------------------
message("Rows in final Supplementary Data: ", nrow(out))
message(
  "Integrated novelty = Novel: ",
  sum(out[["Integrated novelty status"]] == "Novel")
)
message(
  "Integrated novelty = Known in ≥1 evidence source: ",
  sum(out[["Integrated novelty status"]] == "Known in ≥1 evidence source")
)
message(
  "Strict novelty = Known: ",
  sum(out[["Strict novelty across EAS, EUR, and GWAS Catalog"]] == "Known")
)
message(
  "Strict novelty = Strict level: ",
  sum(out[["Strict novelty across EAS, EUR, and GWAS Catalog"]] == "Strict level")
)
message(
  "Strict novelty = Not strict level: ",
  sum(out[["Strict novelty across EAS, EUR, and GWAS Catalog"]] == "Not strict level")
)
message(
  "Unique Phenotype–GWAS Catalog map rows: ",
  nrow(pheno_catalog_map)
)
message(
  "EAS panel meta Bonferroni replicated: ",
  sum(out[["EAS meta Bonferroni replication"]] == "Replicated")
)
message(
  "EUR panel meta Bonferroni replicated: ",
  sum(out[["EUR meta Bonferroni replication"]] == "Replicated")
)

# ------------------------------------------------------------------------------
# 9. Write ONE final workbook
# ------------------------------------------------------------------------------
wb <- createWorkbook()

# Excel sheet names must be <= 31 characters.
# Requested name "Total_summary_Novelty_Replication" is 33 chars, so use
# the closest readable abbreviation that fits the limit.
sheet_merged <- "Total_summary_Novelty_Repl"
sheet_defs <- "Definitions"
sheet_map <- "Phenotype_GWAS_Catalog_map"

addWorksheet(wb, sheet_merged)
addWorksheet(wb, sheet_defs)
addWorksheet(wb, sheet_map)

# Main data table
writeDataTable(
  wb,
  sheet = sheet_merged,
  x = as.data.frame(out),
  startRow = 1,
  startCol = 1,
  tableName = "ExternalValidationSupplementaryData",
  tableStyle = "TableStyleMedium2",
  withFilter = TRUE
)

freezePane(wb, sheet_merged, firstRow = TRUE, firstCol = TRUE)

# Column widths: readable but bounded.
setColWidths(
  wb, sheet_merged,
  cols = seq_len(ncol(out)),
  widths = 16
)

wide_cols <- which(names(out) %in% c(
  "Phenotype",
  "KoGES phenotype code",
  "Novelty evidence pattern",
  "EAS panel contributing cohorts",
  "EUR panel contributing cohorts"
))
if (length(wide_cols) > 0L) {
  setColWidths(wb, sheet_merged, cols = wide_cols, widths = 32)
}

# Compact identifier columns.
for (nm in c("lead_row_id", "Chromosome", "rsID")) {
  j <- match(nm, names(out))
  if (!is.na(j)) setColWidths(wb, sheet_merged, cols = j, widths = 12)
}

# Scientific notation for P values.
p_cols <- grep("P value|Minimum P value", names(out))
if (length(p_cols) > 0L) {
  p_style <- createStyle(numFmt = "0.00E+00")
  addStyle(
    wb, sheet_merged,
    style = p_style,
    rows = 2:(nrow(out) + 1L),
    cols = p_cols,
    gridExpand = TRUE,
    stack = TRUE
  )
}

# Numeric formatting for beta / SE.
effect_cols <- grep("beta|standard error", names(out), ignore.case = TRUE)
if (length(effect_cols) > 0L) {
  effect_style <- createStyle(numFmt = "0.0000")
  addStyle(
    wb, sheet_merged,
    style = effect_style,
    rows = 2:(nrow(out) + 1L),
    cols = effect_cols,
    gridExpand = TRUE,
    stack = TRUE
  )
}

# Definitions sheet
title_style <- createStyle(
  fontSize = 14,
  textDecoration = "bold"
)
wrap_style <- createStyle(
  wrapText = TRUE,
  valign = "top"
)

writeData(wb, sheet_defs, "Field definitions", startRow = 1, startCol = 1)
addStyle(wb, sheet_defs, title_style, rows = 1, cols = 1)

writeDataTable(
  wb, sheet_defs,
  x = as.data.frame(definitions),
  startRow = 3,
  startCol = 1,
  tableName = "FieldDefinitions",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)

threshold_start <- nrow(definitions) + 6L
writeData(
  wb, sheet_defs,
  "Bonferroni thresholds used for replication",
  startRow = threshold_start,
  startCol = 1
)
addStyle(
  wb, sheet_defs,
  title_style,
  rows = threshold_start,
  cols = 1
)

writeDataTable(
  wb, sheet_defs,
  x = as.data.frame(panel_test_info),
  startRow = threshold_start + 2L,
  startCol = 1,
  tableName = "ReplicationThresholds",
  tableStyle = "TableStyleMedium2",
  withFilter = FALSE
)

# Format threshold column.
threshold_col <- 4L
threshold_rows <- (threshold_start + 3L):(threshold_start + 2L + nrow(panel_test_info))
addStyle(
  wb, sheet_defs,
  style = createStyle(numFmt = "0.000E+00"),
  rows = threshold_rows,
  cols = threshold_col,
  gridExpand = TRUE,
  stack = TRUE
)

note_start <- threshold_start + 2L + nrow(panel_test_info) + 3L
writeData(
  wb, sheet_defs,
  "EUR panel primary IVW note",
  startRow = note_start,
  startCol = 1
)
addStyle(wb, sheet_defs, title_style, rows = note_start, cols = 1)
eur_ivw_notes <- c(
  "EUR primary IVW uses the same rule as EAS: all usable cohorts enter the panel meta.",
  "Usable UKB, FinnGen, and MVP are pooled together in primary IVW / Q / LOCO (max k=3).",
  "Not every row is UKB+FinnGen+MVP; composition still depends on which cohorts are usable for that lead.",
  "FinnGen also remains in per-cohort detail sheets and in EUR locus-novelty (min-P within ±1 Mb) assessment.",
  "EUR signed-Z Bonferroni sensitivity uses the same usable-cohort set as primary IVW."
)
writeData(
  wb, sheet_defs,
  data.table(Note = eur_ivw_notes),
  startRow = note_start + 2L,
  startCol = 1,
  colNames = FALSE
)

setColWidths(wb, sheet_defs, cols = 1, widths = 48)
setColWidths(wb, sheet_defs, cols = 2, widths = 95)
setColWidths(wb, sheet_defs, cols = 3:4, widths = 22)
addStyle(
  wb, sheet_defs,
  wrap_style,
  rows = 3:(note_start + 2L + length(eur_ivw_notes)),
  cols = 1:4,
  gridExpand = TRUE,
  stack = TRUE
)
freezePane(wb, sheet_defs, firstRow = TRUE)

# Phenotype ↔ GWAS Catalog mapping sheet (after Definitions)
writeDataTable(
  wb,
  sheet = sheet_map,
  x = as.data.frame(pheno_catalog_map),
  startRow = 1,
  startCol = 1,
  tableName = "PhenotypeGWASCatalogMap",
  tableStyle = "TableStyleMedium2",
  withFilter = TRUE
)
freezePane(wb, sheet_map, firstRow = TRUE)
setColWidths(wb, sheet_map, cols = 1, widths = 28)
setColWidths(wb, sheet_map, cols = 2, widths = 36)

# ------------------------------------------------------------------------------
# 10. Append in-memory summary/detail tables + panel-meta sheets
# ------------------------------------------------------------------------------
cohort_sheet_plan <- list(
  list(src = "Summary_EAS", dest = "Ancestry_summary_EAS", colour = "#1F4E78", kind = "summary"),
  list(src = "BBJ", dest = "Cohort_detail_BBJ", colour = "#2E75B6", kind = "detail"),
  list(src = "TWB", dest = "Cohort_detail_TWB", colour = "#2E75B6", kind = "detail"),
  list(src = "TPMI", dest = "Cohort_detail_TPMI", colour = "#2E75B6", kind = "detail"),
  list(src = "Summary_EUR", dest = "Ancestry_summary_EUR", colour = "#548235", kind = "summary"),
  list(src = "UKB", dest = "Cohort_detail_UKB", colour = "#70AD47", kind = "detail"),
  list(src = "FinnGen", dest = "Cohort_detail_FinnGen", colour = "#70AD47", kind = "detail"),
  list(src = "MVP", dest = "Cohort_detail_MVP", colour = "#70AD47", kind = "detail")
)

for (item in cohort_sheet_plan) {
  dt <- if (identical(item$kind, "summary")) {
    summary_tables[[item$src]]
  } else {
    cohort_detail_tables[[item$src]]
  }
  info <- append_dt_sheet(wb, dt, item$dest, tab_colour = item$colour)
  message(
    "Appended ", item$src, " -> ", info$sheet,
    " (", info$nrow, " x ", info$ncol, ")"
  )
}

panel_sheet_plan <- list(
  list(dt = eas_meta, dest = "Replication_EAS_meta", colour = "#5B9BD5"),
  list(dt = eur_meta, dest = "Replication_EUR_meta", colour = "#A9D08E")
)
for (item in panel_sheet_plan) {
  info <- append_dt_sheet(wb, item$dt, item$dest, tab_colour = item$colour)
  message(
    "Appended panel-meta sheet -> ", info$sheet,
    " (", info$nrow, " x ", info$ncol, ")"
  )
}

dir.create(dirname(output_xlsx), recursive = TRUE, showWarnings = FALSE)
saveWorkbook(wb, output_xlsx, overwrite = TRUE)

message("Wrote: ", normalizePath(output_xlsx))
message("Final sheets: ", paste(wb$sheet_names, collapse = ", "))
