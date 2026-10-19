#!/usr/bin/env Rscript
# ==============================================================================
# Purpose
#   Classify KoGES Cate/Cont lead variants as novel vs known relative to the
#   GWAS Catalog (ontology-annotated associations). A lead is "known" if a
#   Catalog association with P <= GWS_P falls within +/- WINDOW_BP of the
#   KoGES lead on the same chromosome (same GWS threshold used for cohort
#   novelty / COJO). Trait matching uses EFO/MONDO/... IDs from koges_EFOID.xlsx.
#
# Run order (typical manuscript pipeline)
#   1. Lift leads to GRCh38 (liftover_sumstat.R → *_lead_variants_hg38.tsv)
#   2. GROUP=Cate|Cont Rscript gwas_catalog_novelty.R
#   3. {Cate,Cont}_lead_variants_GWAScatalog_novel.tsv feeds Supplementary Data 3 merge
#
# Inputs
#   - Extended_Data/{Cate,Cont}_lead_variants_hg38.tsv (Phenotype, CHR, BP;
#     Cont also needs VarNM for Excel label mapping)
#   - koges_EFOID.xlsx (sheet 1 = Cate, sheet 2 = Cont phenotype ↔ ontology IDs)
#   - gwas-catalog-associations_ontology-annotated-full.zip
#       (internal TSV: gwas-catalog-download-associations-alt-full.tsv)
#
# Outputs
#   - Extended_Data/{Cate|Cont}_lead_variants_GWAScatalog_novel.tsv  (default ±1 Mb)
#   - Non-default WINDOW_BP/MB appends a size tag, e.g. ..._GWAScatalog_novel_500kb.tsv
#   - Columns include novel, gwas_status (novel / known_trait),
#     novel_any_trait / gwas_status_any (any Catalog hit in window),
#     catalog_matched_disease_traits
#
# Environment variables
#   GROUP      Cate (default) or Cont
#   WINDOW_BP  half-window in bp (default 1000000)
#   WINDOW_MB  alternative to WINDOW_BP (e.g. 0.5 → 500000)
#   GWS_P      Catalog P threshold for Known (default 5e-8)
#
# Usage
#   GROUP=Cate Rscript gwas_catalog_novelty.R
#   GROUP=Cont Rscript gwas_catalog_novelty.R
#   WINDOW_MB=0.5 GROUP=Cate Rscript gwas_catalog_novelty.R
#   WINDOW_BP=500000 GROUP=Cont Rscript gwas_catalog_novelty.R
#   GWS_P=5e-8 GROUP=Cont Rscript gwas_catalog_novelty.R
#
# Notes
#   - Default ±1 Mb writes the untagged manuscript filename used by merge / COJO.
#   - Cont phenotypes are remapped VarNM (RIN_B_*, B_FEV1*, ...) → Excel
#     KoGES_Cont labels via CONT_VARNM_BY_EXCEL_ROW (1:1 with Excel sheet 2 rows).
#   - Leads with no Excel phenotype match are dropped; n_efo_ids==0 stays novel.
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

group <- Sys.getenv("GROUP", unset = "Cate")
ext_dir <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
if (group == "Cont") {
  variant_tsv_path <- file.path(ext_dir, "Cont_lead_variants_hg38.tsv")
} else {
  variant_tsv_path <- file.path(ext_dir, "Cate_lead_variants_hg38.tsv")
}
excel_path <- "/home2/sjim/openKoGES/DATA_baseline/forGWAS_annot/koges_EFOID.xlsx"
assoc_zip_path <-  "/home2/sjim/openKoGES/DATA_baseline/forGWAS_annot/gwas-catalog-associations_ontology-annotated-full.zip"
internal_tsv <-  "gwas-catalog-download-associations-alt-full.tsv"

window_bp <- as.integer(Sys.getenv("WINDOW_BP", unset = "1000000"))
wm <- Sys.getenv("WINDOW_MB", unset = "")
if (nzchar(wm)) {
  window_bp <- as.integer(round(as.numeric(wm) * 1e6))
}
if (is.na(window_bp) || window_bp < 1L) {
  stop("Invalid WINDOW_BP / WINDOW_MB (got window_bp=", window_bp, ")")
}
GWS_P <- as.numeric(Sys.getenv("GWS_P", unset = "5e-8"))
if (!is.finite(GWS_P) || GWS_P <= 0 || GWS_P > 1) {
  stop("GWS_P must be in (0, 1].")
}

# Default ±1 Mb: untagged manuscript name; other windows get a size suffix
window_tag <- if (window_bp == 1000000L) {
  ""
} else if (window_bp %% 1000000L == 0L) {
  paste0("_", window_bp %/% 1000000L, "Mb")
} else if (window_bp %% 1000L == 0L) {
  paste0("_", window_bp %/% 1000L, "kb")
} else {
  paste0("_", window_bp, "bp")
}
output_path <- file.path(
  ext_dir, paste0(group, "_lead_variants_GWAScatalog_novel", window_tag, ".tsv")
)
max_rows <- NA
message(
  "group=", group, " window_bp=", window_bp,
  " GWS_P=", format(GWS_P, scientific = TRUE),
  " → ", output_path
)

normalize_ws <- function(x) {
  x <- as.character(x)
  x <- gsub("\u00a0", " ", x, fixed = TRUE)
  gsub("^\\s+|\\s+$", "", x, perl = TRUE)
}

parse_efo_ids <- function(raw) {
  if (length(raw) != 1L || is.na(raw)) return(character(0))
  s <- normalize_ws(raw)
  if (!nzchar(s)) return(character(0))
  parts <- unlist(strsplit(s, "\\s*[/;,]\\s*", perl = TRUE))
  parts <- parts[nzchar(parts)]
  ids <- unique(gsub("^.*[/]", "", parts))
  ids <- gsub("#.*$", "", ids, perl = TRUE)
  ids <- ids[grepl("^(EFO|MONDO|OBA|HP|GO|CHEBI|UBERON|CL)_", ids, perl = TRUE)]
  unique(ids[nzchar(ids)])
}

detect_map_cols <- function(df) {
  nms <- names(df)
  ko <- NA_character_
  if (any(grepl("KoGES", nms, ignore.case = TRUE) & grepl("Cont", nms, ignore.case = TRUE))) {
    ko <- nms[which(grepl("KoGES", nms, ignore.case = TRUE) & grepl("Cont", nms, ignore.case = TRUE))[1]]
  } else if (any(grepl("KoGES", nms, ignore.case = TRUE) & grepl("Cate", nms, ignore.case = TRUE))) {
    ko <- nms[which(grepl("KoGES", nms, ignore.case = TRUE) & grepl("Cate", nms, ignore.case = TRUE))[1]]
  } else if (any(grepl("Cont", nms, ignore.case = TRUE))) {
    ko <- nms[grepl("Cont", nms, ignore.case = TRUE)][1]
  } else if (any(grepl("Cate", nms, ignore.case = TRUE))) {
    ko <- nms[grepl("Cate", nms, ignore.case = TRUE)][1]
  }
  idcol <- NA_character_
  if (any(grepl("EFO", nms, ignore.case = TRUE))) {
    idcol <- nms[grepl("EFO", nms, ignore.case = TRUE)][1]
  } else if (any(grepl("MONDO|OBA", nms, ignore.case = TRUE))) {
    idcol <- nms[grepl("MONDO|OBA", nms, ignore.case = TRUE)][1]
  }
  list(ko_col = ko, id_col = idcol)
}

alias_lead_variant <- function(dt) {
  if ("Lead variant" %in% names(dt)) return(invisible(dt))
  cands <- c("Lead_variant", "rsID", "SNP", "rs_id", "Lead SNP", "variant_id", "Variant", "RSID")
  for (nm in cands) {
    if (nm %in% names(dt)) {
      setnames(dt, nm, "Lead variant")
      return(invisible(dt))
    }
  }
  stop("Variant TSV missing columns: Lead variant (or rsID/SNP/...)")
}

clean_fread_names <- function(dt) {
  old <- names(dt)
  new <- sub("^\ufeff", "", old, perl = TRUE)
  new <- gsub("\u00a0", " ", new, fixed = TRUE)
  new <- gsub("^\\s+|\\s+$", "", new, perl = TRUE)
  if (!identical(old, new)) setnames(dt, old, new)
  invisible(dt)
}

alias_phenotype <- function(dt) {
  if ("Phenotype" %in% names(dt)) return(invisible(dt))
  nms <- names(dt)
  hit <- nms[grepl("^phenotype$", nms, ignore.case = TRUE)]
  if (length(hit) == 1L) {
    setnames(dt, hit, "Phenotype")
    return(invisible(dt))
  }
  cands <- c("Trait", "DISEASE/TRAIT", "DISEASE.TRAIT", "Phenotype name", "KoGES phenotype")
  for (nm in cands) {
    if (nm %in% nms) {
      setnames(dt, nm, "Phenotype")
      return(invisible(dt))
    }
  }
  stop("Variant TSV missing column Phenotype (or Trait / case variant of Phenotype). Found: ", paste(nms, collapse = ", "))
}

resolve_col <- function(names_vec, primary, fallbacks) {
  if (!is.na(primary) && primary %in% names_vec) return(primary)
  for (f in fallbacks) {
    if (f %in% names_vec) return(f)
  }
  NULL
}

# Normalize Chromosome / Base pair position (e.g. KoGES_lead_SNP_table) → CHR / BP
alias_chr_bp <- function(dt) {
  if ("CHR" %in% names(dt) && "BP" %in% names(dt)) return(invisible(dt))
  chr_col <- resolve_col(names(dt), NA_character_, c("CHR", "#CHR", "Chromosome", "chrom", "CHROM"))
  bp_col <- resolve_col(names(dt), NA_character_, c("BP", "POS", "pos", "GENPOS", "Base pair position"))
  if (is.null(chr_col) || is.null(bp_col)) {
    stop("Variant TSV needs CHR & BP (or Chromosome & Base pair position). Found: ", paste(names(dt), collapse = ", "))
  }
  if (chr_col != "CHR") setnames(dt, chr_col, "CHR")
  if (bp_col != "BP") setnames(dt, bp_col, "BP")
  invisible(dt)
}

# Cont Excel sheet 2 (KoGES_Cont) row order ↔ VarNM / B_* codes (NA = no lead table).
# Keep this vector length in sync when Excel rows are added or removed.
CONT_VARNM_BY_EXCEL_ROW <- c(
  "RIN_B_ALT", "RIN_B_ALBUMIN", "RIN_B_ALP", "RIN_B_APOA", "RIN_B_APOB", "RIN_B_AST", "RIN_B_BUN",
  "RIN_B_HEIGHT", "RIN_B_BMI", "RIN_B_WEIGHT", "RIN_B_CPEPTIDE", "RIN_B_CYSTATINC", "RIN_B_DBP",
  "B_FEV1best_noNC", "B_FFRbest_noNC", "B_FVCbest_noNC",
  "RIN_B_FOLATE", "RIN_B_R_GTP", "RIN_B_GLU0", "RIN_B_HDL", "RIN_B_HBA1C", "RIN_B_HCT",
  "RIN_B_HB", "RIN_B_HIP", "RIN_B_HOMOCYST", "RIN_B_INSUL", "RIN_B_IRON", "RIN_B_LDLclean",
  "RIN_B_MCH", "RIN_B_MCHC", "RIN_B_MCV", "RIN_B_PHOSPHATE", "RIN_B_PLAT", "RIN_B_PULSE",
  "RIN_B_RBC_B", "RIN_B_CREATINE", "RIN_B_SBP", "RIN_B_T_BIL", "RIN_B_TCHL",
  "RIN_B_TG", "RIN_B_URICACID", "RIN_B_VITD", "RIN_B_WAIST", "RIN_B_WHR", "RIN_B_WBC_B", "RIN_B_eGFR"
)

build_cont_varnm_lookup <- function(excel_path) {
  raw <- as.data.table(read_excel(excel_path, sheet = 2L))
  nms <- names(raw)
  ko <- nms[grepl("KoGES", nms, ignore.case = TRUE) & grepl("Cont", nms, ignore.case = TRUE)][1L]
  if (is.na(ko)) stop("Excel sheet 2: KoGES_Cont column not found")
  lab <- normalize_ws(raw[[ko]])
  if (length(lab) != length(CONT_VARNM_BY_EXCEL_ROW)) {
    stop(
      "Excel KoGES_Cont rows (", length(lab), ") != CONT_VARNM_BY_EXCEL_ROW (",
      length(CONT_VARNM_BY_EXCEL_ROW), "). Update CONT_VARNM_BY_EXCEL_ROW in script."
    )
  }
  nv <- data.table(varnm = CONT_VARNM_BY_EXCEL_ROW, KoGES_map_label = lab)
  nv[!is.na(varnm)]
}

# ------------------------------------------------------------------------------
# 1. Load GRCh38 lead variants
# ------------------------------------------------------------------------------
variants <- fread(variant_tsv_path, check.names = FALSE)
clean_fread_names(variants)
alias_phenotype(variants)
alias_lead_variant(variants)
alias_chr_bp(variants)
if (!all(c("Phenotype", "CHR", "BP") %in% names(variants))) {
  stop("Variant TSV missing columns: Phenotype, CHR, BP. Found: ", paste(names(variants), collapse = ", "))
}

if (!is.na(max_rows) && nzchar(max_rows)) {
  n <- as.integer(max_rows)
  variants <- variants[seq_len(min(.N, n))]
}

variants[, Phenotype := normalize_ws(Phenotype)]
variants[, row_id := .I]
variants[, CHR := as.character(CHR)]
variants[, BP := suppressWarnings(as.integer(BP))]
if (anyNA(variants$BP)) stop("BP must be integer-coercible.")
if ("VarNM" %in% names(variants)) {
  variants[, VarNM_base := sub("(_BOLT_result|_SAIGE.*)$", "", VarNM, perl = TRUE)]
}
variants[, start := pmax(1L, BP - window_bp)]
variants[, end := BP + window_bp]

# ------------------------------------------------------------------------------
# 2. Phenotype ↔ ontology map (Excel)
# ------------------------------------------------------------------------------
if (group == "Cont") {
  sheet_num <- 2
} else {
  sheet_num <- 1
}
map_df <- as.data.table(read_excel(excel_path, sheet = sheet_num))
mc <- detect_map_cols(map_df)
if (is.na(mc$ko_col) || is.na(mc$id_col)) {
  stop("Could not detect KoGES map column or ID column in Excel. Columns: ", paste(names(map_df), collapse = ", "))
}
map_df <- map_df[, c(mc$ko_col, mc$id_col), with = FALSE]
setnames(map_df, c("KoGES_map_label", "EFO_id_raw"))
map_df[, KoGES_map_label := normalize_ws(KoGES_map_label)]
map_df <- map_df[!is.na(KoGES_map_label) & KoGES_map_label != ""]
map_df <- map_df[, .(EFO_id_raw = paste(unique(as.character(EFO_id_raw)), collapse = " | ")), by = KoGES_map_label]

# Cont: normalize VarNM-style Phenotype → Excel KoGES_map_label, then join
variants[, Phenotype_map := Phenotype]
if (identical(group, "Cont")) {
  nv <- build_cont_varnm_lookup(excel_path)
  if (!"VarNM_base" %in% names(variants)) {
    stop("Cont mapping requires VarNM column to derive VarNM_base.")
  }
  variants[nv, Phenotype_map := i.KoGES_map_label, on = c(VarNM_base = "varnm")]
  variants[is.na(Phenotype_map), Phenotype_map := Phenotype]
  n_hit <- variants[Phenotype != Phenotype_map, .N]
  if (n_hit > 0L) {
    message("[Cont] mapped VarNM-style Phenotype → KoGES label for ", n_hit, " row(s)")
  }
}

variants <- map_df[variants, on = c(KoGES_map_label = "Phenotype_map"), nomatch = 0L]
if ("Phenotype_map" %in% names(variants)) variants[, Phenotype_map := NULL]

if (nrow(variants) == 0L) {
  stop(
    "No rows after joining Phenotype to Excel KoGES map label. Check group (Cate vs Cont), ",
    "VarNM↔Excel mapping (CONT_VARNM_BY_EXCEL_ROW), spelling / NBSP."
  )
}

variants[, efo_ids := lapply(EFO_id_raw, parse_efo_ids)]
variants[, n_efo_ids := lengths(efo_ids)]

# ------------------------------------------------------------------------------
# 3. Load GWAS Catalog associations (GWS filter)
# ------------------------------------------------------------------------------
if (!file.exists(assoc_zip_path)) {
  stop("assoc zip not found: ", assoc_zip_path)
}
message("Loading GWAS associations from zip (heavy)...")
cmd <- sprintf("unzip -p %s %s", shQuote(assoc_zip_path, type = "sh"), shQuote(internal_tsv, type = "sh"))
assoc <- fread(cmd = cmd, quote = "", check.names = FALSE, showProgress = FALSE)
cn <- names(assoc)
dtcol <- cn[cn %in% c("DISEASE/TRAIT", "DISEASE.TRAIT")]
if (length(dtcol) != 1L) stop("Cannot find DISEASE/TRAIT column in association TSV. Names sample: ", paste(head(cn, 10), collapse = ", "))
pcol <- cn[cn %in% c("P-VALUE", "PVALUE", "P_VALUE")]
if (length(pcol) != 1L) {
  stop(
    "Cannot find P-VALUE column in association TSV. Names sample: ",
    paste(head(cn, 20), collapse = ", ")
  )
}
assoc <- assoc[, c("CHR_ID", "CHR_POS", "MAPPED_TRAIT_URI", dtcol, pcol), with = FALSE]
setnames(assoc, dtcol, "catalog_disease_trait")
setnames(assoc, pcol, "catalog_p")

assoc[, CHR := as.character(CHR_ID)]
assoc[, pos := suppressWarnings(as.integer(CHR_POS))]
assoc[, catalog_p := suppressWarnings(as.numeric(catalog_p))]
n_assoc_raw <- nrow(assoc)
assoc <- assoc[
  !is.na(pos) & !is.na(MAPPED_TRAIT_URI) &
    is.finite(catalog_p) & catalog_p <= GWS_P
]
message(
  "Catalog GWS filter: kept ", format(nrow(assoc), big.mark = ","),
  " / ", format(n_assoc_raw, big.mark = ","),
  " rows with finite P <= ", format(GWS_P, scientific = TRUE)
)

# Expand multi-value MAPPED_TRAIT_URI to one ontology ID per row
assoc <- assoc[, {
  uris <- trimws(unlist(strsplit(as.character(MAPPED_TRAIT_URI), ",", fixed = FALSE)))
  uris <- uris[nzchar(uris)]
  ids <- sub("^.*[/]", "", uris)
  ids <- gsub("#.*$", "", ids, perl = TRUE)
  ids <- trimws(ids)
  keep <- !is.na(ids) & nzchar(ids)
  .(efo_id = ids[keep], catalog_disease_trait = rep(catalog_disease_trait[1], sum(keep)))
}, by = .(CHR, pos)]

assoc <- assoc[grepl("^(EFO|MONDO|OBA|HP|GO|CHEBI|UBERON|CL)_", efo_id, perl = TRUE)]

# Point intervals for foverlaps
assoc[, start := pos]
assoc[, end := pos]

assoc_any <- unique(assoc[, .(CHR, start, end)])
setkey(assoc_any, CHR, start, end)

# ------------------------------------------------------------------------------
# 4. Novelty: any Catalog hit in window (trait-agnostic)
# ------------------------------------------------------------------------------
xov <- variants[, .(row_id, CHR, start, end)]
setkey(xov, CHR, start, end)
ov_any <- foverlaps(xov, assoc_any, type = "any", nomatch = 0L)
hit_any <- unique(ov_any$row_id)
variants[, novel_any_trait := !row_id %in% hit_any]
variants[, gwas_status_any := ifelse(novel_any_trait, "novel", "known_other_trait")]

# ------------------------------------------------------------------------------
# 5. Novelty: overlapping Catalog hit with shared ontology ID
# ------------------------------------------------------------------------------
expl <- variants[n_efo_ids > 0, {
  ids <- unique(unlist(efo_ids))
  if (!length(ids)) NULL else data.table(efo_id = ids)
}, by = .(row_id, CHR, start, end)]

if (nrow(expl) == 0L) {
  variants[, novel := TRUE]
  variants[, gwas_status := "novel"]
  variants[, catalog_matched_disease_traits := NA_character_]
} else {
  assoc_sub <- assoc[efo_id %in% unique(expl$efo_id)]
  setkey(assoc_sub, efo_id, CHR, start, end)
  setkey(expl, efo_id, CHR, start, end)
  ov <- foverlaps(expl, assoc_sub, type = "any", nomatch = 0L)
  hit_rows <- unique(ov$row_id)

  trait_by_row <- ov[, .(catalog_matched_disease_traits = paste(sort(unique(catalog_disease_trait)), collapse = " | ")), by = row_id]

  variants[, novel := TRUE]
  variants[, gwas_status := "novel"]
  variants[, catalog_matched_disease_traits := NA_character_]

  variants[row_id %in% hit_rows, novel := FALSE]
  variants[row_id %in% hit_rows, gwas_status := "known_trait"]
  variants[trait_by_row, on = "row_id", catalog_matched_disease_traits := i.catalog_matched_disease_traits]
}

# Phenotypes absent from Excel were already dropped; n_efo_ids==0 remain novel
variants[n_efo_ids == 0L, novel := TRUE]
variants[n_efo_ids == 0L, gwas_status := "novel"]

fwrite(variants, output_path, sep = "\t")
message("Saved: ", output_path)
print(table(variants$gwas_status, useNA = "ifany"))
print(table(variants$novel, useNA = "ifany"))
