#!/usr/bin/env Rscript
# clump_annotate_lead_SAIGE.R
#
# Purpose:
#   Post-GWAS processing for SAIGE cumulative-incidence binary sumstats:
#   (1) PLINK clumping (p1=5e-8, r2=0.1, kb=1000)
#   (2) ANNOVAR gene annotation (genomic REF/ALT from bim_anno when available)
#   (3) Lead-variant table (+ nearest gene)
#
#   Allele convention: BETA / AF_Allele2 are for Allele2 (effect allele; EA).
#   Allele1 = other allele (NEA). Odds ratio = exp(BETA).
#   Clump/avinput read original *.SAIGE.txt; lead-table SE uses *.SAIGE.SEfixed.txt
#   (SE column = SPA-corrected SE from SAIGE.R).
#
# Inputs:
#   .../Result/SAIGE_cumulative_incidence/*_SAIGE_cumulative_inci.SAIGE.txt
#   .../Result/GWAS_sumstat/*_SAIGE_cumulative_inci.SAIGE.SEfixed.txt
#   .../Genotype/KBAv1_...RSupdated.{bed,bim,fam,.frq}
#   DATA_baseline/forGWAS_annot/bim_variantID_annot.rda
#
# Outputs:
#   .../Result/Clumping/*_saige_clump_5e-8_0.1_1000.{clumped,avinput,...}
#   .../Result/ANNOVAR/final_annotation_binary.tsv
#   DATA_manuscript_final/Extended_Data/Cate_lead_variants.tsv
#
# Downstream (GWAS Catalog novelty; not in this script):
#   CODE_Replication/forGWAS_annot/liftover_sumstat.R
#     -> Extended_Data/Cate_lead_variants_hg38.tsv
#   CODE_Replication/forGWAS_annot/gwas_catalog_cate_novelty_from_excel_efo.R
#
# Usage:
#   Rscript clump_annotate_lead_SAIGE.R
#   RUN_STEPS=clump,avinput Rscript clump_annotate_lead_SAIGE.R
#   RUN_STEPS=lead Rscript clump_annotate_lead_SAIGE.R
#
# RUN_STEPS (comma-separated; default = all):
#   clump, avinput, annovar, annotate, lead

library(data.table)
library(parallel)
library(dplyr)

parse_run_steps <- function() {
  raw <- trimws(Sys.getenv("RUN_STEPS", "all"))
  if (!nzchar(raw) || tolower(raw) == "all") {
    return(NULL)
  }
  unique(trimws(strsplit(raw, ",", fixed = TRUE)[[1]]))
}

step_enabled <- function(name, steps = parse_run_steps()) {
  is.null(steps) || name %in% steps
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
saige_dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/SAIGE_cumulative_incidence"
saige_sefixed_dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"
frq_file <- "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated.frq"
geno_file <- "/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated"
clump_dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping"
annot_file <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/ANNOVAR/final_annotation_binary.tsv"
result_dir <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"
bim_anno_file <- "/home2/sjim/openKoGES/DATA_baseline/forGWAS_annot/bim_variantID_annot.rda"

# -----------------------------------------------------------------------------
# Parallel settings
# -----------------------------------------------------------------------------
ncore_clump <- 3
ncore_avinput <- 3
ncore_annovar <- 10
ncore_lead_table <- 10
annovar_threads <- 3

# -----------------------------------------------------------------------------
# Clumping settings (PLINK): p < 5e-8, window 1 Mb, R2 = 0.1
# -----------------------------------------------------------------------------
clump_p1 <- 7.73e-9
clump_r2 <- 0.1
clump_kb <- 1000
clump_tag <- "7.73e-9_0.1_1000"

# -----------------------------------------------------------------------------
# Initial checks / directories
# -----------------------------------------------------------------------------
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annot_file), recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  saige_dir = saige_dir,
  saige_sefixed_dir = saige_sefixed_dir,
  frq_file = frq_file,
  geno_bed = paste0(geno_file, ".bed"),
  geno_bim = paste0(geno_file, ".bim"),
  geno_fam = paste0(geno_file, ".fam"),
  bim_anno_file = bim_anno_file
)
missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Required input path(s) not found:\n",
    paste(names(missing_paths), missing_paths, sep = " = ", collapse = "\n")
  )
}

pheno_bin <- list.files(
  saige_dir,
  pattern = "_SAIGE_cumulative_inci\\.SAIGE\\.txt$",
  full.names = FALSE
)

if (length(pheno_bin) == 0L) {
  stop("No SAIGE result files found in: ", saige_dir)
}

cat("Number of SAIGE phenotype files:", length(pheno_bin), "\n")
cat("SE-fixed sumstat dir:", saige_sefixed_dir, "\n")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
saige_sefixed_name <- function(p) {
  sub("\\.SAIGE\\.txt$", ".SAIGE.SEfixed.txt", p)
}

resolve_saige_sumstat <- function(p, use_sefixed = FALSE) {
  if (!isTRUE(use_sefixed)) {
    return(file.path(saige_dir, p))
  }
  sefixed <- file.path(saige_sefixed_dir, saige_sefixed_name(p))
  if (!file.exists(sefixed)) {
    stop(
      "SE-fixed SAIGE sumstat not found for ", p, ":\n  ", sefixed,
      "\nGenerate with SAIGE.R first."
    )
  }
  sefixed
}

clump_prefix <- function(p) {
  file.path(clump_dir, paste0(p, "_saige_clump_", clump_tag))
}

read_clumped <- function(path, phenotype) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    message("No .clumped file: ", phenotype)
    return(NULL)
  }

  out <- tryCatch(
    fread(path),
    error = function(e) {
      message("Could not read .clumped file for ", phenotype, ": ", e$message)
      NULL
    }
  )

  if (is.null(out) || nrow(out) == 0L || !"SNP" %in% names(out)) {
    message("No clumped index variants: ", phenotype)
    return(NULL)
  }

  out
}

variant_type_from_alleles <- function(a1, a2) {
  a1 <- as.character(a1)
  a2 <- as.character(a2)

  fifelse(
    nchar(a1) == 1L & nchar(a2) == 1L,
    "SNV",
    fifelse(nchar(a1) == nchar(a2), "MNV", "INDEL")
  )
}

load_bim_anno <- function(path) {
  cat("Loading bim_anno:", path, "\n")
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("bim_anno", envir = env, inherits = FALSE)) {
    stop("Object 'bim_anno' not found in: ", path)
  }

  ba <- as.data.table(env$bim_anno)
  required_cols <- c("rsID", "variantID", "BP", "ref", "alt")
  missing_cols <- setdiff(required_cols, names(ba))
  if (length(missing_cols) > 0L) {
    stop("bim_anno missing columns: ", paste(missing_cols, collapse = ", "))
  }

  ba <- ba[
    !is.na(rsID) & nzchar(as.character(rsID)) &
      !is.na(variantID) & nzchar(as.character(variantID)) &
      !is.na(ref) & nzchar(as.character(ref)) &
      !is.na(alt) & nzchar(as.character(alt)),
    .(
      rsID = as.character(rsID),
      variantID = as.character(variantID),
      BP_bim = as.integer(BP),
      Genomic_REF = as.character(ref),
      Genomic_ALT = as.character(alt)
    )
  ]

  if (anyDuplicated(ba$rsID)) {
    stop("bim_anno rsID is not unique after filtering.")
  }

  setkey(ba, rsID)
  cat("bim_anno variants with genomic REF/ALT:", nrow(ba), "\n")
  ba
}

# ANNOVAR interval from genomic REF/ALT (VCF-style left-aligned alleles).
# Falls back to POS-only placeholder when genomic alleles are unavailable.
annovar_coords_from_alleles <- function(pos, g_ref, g_alt, a1, a2, has_genomic) {
  pos <- as.integer(pos)

  if (isTRUE(has_genomic) && !is.na(g_ref) && !is.na(g_alt) &&
        nzchar(g_ref) && nzchar(g_alt)) {
    nref <- nchar(g_ref)
    nalt <- nchar(g_alt)

    if (nref == 1L && nalt == 1L) {
      return(list(
        START = pos, END = pos, AV_REF = g_ref, AV_ALT = g_alt,
        Allele_source = "bim_genomic"
      ))
    }
    if (nref < nalt) {
      return(list(
        START = pos, END = pos, AV_REF = "-",
        AV_ALT = substr(g_alt, nref + 1L, nalt),
        Allele_source = "bim_genomic"
      ))
    }
    if (nref > nalt) {
      return(list(
        START = as.integer(pos + 1L),
        END = as.integer(pos + (nref - nalt)),
        AV_REF = substr(g_ref, nalt + 1L, nref),
        AV_ALT = "-",
        Allele_source = "bim_genomic"
      ))
    }
    return(list(
      START = pos,
      END = as.integer(pos + nref - 1L),
      AV_REF = g_ref,
      AV_ALT = g_alt,
      Allele_source = "bim_genomic"
    ))
  }

  # POS-only fallback: do not treat SAIGE Allele1/Allele2 as genomic REF/ALT.
  if (nchar(a1) == 1L && nchar(a2) == 1L) {
    list(
      START = pos, END = pos, AV_REF = a1, AV_ALT = a2,
      Allele_source = "pos_only_fallback"
    )
  } else {
    list(
      START = pos, END = pos, AV_REF = "A", AV_ALT = "T",
      Allele_source = "pos_only_fallback"
    )
  }
}

bim_anno_dt <- NULL
if (step_enabled("avinput") || step_enabled("lead")) {
  bim_anno_dt <- load_bim_anno(bim_anno_file)
}

# =============================================================================
# 1. PLINK clumping
# =============================================================================
if (step_enabled("clump")) {
  clump_result <- mclapply(
    pheno_bin,
    function(p) {
      cat("CLUMP:", p, "\n")

      pheno_file <- file.path(saige_dir, p)
      out_prefix <- clump_prefix(p)

      args <- c(
        "--bfile", geno_file,
        "--clump", pheno_file,
        "--clump-snp-field", "MarkerID",
        "--clump-field", "p.value",
        "--clump-p1", format(clump_p1, scientific = TRUE),
        "--clump-r2", as.character(clump_r2),
        "--clump-kb", as.character(clump_kb),
        "--out", out_prefix
      )

      status <- system2("/bin2/plink", args = args)
      if (!identical(status, 0L)) {
        stop("PLINK clumping failed for ", p, " (exit status: ", status, ")")
      }

      invisible(out_prefix)
    },
    mc.cores = ncore_clump
  )

  invisible(clump_result)
}

# =============================================================================
# 2. Create ANNOVAR input
#
# Prefer genomic REF/ALT from bim_anno (unique rsID -> variantID/ref/alt).
# When unmatched, use POS-only placeholder alleles (indel/MNV: A/T).
# SAIGE Allele1/Allele2 are retained in Otherinfo as Other/Effect alleles.
# =============================================================================
make_avinput_saige <- function(p, saige_path, bim_anno) {
  cat("AVINPUT:", p, "\n")

  clump_file <- paste0(clump_prefix(p), ".clumped")
  sumstat_file <- file.path(saige_path, p)
  out_file <- paste0(clump_prefix(p), ".avinput")

  clump <- read_clumped(clump_file, p)
  if (is.null(clump)) {
    return(NULL)
  }

  required_cols <- c("MarkerID", "CHR", "POS", "Allele1", "Allele2")
  sumstat <- fread(sumstat_file, select = required_cols)
  setnames(sumstat, "MarkerID", "SNP")

  if (anyDuplicated(sumstat$SNP)) {
    duplicated_ids <- unique(sumstat$SNP[duplicated(sumstat$SNP)])
    stop(
      "Duplicated MarkerID in ", p, ". Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }

  lead <- unique(clump[, .(SNP, CHR_clump = CHR, BP_clump = BP)])
  dt <- merge(lead, sumstat, by = "SNP", all.x = TRUE, sort = FALSE)

  if (anyNA(dt$CHR) || anyNA(dt$POS)) {
    missing_ids <- dt[is.na(CHR) | is.na(POS), SNP]
    stop(
      "Clumped variants missing from SAIGE file ", p, ": ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }

  # Validate that PLINK clump coordinates and SAIGE coordinates agree.
  coord_mismatch <- dt[
    !is.na(CHR_clump) & !is.na(BP_clump) &
      (as.character(CHR_clump) != as.character(CHR) |
         as.integer(BP_clump) != as.integer(POS))
  ]
  if (nrow(coord_mismatch) > 0L) {
    stop(
      "Coordinate mismatch between PLINK clump and SAIGE result for ", p,
      ". First variant: ", coord_mismatch$SNP[1]
    )
  }

  dt[, `:=`(
    A1 = as.character(Allele1),
    A2 = as.character(Allele2),
    BP_original = as.integer(POS)
  )]

  dt <- merge(dt, bim_anno, by.x = "SNP", by.y = "rsID", all.x = TRUE, sort = FALSE)
  dt[, has_genomic := !is.na(Genomic_REF) & !is.na(Genomic_ALT)]

  bp_mismatch <- dt[
    has_genomic == TRUE & !is.na(BP_bim) & as.integer(BP_bim) != BP_original
  ]
  if (nrow(bp_mismatch) > 0L) {
    stop(
      "BP mismatch between SAIGE and bim_anno for ", p,
      ". First variant: ", bp_mismatch$SNP[1]
    )
  }

  allele_mismatch <- dt[
    has_genomic == TRUE &
      !(
        (A1 == Genomic_REF & A2 == Genomic_ALT) |
          (A1 == Genomic_ALT & A2 == Genomic_REF)
      )
  ]
  if (nrow(allele_mismatch) > 0L) {
    stop(
      "SAIGE Allele1/Allele2 do not match bim_anno REF/ALT for ", p,
      ". First variant: ", allele_mismatch$SNP[1]
    )
  }

  if (any(!dt$has_genomic)) {
    warning(
      "bim_anno genomic alleles missing for ",
      sum(!dt$has_genomic), " / ", nrow(dt),
      " clumped variants in ", p, "; using POS-only ANNOVAR fallback."
    )
  }

  dt[, Variant_type := fifelse(
    has_genomic,
    variant_type_from_alleles(Genomic_REF, Genomic_ALT),
    variant_type_from_alleles(A1, A2)
  )]

  coords <- rbindlist(
    Map(
      annovar_coords_from_alleles,
      dt$BP_original,
      dt$Genomic_REF,
      dt$Genomic_ALT,
      dt$A1,
      dt$A2,
      dt$has_genomic
    ),
    use.names = TRUE
  )
  dt[, `:=`(
    START = coords$START,
    END = coords$END,
    AV_REF = coords$AV_REF,
    AV_ALT = coords$AV_ALT,
    Allele_source = coords$Allele_source
  )]

  avinput <- dt[, .(
    CHR = as.character(CHR),
    START = as.integer(START),
    END = as.integer(END),
    AV_REF,
    AV_ALT,
    SNP,
    BP_original,
    A1,
    A2,
    Variant_type,
    Allele_source,
    variantID = fifelse(is.na(variantID), ".", variantID),
    Genomic_REF = fifelse(is.na(Genomic_REF), ".", Genomic_REF),
    Genomic_ALT = fifelse(is.na(Genomic_ALT), ".", Genomic_ALT)
  )]

  fwrite(
    avinput,
    out_file,
    sep = "\t",
    col.names = FALSE,
    quote = FALSE
  )

  cat("Saved:", out_file, "\n")
  invisible(out_file)
}

if (step_enabled("avinput")) {
  if (is.null(bim_anno_dt)) {
    stop("bim_anno_dt is required for the avinput step.")
  }
  avinput_result <- mclapply(
    pheno_bin,
    make_avinput_saige,
    saige_path = saige_dir,
    bim_anno = bim_anno_dt,
    mc.cores = ncore_avinput
  )
  invisible(avinput_result)
}

# =============================================================================
# 3. ANNOVAR gene annotation
# =============================================================================
if (step_enabled("annovar")) {
  annovar_result <- mclapply(
    pheno_bin,
    function(p) {
      cat("ANNOVAR:", p, "\n")

      out_prefix <- clump_prefix(p)
      avinput_file <- paste0(out_prefix, ".avinput")

      if (!file.exists(avinput_file) || file.info(avinput_file)$size == 0L) {
        message("Skipping ANNOVAR (no avinput): ", p)
        return(NULL)
      }

      args <- c(
        avinput_file,
        "/bin2/annovar/humandb/",
        "-buildver", "hg19",
        "-out", out_prefix,
        "-protocol", "refGene",
        "-operation", "g",
        "-nastring", ".",
        "-otherinfo",
        "--thread", as.character(annovar_threads)
      )

      status <- system2("/bin2/annovar/table_annovar.pl", args = args)
      if (!identical(status, 0L)) {
        stop("ANNOVAR failed for ", p, " (exit status: ", status, ")")
      }

      invisible(out_prefix)
    },
    mc.cores = ncore_annovar
  )
  invisible(annovar_result)
}

# =============================================================================
# 4. Post-process nearest-gene annotation
# =============================================================================
get_nearest_gene <- function(gene_str, detail_str) {
  if (is.na(gene_str) || gene_str == "." || gene_str == "") {
    return(NA_character_)
  }

  genes <- trimws(unlist(strsplit(gene_str, ",", fixed = TRUE)))

  if (is.na(detail_str) || detail_str == "." || detail_str == "") {
    return(genes[1])
  }

  dists <- trimws(unlist(strsplit(detail_str, ";", fixed = TRUE)))
  if (length(genes) != length(dists)) {
    return(genes[1])
  }

  dists_num <- suppressWarnings(as.numeric(sub("^dist=", "", dists)))
  if (all(is.na(dists_num))) {
    return(genes[1])
  }

  genes[which.min(abs(dists_num))]
}

clean_gene_non_intergenic <- function(x) {
  if (is.na(x) || x == "." || x == "") {
    return(NA_character_)
  }
  trimws(strsplit(x, "[,;]")[[1]][1])
}

if (step_enabled("annotate")) {
  res_list <- lapply(
    pheno_bin,
    function(p) {
      anno_file_i <- paste0(clump_prefix(p), ".hg19_multianno.txt")

      if (!file.exists(anno_file_i) || file.info(anno_file_i)$size == 0L) {
        message("Skipping annotation post-processing: ", p)
        return(NULL)
      }

      cat("POST-PROCESS ANNOTATION:", p, "\n")
      dt <- fread(anno_file_i)

      required_otherinfo <- paste0("Otherinfo", 1:9)
      if (!all(required_otherinfo %in% names(dt))) {
        stop(
          "Unexpected ANNOVAR Otherinfo columns in ", anno_file_i,
          ". Expected: ", paste(required_otherinfo, collapse = ", "),
          " (SNP, BP_original, Allele1, Allele2, Variant_type, Allele_source,",
          " variantID, Genomic_REF, Genomic_ALT).\n",
          "Re-run: RUN_STEPS=avinput,annovar,annotate Rscript clump_annotate_lead_SAIGE.R"
        )
      }

      dt[, `:=`(
        rsID = as.character(Otherinfo1),
        BP_original = as.integer(Otherinfo2),
        Allele1 = as.character(Otherinfo3),
        Allele2 = as.character(Otherinfo4),
        Variant_type = as.character(Otherinfo5),
        Allele_source = as.character(Otherinfo6),
        variantID = as.character(Otherinfo7),
        Genomic_REF = as.character(Otherinfo8),
        Genomic_ALT = as.character(Otherinfo9)
      )]

      dt[, Gene_clean := fifelse(
        Func.refGene == "intergenic" |
          is.na(Gene.refGene) | Gene.refGene == "." | Gene.refGene == "",
        mapply(get_nearest_gene, Gene.refGene, GeneDetail.refGene),
        vapply(Gene.refGene, clean_gene_non_intergenic, character(1))
      )]

      dt[, `:=`(
        Phenotype = p,
        Annotation_scope = fifelse(
          Allele_source == "bim_genomic",
          "Gene annotation with genomic REF/ALT from bim_anno",
          "POS-only gene annotation; indel/MNV allele consequence not interpreted"
        )
      )]

      dt
    }
  )

  res_list <- Filter(Negate(is.null), res_list)
  if (length(res_list) == 0L) {
    stop("No ANNOVAR annotation outputs were available to combine.")
  }

  final_anno <- rbindlist(res_list, fill = TRUE, use.names = TRUE)
  fwrite(final_anno, annot_file, sep = "\t", quote = FALSE)
  cat("Saved:", annot_file, "\n")
}

# =============================================================================
# 5. Lead-variant table
# =============================================================================
extract_lead_table_saige <- function(p, clump_dir, saige_dir, frq_dt, bim_anno) {
  
  cat("LEAD TABLE:", p, "\n")

  clump_file <- paste0(clump_prefix(p), ".clumped")
  sumstat_file <- resolve_saige_sumstat(p, use_sefixed = TRUE)
  cat("  sumstat (SE-fixed):", sumstat_file, "\n")

  clump <- read_clumped(clump_file, p)
  if (is.null(clump) || !file.exists(sumstat_file)) {
    return(NULL)
  }

  required_cols <- c(
    "MarkerID", "CHR", "POS", "Allele1", "Allele2",
    "BETA", "SE", "p.value", "AF_Allele2"
  )
  sumstat <- fread(sumstat_file, select = required_cols)
  setnames(sumstat, "MarkerID", "SNP")

  if (anyDuplicated(sumstat$SNP)) {
    stop("Duplicated MarkerID in SAIGE result: ", p)
  }

  lead <- unique(clump[, .(SNP, CHR_clump = CHR, BP_clump = BP)])
  dt <- merge(lead, sumstat, by = "SNP", all.x = TRUE, sort = FALSE)

  if (nrow(dt) == 0L) {
    return(NULL)
  }

  coord_mismatch <- dt[
    !is.na(CHR_clump) & !is.na(BP_clump) &
      (as.character(CHR_clump) != as.character(CHR) |
         as.integer(BP_clump) != as.integer(POS))
  ]
  if (nrow(coord_mismatch) > 0L) {
    stop("PLINK/SAIGE coordinate mismatch in lead table for: ", p)
  }

  dt <- merge(dt, frq_dt, by = "SNP", all.x = TRUE, sort = FALSE)
  dt <- merge(dt, bim_anno, by.x = "SNP", by.y = "rsID", all.x = TRUE, sort = FALSE)

  bp_mismatch <- dt[
    !is.na(Genomic_REF) & !is.na(BP_bim) & as.integer(BP_bim) != as.integer(POS)
  ]
  if (nrow(bp_mismatch) > 0L) {
    stop("BP mismatch between SAIGE and bim_anno in lead table for: ", p)
  }

  allele_mismatch <- dt[
    !is.na(Genomic_REF) & !is.na(Genomic_ALT) &
      !(
        (as.character(Allele1) == Genomic_REF & as.character(Allele2) == Genomic_ALT) |
          (as.character(Allele1) == Genomic_ALT & as.character(Allele2) == Genomic_REF)
      )
  ]
  if (nrow(allele_mismatch) > 0L) {
    stop(
      "SAIGE Allele1/Allele2 do not match bim_anno REF/ALT in lead table for: ",
      p, ". First variant: ", allele_mismatch$SNP[1]
    )
  }

  if (any(is.na(dt$Genomic_REF))) {
    warning(
      "bim_anno genomic alleles missing for ",
      sum(is.na(dt$Genomic_REF)), " / ", nrow(dt),
      " lead variants in ", p
    )
  }

  dt[, Variant_type := fifelse(
    !is.na(Genomic_REF) & !is.na(Genomic_ALT),
    variant_type_from_alleles(Genomic_REF, Genomic_ALT),
    variant_type_from_alleles(Allele1, Allele2)
  )]

  dt[, `:=`(
    Phenotype = p,
    `Lead variant` = SNP,
    Chromosome = CHR,
    `Base pair position` = POS,
    `Variant ID` = variantID,
    `Other allele` = Allele1,
    `Effect allele` = Allele2,
    `Reference Allele` = Genomic_REF,
    `Alternative Allele` = Genomic_ALT,
    Beta = BETA,
    `Odds ratio` = exp(BETA),
    `Standard error` = SE,
    `P-value` = p.value,
    `Minor allele frequency` = MAF,
    `Effect allele frequency` = AF_Allele2,
    `Variant type` = Variant_type
  )]

  dt[, .(
    Phenotype,
    `Lead variant`,
    `Variant ID`,
    Chromosome,
    `Base pair position`,
    `Variant type`,
    `Other allele`,
    `Effect allele`,
    `Reference Allele`,
    `Alternative Allele`,
    Beta,
    `Odds ratio`,
    `Standard error`,
    `P-value`,
    `Minor allele frequency`,
    `Effect allele frequency`
  )]
}

if (step_enabled("lead")) {
  if (is.null(bim_anno_dt)) {
    stop("bim_anno_dt is required for the lead-table step.")
  }
  if (!file.exists(annot_file)) {
    stop(
      "The combined annotation file is required for the lead-table step: ",
      annot_file,
      "\nRe-run: RUN_STEPS=avinput,annovar,annotate Rscript clump_annotate_lead_SAIGE.R",
      "\n(or run the complete pipeline with default RUN_STEPS=all)."
    )
  }

  missing_sefixed <- pheno_bin[
    !file.exists(file.path(saige_sefixed_dir, saige_sefixed_name(pheno_bin)))
  ]
  if (length(missing_sefixed) > 0L) {
    stop(
      "SE-fixed sumstat missing for ", length(missing_sefixed), " phenotype(s).\n",
      "Expected under ", saige_sefixed_dir, " as *.SAIGE.SEfixed.txt\n",
      "Examples: ", paste(head(saige_sefixed_name(missing_sefixed), 5L), collapse = ", ")
    )
  }
  cat(
    "Lead table will read corrected SE from ",
    length(pheno_bin), " *.SAIGE.SEfixed.txt files in ",
    saige_sefixed_dir, "\n",
    sep = ""
  )

  frq_dt <- fread(frq_file, select = c("SNP", "MAF"))
  if (anyDuplicated(frq_dt$SNP)) {
    warning("Duplicated SNP IDs found in frequency file; retaining first occurrence.")
    frq_dt <- unique(frq_dt, by = "SNP")
  }

  res <- mclapply(
    pheno_bin,
    extract_lead_table_saige,
    clump_dir = clump_dir,
    saige_dir = saige_dir,
    frq_dt = frq_dt,
    bim_anno = bim_anno_dt,
    mc.cores = ncore_lead_table
  )

  res <- Filter(Negate(is.null), res)
  if (length(res) == 0L) {
    stop("No clumped lead variants were available.")
  }

  final_bin <- rbindlist(res, fill = TRUE, use.names = TRUE)
  cat("Clumped index variants:", nrow(final_bin), "\n")

  pheno_map <- data.table(
    Pheno_code = c(
      "B_ALLER", "B_BRCA", "B_CKD", "B_CAD", "B_DM2", "B_HTN2",
      "B_LIP", "B_GT", "B_FLIV", "B_THY", "B_THYCA", "B_ULCER",
      "B_FRAC", "B_GCA", "B_NOI", "B_OSTE", "B_STOMUL", "B_CATA",
      "B_CLD", "B_CVA", "B_GASTRO", "B_GB", "B_PER", "B_UTCA"
    ),
    Phenotype_label = c(
      "Allergic disease",
      "Breast cancer",
      "Chronic kidney disease (CKD)",
      "Coronary artery disease",
      "Diabetes mellitus",
      "Hypertension",
      "Hyperlipidemia",
      "Gout",
      "Fatty liver",
      "Thyroid disease",
      "Thyroid cancer",
      "Peptic ulcer",
      "Fracture",
      "Gastric cancer",
      "Neurosis",
      "Osteoporosis",
      "Chronic gastritis / Peptic ulcer",
      "Cataract",
      "Chronic lung disease",
      "Cerebrovascular disease",
      "Chronic gastritis",
      "Gallbladder disease",
      "Periodontal Disease",
      "Uterine cancer"
    )
  )
  pheno_map[, Phenotype := paste0(Pheno_code, "_SAIGE_cumulative_inci.SAIGE.txt")]

  final_anno <- fread(annot_file)
  required_anno_cols <- c(
    "rsID", "BP_original", "Gene_clean", "Variant_type", "Annotation_scope"
  )
  if (!all(required_anno_cols %in% names(final_anno))) {
    stop(
      "Annotation file does not contain the expected revised columns: ",
      paste(setdiff(required_anno_cols, names(final_anno)), collapse = ", "),
      "\nRe-run: RUN_STEPS=avinput,annovar,annotate Rscript clump_annotate_lead_SAIGE.R"
    )
  }

  anno_sub_unique <- unique(
    final_anno[, .(
      rsID,
      BP_original,
      Gene_clean,
      Annotation_variant_type = Variant_type,
      Annotation_scope
    )],
    by = c("rsID", "BP_original")
  )

  final_bin <- merge(
    final_bin,
    anno_sub_unique,
    by.x = c("Lead variant", "Base pair position"),
    by.y = c("rsID", "BP_original"),
    all.x = TRUE,
    sort = FALSE
  )
  setnames(final_bin, "Gene_clean", "Nearest gene")

  type_mismatch <- final_bin[
    !is.na(Annotation_variant_type) &
      `Variant type` != Annotation_variant_type
  ]
  if (nrow(type_mismatch) > 0L) {
    warning("Variant-type mismatch found for ", nrow(type_mismatch), " lead variants.")
  }
  final_bin[, Annotation_variant_type := NULL]

  final_bin <- merge(
    final_bin,
    pheno_map[, .(Phenotype, Phenotype_label)],
    by = "Phenotype",
    all.x = TRUE,
    sort = FALSE
  )

  gwascata <- final_bin %>%
    select(
      Phenotype_label,
      Phenotype,
      `Lead variant`,
      `Variant ID`,
      Chromosome,
      `Base pair position`,
      `Variant type`,
      `Nearest gene`,
      `Other allele`,
      `Effect allele`,
      `Reference Allele`,
      `Alternative Allele`,
      Beta,
      `Odds ratio`,
      `Standard error`,
      `P-value`,
      `Minor allele frequency`,
      `Effect allele frequency`
    )

  colnames(gwascata) <- c(
    "Phenotype", "VarNM", "Lead variant", "VariantID", "CHR", "BP",
    "Variant_type", "Nearest_gene", "NEA", "EA", "REF", "ALT",
    "Beta", "OR", "SE", "P", "MAF", "EAF"
  )

  catalogue_file <- file.path(result_dir, "Cate_lead_variants.tsv")
  fwrite(gwascata, catalogue_file, sep = "\t", quote = FALSE)
  cat("Saved:", catalogue_file, "\n")
}
