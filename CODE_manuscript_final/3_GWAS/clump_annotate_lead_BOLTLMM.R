#!/usr/bin/env Rscript
# clump_annotate_lead_BOLTLMM.R
#
# Purpose:
#   Post-GWAS processing for BOLT-LMM continuous-trait sumstats:
#   (1) PLINK clumping (p1=5e-8, r2=0.1, kb=1000)
#   (2) ANNOVAR gene annotation (genomic REF/ALT from bim_anno when available)
#   (3) Lead-variant table (+ nearest gene)
#
#   Allele convention: BETA / A1FREQ are for ALLELE1 (effect allele; EA).
#   ALLELE0 = other allele (NEA). No odds ratio (continuous traits).
#   Association P: P_BOLT_LMM if present, else P_BOLT_LMM_INF.
#
# Inputs:
#   .../Result/GWAS_sumstat/*_BOLT_result
#   .../Genotype/KBAv1_...RSupdated.{bed,bim,fam,.frq}
#   DATA_baseline/forGWAS_annot/bim_variantID_annot.rda
#
# Outputs:
#   .../Result/Clumping_BOLT/*_bolt_clump_5e-8_0.1_1000.{clumped,avinput,...}
#   .../Result/ANNOVAR/final_annotation_continuous.tsv
#   DATA_manuscript_final/Extended_Data/Cont_lead_variants.tsv
#
# Downstream (GWAS Catalog novelty; not in this script):
#   CODE_Replication/forGWAS_annot/liftover_sumstat.R
#     -> Extended_Data/Cont_lead_variants_hg38.tsv
#   CODE_Replication/forGWAS_annot/gwas_catalog_cate_novelty_from_excel_efo.R
#
# Usage:
#   Rscript clump_annotate_lead_BOLTLMM.R
#   RUN_STEPS=clump,avinput Rscript clump_annotate_lead_BOLTLMM.R
#   RUN_STEPS=lead Rscript clump_annotate_lead_BOLTLMM.R
#
# RUN_STEPS (comma-separated; default = all):
#   clump, avinput, annovar, annotate, lead

library(data.table)
library(parallel)
library(dplyr)
pheno_map <- data.table(
  Pheno_code = c(
    "B_FEV1best_noNC_BOLT_result", "B_FFRbest_noNC_BOLT_result",  "B_FVCbest_noNC_BOLT_result", 
    "RIN_B_ALBUMIN_BOLT_result",   "RIN_B_ALP_BOLT_result",       "RIN_B_ALT_BOLT_result",      
    "RIN_B_APOA_BOLT_result",      "RIN_B_APOB_BOLT_result",      "RIN_B_AST_BOLT_result",      
    "RIN_B_BMI_BOLT_result",       "RIN_B_BUN_BOLT_result",       "RIN_B_CPEPTIDE_BOLT_result", 
    "RIN_B_CREATINE_BOLT_result",  "RIN_B_CYSTATINC_BOLT_result", "RIN_B_DBP_BOLT_result",      
    "RIN_B_eGFR_BOLT_result",      "RIN_B_FOLATE_BOLT_result",    "RIN_B_GLU0_BOLT_result",     
    "RIN_B_HB_BOLT_result",        "RIN_B_HBA1C_BOLT_result",     "RIN_B_HCT_BOLT_result",      
    "RIN_B_HDL_BOLT_result",       "RIN_B_HEIGHT_BOLT_result",    "RIN_B_HIP_BOLT_result",      
    "RIN_B_HOMOCYST_BOLT_result",  "RIN_B_INSUL_BOLT_result",     "RIN_B_IRON_BOLT_result",     
    "RIN_B_LDLclean_BOLT_result",  "RIN_B_MCH_BOLT_result",       "RIN_B_MCHC_BOLT_result",     
    "RIN_B_MCV_BOLT_result",       "RIN_B_PHOSPHATE_BOLT_result", "RIN_B_PLAT_BOLT_result",     
    "RIN_B_PULSE_BOLT_result",     "RIN_B_R_GTP_BOLT_result",     "RIN_B_RBC_B_BOLT_result",    
    "RIN_B_SBP_BOLT_result",       "RIN_B_T_BIL_BOLT_result",    "RIN_B_TCHL_BOLT_result",     
    "RIN_B_TG_BOLT_result",        "RIN_B_URICACID_BOLT_result",  "RIN_B_VITD_BOLT_result",     
    "RIN_B_WAIST_BOLT_result",     "RIN_B_WBC_B_BOLT_result",     "RIN_B_WEIGHT_BOLT_result",   
    "RIN_B_WHR_BOLT_result"
  ),
  Phenotype_label = c(
    "FEV1",
    "FEV1/FVC", 
    "FVC", 
    "Albumin (Alb)",   
    "Alkaline phosphatase (ALP)",      
    "Alanine aminotransferase (ALT)",      
    "Apolipoprotein A",      
    "Apolipoprotein B",      
    "Aspartate transaminase (AST)",      
    "Body mass index (Body mass index (BMI))",
    "Blood urea nitrogen (BUN)",
    "C-peptide", 
    "Serum creatinine (sCr)", 
    "Cystatin C", 
    "Diastolic blood pressure (DBP)",      
    "eGFR (CKD-EPI)",      
    "Folate",    
    "Glucose",     
    "Hemoglobin (Hb)",        
    "HbA1C",     
    "Hematocrit (HCT)",      
    "HDL-cholesterol (HDL-C)",
    "Body height",   
    "Hip",      
    "Homocysteine",  
    "Insulin",     
    "Iron",     
    "LDL-cholesterol (LDL-C)",
    "Mean corpuscular hemoglobin (MCH)",
    "Mean corpuscular hemoglobin concentration (MCHC)",     
    "Mean corpuscular volume (MCV)",
    "Phosphate", 
    "Platelet count (Plt)",     
    "Pulse", 
    "r_Gtp (Gamma glutamyl transferase)",    
    "Red blood cell count (RBC)",    
    "Systolic blood pressure (SBP)", 
    "Total Bilirubin (TBil)",    
    "Total cholesterol (TC)",     
    "Triglycerides (TG)",        
    "Uric acids (UA)",  
    "Vit.D",     
    "Waist",     "White blood cell count (WBC)",     "Body weight (BW)",   
    "Waist-to-hip ratio (WHR)"
  )
)

parse_run_steps <- function() {
  raw <- trimws(Sys.getenv("RUN_STEPS", "all"))
  if (!nzchar(raw) || tolower(raw) == "all") {
    return(NULL)
  }
  
  steps <- unique(trimws(strsplit(raw, ",", fixed = TRUE)[[1]]))
  valid_steps <- c("clump", "avinput", "annovar", "annotate", "lead")
  invalid_steps <- setdiff(steps, valid_steps)
  
  if (length(invalid_steps) > 0L) {
    stop(
      "Unknown RUN_STEPS value(s): ",
      paste(invalid_steps, collapse = ", "),
      "\nAllowed values: ", paste(valid_steps, collapse = ", ")
    )
  }
  
  steps
}

run_steps <- parse_run_steps()

step_enabled <- function(name) {
  is.null(run_steps) || name %in% run_steps
}

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
bolt_dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/GWAS_sumstat"

frq_file <- paste0(
  "/data6/OpenGWAS/KBA160K/Genotype/",
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated.frq"
)

geno_file <- paste0(
  "/data6/OpenGWAS/KBA160K/Genotype/",
  "KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated"
)

clump_dir <- "/data6/OpenGWAS/KBA160K_longitudinal/Result/Clumping_BOLT"
annot_file <- paste0(
  "/data6/OpenGWAS/KBA160K_longitudinal/Result/ANNOVAR/",
  "final_annotation_bolt.tsv"
)
result_dir <- "/home2/sjim/openKoGES/DATA_manuscript_final/Extended_Data"

bim_anno_file <- paste0(
  "/home2/sjim/openKoGES/DATA_baseline/forGWAS_annot/",
  "bim_variantID_annot.rda"
)

plink_bin <- "/bin2/plink"
annovar_bin <- "/bin2/annovar/table_annovar.pl"
annovar_db <- "/bin2/annovar/humandb/"

# -----------------------------------------------------------------------------
# Parallel settings
# -----------------------------------------------------------------------------
ncore_clump <- 3
ncore_avinput <- 3
ncore_annovar <- 10
ncore_lead_table <- 10
annovar_threads <- 3

# -----------------------------------------------------------------------------
# Clumping settings: p < 5e-8, window 1 Mb, r2 = 0.1
# -----------------------------------------------------------------------------
clump_p1 <- 7.73e-9
clump_r2 <- 0.1
clump_kb <- 1000
clump_tag <- "7.73e-9_0.1_1000"

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
phenotype_code <- function(p) {
  sub("_BOLT_result$", "", basename(p))
}

read_bolt_header <- function(path) {
  line <- readLines(path, n = 1L, warn = FALSE)
  if (length(line) == 0L || !nzchar(line)) {
    stop("Could not read header from: ", path)
  }
  strsplit(line, "\t", fixed = TRUE)[[1]]
}

# Phenotypes whose BOLT result files lack an N column.
# Keys are full filenames under bolt_dir (*_BOLT_result).
bolt_known_n <- c(
  "B_FEV1best_noNC_BOLT_result" = 7714,
  "B_FEV1pp_noNC_BOLT_result" = 7714,
  "B_FFRbest_noNC_BOLT_result" = 7714,
  "B_FVCbest_noNC_BOLT_result" = 7714,
  "B_FVCpp_noNC_BOLT_result" = 7714,
  "RIN_B_LDLclean_BOLT_result" = 146583,
  "RIN_B_HOMOCYST_BOLT_result" = 43534
)

resolve_bolt_n <- function(phenotype, has_n_column) {
  if (isTRUE(has_n_column)) {
    return(list(n = NA_real_, source = "BOLT_file"))
  }

  pheno_file <- basename(phenotype)
  if (pheno_file %in% names(bolt_known_n)) {
    return(list(n = unname(bolt_known_n[[pheno_file]]), source = "manual_known_n"))
  }

  stop(
    "BOLT file ", pheno_file, " has no N column and no entry in bolt_known_n. ",
    "Add the sample size to bolt_known_n or regenerate the BOLT output with N."
  )
}

# Prefer non-infinitesimal P when BOLT wrote it; else infinitesimal P.
resolve_bolt_p_field <- function(cols, phenotype = NULL) {
  cols <- as.character(cols)
  if ("P_BOLT_LMM" %in% cols) {
    return("P_BOLT_LMM")
  }
  if ("P_BOLT_LMM_INF" %in% cols) {
    return("P_BOLT_LMM_INF")
  }
  stop(
    "No usable BOLT P column found",
    if (!is.null(phenotype)) paste0(" in ", phenotype) else "",
    ". Expected P_BOLT_LMM or P_BOLT_LMM_INF."
  )
}

read_bolt_sumstat_for_lead <- function(path, phenotype) {
  header <- read_bolt_header(path)
  p_field <- resolve_bolt_p_field(header, phenotype)

  base_cols <- c(
    "SNP", "CHR", "BP", "ALLELE1", "ALLELE0",
    "BETA", "SE", "A1FREQ", p_field
  )
  missing_base <- setdiff(base_cols, header)
  if (length(missing_base) > 0L) {
    stop(
      "BOLT file ", phenotype, " missing required columns: ",
      paste(missing_base, collapse = ", ")
    )
  }

  select_cols <- base_cols
  has_n <- "N" %in% header
  if (has_n) {
    select_cols <- c(select_cols, "N")
  }

  sumstat <- fread(path, select = select_cols)
  setnames(sumstat, p_field, "ASSOC_P")

  n_info <- resolve_bolt_n(phenotype, has_n)
  if (!has_n) {
    sumstat[, N := as.numeric(n_info$n)]
  }

  list(
    sumstat = sumstat,
    p_field = p_field,
    has_n = has_n,
    n_source = n_info$source,
    n_value = if (has_n) NA_real_ else as.numeric(n_info$n)
  )
}

clump_prefix <- function(p) {
  file.path(
    clump_dir,
    paste0(phenotype_code(p), "_bolt_clump_", clump_tag)
  )
}

read_clumped <- function(path, phenotype) {
  if (!file.exists(path) || file.info(path)$size == 0L) {
    message("No .clumped file: ", phenotype)
    return(NULL)
  }
  
  out <- tryCatch(
    fread(path),
    error = function(e) {
      message(
        "Could not read .clumped file for ", phenotype, ": ", e$message
      )
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
    duplicated_ids <- unique(ba$rsID[duplicated(ba$rsID)])
    stop(
      "bim_anno rsID is not unique after filtering. Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }
  
  setkey(ba, rsID)
  cat("bim_anno variants with genomic REF/ALT:", nrow(ba), "\n")
  ba
}

# ANNOVAR interval from genomic REF/ALT (VCF-style alleles).
# Falls back to position-only placeholder alleles when genomic alleles are absent.
annovar_coords_from_alleles <- function(
  pos, g_ref, g_alt, other_allele, effect_allele, has_genomic) {
  pos <- as.integer(pos)
  
  if (
    isTRUE(has_genomic) &&
    !is.na(g_ref) && !is.na(g_alt) &&
    nzchar(g_ref) && nzchar(g_alt)
  ) {
    nref <- nchar(g_ref)
    nalt <- nchar(g_alt)
    
    if (nref == 1L && nalt == 1L) {
      return(list(
        START = pos,
        END = pos,
        AV_REF = g_ref,
        AV_ALT = g_alt,
        Allele_source = "bim_genomic"
      ))
    }
    
    if (nref < nalt) {
      return(list(
        START = pos,
        END = pos,
        AV_REF = "-",
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
  
  # Position-only fallback. Do not interpret ALLELE0/ALLELE1 as genomic REF/ALT.
  if (
    !is.na(other_allele) && !is.na(effect_allele) &&
    nchar(other_allele) == 1L && nchar(effect_allele) == 1L
  ) {
    return(list(
      START = pos,
      END = pos,
      AV_REF = other_allele,
      AV_ALT = effect_allele,
      Allele_source = "pos_only_fallback"
    ))
  }
  
  list(
    START = pos,
    END = pos,
    AV_REF = "A",
    AV_ALT = "T",
    Allele_source = "pos_only_fallback"
  )
}

check_bolt_alleles_against_bim <- function(dt, phenotype, context) {
  allele_mismatch <- dt[
    !is.na(Genomic_REF) & !is.na(Genomic_ALT) &
      !(
        (Other_allele == Genomic_REF & Effect_allele == Genomic_ALT) |
          (Other_allele == Genomic_ALT & Effect_allele == Genomic_REF)
      )
    ]
  
  if (nrow(allele_mismatch) > 0L) {
    stop(
      "BOLT ALLELE0/ALLELE1 do not match bim_anno REF/ALT during ",
      context, " for ", phenotype,
      ". First variant: ", allele_mismatch$SNP[1],
      " (ALLELE0=", allele_mismatch$Other_allele[1],
      ", ALLELE1=", allele_mismatch$Effect_allele[1],
      ", REF=", allele_mismatch$Genomic_REF[1],
      ", ALT=", allele_mismatch$Genomic_ALT[1], ")"
    )
  }
  
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# Initial checks / directories
# -----------------------------------------------------------------------------
dir.create(clump_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(annot_file), recursive = TRUE, showWarnings = FALSE)
dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)

required_paths <- c(
  bolt_dir = bolt_dir,
  frq_file = frq_file,
  geno_bed = paste0(geno_file, ".bed"),
  geno_bim = paste0(geno_file, ".bim"),
  geno_fam = paste0(geno_file, ".fam"),
  bim_anno_file = bim_anno_file,
  plink_bin = plink_bin,
  annovar_bin = annovar_bin,
  annovar_db = annovar_db
)

missing_paths <- required_paths[!file.exists(required_paths)]
if (length(missing_paths) > 0L) {
  stop(
    "Required input path(s) not found:\n",
    paste(names(missing_paths), missing_paths, sep = " = ", collapse = "\n")
  )
}

pheno_bolt <- list.files(
  bolt_dir,
  pattern = "_BOLT_result$",
  full.names = FALSE
)

if (length(pheno_bolt) == 0L) {
  stop("No *_BOLT_result files found in: ", bolt_dir)
}

cat("Number of BOLT phenotype files:", length(pheno_bolt), "\n")
cat("First files:\n")
print(head(pheno_bolt, 10L))

bim_anno_dt <- NULL
if (step_enabled("avinput") || step_enabled("lead")) {
  bim_anno_dt <- load_bim_anno(bim_anno_file)
}

# =============================================================================
# 1. PLINK clumping
#
# BOLT file fields:
#   SNP field = SNP
#   association P field = P_BOLT_LMM if present, else P_BOLT_LMM_INF
# =============================================================================
if (step_enabled("clump")) {
  clump_result <- mclapply(
    pheno_bolt,
    function(p) {
      cat("CLUMP:", p, "\n")

      pheno_file <- file.path(bolt_dir, p)
      out_prefix <- clump_prefix(p)
      p_field <- resolve_bolt_p_field(read_bolt_header(pheno_file), p)
      cat("  using association P field:", p_field, "\n")

      args <- c(
        "--bfile", geno_file,
        "--clump", pheno_file,
        "--clump-snp-field", "SNP",
        "--clump-field", p_field,
        "--clump-p1", format(clump_p1, scientific = TRUE),
        "--clump-r2", as.character(clump_r2),
        "--clump-kb", as.character(clump_kb),
        "--out", out_prefix
      )

      status <- system2(plink_bin, args = args)
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
# BOLT convention used here:
# - ALLELE1 = effect allele (EA)
# - ALLELE0 = other/non-effect allele (NEA)
# - BETA and A1FREQ refer to ALLELE1
# =============================================================================
make_avinput_bolt <- function(p, bolt_path, bim_anno) {
  cat("AVINPUT:", p, "\n")
  
  clump_file <- paste0(clump_prefix(p), ".clumped")
  sumstat_file <- file.path(bolt_path, p)
  out_file <- paste0(clump_prefix(p), ".avinput")
  
  clump <- read_clumped(clump_file, p)
  if (is.null(clump)) {
    return(NULL)
  }
  
  required_cols <- c("SNP", "CHR", "BP", "ALLELE1", "ALLELE0")
  sumstat <- fread(sumstat_file, select = required_cols)
  
  if (anyDuplicated(sumstat$SNP)) {
    duplicated_ids <- unique(sumstat$SNP[duplicated(sumstat$SNP)])
    stop(
      "Duplicated SNP IDs in ", p, ". Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }
  
  lead <- unique(clump[, .(SNP, CHR_clump = CHR, BP_clump = BP)])
  dt <- merge(lead, sumstat, by = "SNP", all.x = TRUE, sort = FALSE)
  
  if (anyNA(dt$CHR) || anyNA(dt$BP)) {
    missing_ids <- dt[is.na(CHR) | is.na(BP), SNP]
    stop(
      "Clumped variants missing from BOLT file ", p, ": ",
      paste(head(missing_ids, 10L), collapse = ", ")
    )
  }
  
  coord_mismatch <- dt[
    !is.na(CHR_clump) & !is.na(BP_clump) &
      (
        as.character(CHR_clump) != as.character(CHR) |
          as.integer(BP_clump) != as.integer(BP)
      )
    ]
  
  if (nrow(coord_mismatch) > 0L) {
    stop(
      "Coordinate mismatch between PLINK clump and BOLT result for ", p,
      ". First variant: ", coord_mismatch$SNP[1]
    )
  }
  
  dt[, `:=`(
    Other_allele = as.character(ALLELE0),
    Effect_allele = as.character(ALLELE1),
    BP_original = as.integer(BP)
  )]
  
  dt <- merge(
    dt,
    bim_anno,
    by.x = "SNP",
    by.y = "rsID",
    all.x = TRUE,
    sort = FALSE
  )
  
  dt[, has_genomic := !is.na(Genomic_REF) & !is.na(Genomic_ALT)]
  
  bp_mismatch <- dt[
    has_genomic == TRUE &
      !is.na(BP_bim) &
      as.integer(BP_bim) != BP_original
    ]
  
  if (nrow(bp_mismatch) > 0L) {
    stop(
      "BP mismatch between BOLT and bim_anno for ", p,
      ". First variant: ", bp_mismatch$SNP[1]
    )
  }
  
  check_bolt_alleles_against_bim(dt, p, context = "avinput creation")
  
  if (any(!dt$has_genomic)) {
    warning(
      "bim_anno genomic alleles missing for ",
      sum(!dt$has_genomic), " / ", nrow(dt),
      " clumped variants in ", p,
      "; using position-only ANNOVAR fallback."
    )
  }
  
  dt[, Variant_type := fifelse(
    has_genomic,
    variant_type_from_alleles(Genomic_REF, Genomic_ALT),
    variant_type_from_alleles(Other_allele, Effect_allele)
  )]
  
  coords <- rbindlist(
    Map(
      annovar_coords_from_alleles,
      dt$BP_original,
      dt$Genomic_REF,
      dt$Genomic_ALT,
      dt$Other_allele,
      dt$Effect_allele,
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
    Other_allele,
    Effect_allele,
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
    pheno_bolt,
    make_avinput_bolt,
    bolt_path = bolt_dir,
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
    pheno_bolt,
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
        annovar_db,
        "-buildver", "hg19",
        "-out", out_prefix,
        "-protocol", "refGene",
        "-operation", "g",
        "-nastring", ".",
        "-otherinfo",
        "--thread", as.character(annovar_threads)
      )
      
      status <- system2(annovar_bin, args = args)
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
    pheno_bolt,
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
          " (SNP, BP_original, Other_allele, Effect_allele, Variant_type, ",
          "Allele_source, variantID, Genomic_REF, Genomic_ALT).\n",
          "Re-run: RUN_STEPS=avinput,annovar,annotate ",
          "Rscript clump_annotate_lead_BOLTLMM.R"
        )
      }
      
      dt[, `:=`(
        rsID = as.character(Otherinfo1),
        BP_original = as.integer(Otherinfo2),
        Other_allele = as.character(Otherinfo3),
        Effect_allele = as.character(Otherinfo4),
        Variant_type = as.character(Otherinfo5),
        Allele_source = as.character(Otherinfo6),
        variantID = as.character(Otherinfo7),
        Genomic_REF = as.character(Otherinfo8),
        Genomic_ALT = as.character(Otherinfo9)
      )]
      
      dt[, Gene_clean := fifelse(
        Func.refGene == "intergenic" |
          is.na(Gene.refGene) |
          Gene.refGene == "." |
          Gene.refGene == "",
        mapply(get_nearest_gene, Gene.refGene, GeneDetail.refGene),
        vapply(Gene.refGene, clean_gene_non_intergenic, character(1))
      )]
      
      dt[, `:=`(
        Phenotype = phenotype_code(p),
        Phenotype_file = p,
        Annotation_scope = fifelse(
          Allele_source == "bim_genomic",
          "Gene annotation with genomic REF/ALT from bim_anno",
          paste0(
            "Position-only gene annotation; indel/MNV allele consequence ",
            "not interpreted"
          )
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
extract_lead_table_bolt <- function(p, bolt_path, frq_dt, bim_anno) {
  cat("LEAD TABLE:", p, "\n")

  clump_file <- paste0(clump_prefix(p), ".clumped")
  sumstat_file <- file.path(bolt_path, p)

  clump <- read_clumped(clump_file, p)
  if (is.null(clump) || !file.exists(sumstat_file)) {
    return(NULL)
  }

  bolt_read <- read_bolt_sumstat_for_lead(sumstat_file, p)
  sumstat <- bolt_read$sumstat
  cat("  using association P field:", bolt_read$p_field, "\n")
  if (!bolt_read$has_n) {
    message(
      "  N column absent in ", p,
      "; using bolt_known_n = ", bolt_read$n_value,
      " (source: ", bolt_read$n_source, ")."
    )
  }

  if (anyDuplicated(sumstat$SNP)) {
    duplicated_ids <- unique(sumstat$SNP[duplicated(sumstat$SNP)])
    stop(
      "Duplicated SNP IDs in BOLT result ", p, ". Examples: ",
      paste(head(duplicated_ids, 10L), collapse = ", ")
    )
  }

  lead <- unique(clump[, .(SNP, CHR_clump = CHR, BP_clump = BP)])
  dt <- merge(lead, sumstat, by = "SNP", all.x = TRUE, sort = FALSE)

  if (nrow(dt) == 0L) {
    return(NULL)
  }

  missing_from_sumstat <- dt[is.na(CHR) | is.na(BP), SNP]
  if (length(missing_from_sumstat) > 0L) {
    stop(
      "Clumped variants missing from BOLT result ", p, ": ",
      paste(head(missing_from_sumstat, 10L), collapse = ", ")
    )
  }

  coord_mismatch <- dt[
    !is.na(CHR_clump) & !is.na(BP_clump) &
      (
        as.character(CHR_clump) != as.character(CHR) |
          as.integer(BP_clump) != as.integer(BP)
      )
    ]

  if (nrow(coord_mismatch) > 0L) {
    stop(
      "PLINK/BOLT coordinate mismatch in lead table for ", p,
      ". First variant: ", coord_mismatch$SNP[1]
    )
  }

  dt[, `:=`(
    Other_allele = as.character(ALLELE0),
    Effect_allele = as.character(ALLELE1)
  )]

  dt <- merge(dt, frq_dt, by = "SNP", all.x = TRUE, sort = FALSE)
  dt <- merge(
    dt,
    bim_anno,
    by.x = "SNP",
    by.y = "rsID",
    all.x = TRUE,
    sort = FALSE
  )

  bp_mismatch <- dt[
    !is.na(Genomic_REF) &
      !is.na(BP_bim) &
      as.integer(BP_bim) != as.integer(BP)
    ]

  if (nrow(bp_mismatch) > 0L) {
    stop(
      "BP mismatch between BOLT and bim_anno in lead table for ", p,
      ". First variant: ", bp_mismatch$SNP[1]
    )
  }

  check_bolt_alleles_against_bim(dt, p, context = "lead-table creation")

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
    variant_type_from_alleles(Other_allele, Effect_allele)
  )]

  dt[, `:=`(
    Phenotype = phenotype_code(p),
    Phenotype_file = p,
    `Lead variant` = SNP,
    Chromosome = CHR,
    `Base pair position` = BP,
    `Variant ID` = variantID,
    `Other allele` = Other_allele,
    `Effect allele` = Effect_allele,
    `Reference allele` = Genomic_REF,
    `Alternative allele` = Genomic_ALT,
    Beta = BETA,
    `Standard error` = SE,
    `P-value` = ASSOC_P,
    `Association P field` = bolt_read$p_field,
    `Minor allele frequency` = MAF,
    `Effect allele frequency` = A1FREQ,
    `Sample size` = N,
    `Variant type` = Variant_type
  )]

  dt[, .(
    Phenotype,
    Phenotype_file,
    `Lead variant`,
    `Variant ID`,
    Chromosome,
    `Base pair position`,
    `Variant type`,
    `Other allele`,
    `Effect allele`,
    `Reference allele`,
    `Alternative allele`,
    Beta,
    `Standard error`,
    `P-value`,
    `Association P field`,
    `Minor allele frequency`,
    `Effect allele frequency`,
    `Sample size`
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
      "\nRe-run: RUN_STEPS=avinput,annovar,annotate ",
      "Rscript clump_annotate_lead_BOLTLMM.R",
      "\nOr run the complete pipeline with RUN_STEPS=all."
    )
  }
  
  frq_dt <- fread(frq_file, select = c("SNP", "MAF"))
  if (anyDuplicated(frq_dt$SNP)) {
    warning("Duplicated SNP IDs found in frequency file; retaining first occurrence.")
    frq_dt <- unique(frq_dt, by = "SNP")
  }
  
  res <- mclapply(
    pheno_bolt,
    extract_lead_table_bolt,
    bolt_path = bolt_dir,
    frq_dt = frq_dt,
    bim_anno = bim_anno_dt,
    mc.cores = ncore_lead_table
  )
  
  res <- Filter(Negate(is.null), res)
  if (length(res) == 0L) {
    stop("No clumped BOLT lead variants were available.")
  }
  
  final_bolt <- rbindlist(res, fill = TRUE, use.names = TRUE)
  cat("Clumped index variants:", nrow(final_bolt), "\n")
  
  final_anno <- fread(annot_file)
  required_anno_cols <- c(
    "Phenotype", "rsID", "BP_original", "Gene_clean",
    "Variant_type", "Annotation_scope"
  )
  
  if (!all(required_anno_cols %in% names(final_anno))) {
    stop(
      "Annotation file does not contain the expected revised columns: ",
      paste(setdiff(required_anno_cols, names(final_anno)), collapse = ", "),
      "\nRe-run: RUN_STEPS=avinput,annovar,annotate ",
      "Rscript clump_annotate_lead_BOLTLMM.R"
    )
  }
  
  anno_sub_unique <- unique(
    final_anno[, .(
      Phenotype,
      rsID,
      BP_original,
      Gene_clean,
      Annotation_variant_type = Variant_type,
      Annotation_scope
    )],
    by = c("Phenotype", "rsID", "BP_original")
  )
  
  final_bolt <- merge(
    final_bolt,
    anno_sub_unique,
    by.x = c("Phenotype", "Lead variant", "Base pair position"),
    by.y = c("Phenotype", "rsID", "BP_original"),
    all.x = TRUE,
    sort = FALSE
  )
  
  setnames(final_bolt, "Gene_clean", "Nearest gene")
  
  type_mismatch <- final_bolt[
    !is.na(Annotation_variant_type) &
      `Variant type` != Annotation_variant_type
    ]
  
  if (nrow(type_mismatch) > 0L) {
    warning(
      "Variant-type mismatch found for ",
      nrow(type_mismatch), " lead variants."
    )
  }
  
  final_bolt[, Annotation_variant_type := NULL]
  
  final_bolt <- final_bolt[, .(
    Phenotype,
    Phenotype_file,
    `Lead variant`,
    `Variant ID`,
    Chromosome,
    `Base pair position`,
    `Variant type`,
    `Nearest gene`,
    `Other allele`,
    `Effect allele`,
    `Reference allele`,
    `Alternative allele`,
    Beta,
    `Standard error`,
    `P-value`,
    `Association P field`,
    `Minor allele frequency`,
    `Effect allele frequency`,
    `Sample size`
  )]
  
  colnames(final_bolt)[2] <- "Pheno_code"
  final_bolt <- left_join(final_bolt, pheno_map, by = "Pheno_code")
  final_bolt$Phenotype <- final_bolt$Phenotype_label
  final_bolt <- final_bolt %>% select(-Phenotype_label)
  
  colnames(final_bolt)[2] <- "Phenotype_file"
  bolt_catalogue <- final_bolt %>%
    select(
      Phenotype,
      Phenotype_file,
      `Lead variant`,
      `Variant ID`,
      Chromosome,
      `Base pair position`,
      `Variant type`,
      `Nearest gene`,
      `Other allele`,
      `Effect allele`,
      `Reference allele`,
      `Alternative allele`,
      Beta,
      `Standard error`,
      `P-value`,
      `Association P field`,
      `Minor allele frequency`,
      `Effect allele frequency`,
      `Sample size`
    )

  colnames(bolt_catalogue) <- c(
    "Phenotype", "VarNM", "Lead_variant", "VariantID", "CHR", "BP",
    "Variant_type", "Nearest_gene", "NEA", "EA", "REF", "ALT",
    "Beta", "SE", "P", "P_field", "MAF", "EAF", "N"
  )
  
  catalogue_file <- file.path(result_dir, "Cont_lead_variants.tsv")
  
  fwrite(bolt_catalogue, catalogue_file, sep = "\t", quote = FALSE)
  cat("Saved:", catalogue_file, "\n")
}

