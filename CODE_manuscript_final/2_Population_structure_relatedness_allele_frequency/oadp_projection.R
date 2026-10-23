#!/usr/bin/env Rscript
# =============================================================================
# oadp_projection.R
# =============================================================================
#
# Companion to 1000G_projection_OADP.sh (Figure 2A / Supplementary Figure 4).
#
# Purpose:
#   Project KoGES samples onto a fixed 1000 Genomes reference PCA using
#   bigsnpr OADP correction. Supports EAS (Figure 2A) and ALL
#   (Supplementary Figure 4) reference modes.
#
# Reference:
#   Privé et al. (2020). Bioinformatics, 36(16), 4482-4484.
#   https://doi.org/10.1093/bioinformatics/btaa520
#
# Arguments:
#   --ref_bed      Reference BED (PLINK format; LD-pruned)
#   --new_bed      KoGES BED (same variant set as reference)
#   --panel        1000G sample panel (sample, pop, super_pop)
#   --pca_n        Number of PCs (default: 20)
#   --ncores       Number of cores (default: 1)
#   --source_label Reference SOURCE label ("1KG_EAS" or "1KG")
#   --out          OADP coordinate TSV
#   --out_simple   Simple-projection coordinate TSV (for comparison)
#   --out_scaling  Shrinkage diagnostic TSV
#
# Usage:
#   # EAS mode
#   Rscript oadp_projection.R \
#     --ref_bed g1k.reference.final.bed \
#     --new_bed koges.projection.final.bed \
#     --panel 1000Gsamples.ALL.panel \
#     --pca_n 20 --ncores 30 \
#     --source_label 1KG_EAS \
#     --out pca_1kg_eas_reference.oadp.coordinates.tsv \
#     --out_simple pca_1kg_eas_reference.simple.coordinates.tsv \
#     --out_scaling pca_1kg_eas_reference.oadp.scaling.tsv
#
#   # ALL mode
#   Rscript oadp_projection.R \
#     --ref_bed g1k.reference.final.bed \
#     --new_bed koges.projection.final.bed \
#     --panel 1000Gsamples.ALL.panel \
#     --pca_n 20 --ncores 30 \
#     --source_label 1KG \
#     --out pca_1kg_reference.oadp.coordinates.tsv \
#     --out_simple pca_1kg_reference.simple.coordinates.tsv \
#     --out_scaling pca_1kg_reference.oadp.scaling.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(bigsnpr)
  library(data.table)
})

# Check bigsnpr version
pkg_ver <- tryCatch(packageVersion("bigsnpr"), error = function(e) NULL)
if (is.null(pkg_ver)) {
  stop(
    "[ERROR] bigsnpr is not installed.\n",
    "  install.packages('bigsnpr')  or\n",
    "  remotes::install_github('privefl/bigsnpr')"
  )
}
cat("[INFO] bigsnpr version:", as.character(pkg_ver), "\n")
if (pkg_ver < "1.10.0") {
  stop(
    "[ERROR] bigsnpr >= 1.10.0 required (found: ", as.character(pkg_ver), ").\n",
    "  remotes::install_github('privefl/bigsnpr')"
  )
}

# =============================================================================
# Argument parsing
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0L) {
    if (is.null(default)) stop("[ERROR] Missing required argument: ", flag)
    return(default)
  }
  args[idx[1L] + 1L]
}

ref_bed_path <- get_arg("--ref_bed")
new_bed_path <- get_arg("--new_bed")
panel_path   <- get_arg("--panel")
pca_n        <- as.integer(get_arg("--pca_n",        "20"))
ncores       <- as.integer(get_arg("--ncores",       "1"))
source_lbl   <- get_arg("--source_label", "1KG_EAS")
out_oadp     <- get_arg("--out",          "pca_oadp.coordinates.tsv")
out_simple   <- get_arg("--out_simple",   "pca_simple.coordinates.tsv")
out_scaling  <- get_arg("--out_scaling",  "pca_oadp.scaling.tsv")

# Validate source_label
if (!source_lbl %in% c("1KG_EAS", "1KG")) {
  warning("[WARN] Non-standard source_label '", source_lbl,
          "'. Expected '1KG_EAS' (EAS) or '1KG' (ALL).")
}

cat("[INFO] ref_bed     :", ref_bed_path, "\n")
cat("[INFO] new_bed     :", new_bed_path, "\n")
cat("[INFO] panel       :", panel_path,   "\n")
cat("[INFO] pca_n       :", pca_n,        "\n")
cat("[INFO] ncores      :", ncores,       "\n")
cat("[INFO] source_label:", source_lbl,   "\n")
cat("[INFO] mode        :", ifelse(source_lbl == "1KG_EAS", "EAS", "ALL"), "\n")
cat("[INFO] out (OADP)  :", out_oadp,     "\n")
cat("[INFO] out (simple):", out_simple,   "\n")
cat("[INFO] out_scaling :", out_scaling,  "\n")

for (f in c(ref_bed_path, new_bed_path, panel_path)) {
  if (!file.exists(f)) stop("[ERROR] File not found: ", f)
}

# =============================================================================
# STEP 1. Open BED files
# =============================================================================
cat("\n[STEP 1] Open BED files (bigsnpr::bed)\n")

obj_ref <- bed(ref_bed_path)
obj_new <- bed(new_bed_path)

n_ref <- nrow(obj_ref)
n_new <- nrow(obj_new)
n_var <- ncol(obj_ref)

cat("[INFO] Reference (", source_lbl, "):", n_ref, "samples,", n_var, "variants\n")
cat("[INFO] New       (KoGES)         :", n_new, "samples,", ncol(obj_new), "variants\n")
cat("[INFO] p/n ratio                 :", round(n_var / n_ref, 1),
    "(p >> n; OADP correction expected)\n")

if (ncol(obj_new) != n_var) {
  stop("[ERROR] Variant count mismatch: ref=", n_var, " new=", ncol(obj_new))
}

# Verify variant order
ref_map <- obj_ref$map
new_map <- obj_new$map

if (identical(ref_map$marker.ID, new_map$marker.ID)) {
  cat("[INFO] Variant IDs match\n")
} else {
  cat("[WARN] Marker IDs differ; checking positions\n")
  ref_pos <- paste(ref_map$chromosome, ref_map$physical.pos, sep = ":")
  new_pos <- paste(new_map$chromosome, new_map$physical.pos, sep = ":")
  if (!identical(ref_pos, new_pos)) {
    stop("[ERROR] Variant positions differ. Check BED variant order.")
  }
  cat("[INFO] Positions match\n")
}

# =============================================================================
# STEP 2-3. Reference PCA + KoGES OADP projection (bed_projectPCA)
# =============================================================================
# In bigsnpr 1.12.x, bed_projectPCA does not accept a precomputed SVD.
# obj.svd.ref is a return-list name; extra args (...) are passed to bed_autoSVD.
# Returns: $obj.svd.ref, $simple_proj, $OADP_proj  (not $u / $OADP)
#
# thr.r2 = Inf  : skip LD pruning (already done in PLINK2)
# max.iter = 0  : no iterative outlier removal (already QC'd)
# =============================================================================
cat("\n[STEP 2-3] Reference PCA + KoGES OADP projection (bed_projectPCA)\n")

proj_result <- bed_projectPCA(
  obj_ref,
  obj_new,
  k         = pca_n,
  thr.r2    = Inf,
  roll.size = 50L,
  min.mac   = 1L,
  max.iter  = 0L,
  ncores    = ncores
)

obj_svd <- proj_result$obj.svd.ref
k_actual <- length(obj_svd$d)
cat("[INFO] Number of PCs:", k_actual, "\n")
cat("[INFO] Singular values (top 5):",
    paste(round(obj_svd$d[seq_len(min(5L, k_actual))], 3), collapse = ", "), "\n")

# predict() = U * D (standard PCA coordinates)
ref_pcs <- predict(obj_svd)    # n_ref × k_actual
cat("[INFO] Reference PC matrix:", paste(dim(ref_pcs), collapse = " x "), "\n")

as_pc_matrix <- function(x) {
  if (inherits(x, "FBM")) x[] else as.matrix(x)
}
new_simple <- as_pc_matrix(proj_result$simple_proj)  # simple projection (with shrinkage)
new_oadp   <- as_pc_matrix(proj_result$OADP_proj)    # OADP-corrected projection

cat("[INFO] Simple dim:", paste(dim(new_simple), collapse = " x "), "\n")
cat("[INFO] OADP   dim:", paste(dim(new_oadp),   collapse = " x "), "\n")

# =============================================================================
# STEP 4. Shrinkage diagnostics
# =============================================================================
cat("\n[STEP 4] Shrinkage diagnostics\n")

sd_simple      <- apply(new_simple, 2L, sd)
sd_oadp        <- apply(new_oadp,   2L, sd)
correction_obs <- sd_oadp / sd_simple

scaling_dt <- data.table(
  PC                   = seq_len(k_actual),
  n_ref                = n_ref,
  n_var                = n_var,
  c_ratio_p_over_n     = round(n_var / (n_ref - 1L), 3),
  eigenvalue           = round(obj_svd$d^2 / (n_ref - 1L), 6),
  singular_value       = round(obj_svd$d, 4),
  sd_ref_direct        = round(apply(ref_pcs,    2L, sd), 6),
  sd_simple_projection = round(sd_simple, 6),
  sd_oadp_projection   = round(sd_oadp,   6),
  observed_correction  = round(correction_obs, 4)
)

cat("\n--- Shrinkage summary (first 10 PCs) ---\n")
cat(sprintf("%-4s | %-10s | %-10s | %-10s | %-10s\n",
  "PC", "eigenvalue", "SD_simple", "SD_OADP", "correction"))
cat(strrep("-", 55), "\n")
for (i in seq_len(min(10L, k_actual))) {
  cat(sprintf("%-4d | %-10.5f | %-10.5f | %-10.5f | x%-9.4f\n",
    i, scaling_dt$eigenvalue[i],
    sd_simple[i], sd_oadp[i], correction_obs[i]))
}
if (k_actual > 10L) cat("  ... (", k_actual - 10L, " PCs omitted)\n")

cat("[INFO]", sum(correction_obs > 1.01), "/", k_actual,
    "PCs with >=1% OADP correction\n")

fwrite(scaling_dt, out_scaling, sep = "\t")
cat("[INFO] Wrote shrinkage diagnostics:", out_scaling, "\n")

# =============================================================================
# STEP 5. Load 1000G panel
# =============================================================================
cat("\n[STEP 5] Load 1000G panel\n")

# Official panel header ends with empty trailing tabs, so column count exceeds
# the data. Without fill=TRUE, fread drops the header and uses the first sample
# as column names.
panel <- fread(panel_path, header = TRUE, sep = "\t", fill = TRUE)
setnames(panel, tolower(names(panel)))
empty_cols <- names(panel) == "" | grepl("^v[0-9]+$", names(panel))
if (any(empty_cols)) panel[, (names(panel)[empty_cols]) := NULL]

# Standardize column names
if (!"sample" %in% names(panel)) {
  alt <- c("sample_id", "id", "iid", "#sample_id", "#sample")
  found <- alt[alt %in% names(panel)]
  if (length(found)) setnames(panel, found[1L], "sample")
  else stop("[ERROR] No 'sample' column in panel. Columns: ",
            paste(names(panel), collapse = ", "))
}
if (!"super_pop" %in% names(panel)) {
  sp <- grep("super", names(panel), value = TRUE, ignore.case = TRUE)[1L]
  if (!is.na(sp)) setnames(panel, sp, "super_pop")
  else stop("[ERROR] No 'super_pop' column in panel.")
}
stopifnot("pop" %in% names(panel))

panel_key <- unique(panel[, .(sample, pop, super_pop)])
setkeyv(panel_key, "sample")
cat("[INFO] Panel loaded:", nrow(panel_key), "samples\n")

# =============================================================================
# STEP 6. Assemble and write coordinate tables
# =============================================================================
cat("\n[STEP 6] Assemble coordinate tables\n")

ref_fam <- as.data.table(obj_ref$fam)
new_fam <- as.data.table(obj_new$fam)

# Panel look-up (key join)
ref_info  <- panel_key[.(ref_fam$sample.ID)]
n_miss    <- sum(is.na(ref_info$pop))
if (n_miss > 0L) {
  warning("[WARN] ", n_miss, " reference samples missing from panel; set to 'NA'.")
  ref_info[is.na(pop),   pop       := "NA"]
  ref_info[is.na(super_pop), super_pop := "NA"]
}

# Reference block
ref_pc_dt <- as.data.table(ref_pcs[, seq_len(k_actual), drop = FALSE])
setnames(ref_pc_dt, paste0("PC", seq_len(k_actual)))

ref_block <- cbind(
  data.table(
    FID      = ref_fam$family.ID,
    IID      = ref_fam$sample.ID,
    SOURCE   = source_lbl,
    POP      = ref_info$pop,
    SUPERPOP = ref_info$super_pop
  ),
  ref_pc_dt
)

# KoGES block — OADP
oadp_pc_dt <- as.data.table(new_oadp[, seq_len(k_actual), drop = FALSE])
setnames(oadp_pc_dt, paste0("PC", seq_len(k_actual)))

new_oadp_block <- cbind(
  data.table(
    FID      = new_fam$family.ID,
    IID      = new_fam$sample.ID,
    SOURCE   = "KoGES",
    POP      = "KoGES",
    SUPERPOP = "EAS"
  ),
  oadp_pc_dt
)

# KoGES block — simple projection (comparison)
simple_pc_dt <- as.data.table(new_simple[, seq_len(k_actual), drop = FALSE])
setnames(simple_pc_dt, paste0("PC", seq_len(k_actual)))

new_simple_block <- cbind(
  data.table(
    FID      = new_fam$family.ID,
    IID      = new_fam$sample.ID,
    SOURCE   = "KoGES",
    POP      = "KoGES",
    SUPERPOP = "EAS"
  ),
  simple_pc_dt
)

# Write outputs
fwrite(rbindlist(list(ref_block, new_oadp_block)),   out_oadp,   sep = "\t")
fwrite(rbindlist(list(ref_block, new_simple_block)), out_simple, sep = "\t")

cat("[DONE] OADP coordinates :", out_oadp,   "\n")
cat("[DONE] Simple coordinates:", out_simple, "\n")
cat("[INFO] Reference:", nrow(ref_block),
    " | KoGES:", nrow(new_oadp_block), "\n")

# =============================================================================
# STEP 7. Simple vs OADP distribution summary
# =============================================================================
cat("\n[STEP 7] KoGES distribution (Simple -> OADP, first 3 PCs)\n")
cat(sprintf("%-6s | %-22s | %-22s | %s\n",
    "PC", "Simple (mean ± SD)", "OADP (mean ± SD)", "correction"))
cat(strrep("-", 75), "\n")
for (i in seq_len(min(3L, k_actual))) {
  cat(sprintf("PC%-4d | %8.4f ± %-12.4f | %8.4f ± %-12.4f | x%.3f\n",
    i,
    mean(new_simple[, i]), sd(new_simple[, i]),
    mean(new_oadp[, i]),   sd(new_oadp[, i]),
    correction_obs[i]))
}

cat("\n[ALL DONE] oadp_projection.R finished\n")
