#!/bin/bash
# Kinship.sh
#
# Purpose:
#   KING-based relatedness pipeline for KoGES (LD-pruned SNPs).
#   Produces pairwise kinship tables and unrelated-sample keep/remove lists
#   used by Figure 2B.
#
# Input:
#   KoGES QC-filtered imputed BGEN + sample file
# Output (private; sample IDs) under DATA_manuscript_final/Figure2/raw/:
#   (1) *.prune.in / *.prune.out
#   (2) *_pr0.2_king.kin0          (pairs with kinship >= 0.0884)
#   (3) *_1stdeg_2.king.cutoff.*.id (remove one of each 1st-degree pair)
#   (4) *_2nddeg_2.king.cutoff.*.id (remove one of each 2nd-degree pair)
#
# Downstream:
#   Rscript Figure2B.r  (builds publishable frequency tables + Fig. 2B plots)
#
# How the steps differ
# --------------------
# (1) LD prune:
#     Independent SNP set (--indep-pairwise 50 5 0.2) for downstream KING.
# (2) KING relatedness table:
#     --make-king-table + --king-table-filter 0.0884 writes pairwise kinship.
# (3) KING cutoff (1st-degree, 0.177):
#     Builds an unrelated keep/remove list (not a full kin0 table).
# (4) KING cutoff (2nd-degree, 0.0884):
#     Same as (3) with a stricter relatedness threshold.

set -euo pipefail

KING_DIR="/home2/sjim/openKoGES/DATA_manuscript_final/Figure2/raw"
BGEN="/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.bgen"
SAMP="/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.sample"
PRUNE_IN="${KING_DIR}/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.prune.in"
mkdir -p "${KING_DIR}"

# (1) LD pruning → .prune.in
plink2 --bgen "${BGEN}" ref-first \
  --sample "${SAMP}" \
  --indep-pairwise 50 5 0.2 \
  --out "${KING_DIR}/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8"

# (2) Pairwise KING table (filter 0.0884) → .kin0 (private)
plink2 --bgen "${BGEN}" ref-first \
  --sample "${SAMP}" \
  --extract "${PRUNE_IN}" \
  --make-king-table --king-table-filter 0.0884 \
  --out "${KING_DIR}/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king" \
  --threads 20

# (3) Unrelated sample set by 1st-degree cutoff (0.177)
plink2 --bgen "${BGEN}" ref-first \
  --sample "${SAMP}" \
  --extract "${PRUNE_IN}" \
  --king-cutoff 0.177 \
  --out "${KING_DIR}/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_1stdeg_2" \
  --threads 20

# (4) Unrelated sample set by 2nd-degree cutoff (0.0884)
plink2 --bgen "${BGEN}" ref-first \
  --sample "${SAMP}" \
  --extract "${PRUNE_IN}" \
  --king-cutoff 0.0884 \
  --out "${KING_DIR}/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8_pr0.2_king_2nddeg_2" \
  --threads 20
