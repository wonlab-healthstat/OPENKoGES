#!/bin/bash
# Allele_frequency.sh
#
# Purpose:
#   Compute autosomal allele frequencies for KoGES imputed genotypes
#   (PLINK2 --freq). Output is used by Figure2C.r (MAF / MAC spectrum).
#
# Input:
#   KoGES QC-filtered imputed BGEN + sample file
# Output:
#   DATA_manuscript_final/Figure2/KoGES_autosome_MAF.afreq
#
# Downstream:
#   Rscript Figure2C.r

set -euo pipefail

plink2 \
  --bgen /data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.bgen ref-first \
  --sample /data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.sample \
  --autosome \
  --freq \
  --out /home2/sjim/openKoGES/DATA_manuscript_final/Figure2/KoGES_autosome_MAF
