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

# Resolve openKoGES repo root without hard-coding host absolute paths.
_find_repo_root() {
  local d="$1"
  local i
  for i in $(seq 1 12); do
    if [[ -d "$d/CODE_manuscript_final" && -d "$d/DATA_manuscript_final" ]]; then
      printf '%s' "$d"
      return 0
    fi
    local parent
    parent="$(dirname "$d")"
    [[ "$parent" == "$d" ]] && break
    d="$parent"
  done
  echo "Could not locate openKoGES root; set OPENKOGES_ROOT." >&2
  return 1
}
if [[ -z "${SCRIPT_DIR:-}" ]]; then
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
fi
REPO_ROOT="${OPENKOGES_ROOT:-$(_find_repo_root "$SCRIPT_DIR")}"

: "${BGEN_DIR:?Set BGEN_DIR to the KoGES BGEN directory}"


plink2 \
  --bgen "${BGEN_DIR}/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.bgen" ref-first \
  --sample "${BGEN_DIR}/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.sample" \
  --autosome \
  --freq \
  --out "${REPO_ROOT}/DATA_manuscript_final/Figure2/KoGES_autosome_MAF"
