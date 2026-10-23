#!/usr/bin/env bash
# =============================================================================
# 1000G_projection_OADP.sh
# =============================================================================
#
# Purpose:
#   Project KoGES onto 1000G reference PCs with bigsnpr OADP correction.
#   Manuscript analysis backend for Figure 2A (EAS) and Supplementary Figure 4
#   (ALL).
#
# Usage:
#   bash 1000G_projection_OADP.sh [eas|all|both]
#     eas  : 1KG EAS reference (Figure 2A)
#     all  : 1KG ALL reference (Supplementary Figure 4)
#     both : eas then all (default)
#
# Upstream genotype prep (non-OADP PLINK --score archive):
#   bash CODE_openkoges_additional/forPCA_projection/1000G_projection_PCA.sh eas|all
#   Produces the *.reference.final / koges.projection.final PGEN filesets
#   reused below.
#
# Outputs:
#   [eas]  pca_1kg_eas_reference.oadp.coordinates.tsv   <- Figure2A.R default
#          pca_1kg_eas_reference.simple.coordinates.tsv
#          pca_1kg_eas_reference.oadp.scaling.tsv
#   [all]  pca_1kg_reference.oadp.coordinates.tsv       <- SuppleFigure4.R default
#          pca_1kg_reference.simple.coordinates.tsv
#          pca_1kg_reference.oadp.scaling.tsv
#
# Downstream plots:
#   Rscript Figure2A.R
#   Rscript SuppleFigure4.R
# =============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REF_ARG="${1:-both}"

case "${REF_ARG}" in
  eas|all|both) ;;
  -h|--help)
    echo "Usage: $0 [eas|all|both]"
    exit 0
    ;;
  *)
    echo "[ERROR] Unknown mode '${REF_ARG}'. Use: eas | all | both"
    exit 1
    ;;
esac

# both: run eas then all
if [[ "${REF_ARG}" == "both" ]]; then
  bash "${BASH_SOURCE[0]}" eas
  bash "${BASH_SOURCE[0]}" all
  exit 0
fi

# =============================================================================
# Shared settings
# =============================================================================
THREADS=30
PLINK2=/bin2/plink2
BASE_WD=/data4/sjim/openKoGES_revision/DATA_projectionPCA
PCA_N=20

# Companion scripts
OADP_RSCRIPT="${SCRIPT_DIR}/oadp_projection.R"
COMPARE_RSCRIPT="/home2/sjim/openKoGES/CODE_openkoges_additional/forPCA_projection/Figure_compare_oadp.R"

# =============================================================================
# Mode-specific settings
# =============================================================================
if [[ "${REF_ARG}" == "eas" ]]; then

  WD="${BASE_WD}/1kg_eas_reference_projection_all"
  REF_BED_PREFIX="g1k.eas.reference.final"   # EAS pipeline fileset prefix
  KOGES_BED_PREFIX="koges.projection.final"
  SOURCE_LABEL="1KG_EAS"
  OUT_PREFIX="pca_1kg_eas_reference"
  G1K_PANEL=/data/1000G/Phase3_v5a/integrated_call_samples_v3.20130502.ALL.panel
  EXPECTED_REF_N=504
  FIGURE_RSCRIPT="${SCRIPT_DIR}/Figure2A.R"
  LOG_NAME="oadp_projection_eas.log"
  MODE_DESC="EAS reference (Figure 2A)"

else  # all

  WD="${BASE_WD}/1kg_reference_projection_all"
  REF_BED_PREFIX="g1k.reference.final"       # ALL pipeline fileset prefix
  KOGES_BED_PREFIX="koges.projection.final"
  SOURCE_LABEL="1KG"
  OUT_PREFIX="pca_1kg_reference"
  G1K_PANEL=/data/1000G/Phase3_v5a/integrated_call_samples_v3.20130502.ALL.panel
  EXPECTED_REF_N=2504
  FIGURE_RSCRIPT="${SCRIPT_DIR}/SuppleFigure4.R"
  LOG_NAME="oadp_projection_all.log"
  MODE_DESC="ALL reference (Supplementary Figure 4)"

fi

LOG="${WD}/${LOG_NAME}"

# =============================================================================
# Helpers
# =============================================================================
run() {
  echo -e "\n[CMD] $*" | tee -a "${LOG}"
  "$@" 2>&1 | tee -a "${LOG}"
}

need_file() {
  [[ -f "$1" ]] || {
    echo "[ERROR] File not found: $1" | tee -a "${LOG}"
    exit 1
  }
}

skip_if_file() {
  if [[ -f "$1" ]]; then
    echo "[SKIP] Already exists: $1" | tee -a "${LOG}"
    return 0
  fi
  return 1
}

# =============================================================================
# Start
# =============================================================================
cd "${WD}"
: > "${LOG}"

echo "============================================================" | tee -a "${LOG}"
echo "[START] OADP Projection Patch: ${MODE_DESC}" | tee -a "${LOG}"
echo "[INFO] WD      : ${WD}"           | tee -a "${LOG}"
echo "[INFO] MODE    : ${REF_ARG}"      | tee -a "${LOG}"
echo "[INFO] SOURCE  : ${SOURCE_LABEL}" | tee -a "${LOG}"
echo "[INFO] THREADS : ${THREADS}"      | tee -a "${LOG}"
echo "[INFO] PCA_N   : ${PCA_N}"        | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"

# =============================================================================
# Pre-flight checks
# =============================================================================
echo "[CHECK] Required inputs" | tee -a "${LOG}"
echo "[INFO] Reference prefix : ${REF_BED_PREFIX}" | tee -a "${LOG}"
echo "[INFO] KoGES prefix     : ${KOGES_BED_PREFIX}" | tee -a "${LOG}"

need_file "${REF_BED_PREFIX}.pgen"
need_file "${REF_BED_PREFIX}.pvar"
need_file "${REF_BED_PREFIX}.psam"
need_file "${KOGES_BED_PREFIX}.pgen"
need_file "${KOGES_BED_PREFIX}.pvar"
need_file "${KOGES_BED_PREFIX}.psam"
need_file "${G1K_PANEL}"
need_file "${OADP_RSCRIPT}"

# Check reference sample count (EAS: 504, ALL: 2504)
REF_ACTUAL=$(awk '$0 !~ /^#/ {n++} END{print n+0}' "${REF_BED_PREFIX}.psam")
if [[ "${REF_ACTUAL}" -ne "${EXPECTED_REF_N}" ]]; then
  echo "[ERROR] Reference sample count mismatch: found=${REF_ACTUAL}, expected=${EXPECTED_REF_N}" \
    | tee -a "${LOG}"
  exit 1
fi
echo "[OK] Reference samples: ${REF_ACTUAL} (${SOURCE_LABEL})" | tee -a "${LOG}"

# =============================================================================
# STEP A. Convert PGEN to BED (bigsnpr input)
# =============================================================================
echo -e "\n[STEP A] Convert PGEN to BED" | tee -a "${LOG}"

if skip_if_file "${REF_BED_PREFIX}.bed"; then
  : # skip
else
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile "${REF_BED_PREFIX}" \
    --make-bed \
    --out "${REF_BED_PREFIX}"
fi

if skip_if_file "${KOGES_BED_PREFIX}.bed"; then
  : # skip
else
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile "${KOGES_BED_PREFIX}" \
    --make-bed \
    --out "${KOGES_BED_PREFIX}"
fi

need_file "${REF_BED_PREFIX}.bed"
need_file "${KOGES_BED_PREFIX}.bed"

# Log variant/sample counts
REF_VAR=$(awk 'END{print NR}' "${REF_BED_PREFIX}.bim")
NEW_VAR=$(awk 'END{print NR}' "${KOGES_BED_PREFIX}.bim")
NEW_SAMP=$(awk 'END{print NR}' "${KOGES_BED_PREFIX}.fam")

echo "[INFO] Reference : ${REF_ACTUAL} samples, ${REF_VAR} variants" | tee -a "${LOG}"
echo "[INFO] KoGES     : ${NEW_SAMP} samples, ${NEW_VAR} variants"   | tee -a "${LOG}"
echo "[INFO] p/n ratio : $(echo "scale=1; ${REF_VAR}/${REF_ACTUAL}" | bc)" | tee -a "${LOG}"

if [[ "${REF_VAR}" -ne "${NEW_VAR}" ]]; then
  echo "[ERROR] Variant count mismatch: ref=${REF_VAR}, new=${NEW_VAR}" | tee -a "${LOG}"
  exit 1
fi

# =============================================================================
# STEP B. bigsnpr OADP projection
# =============================================================================
echo -e "\n[STEP B] Run bigsnpr OADP projection (mode=${REF_ARG})" | tee -a "${LOG}"

OADP_COORD="${OUT_PREFIX}.oadp.coordinates.tsv"
SIMPLE_COORD="${OUT_PREFIX}.simple.coordinates.tsv"
SCALING_FILE="${OUT_PREFIX}.oadp.scaling.tsv"

if skip_if_file "${OADP_COORD}"; then
  : # skip
else
  run Rscript "${OADP_RSCRIPT}" \
    --ref_bed      "${REF_BED_PREFIX}.bed" \
    --new_bed      "${KOGES_BED_PREFIX}.bed" \
    --panel        "${G1K_PANEL}" \
    --pca_n        "${PCA_N}" \
    --ncores       "${THREADS}" \
    --source_label "${SOURCE_LABEL}" \
    --out          "${OADP_COORD}" \
    --out_simple   "${SIMPLE_COORD}" \
    --out_scaling  "${SCALING_FILE}"
fi

need_file "${OADP_COORD}"
need_file "${SIMPLE_COORD}"
need_file "${SCALING_FILE}"

echo "[INFO] OADP coordinate rows: $(wc -l < "${OADP_COORD}")" | tee -a "${LOG}"

# =============================================================================
# STEP C. Generate figures from OADP coordinates
# =============================================================================
# Figure2A.R expects SOURCE = "1KG_EAS" | "KoGES" (EAS mode).
# SuppleFigure4.R expects SOURCE = "1KG" | "KoGES" (ALL mode).
# OADP outputs use matching SOURCE labels and can be passed directly.
# =============================================================================
echo -e "\n[STEP C] Generate figures (OADP)" | tee -a "${LOG}"

if [[ -f "${FIGURE_RSCRIPT}" ]]; then
  run Rscript "${FIGURE_RSCRIPT}" all \
    "${WD}/${OADP_COORD}"
  echo "[INFO] Figure generation complete (OADP)" | tee -a "${LOG}"
else
  echo "[WARN] Figure script missing; skip: ${FIGURE_RSCRIPT}" | tee -a "${LOG}"
fi

# =============================================================================
# STEP D. Simple vs OADP comparison plots
# =============================================================================
echo -e "\n[STEP D] Simple vs OADP comparison plots" | tee -a "${LOG}"

if [[ -f "${COMPARE_RSCRIPT}" ]]; then
  run Rscript "${COMPARE_RSCRIPT}" \
    --mode         "${REF_ARG}" \
    --source_label "${SOURCE_LABEL}" \
    --oadp_coord   "${OADP_COORD}" \
    --simple_coord "${SIMPLE_COORD}" \
    --scaling      "${SCALING_FILE}" \
    --out_dir      "${WD}"
  echo "[INFO] Comparison plots complete" | tee -a "${LOG}"
else
  echo "[WARN] Comparison script missing; skip: ${COMPARE_RSCRIPT}" | tee -a "${LOG}"
fi

# =============================================================================
# FINAL SUMMARY
# =============================================================================
echo -e "\n============================================================" | tee -a "${LOG}"
echo "[FINAL SUMMARY] ${MODE_DESC}" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/${OADP_COORD}"   | tee -a "${LOG}"
echo "  -> Direct input for ${FIGURE_RSCRIPT} (SOURCE: ${SOURCE_LABEL} / KoGES)" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/${SIMPLE_COORD}" | tee -a "${LOG}"
echo "  -> Simple projection for comparison with --score" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/${SCALING_FILE}" | tee -a "${LOG}"
echo "  -> Per-PC shrinkage correction diagnostics" | tee -a "${LOG}"
echo "[LOG]    ${LOG}" | tee -a "${LOG}"
echo "[ALL DONE] mode=${REF_ARG}" | tee -a "${LOG}"
