#!/usr/bin/env bash
# 1000G_projection_PCA.sh
#
# Purpose:
#   Build a 1000 Genomes Phase 3 reference PCA and project KoGES samples
#   onto those fixed PC axes (genome-wide imputed markers; no K-chip filter).
#   Analysis backend for Figure 2A (EAS) and Supplementary Figure 4 (ALL).
#
# Usage:
#   bash 1000G_projection_PCA.sh [eas|all|both]
#     eas  : 1KG EAS-only reference PCA -> KoGES projection (Figure 2A)
#            Output dir: .../DATA_projectionPCA/1kg_eas_reference_projection_all
#     all  : 1KG ALL reference PCA -> KoGES projection (Supplementary Figure 4)
#            Output dir: .../DATA_projectionPCA/1kg_reference_projection_all
#     both : run eas then all (default)
#
# Downstream plots:
#   Rscript Figure2A.R
#   Rscript SuppleFigure4.R

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

if [[ "${REF_ARG}" == "both" ]]; then
  bash "${BASH_SOURCE[0]}" eas
  bash "${BASH_SOURCE[0]}" all
  exit 0
fi

_run_embedded() {
  local label="$1"
  local tmp rc
  tmp="$(mktemp "${TMPDIR:-/tmp}/1000G_projection_PCA.${label}.XXXXXX.sh")"
  cat > "${tmp}"
  set +e
  bash "${tmp}"
  rc=$?
  set -e
  rm -f "${tmp}"
  return "${rc}"
}

if [[ "${REF_ARG}" == "eas" ]]; then
  _run_embedded eas <<'__EAS_BODY__'
set -euo pipefail

# ============================================================================
# Figure 2A: KoGES projected onto 1000G EAS-only PCs (genome-wide imputed)
# Usage: bash Figure2A.sh
# ============================================================================
#
# PURPOSE
# -------
# This script creates a 1000 Genomes Phase 3 East Asian (1KG EAS) reference PCA
# and projects all KoGES samples onto those fixed EAS-specific PC axes
# (1kg_eas_reference_projection_all).
#
# Contrast with joint_pca_kchip_variant.sh
# ---------------------------------------
# joint_pca_kchip_variant.sh:
#   KoGES + 1KG ALL joint PCA (KoGES influences the axes).
#
# This script:
#   PCA axes from 1KG EAS only; KoGES is projected onto those axes.
#
# Key design choices
# ------------------
# - The available KoGES input is imputed BGEN.  K-chip manifests are used only
#   to select positions targeted by BOTH KoreanChip v1.0 and v1.1.  This is not
#   a reconstruction of directly genotyped calls.
# - PCA axes, reference allele frequencies, MAF filtering, and LD pruning are
#   determined by 1KG EAS (504 individuals). KoGES therefore cannot dominate the PCA axes.
# - KoGES is used only to retain technically available shared markers
#   (--geno 0.02 dosage) and for projection.  KoGES MAF is NOT used to select
#   reference PCA markers.
# - All KoGES samples, including related samples, are projected.  They do not
#   influence the 1KG PCA axes.  For a formal KoGES centroid/distance analysis,
#   use an unrelated subset separately.
# - The final coordinates file contains direct 1KG PC coordinates and
#   scale-calibrated KoGES projected coordinates in the same PC scale.
#
# Shared intermediates
# --------------------
# Steps 1-7 match joint_pca_kchip_variant.sh (same inputs and allele
# harmonization through koges.common / g1k.common).  When those files already
# exist under SHARED_WD, this script reuses them via symlink instead of
# recomputing.  Projection-specific QC, PCA, and outputs live under WD.
#
# Genome build
# ------------
# K-chip v1.0/v1.1, KoGES BGEN, and 1KG Phase 3 VCFs must all be GRCh37/hg19.
# Verify this once from the supplied manifests before running the full pipeline.
#
# Outputs
# -------
# pca_1kg_eas_reference.coordinates.tsv
#   Ready-to-plot table: direct 1KG PCs + scaled projected KoGES PCs.
# pca_1kg_eas_reference.scaling.tsv
#   Self-projection calibration diagnostics for each PC.
#
# ============================================================================


# ============================================================================
# SETTINGS
# ============================================================================

THREADS=30
PLINK2=/bin2/plink2

# Reference-marker QC
MAF=0.05
GENO=0.02

# PCA / projection settings
PCA_N=20
MIN_SELFPROJ_COR=0.999

# Allele-frequency flip check: KoGES vs 1KG EAS
# Few residual flips are excluded; many flips abort as systematic.
AF_DIFF_WARN=0.20
AF_FLIP_MATCH=0.05
MAX_SUSPECT_FLIPS_ABORT=200
EXPECTED_1KG_EAS_N=504

# Engineering guards, not biological cutoffs
MIN_COMMON_SNPS=10000
MIN_PRUNED_SNPS=5000

# Genome-wide imputed markers only (1kg_eas_reference_projection_all)
MODE_ARG=all
USE_KCHIP=0

BASE_WD=/data4/sjim/openKoGES_revision/DATA_projectionPCA
SHARED_RAW_WD="${BASE_WD}"
SHARED_WD="${BASE_WD}/all_variants"
WD="${BASE_WD}/1kg_eas_reference_projection_all"
LOG=koges_1kg_eas_reference_projection_all.log
MARKER_DESC="genome-wide imputed (no K-chip filter)"
PLOT_OUTPREF=Figure2A_eas_projection_all
PLOT_RSCRIPT=/home2/sjim/openKoGES/CODE_manuscript_final/2_Population_structure_relatedness_allele_frequency/Figure2A.R

mkdir -p "${SHARED_WD}"
mkdir -p "${WD}"
cd "${WD}"

: > "${LOG}"


# ============================================================================
# INPUTS
# ============================================================================

KOGES_BGEN=/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.bgen
KOGES_SAMP=/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.sample

G1K_DIR=/data/1000G/Phase3_v5a
G1K_VCF_PATTERN="${G1K_DIR}/chr%s.1kg.phase3.v5a.vcf.gz"
G1K_PANEL=/home2/sjim/openKoGES/DATA_baseline/forPCA/1000Gsamples.ALL.panel


# ============================================================================
# HELPERS
# ============================================================================

# Skip only the product type implied by the PLINK flags.
_run_should_skip() {
  local out="$1"
  shift
  local make_bed=0 make_pgen=0 do_pca=0 do_score=0 do_prune=0 write_snplist=0 do_freq=0
  local a
  for a in "$@"; do
    case "${a}" in
      --make-bed) make_bed=1 ;;
      --make-pgen) make_pgen=1 ;;
      --pca) do_pca=1 ;;
      --score) do_score=1 ;;
      --indep-pairwise) do_prune=1 ;;
      --write-snplist) write_snplist=1 ;;
      --freq) do_freq=1 ;;
    esac
  done

  if (( make_bed == 1 )); then
    [[ -f "${out}.bed" && -f "${out}.bim" && -f "${out}.fam" ]] || return 1
    echo "[INFO] Skip (exists bed): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( make_pgen == 1 )); then
    [[ -f "${out}.pgen" && -f "${out}.pvar" && -f "${out}.psam" ]] || return 1
    echo "[INFO] Skip (exists pfile): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_pca == 1 )); then
    [[ -f "${out}.eigenvec" ]] || return 1
    echo "[INFO] Skip (exists eigenvec): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_score == 1 )); then
    [[ -f "${out}.sscore" ]] || return 1
    echo "[INFO] Skip (exists sscore): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_prune == 1 )); then
    [[ -f "${out}.prune.in" ]] || return 1
    echo "[INFO] Skip (exists prune): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( write_snplist == 1 )); then
    [[ -f "${out}.snplist" ]] || return 1
    echo "[INFO] Skip (exists snplist): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_freq == 1 )); then
    [[ -f "${out}.acount" || -f "${out}.afreq" ]] || return 1
    echo "[INFO] Skip (exists freq): ${out}" | tee -a "${LOG}"; return 0
  fi
  return 1
}

skip_if_file() {
  local f="$1"
  if [[ -f "${f}" ]]; then
    echo "[INFO] Skip (exists): ${f}" | tee -a "${LOG}"
    return 0
  fi
  return 1
}

run() {
  local out="" prev=""
  for a in "$@"; do
    if [[ "${prev}" == "--out" ]]; then out="${a}"; fi
    prev="${a}"
  done
  if [[ -n "${out}" ]] && _run_should_skip "${out}" "$@"; then
    return 0
  fi
  echo -n "[CMD]" | tee -a "${LOG}"
  for a in "$@"; do
    echo -n " ${a}" | tee -a "${LOG}"
  done
  echo "" | tee -a "${LOG}"
  "$@" 2>&1 | tee -a "${LOG}"
}

need_file() {
  local f="$1"
  [[ -f "${f}" ]] || {
    echo "[ERROR] Missing required file: ${f}" | tee -a "${LOG}"
    exit 1
  }
}

need_exec() {
  local x="$1"
  [[ -x "${x}" ]] || {
    echo "[ERROR] Not executable: ${x}" | tee -a "${LOG}"
    exit 1
  }
}

count_pvar() {
  awk '$0 !~ /^#/ {n++} END {print n+0}' "$1"
}

count_psam() {
  awk '$0 !~ /^#/ {n++} END {print n+0}' "$1"
}

assert_min_lines() {
  local f="$1"
  local min="$2"
  local label="$3"
  local n
  n=$(wc -l < "${f}" | tr -d ' ')

  if (( n < min )); then
    echo "[ERROR] ${label}: ${n} < minimum ${min}" | tee -a "${LOG}"
    exit 1
  fi
  echo "[INFO] ${label}: ${n} (min ${min} OK)" | tee -a "${LOG}"
}

# Reuse an identical intermediate from SHARED_WD via symlink when present.
# Usage: reuse_shared <local_name> [shared_name]
# Returns 0 if reused, 1 if missing in SHARED_WD (caller should compute).
reuse_shared() {
  local local_name="$1"
  local shared_name="${2:-$1}"
  local src="${SHARED_WD}/${shared_name}"
  local dst="${WD}/${local_name}"

  if [[ -e "${dst}" || -L "${dst}" ]]; then
    echo "[INFO] Reusing existing local file: ${local_name}" | tee -a "${LOG}"
    return 0
  fi

  if [[ -f "${src}" ]]; then
    ln -s "${src}" "${dst}"
    echo "[INFO] Linked shared intermediate: ${shared_name} -> ${local_name}" \
      | tee -a "${LOG}"
    return 0
  fi

  return 1
}

# Drop a local path so a later write cannot clobber SHARED_WD through a symlink.
clear_local_target() {
  local path="$1"
  if [[ -L "${path}" ]]; then
    rm -f "${path}"
  elif [[ -e "${path}" ]]; then
    rm -f "${path}"
  fi
}

# Reuse a PLINK2 prefix trio (.pgen/.pvar/.psam) from SHARED_WD.
# Links only when all three shared files exist (avoids partial reuse).
reuse_shared_pfile() {
  local local_prefix="$1"
  local shared_prefix="${2:-$1}"
  local ext
  local all_local=1
  local all_shared=1

  for ext in pgen pvar psam; do
    if [[ ! -e "${WD}/${local_prefix}.${ext}" && ! -L "${WD}/${local_prefix}.${ext}" ]]; then
      all_local=0
    fi
    if [[ ! -f "${SHARED_WD}/${shared_prefix}.${ext}" ]]; then
      all_shared=0
    fi
  done

  if (( all_local == 1 )); then
    echo "[INFO] Reusing existing local pfile: ${local_prefix}" | tee -a "${LOG}"
    return 0
  fi

  if (( all_shared == 1 )); then
    for ext in pgen pvar psam; do
      clear_local_target "${WD}/${local_prefix}.${ext}"
      ln -s "${SHARED_WD}/${shared_prefix}.${ext}" "${WD}/${local_prefix}.${ext}"
    done
    echo "[INFO] Linked shared pfile: ${shared_prefix} -> ${local_prefix}" \
      | tee -a "${LOG}"
    return 0
  fi

  return 1
}


# ============================================================================
# PRE-FLIGHT
# ============================================================================

echo "============================================================" | tee -a "${LOG}"
echo "[START] KoGES projection onto 1KG EAS-only (${MARKER_DESC}, mode=${MODE_ARG})" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
echo "[INFO] SHARED_WD=${SHARED_WD}" | tee -a "${LOG}"
echo "[INFO] WD=${WD}" | tee -a "${LOG}"
echo "[INFO] THREADS=${THREADS} MAF=${MAF} GENO=${GENO} PCA_N=${PCA_N}" | tee -a "${LOG}"

need_exec "${PLINK2}"
need_file "${KOGES_BGEN}"
need_file "${KOGES_SAMP}"
need_file "${G1K_PANEL}"
if (( USE_KCHIP == 1 )); then
  need_file "${KCHIP10}"
  need_file "${KCHIP11}"
fi

if [[ -f "pca_1kg_eas_reference.coordinates.tsv" ]]; then
  echo "[INFO] Final coordinates already present for mode=${MODE_ARG}; skip PCA steps." | tee -a "${LOG}"
  echo "[OUTPUT] ${WD}/pca_1kg_eas_reference.coordinates.tsv" | tee -a "${LOG}"
  echo "[PLOT] Rscript ${PLOT_RSCRIPT} all" | tee -a "${LOG}"
  Rscript "${PLOT_RSCRIPT}" all \
    "${WD}/pca_1kg_eas_reference.coordinates.tsv"
  echo "[ALL DONE]" | tee -a "${LOG}"
  exit 0
fi


# ============================================================================
# STEP 1. K-chip v1.0 ∩ v1.1 autosomal positions
# ============================================================================

if (( USE_KCHIP == 1 )); then
echo "[STEP 1] Extract K-chip v1.0 ∩ v1.1 positions" | tee -a "${LOG}"

if ! reuse_shared kchip.hg19.intersection.positions kchip.hg19.positions; then
  clear_local_target kchip.hg19.intersection.positions
  python3 - "${KCHIP10}" "${KCHIP11}" \
  > kchip.hg19.intersection.positions \
  2> >(tee -a "${LOG}" >&2) <<'PY'
import csv
import sys

kchip10, kchip11 = sys.argv[1:]
autosomes = {str(i) for i in range(1, 23)}


def read_positions(path, chr_col, pos_col):
    positions = set()
    n_total = 0
    n_valid = 0
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(
            line for line in handle if line.strip() and not line.startswith("#")
        )
        if reader.fieldnames is None:
            raise RuntimeError(f"No table header detected in {path}")
        if chr_col not in reader.fieldnames or pos_col not in reader.fieldnames:
            raise RuntimeError(
                f"Expected columns {chr_col!r}, {pos_col!r}; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            n_total += 1
            chrom = (row.get(chr_col) or "").strip().replace("chr", "").replace("CHR", "")
            pos = (row.get(pos_col) or "").strip()
            if chrom not in autosomes:
                continue
            try:
                pos = int(pos)
            except ValueError:
                continue
            if pos <= 0:
                continue
            positions.add((int(chrom), pos))
            n_valid += 1
    print(
        f"[INFO] {path}: total_rows={n_total}, valid_autosomal_rows={n_valid}, "
        f"unique_positions={len(positions)}",
        file=sys.stderr,
    )
    return positions


pos10 = read_positions(kchip10, "Chr_id", "Start")
pos11 = read_positions(kchip11, "Chromosome", "Physical Position")
common = pos10 & pos11

if not common:
    raise RuntimeError("K-chip v1.0 ∩ v1.1 is empty; check genome builds/columns.")

print(f"[INFO] Common autosomal K-chip positions={len(common)}", file=sys.stderr)
for chrom, pos in sorted(common):
    print(f"{chrom}:{pos}")
PY
fi

need_file kchip.hg19.intersection.positions
echo "[INFO] K-chip intersection positions=$(wc -l < kchip.hg19.intersection.positions)" | tee -a "${LOG}"
else
  echo "[STEP 1] Skip K-chip position extraction (mode=all)" | tee -a "${LOG}"
fi


# ============================================================================
# STEP 2. KoGES BGEN -> PGEN (cached / shared)
# ============================================================================

echo "[STEP 2] Convert KoGES BGEN to PGEN" | tee -a "${LOG}"

# koges.raw is mode-independent; prefer SHARED_RAW_WD then SHARED_WD
if ! reuse_shared_pfile koges.raw; then
  if [[ ! -e koges.raw.pgen && -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
    for ext in pgen pvar psam; do
      ln -s "${SHARED_RAW_WD}/koges.raw.${ext}" "koges.raw.${ext}"
    done
    echo "[INFO] Linked koges.raw from SHARED_RAW_WD=${SHARED_RAW_WD}" | tee -a "${LOG}"
  elif [[ -f koges.raw.pgen && -f koges.raw.pvar && -f koges.raw.psam ]]; then
    echo "[INFO] Using cached koges.raw PGEN." | tee -a "${LOG}"
  else
    # Write into SHARED_RAW_WD for cross-mode reuse
    if [[ ! -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
      (
        cd "${SHARED_RAW_WD}"
        "${PLINK2}" \
          --threads "${THREADS}" \
          --bgen "${KOGES_BGEN}" ref-first \
          --sample "${KOGES_SAMP}" \
          --make-pgen \
          --out koges.raw 2>&1 | tee -a "${WD}/${LOG}"
      )
    fi
    if [[ -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
      for ext in pgen pvar psam; do
        [[ -e "koges.raw.${ext}" ]] || ln -s "${SHARED_RAW_WD}/koges.raw.${ext}" "koges.raw.${ext}"
      done
    else
      run "${PLINK2}" \
        --threads "${THREADS}" \
        --bgen "${KOGES_BGEN}" ref-first \
        --sample "${KOGES_SAMP}" \
        --make-pgen \
        --out koges.raw
    fi
  fi
fi

need_file koges.raw.pgen
need_file koges.raw.pvar
need_file koges.raw.psam
echo "[INFO] KoGES raw samples=$(count_psam koges.raw.psam), variants=$(count_pvar koges.raw.pvar)" | tee -a "${LOG}"


# ============================================================================
# STEP 3. KoGES K-chip intersection SNPs + basic marker QC
# ============================================================================

echo "[STEP 3] Select ${MARKER_DESC} KoGES SNPs" | tee -a "${LOG}"

if ! reuse_shared koges.kchip.qc.ids; then
  clear_local_target koges.kchip.qc.ids
  if (( USE_KCHIP == 1 )); then
    awk '
    BEGIN {
      FS="\t"
      while ((getline x < "kchip.hg19.intersection.positions") > 0) wanted[x]=1
    }
    $0 ~ /^#/ {next}
    {
      chr=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5)
      c=chr; sub(/^chr/, "", c)
      if (c !~ /^[0-9]+$/ || (c+0)<1 || (c+0)>22) next
      key=(c+0) ":" pos
      if (!(key in wanted)) next
      if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
      pair=ref alt
      if (pair=="AT" || pair=="TA" || pair=="CG" || pair=="GC") next
      if ((c+0)==6 && pos>=25000000 && pos<=34000000) next
      print id
    }
    ' koges.raw.pvar > koges.kchip.qc.ids
  else
    awk '
    BEGIN { FS="\t" }
    $0 ~ /^#/ {next}
    {
      chr=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5)
      c=chr; sub(/^chr/, "", c)
      if (c !~ /^[0-9]+$/ || (c+0)<1 || (c+0)>22) next
      if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
      pair=ref alt
      if (pair=="AT" || pair=="TA" || pair=="CG" || pair=="GC") next
      if ((c+0)==6 && pos>=25000000 && pos<=34000000) next
      print id
    }
    ' koges.raw.pvar > koges.kchip.qc.ids
  fi
fi

need_file koges.kchip.qc.ids
echo "[INFO] KoGES ${MARKER_DESC} SNPs after basic QC=$(wc -l < koges.kchip.qc.ids)" | tee -a "${LOG}"

if ! reuse_shared_pfile koges.kchip.tmp; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.raw \
    --extract koges.kchip.qc.ids \
    --make-pgen \
    --out koges.kchip.tmp
fi

if ! reuse_shared_pfile koges.kchip.ra; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.kchip.tmp \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 1 \
    --rm-dup exclude-all \
    --make-pgen \
    --out koges.kchip.ra
fi

need_file koges.kchip.ra.pgen
need_file koges.kchip.ra.pvar
need_file koges.kchip.ra.psam


# ============================================================================
# STEP 4. KoGES allele-aware candidate IDs for 1KG extraction
# ============================================================================

echo "[STEP 4] Build KoGES allele-aware candidate IDs" | tee -a "${LOG}"

if ! reuse_shared koges.kchip.alleles.both; then
  clear_local_target koges.kchip.alleles.both
  awk '
  $0 !~ /^#/ {
    split($3, a, ":")
    if (length(a[3]) != 1 || length(a[4]) != 1) next
    print a[1] ":" a[2] ":" a[3] ":" a[4]
    print a[1] ":" a[2] ":" a[4] ":" a[3]
  }
  ' koges.kchip.ra.pvar | sort -u > koges.kchip.alleles.both
fi

need_file koges.kchip.alleles.both


# ============================================================================
# STEP 5. Extract candidate SNPs from 1KG ALL
# ============================================================================

echo "[STEP 5] Extract candidate SNPs from 1KG ALL" | tee -a "${LOG}"

if ! reuse_shared_pfile g1k.kchip.ra; then
  rm -f g1k.merge_list

  for chr in $(seq 1 22); do
    if ! reuse_shared_pfile "g1k.chr${chr}"; then
      vcf=$(printf "${G1K_VCF_PATTERN}" "${chr}")
      need_file "${vcf}"
      run "${PLINK2}" \
        --threads "${THREADS}" \
        --vcf "${vcf}" \
        --snps-only just-acgt \
        --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1 \
        --rm-dup exclude-all \
        --extract koges.kchip.alleles.both \
        --make-pgen \
        --out "g1k.chr${chr}"
    fi
    echo "g1k.chr${chr}" >> g1k.merge_list
  done

  head -n 1 g1k.merge_list > g1k.first
  tail -n +2 g1k.merge_list > g1k.rest

  # Chromosome concatenation (same samples, different loci) is supported by
  # PLINK2 --pmerge-list.  Sample-wise merges (different samples, same loci)
  # are NOT supported in v2.00a3.7 and must use PLINK1.9 --bmerge instead.
  # This script never merges KoGES with 1KG; projection uses --score.
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile "$(cat g1k.first)" \
    --pmerge-list g1k.rest \
    --make-pgen \
    --out g1k.kchip.ra
fi

need_file g1k.kchip.ra.pgen
need_file g1k.kchip.ra.pvar
need_file g1k.kchip.ra.psam


# ============================================================================
# STEP 6. Exact unordered-allele match between KoGES and 1KG
# ============================================================================

echo "[STEP 6] Match exact KoGES-1KG allele pairs" | tee -a "${LOG}"

if [[ -f "${SHARED_WD}/koges.common.ids" && -f "${SHARED_WD}/g1k.common.ids" ]]; then
  reuse_shared koges.common.ids
  reuse_shared g1k.common.ids
else
  clear_local_target koges.common.ids
  clear_local_target g1k.common.ids
  python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def read_pvar(path):
    result = {}
    duplicates = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) < 5:
                continue
            chrom = x[0].replace("chr", "")
            ref, alt = x[3].upper(), x[4].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            a1, a2 = sorted((ref, alt))
            key = f"{int(chrom)}:{x[1]}:{a1}:{a2}"
            if key in result:
                duplicates.add(key)
            else:
                result[key] = x[2]
    for key in duplicates:
        result.pop(key, None)
    return result, duplicates


koges, koges_dup = read_pvar("koges.kchip.ra.pvar")
g1k, g1k_dup = read_pvar("g1k.kchip.ra.pvar")
common = sorted(set(koges) & set(g1k), key=lambda z: (int(z.split(":")[0]), int(z.split(":")[1]), z))

with open("koges.common.ids", "w") as fk, open("g1k.common.ids", "w") as fg:
    for key in common:
        fk.write(koges[key] + "\n")
        fg.write(g1k[key] + "\n")

print(f"[INFO] matched common SNPs={len(common)}", file=__import__("sys").stderr)
print(f"[INFO] duplicate canonical keys excluded: KoGES={len(koges_dup)}, 1KG={len(g1k_dup)}", file=__import__("sys").stderr)
PY
fi

assert_min_lines koges.common.ids "${MIN_COMMON_SNPS}" "Allele-matched KoGES-1KG SNPs"


# ============================================================================
# STEP 7. Extract first, then rename to chr:pos, then harmonize REF alleles
# ============================================================================

echo "[STEP 7] Extract common markers and harmonize KoGES REF" | tee -a "${LOG}"

# --set-all-var-ids is applied before --extract by PLINK2.  Therefore extract
# chr:pos:REF:ALT IDs first, and only then rename them to chr:pos.
if ! reuse_shared_pfile koges.common.extract; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.kchip.ra \
    --extract koges.common.ids \
    --make-pgen \
    --out koges.common.extract
fi

if ! reuse_shared_pfile koges.common.pre; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.common.extract \
    --set-all-var-ids '@:#' \
    --rm-dup exclude-all \
    --make-pgen \
    --out koges.common.pre
fi

if ! reuse_shared_pfile g1k.common.extract; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.kchip.ra \
    --extract g1k.common.ids \
    --make-pgen \
    --out g1k.common.extract
fi

if ! reuse_shared_pfile g1k.common; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.common.extract \
    --set-all-var-ids '@:#' \
    --rm-dup exclude-all \
    --make-pgen \
    --out g1k.common
fi

need_file koges.common.pre.pvar
need_file g1k.common.pvar

if ! reuse_shared g1k.ref.alleles; then
  clear_local_target g1k.ref.alleles
  awk '$0 !~ /^#/ {print $3 "\t" $4}' g1k.common.pvar > g1k.ref.alleles
fi

if ! reuse_shared_pfile koges.common; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.common.pre \
    --ref-allele force g1k.ref.alleles 2 1 \
    --make-pgen \
    --out koges.common
fi

python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def marker_map(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            out[x[2]] = (x[0].replace("chr", ""), x[1], x[3].upper(), x[4].upper())
    return out


koges = marker_map("koges.common.pvar")
g1k = marker_map("g1k.common.pvar")
if set(koges) != set(g1k):
    raise SystemExit(f"[ERROR] Common marker ID sets differ: KoGES={len(koges)}, 1KG={len(g1k)}")
mismatch = [vid for vid in koges if koges[vid] != g1k[vid]]
if mismatch:
    raise SystemExit(f"[ERROR] {len(mismatch)} CHR/POS/REF/ALT mismatches remain after harmonization; first={mismatch[:5]}")
print(f"[INFO] REF/ALT harmonization verified for {len(koges)} shared markers", file=__import__("sys").stderr)
PY


# ============================================================================
# STEP 7b. Define 1KG EAS reference samples (expected N=504)
# ============================================================================

echo "[STEP 7b] Build 1KG EAS keep list" | tee -a "${LOG}"

python3 - "${G1K_PANEL}" "${EXPECTED_1KG_EAS_N}" > g1k.eas.iids <<'PY'
import csv
import sys

path, expected = sys.argv[1], int(sys.argv[2])
with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise RuntimeError("1KG panel header not detected")
    lower = {x.lower().strip(): x for x in reader.fieldnames}
    if "sample" not in lower or "super_pop" not in lower:
        raise RuntimeError(f"Panel requires sample/super_pop columns; found {reader.fieldnames}")
    ids = [
        row[lower["sample"]].strip()
        for row in reader
        if row[lower["super_pop"]].strip() == "EAS" and row[lower["sample"]].strip()
    ]

if len(ids) != expected:
    raise RuntimeError(f"Expected {expected} EAS panel samples, found {len(ids)}")
for iid in ids:
    print(iid)
print(f"[INFO] 1KG EAS panel samples={len(ids)}", file=sys.stderr)
PY

python3 - <<'PY'
wanted = set(x.strip() for x in open("g1k.eas.iids") if x.strip())

with open("g1k.common.psam") as handle:
    for line in handle:
        if line.startswith("#"):
            header = line.rstrip().lstrip("#").split()
            break
    else:
        raise RuntimeError("PSAM header not found")

cols = {name: i for i, name in enumerate(header)}
if "IID" not in cols:
    raise RuntimeError(f"IID column not found: {header}")

found = 0
with open("g1k.common.psam") as source, open("g1k.eas.keep", "w") as out:
    for line in source:
        if line.startswith("#"):
            continue
        x = line.rstrip().split()
        if not x or x[cols["IID"]] not in wanted:
            continue
        if "FID" in cols:
            out.write(f"{x[cols['FID']]}\t{x[cols['IID']]}\n")
        else:
            out.write(f"{x[cols['IID']]}\n")
        found += 1

if found != len(wanted):
    raise RuntimeError(f"Only {found}/{len(wanted)} EAS IDs matched the 1KG PSAM")
print(f"[INFO] EAS samples matched to 1KG PGEN={found}")
PY

# Materialize the EAS-only reference once; all downstream reference QC, LD pruning,
# allele frequencies, and PCA are performed on these 504 EAS individuals only.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.common \
  --keep g1k.eas.keep \
  --make-pgen \
  --out g1k.eas.common

need_file g1k.eas.common.pgen
need_file g1k.eas.common.pvar
need_file g1k.eas.common.psam

EAS_N=$(count_psam g1k.eas.common.psam)
if (( EAS_N != EXPECTED_1KG_EAS_N )); then
  echo "[ERROR] 1KG EAS reference sample count=${EAS_N}, expected=${EXPECTED_1KG_EAS_N}" | tee -a "${LOG}"
  exit 1
fi
echo "[INFO] 1KG EAS reference samples=${EAS_N}" | tee -a "${LOG}"


# ============================================================================
# STEP 8. Reference-defined marker QC
# (Projection-specific from here: do NOT reuse joint PCA SNP lists)
# ============================================================================

echo "[STEP 8] Reference-defined MAF/missingness QC" | tee -a "${LOG}"

# 1KG EAS determines MAF and reference marker eligibility.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.common \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --write-snplist \
  --out g1k.eas.ref.maf05.geno02

# KoGES contributes only a technical-availability check.  Do not impose a
# KoGES MAF threshold, since KoGES must not define the reference PCA markers.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.common \
  --geno "${GENO}" dosage \
  --write-snplist \
  --out koges.geno02

sort -u g1k.eas.ref.maf05.geno02.snplist > g1k.eas.ref.maf05.geno02.sorted
sort -u koges.geno02.snplist > koges.geno02.sorted
comm -12 g1k.eas.ref.maf05.geno02.sorted koges.geno02.sorted > common.reference_qc.snps

assert_min_lines common.reference_qc.snps "${MIN_COMMON_SNPS}" "Reference-QC shared SNPs"


# ============================================================================
# STEP 9. LD pruning only in 1KG EAS
# ============================================================================

echo "[STEP 9] LD prune reference SNPs in 1KG EAS" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.common \
  --extract common.reference_qc.snps \
  --indep-pairwise 50 5 0.2 \
  --out g1k.eas.reference.prune

need_file g1k.eas.reference.prune.prune.in
assert_min_lines g1k.eas.reference.prune.prune.in "${MIN_PRUNED_SNPS}" "LD-pruned reference SNPs"


# ============================================================================
# STEP 10. Apply identical final markers to 1KG and KoGES
# ============================================================================

echo "[STEP 10] Create final reference and target PGENs" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.common \
  --extract g1k.eas.reference.prune.prune.in \
  --make-pgen \
  --out g1k.eas.reference.final

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.common \
  --extract g1k.eas.reference.prune.prune.in \
  --make-pgen \
  --out koges.projection.final

need_file g1k.eas.reference.final.pgen
need_file g1k.eas.reference.final.pvar
need_file g1k.eas.reference.final.psam
need_file koges.projection.final.pgen
need_file koges.projection.final.pvar
need_file koges.projection.final.psam

python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def marker_set(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            rows.append((x[2], x[0].replace("chr", ""), x[1], x[3].upper(), x[4].upper()))
    if len(rows) != len(set(rows)):
        raise RuntimeError(f"Duplicate final markers found in {path}")
    return set(rows)


ref = marker_set("g1k.eas.reference.final.pvar")
target = marker_set("koges.projection.final.pvar")
if ref != target:
    only_ref = sorted(ref - target)[:5]
    only_target = sorted(target - ref)[:5]
    raise SystemExit(
        f"[ERROR] Final marker sets differ: reference_only={only_ref}, target_only={only_target}"
    )
print(f"[INFO] Final marker set identity verified: n={len(ref)}", file=__import__("sys").stderr)
PY


# ============================================================================
# STEP 10b. Allele-frequency flip diagnostic against 1KG EAS
# ============================================================================

echo "[STEP 10b] KoGES-vs-1KG EAS allele-frequency check" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.projection.final \
  --freq \
  --out koges.projection.af

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.reference.final \
  --freq \
  --out g1k.eas.reference.af

python3 - "${AF_DIFF_WARN}" "${AF_FLIP_MATCH}" "${MAX_SUSPECT_FLIPS_ABORT}" \
  2> >(tee -a "${LOG}" >&2) <<'PY'
import sys

diff_warn, flip_match, max_abort = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])

def read_af(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) < 5:
                continue
            out[x[1]] = (x[2].upper(), x[3].upper(), float(x[4].split(",")[0]))
    return out

koges = read_af("koges.projection.af.afreq")
eas = read_af("g1k.eas.reference.af.afreq")
suspect = []
allele_errors = []

for vid in sorted(set(koges) & set(eas)):
    kr, ka, kaf = koges[vid]
    er, ea, eaf = eas[vid]
    if (kr, ka) != (er, ea):
        allele_errors.append((vid, kr, ka, er, ea))
        continue
    d = abs(kaf - eaf)
    d_flip = abs(kaf - (1.0 - eaf))
    if d >= diff_warn and d_flip <= flip_match:
        suspect.append((vid, kaf, eaf, d, d_flip))

with open("af_flip_check.suspect.tsv", "w") as out:
    out.write("ID\tKoGES_ALT_AF\tEAS_ALT_AF\tABS_DIFF\tABS_DIFF_TO_1_MINUS_EAS\n")
    for row in suspect:
        out.write("\t".join(map(str, row)) + "\n")

with open("af_flip_check.exclude.ids", "w") as out:
    for vid, *_ in suspect:
        out.write(vid + "\n")

if allele_errors:
    raise SystemExit(f"[ERROR] REF/ALT mismatch in AF files; first={allele_errors[:5]}")

if len(suspect) > max_abort:
    raise SystemExit(
        f"[ERROR] {len(suspect)} suspect allele flips "
        f"(> abort threshold {max_abort}); see af_flip_check.suspect.tsv"
    )

if suspect:
    for vid, kaf, eaf, d, d_flip in suspect[:20]:
        print(
            f"[WARN] Suspect allele flip (will exclude): {vid} "
            f"KoGES_AF={kaf:.4f} EAS_AF={eaf:.4f} "
            f"|diff|={d:.4f} |diff_to_1-AF|={d_flip:.4f}",
            file=sys.stderr,
        )
    print(
        f"[INFO] Excluding {len(suspect)} suspect allele-flip SNP(s). "
        f"See af_flip_check.exclude.ids",
        file=sys.stderr,
    )
else:
    print("[INFO] AF flip check passed: suspect_flips=0", file=sys.stderr)
PY

# Exclude rare residual allele-flip SNPs from both final datasets.
if [[ -s af_flip_check.exclude.ids ]]; then
  NEX=$(wc -l < af_flip_check.exclude.ids | tr -d ' ')
  echo "[INFO] Removing ${NEX} suspect-flip SNP(s) from EAS reference/KoGES target finals" \
    | tee -a "${LOG}"

  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.eas.reference.final \
    --exclude af_flip_check.exclude.ids \
    --make-pgen \
    --out g1k.eas.reference.final.tmp

  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.projection.final \
    --exclude af_flip_check.exclude.ids \
    --make-pgen \
    --out koges.projection.final.tmp

  for ext in pgen pvar psam; do
    mv -f "g1k.eas.reference.final.tmp.${ext}" "g1k.eas.reference.final.${ext}"
    mv -f "koges.projection.final.tmp.${ext}" "koges.projection.final.${ext}"
  done
  rm -f g1k.eas.reference.final.tmp.log koges.projection.final.tmp.log

  echo "[INFO] After flip exclusion: 1KG_EAS=$(count_pvar g1k.eas.reference.final.pvar) KoGES=$(count_pvar koges.projection.final.pvar)" \
    | tee -a "${LOG}"
fi


# ============================================================================
# STEP 11. 1KG EAS reference allele frequencies and PCA
# ============================================================================

echo "[STEP 11] Compute 1KG EAS-only reference PCA" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.reference.final \
  --freq counts \
  --out g1k.eas.reference

need_file g1k.eas.reference.acount

# Exact PCA is practical for 504 1KG EAS samples.  meanimpute makes the
# reference PCA standardization consistent with the projection score command.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.reference.final \
  --read-freq g1k.eas.reference.acount \
  --pca allele-wts "${PCA_N}" meanimpute vcols=chrom,ref,alt \
  --out g1k.eas.reference.pca

need_file g1k.eas.reference.pca.eigenvec
need_file g1k.eas.reference.pca.eigenval
need_file g1k.eas.reference.pca.eigenvec.allele


# ============================================================================
# STEP 12. Reference self-projection and KoGES projection
# ============================================================================

echo "[STEP 12] Project 1KG (calibration) and KoGES onto reference PCs" | tee -a "${LOG}"

SCORE_COLS="6-$((PCA_N + 5))"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.eas.reference.final \
  --read-freq g1k.eas.reference.acount \
  --score g1k.eas.reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums "${SCORE_COLS}" \
  --out g1k.eas.reference.selfproject

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.projection.final \
  --read-freq g1k.eas.reference.acount \
  --score g1k.eas.reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums "${SCORE_COLS}" \
  --out koges.projected

need_file g1k.eas.reference.selfproject.sscore
need_file koges.projected.sscore


# ============================================================================
# STEP 13. Scale-calibrate projected scores and write plot-ready coordinates
# ============================================================================

echo "[STEP 13] Calibrate projected PCs and write coordinate table" | tee -a "${LOG}"

if skip_if_file pca_1kg_eas_reference.coordinates.tsv; then
  echo "[INFO] Coordinates already present; skip calibration" | tee -a "${LOG}"
else

python3 - "${PCA_N}" "${MIN_SELFPROJ_COR}" "${G1K_PANEL}" <<'PY' 2> >(tee -a "${LOG}" >&2)
import csv
import math
import sys

n_pc = int(sys.argv[1])
min_corr = float(sys.argv[2])
panel_path = sys.argv[3]


def read_ws_table(path):
    header = None
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if header is None:
                header = line.rstrip().lstrip("#").split()
                continue
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) != len(header):
                raise RuntimeError(f"Column count mismatch in {path}: {x[:5]}")
            rows.append(x)
    if header is None:
        raise RuntimeError(f"No header found in {path}")
    return header, rows


def column_index(header, name):
    try:
        return header.index(name)
    except ValueError as exc:
        raise RuntimeError(f"Column {name!r} absent. Header={header}") from exc


def row_key(header, row):
    iid = row[column_index(header, "IID")]
    fid = row[column_index(header, "FID")] if "FID" in header else "0"
    return (fid, iid)


def pc_columns(header, count, score=False):
    result = []
    for k in range(1, count + 1):
        choices = [f"PC{k}"]
        if score:
            choices = [f"PC{k}_SUM", f"PC{k}_AVG", f"PC{k}"]
        for choice in choices:
            if choice in header:
                result.append(column_index(header, choice))
                break
        else:
            raise RuntimeError(f"Cannot find PC{k} column in header={header}")
    return result


eig_header, eig_rows = read_ws_table("g1k.eas.reference.pca.eigenvec")
self_header, self_rows = read_ws_table("g1k.eas.reference.selfproject.sscore")
koges_header, koges_rows = read_ws_table("koges.projected.sscore")

eig_pc_idx = pc_columns(eig_header, n_pc, score=False)
self_pc_idx = pc_columns(self_header, n_pc, score=True)
koges_pc_idx = pc_columns(koges_header, n_pc, score=True)

eig = {row_key(eig_header, row): [float(row[i]) for i in eig_pc_idx] for row in eig_rows}
self_score = {row_key(self_header, row): [float(row[i]) for i in self_pc_idx] for row in self_rows}
koges_score = {row_key(koges_header, row): [float(row[i]) for i in koges_pc_idx] for row in koges_rows}

shared = sorted(set(eig) & set(self_score))
if len(shared) != len(eig):
    raise RuntimeError(f"Reference self-projection matched {len(shared)}/{len(eig)} samples")

scales = []
diagnostics = []
for pc in range(n_pc):
    x = [self_score[key][pc] for key in shared]
    y = [eig[key][pc] for key in shared]
    sxx = sum(value * value for value in x)
    syy = sum(value * value for value in y)
    sxy = sum(a * b for a, b in zip(x, y))
    if sxx <= 0 or syy <= 0:
        raise RuntimeError(f"Degenerate self-projection for PC{pc + 1}")
    scale = sxy / sxx
    corr = sxy / math.sqrt(sxx * syy)
    if abs(corr) < min_corr:
        raise RuntimeError(
            f"Self-projection calibration failed for PC{pc + 1}: r={corr:.8f} < {min_corr}"
        )
    scales.append(scale)
    diagnostics.append((pc + 1, len(shared), scale, corr))


def read_panel(path):
    out = {}
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("1KG panel header not detected")
        lower = {x.lower().strip(): x for x in reader.fieldnames}
        for required in ("sample", "pop", "super_pop"):
            if required not in lower:
                raise RuntimeError(f"Panel lacks {required!r}: {reader.fieldnames}")
        for row in reader:
            iid = row[lower["sample"]].strip()
            if iid:
                out[iid] = (row[lower["pop"]].strip(), row[lower["super_pop"]].strip())
    return out


panel = read_panel(panel_path)

with open("pca_1kg_eas_reference.scaling.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["PC", "N_REFERENCE", "SCALE_TO_DIRECT_1KG_PC", "SELF_PROJECTION_CORRELATION"])
    for row in diagnostics:
        writer.writerow(row)

with open("pca_1kg_eas_reference.coordinates.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["FID", "IID", "SOURCE", "POP", "SUPERPOP"] + [f"PC{k}" for k in range(1, n_pc + 1)])

    for fid, iid in sorted(eig):
        pop, superpop = panel.get(iid, ("NA", "NA"))
        writer.writerow([fid, iid, "1KG_EAS", pop, superpop] + eig[(fid, iid)])

    for fid, iid in sorted(koges_score):
        values = [koges_score[(fid, iid)][pc] * scales[pc] for pc in range(n_pc)]
        writer.writerow([fid, iid, "KoGES", "KoGES", "EAS"] + values)

print(
    f"[INFO] Wrote coordinates: 1KG_EAS={len(eig)}, KoGES={len(koges_score)}, PCs={n_pc}",
    file=sys.stderr,
)
PY
fi

need_file pca_1kg_eas_reference.coordinates.tsv
need_file pca_1kg_eas_reference.scaling.tsv


# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo "============================================================" | tee -a "${LOG}"
echo "[FINAL SUMMARY]" | tee -a "${LOG}"
if (( USE_KCHIP == 1 )); then
  echo "[RESULT] K-chip v1.0 ∩ v1.1 positions=$(wc -l < kchip.hg19.intersection.positions)" | tee -a "${LOG}"
else
  echo "[RESULT] K-chip filter skipped (mode=all)" | tee -a "${LOG}"
fi
echo "[RESULT] Final LD-pruned reference SNPs=$(count_pvar g1k.eas.reference.final.pvar)" | tee -a "${LOG}"
echo "[RESULT] 1KG EAS reference samples=$(count_psam g1k.eas.reference.final.psam)" | tee -a "${LOG}"
echo "[RESULT] KoGES projected samples=$(count_psam koges.projection.final.psam)" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/pca_1kg_eas_reference.coordinates.tsv" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/pca_1kg_eas_reference.scaling.tsv" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/${LOG}" | tee -a "${LOG}"
echo "[PLOT] Rscript ${PLOT_RSCRIPT} all" | tee -a "${LOG}"
Rscript "${PLOT_RSCRIPT}" all \
  "${WD}/pca_1kg_eas_reference.coordinates.tsv"
echo "[ALL DONE]" | tee -a "${LOG}"
__EAS_BODY__
  exit $?
fi

_run_embedded all <<'__ALL_BODY__'
set -euo pipefail

# ============================================================================
# Supplementary Figure 4: KoGES projected onto 1000G ALL PCs (genome-wide imputed)
# Usage: bash SuppleFigure4.sh
# ============================================================================
#
# PURPOSE
# -------
# This script creates a 1000 Genomes Phase 3 (1KG ALL) reference PCA and
# projects all KoGES samples onto those fixed PC axes
# (1kg_reference_projection_all).
#
# Contrast with joint_pca_kchip_variant.sh
# ---------------------------------------
# joint_pca_kchip_variant.sh:
#   KoGES + 1KG ALL joint PCA (KoGES influences the axes).
#
# This script:
#   PCA axes from 1KG ALL only; KoGES is projected onto those axes.
#
# Key design choices
# ------------------
# - The available KoGES input is imputed BGEN.  K-chip manifests are used only
#   to select positions targeted by BOTH KoreanChip v1.0 and v1.1.  This is not
#   a reconstruction of directly genotyped calls.
# - PCA axes, reference allele frequencies, MAF filtering, and LD pruning are
#   determined by 1KG ALL.  KoGES therefore cannot dominate the PCA axes.
# - KoGES is used only to retain technically available shared markers
#   (--geno 0.02 dosage) and for projection.  KoGES MAF is NOT used to select
#   reference PCA markers.
# - All KoGES samples, including related samples, are projected.  They do not
#   influence the 1KG PCA axes.  For a formal KoGES centroid/distance analysis,
#   use an unrelated subset separately.
# - The final coordinates file contains direct 1KG PC coordinates and
#   scale-calibrated KoGES projected coordinates in the same PC scale.
#
# Shared intermediates
# --------------------
# Steps 1-7 match joint_pca_kchip_variant.sh (same inputs and allele
# harmonization through koges.common / g1k.common).  When those files already
# exist under SHARED_WD, this script reuses them via symlink instead of
# recomputing.  Projection-specific QC, PCA, and outputs live under WD.
#
# Genome build
# ------------
# K-chip v1.0/v1.1, KoGES BGEN, and 1KG Phase 3 VCFs must all be GRCh37/hg19.
# Verify this once from the supplied manifests before running the full pipeline.
#
# Outputs
# -------
# pca_1kg_reference.coordinates.tsv
#   Ready-to-plot table: direct 1KG PCs + scaled projected KoGES PCs.
# pca_1kg_reference.scaling.tsv
#   Self-projection calibration diagnostics for each PC.
#
# ============================================================================


# ============================================================================
# SETTINGS
# ============================================================================

THREADS=30
PLINK2=/bin2/plink2

# Reference-marker QC
MAF=0.05
GENO=0.02

# PCA / projection settings
PCA_N=20
MIN_SELFPROJ_COR=0.999

# Allele-frequency flip check: KoGES vs 1KG EAS
# Few residual flips are excluded; many flips abort as systematic.
AF_DIFF_WARN=0.20
AF_FLIP_MATCH=0.05
MAX_SUSPECT_FLIPS_ABORT=200
EXPECTED_1KG_EAS_N=504

# Engineering guards, not biological cutoffs
MIN_COMMON_SNPS=10000
MIN_PRUNED_SNPS=5000

# Genome-wide imputed markers only (1kg_reference_projection_all)
MODE_ARG=all
USE_KCHIP=0

BASE_WD=/data4/sjim/openKoGES_revision/DATA_projectionPCA
SHARED_RAW_WD="${BASE_WD}"
SHARED_WD="${BASE_WD}/all_variants"
WD="${BASE_WD}/1kg_reference_projection_all"
LOG=koges_1kg_reference_projection_all.log
MARKER_DESC="genome-wide imputed (no K-chip filter)"
PLOT_OUTPREF=Figure2A_projection_all
PLOT_RSCRIPT=/home2/sjim/openKoGES/CODE_manuscript_final/2_Population_structure_relatedness_allele_frequency/SuppleFigure4.R

mkdir -p "${SHARED_WD}"
mkdir -p "${WD}"
cd "${WD}"

: > "${LOG}"


# ============================================================================
# INPUTS
# ============================================================================

KOGES_BGEN=/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.bgen
KOGES_SAMP=/data6/koges/bgen/CODA_ID_KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_Filtered_HWE1e6_R0.8.sample

G1K_DIR=/data/1000G/Phase3_v5a
G1K_VCF_PATTERN="${G1K_DIR}/chr%s.1kg.phase3.v5a.vcf.gz"
G1K_PANEL=/home2/sjim/openKoGES/DATA_baseline/forPCA/1000Gsamples.ALL.panel


# ============================================================================
# HELPERS
# ============================================================================

# Skip only the product type implied by the PLINK flags.
_run_should_skip() {
  local out="$1"
  shift
  local make_bed=0 make_pgen=0 do_pca=0 do_score=0 do_prune=0 write_snplist=0 do_freq=0
  local a
  for a in "$@"; do
    case "${a}" in
      --make-bed) make_bed=1 ;;
      --make-pgen) make_pgen=1 ;;
      --pca) do_pca=1 ;;
      --score) do_score=1 ;;
      --indep-pairwise) do_prune=1 ;;
      --write-snplist) write_snplist=1 ;;
      --freq) do_freq=1 ;;
    esac
  done

  if (( make_bed == 1 )); then
    [[ -f "${out}.bed" && -f "${out}.bim" && -f "${out}.fam" ]] || return 1
    echo "[INFO] Skip (exists bed): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( make_pgen == 1 )); then
    [[ -f "${out}.pgen" && -f "${out}.pvar" && -f "${out}.psam" ]] || return 1
    echo "[INFO] Skip (exists pfile): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_pca == 1 )); then
    [[ -f "${out}.eigenvec" ]] || return 1
    echo "[INFO] Skip (exists eigenvec): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_score == 1 )); then
    [[ -f "${out}.sscore" ]] || return 1
    echo "[INFO] Skip (exists sscore): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_prune == 1 )); then
    [[ -f "${out}.prune.in" ]] || return 1
    echo "[INFO] Skip (exists prune): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( write_snplist == 1 )); then
    [[ -f "${out}.snplist" ]] || return 1
    echo "[INFO] Skip (exists snplist): ${out}" | tee -a "${LOG}"; return 0
  fi
  if (( do_freq == 1 )); then
    [[ -f "${out}.acount" || -f "${out}.afreq" ]] || return 1
    echo "[INFO] Skip (exists freq): ${out}" | tee -a "${LOG}"; return 0
  fi
  return 1
}

skip_if_file() {
  local f="$1"
  if [[ -f "${f}" ]]; then
    echo "[INFO] Skip (exists): ${f}" | tee -a "${LOG}"
    return 0
  fi
  return 1
}

run() {
  local out="" prev=""
  for a in "$@"; do
    if [[ "${prev}" == "--out" ]]; then out="${a}"; fi
    prev="${a}"
  done
  if [[ -n "${out}" ]] && _run_should_skip "${out}" "$@"; then
    return 0
  fi
  echo -n "[CMD]" | tee -a "${LOG}"
  for a in "$@"; do
    echo -n " ${a}" | tee -a "${LOG}"
  done
  echo "" | tee -a "${LOG}"
  "$@" 2>&1 | tee -a "${LOG}"
}

need_file() {
  local f="$1"
  [[ -f "${f}" ]] || {
    echo "[ERROR] Missing required file: ${f}" | tee -a "${LOG}"
    exit 1
  }
}

need_exec() {
  local x="$1"
  [[ -x "${x}" ]] || {
    echo "[ERROR] Not executable: ${x}" | tee -a "${LOG}"
    exit 1
  }
}

count_pvar() {
  awk '$0 !~ /^#/ {n++} END {print n+0}' "$1"
}

count_psam() {
  awk '$0 !~ /^#/ {n++} END {print n+0}' "$1"
}

assert_min_lines() {
  local f="$1"
  local min="$2"
  local label="$3"
  local n
  n=$(wc -l < "${f}" | tr -d ' ')

  if (( n < min )); then
    echo "[ERROR] ${label}: ${n} < minimum ${min}" | tee -a "${LOG}"
    exit 1
  fi
  echo "[INFO] ${label}: ${n} (min ${min} OK)" | tee -a "${LOG}"
}

# Reuse an identical intermediate from SHARED_WD via symlink when present.
# Usage: reuse_shared <local_name> [shared_name]
# Returns 0 if reused, 1 if missing in SHARED_WD (caller should compute).
reuse_shared() {
  local local_name="$1"
  local shared_name="${2:-$1}"
  local src="${SHARED_WD}/${shared_name}"
  local dst="${WD}/${local_name}"

  if [[ -e "${dst}" || -L "${dst}" ]]; then
    echo "[INFO] Reusing existing local file: ${local_name}" | tee -a "${LOG}"
    return 0
  fi

  if [[ -f "${src}" ]]; then
    ln -s "${src}" "${dst}"
    echo "[INFO] Linked shared intermediate: ${shared_name} -> ${local_name}" \
      | tee -a "${LOG}"
    return 0
  fi

  return 1
}

# Drop a local path so a later write cannot clobber SHARED_WD through a symlink.
clear_local_target() {
  local path="$1"
  if [[ -L "${path}" ]]; then
    rm -f "${path}"
  elif [[ -e "${path}" ]]; then
    rm -f "${path}"
  fi
}

# Reuse a PLINK2 prefix trio (.pgen/.pvar/.psam) from SHARED_WD.
# Links only when all three shared files exist (avoids partial reuse).
reuse_shared_pfile() {
  local local_prefix="$1"
  local shared_prefix="${2:-$1}"
  local ext
  local all_local=1
  local all_shared=1

  for ext in pgen pvar psam; do
    if [[ ! -e "${WD}/${local_prefix}.${ext}" && ! -L "${WD}/${local_prefix}.${ext}" ]]; then
      all_local=0
    fi
    if [[ ! -f "${SHARED_WD}/${shared_prefix}.${ext}" ]]; then
      all_shared=0
    fi
  done

  if (( all_local == 1 )); then
    echo "[INFO] Reusing existing local pfile: ${local_prefix}" | tee -a "${LOG}"
    return 0
  fi

  if (( all_shared == 1 )); then
    for ext in pgen pvar psam; do
      clear_local_target "${WD}/${local_prefix}.${ext}"
      ln -s "${SHARED_WD}/${shared_prefix}.${ext}" "${WD}/${local_prefix}.${ext}"
    done
    echo "[INFO] Linked shared pfile: ${shared_prefix} -> ${local_prefix}" \
      | tee -a "${LOG}"
    return 0
  fi

  return 1
}


# ============================================================================
# PRE-FLIGHT
# ============================================================================

echo "============================================================" | tee -a "${LOG}"
echo "[START] KoGES projection onto 1KG ALL (${MARKER_DESC}, mode=${MODE_ARG})" | tee -a "${LOG}"
echo "============================================================" | tee -a "${LOG}"
echo "[INFO] SHARED_WD=${SHARED_WD}" | tee -a "${LOG}"
echo "[INFO] WD=${WD}" | tee -a "${LOG}"
echo "[INFO] THREADS=${THREADS} MAF=${MAF} GENO=${GENO} PCA_N=${PCA_N}" | tee -a "${LOG}"

need_exec "${PLINK2}"
need_file "${KOGES_BGEN}"
need_file "${KOGES_SAMP}"
need_file "${G1K_PANEL}"
if (( USE_KCHIP == 1 )); then
  need_file "${KCHIP10}"
  need_file "${KCHIP11}"
fi

if [[ -f "pca_1kg_reference.coordinates.tsv" ]]; then
  echo "[INFO] Final coordinates already present for mode=${MODE_ARG}; skip PCA steps." | tee -a "${LOG}"
  echo "[OUTPUT] ${WD}/pca_1kg_reference.coordinates.tsv" | tee -a "${LOG}"
  echo "[PLOT] Rscript ${PLOT_RSCRIPT} all" | tee -a "${LOG}"
  Rscript "${PLOT_RSCRIPT}" all \
    "${WD}/pca_1kg_reference.coordinates.tsv"
  echo "[ALL DONE]" | tee -a "${LOG}"
  exit 0
fi


# ============================================================================
# STEP 1. K-chip v1.0 ∩ v1.1 autosomal positions
# ============================================================================

if (( USE_KCHIP == 1 )); then
echo "[STEP 1] Extract K-chip v1.0 ∩ v1.1 positions" | tee -a "${LOG}"

if ! reuse_shared kchip.hg19.intersection.positions kchip.hg19.positions; then
  clear_local_target kchip.hg19.intersection.positions
  python3 - "${KCHIP10}" "${KCHIP11}" \
  > kchip.hg19.intersection.positions \
  2> >(tee -a "${LOG}" >&2) <<'PY'
import csv
import sys

kchip10, kchip11 = sys.argv[1:]
autosomes = {str(i) for i in range(1, 23)}


def read_positions(path, chr_col, pos_col):
    positions = set()
    n_total = 0
    n_valid = 0
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(
            line for line in handle if line.strip() and not line.startswith("#")
        )
        if reader.fieldnames is None:
            raise RuntimeError(f"No table header detected in {path}")
        if chr_col not in reader.fieldnames or pos_col not in reader.fieldnames:
            raise RuntimeError(
                f"Expected columns {chr_col!r}, {pos_col!r}; "
                f"found {reader.fieldnames}"
            )
        for row in reader:
            n_total += 1
            chrom = (row.get(chr_col) or "").strip().replace("chr", "").replace("CHR", "")
            pos = (row.get(pos_col) or "").strip()
            if chrom not in autosomes:
                continue
            try:
                pos = int(pos)
            except ValueError:
                continue
            if pos <= 0:
                continue
            positions.add((int(chrom), pos))
            n_valid += 1
    print(
        f"[INFO] {path}: total_rows={n_total}, valid_autosomal_rows={n_valid}, "
        f"unique_positions={len(positions)}",
        file=sys.stderr,
    )
    return positions


pos10 = read_positions(kchip10, "Chr_id", "Start")
pos11 = read_positions(kchip11, "Chromosome", "Physical Position")
common = pos10 & pos11

if not common:
    raise RuntimeError("K-chip v1.0 ∩ v1.1 is empty; check genome builds/columns.")

print(f"[INFO] Common autosomal K-chip positions={len(common)}", file=sys.stderr)
for chrom, pos in sorted(common):
    print(f"{chrom}:{pos}")
PY
fi

need_file kchip.hg19.intersection.positions
echo "[INFO] K-chip intersection positions=$(wc -l < kchip.hg19.intersection.positions)" | tee -a "${LOG}"
else
  echo "[STEP 1] Skip K-chip position extraction (mode=all)" | tee -a "${LOG}"
fi


# ============================================================================
# STEP 2. KoGES BGEN -> PGEN (cached / shared)
# ============================================================================

echo "[STEP 2] Convert KoGES BGEN to PGEN" | tee -a "${LOG}"

# koges.raw is mode-independent; prefer SHARED_RAW_WD then SHARED_WD
if ! reuse_shared_pfile koges.raw; then
  if [[ ! -e koges.raw.pgen && -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
    for ext in pgen pvar psam; do
      ln -s "${SHARED_RAW_WD}/koges.raw.${ext}" "koges.raw.${ext}"
    done
    echo "[INFO] Linked koges.raw from SHARED_RAW_WD=${SHARED_RAW_WD}" | tee -a "${LOG}"
  elif [[ -f koges.raw.pgen && -f koges.raw.pvar && -f koges.raw.psam ]]; then
    echo "[INFO] Using cached koges.raw PGEN." | tee -a "${LOG}"
  else
    # Write into SHARED_RAW_WD for cross-mode reuse
    if [[ ! -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
      (
        cd "${SHARED_RAW_WD}"
        "${PLINK2}" \
          --threads "${THREADS}" \
          --bgen "${KOGES_BGEN}" ref-first \
          --sample "${KOGES_SAMP}" \
          --make-pgen \
          --out koges.raw 2>&1 | tee -a "${WD}/${LOG}"
      )
    fi
    if [[ -f "${SHARED_RAW_WD}/koges.raw.pgen" ]]; then
      for ext in pgen pvar psam; do
        [[ -e "koges.raw.${ext}" ]] || ln -s "${SHARED_RAW_WD}/koges.raw.${ext}" "koges.raw.${ext}"
      done
    else
      run "${PLINK2}" \
        --threads "${THREADS}" \
        --bgen "${KOGES_BGEN}" ref-first \
        --sample "${KOGES_SAMP}" \
        --make-pgen \
        --out koges.raw
    fi
  fi
fi

need_file koges.raw.pgen
need_file koges.raw.pvar
need_file koges.raw.psam
echo "[INFO] KoGES raw samples=$(count_psam koges.raw.psam), variants=$(count_pvar koges.raw.pvar)" | tee -a "${LOG}"


# ============================================================================
# STEP 3. KoGES K-chip intersection SNPs + basic marker QC
# ============================================================================

echo "[STEP 3] Select ${MARKER_DESC} KoGES SNPs" | tee -a "${LOG}"

if ! reuse_shared koges.kchip.qc.ids; then
  clear_local_target koges.kchip.qc.ids
  if (( USE_KCHIP == 1 )); then
    awk '
    BEGIN {
      FS="\t"
      while ((getline x < "kchip.hg19.intersection.positions") > 0) wanted[x]=1
    }
    $0 ~ /^#/ {next}
    {
      chr=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5)
      c=chr; sub(/^chr/, "", c)
      if (c !~ /^[0-9]+$/ || (c+0)<1 || (c+0)>22) next
      key=(c+0) ":" pos
      if (!(key in wanted)) next
      if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
      pair=ref alt
      if (pair=="AT" || pair=="TA" || pair=="CG" || pair=="GC") next
      if ((c+0)==6 && pos>=25000000 && pos<=34000000) next
      print id
    }
    ' koges.raw.pvar > koges.kchip.qc.ids
  else
    awk '
    BEGIN { FS="\t" }
    $0 ~ /^#/ {next}
    {
      chr=$1; pos=$2; id=$3; ref=toupper($4); alt=toupper($5)
      c=chr; sub(/^chr/, "", c)
      if (c !~ /^[0-9]+$/ || (c+0)<1 || (c+0)>22) next
      if (ref !~ /^[ACGT]$/ || alt !~ /^[ACGT]$/) next
      pair=ref alt
      if (pair=="AT" || pair=="TA" || pair=="CG" || pair=="GC") next
      if ((c+0)==6 && pos>=25000000 && pos<=34000000) next
      print id
    }
    ' koges.raw.pvar > koges.kchip.qc.ids
  fi
fi

need_file koges.kchip.qc.ids
echo "[INFO] KoGES ${MARKER_DESC} SNPs after basic QC=$(wc -l < koges.kchip.qc.ids)" | tee -a "${LOG}"

if ! reuse_shared_pfile koges.kchip.tmp; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.raw \
    --extract koges.kchip.qc.ids \
    --make-pgen \
    --out koges.kchip.tmp
fi

if ! reuse_shared_pfile koges.kchip.ra; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.kchip.tmp \
    --set-all-var-ids '@:#:$r:$a' \
    --new-id-max-allele-len 1 \
    --rm-dup exclude-all \
    --make-pgen \
    --out koges.kchip.ra
fi

need_file koges.kchip.ra.pgen
need_file koges.kchip.ra.pvar
need_file koges.kchip.ra.psam


# ============================================================================
# STEP 4. KoGES allele-aware candidate IDs for 1KG extraction
# ============================================================================

echo "[STEP 4] Build KoGES allele-aware candidate IDs" | tee -a "${LOG}"

if ! reuse_shared koges.kchip.alleles.both; then
  clear_local_target koges.kchip.alleles.both
  awk '
  $0 !~ /^#/ {
    split($3, a, ":")
    if (length(a[3]) != 1 || length(a[4]) != 1) next
    print a[1] ":" a[2] ":" a[3] ":" a[4]
    print a[1] ":" a[2] ":" a[4] ":" a[3]
  }
  ' koges.kchip.ra.pvar | sort -u > koges.kchip.alleles.both
fi

need_file koges.kchip.alleles.both


# ============================================================================
# STEP 5. Extract candidate SNPs from 1KG ALL
# ============================================================================

echo "[STEP 5] Extract candidate SNPs from 1KG ALL" | tee -a "${LOG}"

if ! reuse_shared_pfile g1k.kchip.ra; then
  rm -f g1k.merge_list

  for chr in $(seq 1 22); do
    if ! reuse_shared_pfile "g1k.chr${chr}"; then
      vcf=$(printf "${G1K_VCF_PATTERN}" "${chr}")
      need_file "${vcf}"
      run "${PLINK2}" \
        --threads "${THREADS}" \
        --vcf "${vcf}" \
        --snps-only just-acgt \
        --max-alleles 2 \
        --set-all-var-ids '@:#:$r:$a' \
        --new-id-max-allele-len 1 \
        --rm-dup exclude-all \
        --extract koges.kchip.alleles.both \
        --make-pgen \
        --out "g1k.chr${chr}"
    fi
    echo "g1k.chr${chr}" >> g1k.merge_list
  done

  head -n 1 g1k.merge_list > g1k.first
  tail -n +2 g1k.merge_list > g1k.rest

  # Chromosome concatenation (same samples, different loci) is supported by
  # PLINK2 --pmerge-list.  Sample-wise merges (different samples, same loci)
  # are NOT supported in v2.00a3.7 and must use PLINK1.9 --bmerge instead.
  # This script never merges KoGES with 1KG; projection uses --score.
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile "$(cat g1k.first)" \
    --pmerge-list g1k.rest \
    --make-pgen \
    --out g1k.kchip.ra
fi

need_file g1k.kchip.ra.pgen
need_file g1k.kchip.ra.pvar
need_file g1k.kchip.ra.psam


# ============================================================================
# STEP 6. Exact unordered-allele match between KoGES and 1KG
# ============================================================================

echo "[STEP 6] Match exact KoGES-1KG allele pairs" | tee -a "${LOG}"

if [[ -f "${SHARED_WD}/koges.common.ids" && -f "${SHARED_WD}/g1k.common.ids" ]]; then
  reuse_shared koges.common.ids
  reuse_shared g1k.common.ids
else
  clear_local_target koges.common.ids
  clear_local_target g1k.common.ids
  python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def read_pvar(path):
    result = {}
    duplicates = set()
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) < 5:
                continue
            chrom = x[0].replace("chr", "")
            ref, alt = x[3].upper(), x[4].upper()
            if len(ref) != 1 or len(alt) != 1:
                continue
            a1, a2 = sorted((ref, alt))
            key = f"{int(chrom)}:{x[1]}:{a1}:{a2}"
            if key in result:
                duplicates.add(key)
            else:
                result[key] = x[2]
    for key in duplicates:
        result.pop(key, None)
    return result, duplicates


koges, koges_dup = read_pvar("koges.kchip.ra.pvar")
g1k, g1k_dup = read_pvar("g1k.kchip.ra.pvar")
common = sorted(set(koges) & set(g1k), key=lambda z: (int(z.split(":")[0]), int(z.split(":")[1]), z))

with open("koges.common.ids", "w") as fk, open("g1k.common.ids", "w") as fg:
    for key in common:
        fk.write(koges[key] + "\n")
        fg.write(g1k[key] + "\n")

print(f"[INFO] matched common SNPs={len(common)}", file=__import__("sys").stderr)
print(f"[INFO] duplicate canonical keys excluded: KoGES={len(koges_dup)}, 1KG={len(g1k_dup)}", file=__import__("sys").stderr)
PY
fi

assert_min_lines koges.common.ids "${MIN_COMMON_SNPS}" "Allele-matched KoGES-1KG SNPs"


# ============================================================================
# STEP 7. Extract first, then rename to chr:pos, then harmonize REF alleles
# ============================================================================

echo "[STEP 7] Extract common markers and harmonize KoGES REF" | tee -a "${LOG}"

# --set-all-var-ids is applied before --extract by PLINK2.  Therefore extract
# chr:pos:REF:ALT IDs first, and only then rename them to chr:pos.
if ! reuse_shared_pfile koges.common.extract; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.kchip.ra \
    --extract koges.common.ids \
    --make-pgen \
    --out koges.common.extract
fi

if ! reuse_shared_pfile koges.common.pre; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.common.extract \
    --set-all-var-ids '@:#' \
    --rm-dup exclude-all \
    --make-pgen \
    --out koges.common.pre
fi

if ! reuse_shared_pfile g1k.common.extract; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.kchip.ra \
    --extract g1k.common.ids \
    --make-pgen \
    --out g1k.common.extract
fi

if ! reuse_shared_pfile g1k.common; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.common.extract \
    --set-all-var-ids '@:#' \
    --rm-dup exclude-all \
    --make-pgen \
    --out g1k.common
fi

need_file koges.common.pre.pvar
need_file g1k.common.pvar

if ! reuse_shared g1k.ref.alleles; then
  clear_local_target g1k.ref.alleles
  awk '$0 !~ /^#/ {print $3 "\t" $4}' g1k.common.pvar > g1k.ref.alleles
fi

if ! reuse_shared_pfile koges.common; then
  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.common.pre \
    --ref-allele force g1k.ref.alleles 2 1 \
    --make-pgen \
    --out koges.common
fi

python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def marker_map(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            out[x[2]] = (x[0].replace("chr", ""), x[1], x[3].upper(), x[4].upper())
    return out


koges = marker_map("koges.common.pvar")
g1k = marker_map("g1k.common.pvar")
if set(koges) != set(g1k):
    raise SystemExit(f"[ERROR] Common marker ID sets differ: KoGES={len(koges)}, 1KG={len(g1k)}")
mismatch = [vid for vid in koges if koges[vid] != g1k[vid]]
if mismatch:
    raise SystemExit(f"[ERROR] {len(mismatch)} CHR/POS/REF/ALT mismatches remain after harmonization; first={mismatch[:5]}")
print(f"[INFO] REF/ALT harmonization verified for {len(koges)} shared markers", file=__import__("sys").stderr)
PY


# ============================================================================
# STEP 8. Reference-defined marker QC
# (Projection-specific from here: do NOT reuse joint PCA SNP lists)
# ============================================================================

echo "[STEP 8] Reference-defined MAF/missingness QC" | tee -a "${LOG}"

# 1KG determines MAF and reference marker eligibility.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.common \
  --maf "${MAF}" \
  --geno "${GENO}" \
  --write-snplist \
  --out g1k.ref.maf05.geno02

# KoGES contributes only a technical-availability check.  Do not impose a
# KoGES MAF threshold, since KoGES must not define the reference PCA markers.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.common \
  --geno "${GENO}" dosage \
  --write-snplist \
  --out koges.geno02

sort -u g1k.ref.maf05.geno02.snplist > g1k.ref.maf05.geno02.sorted
sort -u koges.geno02.snplist > koges.geno02.sorted
comm -12 g1k.ref.maf05.geno02.sorted koges.geno02.sorted > common.reference_qc.snps

assert_min_lines common.reference_qc.snps "${MIN_COMMON_SNPS}" "Reference-QC shared SNPs"


# ============================================================================
# STEP 9. LD pruning only in 1KG ALL
# ============================================================================

echo "[STEP 9] LD prune reference SNPs in 1KG ALL" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.common \
  --extract common.reference_qc.snps \
  --indep-pairwise 50 5 0.2 \
  --out g1k.reference.prune

need_file g1k.reference.prune.prune.in
assert_min_lines g1k.reference.prune.prune.in "${MIN_PRUNED_SNPS}" "LD-pruned reference SNPs"


# ============================================================================
# STEP 10. Apply identical final markers to 1KG and KoGES
# ============================================================================

echo "[STEP 10] Create final reference and target PGENs" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.common \
  --extract g1k.reference.prune.prune.in \
  --make-pgen \
  --out g1k.reference.final

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.common \
  --extract g1k.reference.prune.prune.in \
  --make-pgen \
  --out koges.projection.final

need_file g1k.reference.final.pgen
need_file g1k.reference.final.pvar
need_file g1k.reference.final.psam
need_file koges.projection.final.pgen
need_file koges.projection.final.pvar
need_file koges.projection.final.psam

python3 - <<'PY' 2> >(tee -a "${LOG}" >&2)
def marker_set(path):
    rows = []
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            rows.append((x[2], x[0].replace("chr", ""), x[1], x[3].upper(), x[4].upper()))
    if len(rows) != len(set(rows)):
        raise RuntimeError(f"Duplicate final markers found in {path}")
    return set(rows)


ref = marker_set("g1k.reference.final.pvar")
target = marker_set("koges.projection.final.pvar")
if ref != target:
    only_ref = sorted(ref - target)[:5]
    only_target = sorted(target - ref)[:5]
    raise SystemExit(
        f"[ERROR] Final marker sets differ: reference_only={only_ref}, target_only={only_target}"
    )
print(f"[INFO] Final marker set identity verified: n={len(ref)}", file=__import__("sys").stderr)
PY


# ============================================================================
# STEP 10b. Allele-frequency flip diagnostic against 1KG EAS
# ============================================================================

echo "[STEP 10b] KoGES-vs-1KG EAS allele-frequency check" | tee -a "${LOG}"

python3 - "${G1K_PANEL}" "${EXPECTED_1KG_EAS_N}" > g1k.eas.iids <<'PY'
import csv
import sys

path, expected = sys.argv[1], int(sys.argv[2])
with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    if reader.fieldnames is None:
        raise RuntimeError("1KG panel header not detected")
    lower = {x.lower().strip(): x for x in reader.fieldnames}
    if "sample" not in lower or "super_pop" not in lower:
        raise RuntimeError(f"Panel requires sample/super_pop columns; found {reader.fieldnames}")
    ids = [
        row[lower["sample"]].strip()
        for row in reader
        if row[lower["super_pop"]].strip() == "EAS" and row[lower["sample"]].strip()
    ]

if len(ids) != expected:
    raise RuntimeError(f"Expected {expected} EAS panel samples, found {len(ids)}")
for iid in ids:
    print(iid)
print(f"[INFO] 1KG EAS panel samples={len(ids)}", file=sys.stderr)
PY

python3 - <<'PY'
wanted = set(x.strip() for x in open("g1k.eas.iids") if x.strip())

with open("g1k.reference.final.psam") as handle:
    for line in handle:
        if line.startswith("#"):
            header = line.rstrip().lstrip("#").split()
            break
    else:
        raise RuntimeError("PSAM header not found")

cols = {name: i for i, name in enumerate(header)}
if "IID" not in cols:
    raise RuntimeError(f"IID column not found: {header}")

found = 0
with open("g1k.reference.final.psam") as source, open("g1k.eas.keep", "w") as out:
    for line in source:
        if line.startswith("#"):
            continue
        x = line.rstrip().split()
        if not x or x[cols["IID"]] not in wanted:
            continue
        if "FID" in cols:
            out.write(f"{x[cols['FID']]}\t{x[cols['IID']]}\n")
        else:
            out.write(f"{x[cols['IID']]}\n")
        found += 1

if found != len(wanted):
    raise RuntimeError(f"Only {found}/{len(wanted)} EAS IDs matched the 1KG PSAM")
print(f"[INFO] EAS samples matched to 1KG PGEN={found}")
PY

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.projection.final \
  --freq \
  --out koges.projection.af

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.reference.final \
  --keep g1k.eas.keep \
  --freq \
  --out g1k.eas.reference.af

python3 - "${AF_DIFF_WARN}" "${AF_FLIP_MATCH}" "${MAX_SUSPECT_FLIPS_ABORT}" \
  2> >(tee -a "${LOG}" >&2) <<'PY'
import sys

diff_warn, flip_match, max_abort = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])


def read_af(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) < 5:
                continue
            out[x[1]] = (x[2].upper(), x[3].upper(), float(x[4].split(",")[0]))
    return out


koges = read_af("koges.projection.af.afreq")
eas = read_af("g1k.eas.reference.af.afreq")
suspect = []
allele_errors = []

for vid in sorted(set(koges) & set(eas)):
    kr, ka, kaf = koges[vid]
    er, ea, eaf = eas[vid]
    if (kr, ka) != (er, ea):
        allele_errors.append((vid, kr, ka, er, ea))
        continue
    d = abs(kaf - eaf)
    d_flip = abs(kaf - (1.0 - eaf))
    if d >= diff_warn and d_flip <= flip_match:
        suspect.append((vid, kaf, eaf, d, d_flip))

with open("af_flip_check.suspect.tsv", "w") as out:
    out.write("ID\tKoGES_ALT_AF\tEAS_ALT_AF\tABS_DIFF\tABS_DIFF_TO_1_MINUS_EAS\n")
    for row in suspect:
        out.write("\t".join(map(str, row)) + "\n")

with open("af_flip_check.exclude.ids", "w") as out:
    for vid, *_ in suspect:
        out.write(vid + "\n")

if allele_errors:
    raise SystemExit(f"[ERROR] REF/ALT mismatch in AF files; first={allele_errors[:5]}")

if len(suspect) > max_abort:
    raise SystemExit(
        f"[ERROR] {len(suspect)} suspect allele flips "
        f"(> abort threshold {max_abort}); see af_flip_check.suspect.tsv"
    )

if suspect:
    for vid, kaf, eaf, d, d_flip in suspect[:20]:
        print(
            f"[WARN] Suspect allele flip (will exclude): {vid} "
            f"KoGES_AF={kaf:.4f} EAS_AF={eaf:.4f} "
            f"|diff|={d:.4f} |diff_to_1-AF|={d_flip:.4f}",
            file=sys.stderr,
        )
    print(
        f"[INFO] Excluding {len(suspect)} suspect allele-flip SNP(s). "
        f"See af_flip_check.exclude.ids",
        file=sys.stderr,
    )
else:
    print(
        f"[INFO] AF flip check passed: suspect_flips=0",
        file=sys.stderr,
    )
PY

# Exclude rare residual allele-flip SNPs from both final datasets
if [[ -s af_flip_check.exclude.ids ]]; then
  NEX=$(wc -l < af_flip_check.exclude.ids | tr -d ' ')
  echo "[INFO] Removing ${NEX} suspect-flip SNP(s) from reference/target finals" \
    | tee -a "${LOG}"

  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile g1k.reference.final \
    --exclude af_flip_check.exclude.ids \
    --make-pgen \
    --out g1k.reference.final.tmp

  run "${PLINK2}" \
    --threads "${THREADS}" \
    --pfile koges.projection.final \
    --exclude af_flip_check.exclude.ids \
    --make-pgen \
    --out koges.projection.final.tmp

  for ext in pgen pvar psam; do
    mv -f "g1k.reference.final.tmp.${ext}" "g1k.reference.final.${ext}"
    mv -f "koges.projection.final.tmp.${ext}" "koges.projection.final.${ext}"
  done
  rm -f g1k.reference.final.tmp.log koges.projection.final.tmp.log

  echo "[INFO] After flip exclusion: 1KG=$(count_pvar g1k.reference.final.pvar) KoGES=$(count_pvar koges.projection.final.pvar)" \
    | tee -a "${LOG}"
fi


# ============================================================================
# STEP 11. 1KG ALL reference allele frequencies and PCA
# ============================================================================

echo "[STEP 11] Compute 1KG ALL reference PCA" | tee -a "${LOG}"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.reference.final \
  --freq counts \
  --out g1k.reference

need_file g1k.reference.acount

# Exact PCA is practical for 2,504 1KG samples.  meanimpute makes the
# reference PCA standardization consistent with the projection score command.
run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.reference.final \
  --read-freq g1k.reference.acount \
  --pca allele-wts "${PCA_N}" meanimpute vcols=chrom,ref,alt \
  --out g1k.reference.pca

need_file g1k.reference.pca.eigenvec
need_file g1k.reference.pca.eigenval
need_file g1k.reference.pca.eigenvec.allele


# ============================================================================
# STEP 12. Reference self-projection and KoGES projection
# ============================================================================

echo "[STEP 12] Project 1KG (calibration) and KoGES onto reference PCs" | tee -a "${LOG}"

SCORE_COLS="6-$((PCA_N + 5))"

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile g1k.reference.final \
  --read-freq g1k.reference.acount \
  --score g1k.reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums "${SCORE_COLS}" \
  --out g1k.reference.selfproject

run "${PLINK2}" \
  --threads "${THREADS}" \
  --pfile koges.projection.final \
  --read-freq g1k.reference.acount \
  --score g1k.reference.pca.eigenvec.allele 2 5 header-read no-mean-imputation variance-standardize \
  --score-col-nums "${SCORE_COLS}" \
  --out koges.projected

need_file g1k.reference.selfproject.sscore
need_file koges.projected.sscore


# ============================================================================
# STEP 13. Scale-calibrate projected scores and write plot-ready coordinates
# ============================================================================

echo "[STEP 13] Calibrate projected PCs and write coordinate table" | tee -a "${LOG}"

if skip_if_file pca_1kg_reference.coordinates.tsv; then
  echo "[INFO] Coordinates already present; skip calibration" | tee -a "${LOG}"
else

python3 - "${PCA_N}" "${MIN_SELFPROJ_COR}" "${G1K_PANEL}" <<'PY' 2> >(tee -a "${LOG}" >&2)
import csv
import math
import sys

n_pc = int(sys.argv[1])
min_corr = float(sys.argv[2])
panel_path = sys.argv[3]


def read_ws_table(path):
    header = None
    rows = []
    with open(path) as handle:
        for line in handle:
            if not line.strip():
                continue
            if header is None:
                header = line.rstrip().lstrip("#").split()
                continue
            if line.startswith("#"):
                continue
            x = line.rstrip().split()
            if len(x) != len(header):
                raise RuntimeError(f"Column count mismatch in {path}: {x[:5]}")
            rows.append(x)
    if header is None:
        raise RuntimeError(f"No header found in {path}")
    return header, rows


def column_index(header, name):
    try:
        return header.index(name)
    except ValueError as exc:
        raise RuntimeError(f"Column {name!r} absent. Header={header}") from exc


def row_key(header, row):
    iid = row[column_index(header, "IID")]
    fid = row[column_index(header, "FID")] if "FID" in header else "0"
    return (fid, iid)


def pc_columns(header, count, score=False):
    result = []
    for k in range(1, count + 1):
        choices = [f"PC{k}"]
        if score:
            choices = [f"PC{k}_SUM", f"PC{k}_AVG", f"PC{k}"]
        for choice in choices:
            if choice in header:
                result.append(column_index(header, choice))
                break
        else:
            raise RuntimeError(f"Cannot find PC{k} column in header={header}")
    return result


eig_header, eig_rows = read_ws_table("g1k.reference.pca.eigenvec")
self_header, self_rows = read_ws_table("g1k.reference.selfproject.sscore")
koges_header, koges_rows = read_ws_table("koges.projected.sscore")

eig_pc_idx = pc_columns(eig_header, n_pc, score=False)
self_pc_idx = pc_columns(self_header, n_pc, score=True)
koges_pc_idx = pc_columns(koges_header, n_pc, score=True)

eig = {row_key(eig_header, row): [float(row[i]) for i in eig_pc_idx] for row in eig_rows}
self_score = {row_key(self_header, row): [float(row[i]) for i in self_pc_idx] for row in self_rows}
koges_score = {row_key(koges_header, row): [float(row[i]) for i in koges_pc_idx] for row in koges_rows}

shared = sorted(set(eig) & set(self_score))
if len(shared) != len(eig):
    raise RuntimeError(f"Reference self-projection matched {len(shared)}/{len(eig)} samples")

scales = []
diagnostics = []
for pc in range(n_pc):
    x = [self_score[key][pc] for key in shared]
    y = [eig[key][pc] for key in shared]
    sxx = sum(value * value for value in x)
    syy = sum(value * value for value in y)
    sxy = sum(a * b for a, b in zip(x, y))
    if sxx <= 0 or syy <= 0:
        raise RuntimeError(f"Degenerate self-projection for PC{pc + 1}")
    scale = sxy / sxx
    corr = sxy / math.sqrt(sxx * syy)
    if abs(corr) < min_corr:
        raise RuntimeError(
            f"Self-projection calibration failed for PC{pc + 1}: r={corr:.8f} < {min_corr}"
        )
    scales.append(scale)
    diagnostics.append((pc + 1, len(shared), scale, corr))


def read_panel(path):
    out = {}
    with open(path, newline="", encoding="utf-8-sig", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("1KG panel header not detected")
        lower = {x.lower().strip(): x for x in reader.fieldnames}
        for required in ("sample", "pop", "super_pop"):
            if required not in lower:
                raise RuntimeError(f"Panel lacks {required!r}: {reader.fieldnames}")
        for row in reader:
            iid = row[lower["sample"]].strip()
            if iid:
                out[iid] = (row[lower["pop"]].strip(), row[lower["super_pop"]].strip())
    return out


panel = read_panel(panel_path)

with open("pca_1kg_reference.scaling.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["PC", "N_REFERENCE", "SCALE_TO_DIRECT_1KG_PC", "SELF_PROJECTION_CORRELATION"])
    for row in diagnostics:
        writer.writerow(row)

with open("pca_1kg_reference.coordinates.tsv", "w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["FID", "IID", "SOURCE", "POP", "SUPERPOP"] + [f"PC{k}" for k in range(1, n_pc + 1)])

    for fid, iid in sorted(eig):
        pop, superpop = panel.get(iid, ("NA", "NA"))
        writer.writerow([fid, iid, "1KG", pop, superpop] + eig[(fid, iid)])

    for fid, iid in sorted(koges_score):
        values = [koges_score[(fid, iid)][pc] * scales[pc] for pc in range(n_pc)]
        writer.writerow([fid, iid, "KoGES", "KoGES", "EAS"] + values)

print(
    f"[INFO] Wrote coordinates: 1KG={len(eig)}, KoGES={len(koges_score)}, PCs={n_pc}",
    file=sys.stderr,
)
PY
fi

need_file pca_1kg_reference.coordinates.tsv
need_file pca_1kg_reference.scaling.tsv


# ============================================================================
# FINAL SUMMARY
# ============================================================================

echo "============================================================" | tee -a "${LOG}"
echo "[FINAL SUMMARY]" | tee -a "${LOG}"
if (( USE_KCHIP == 1 )); then
  echo "[RESULT] K-chip v1.0 ∩ v1.1 positions=$(wc -l < kchip.hg19.intersection.positions)" | tee -a "${LOG}"
else
  echo "[RESULT] K-chip filter skipped (mode=all)" | tee -a "${LOG}"
fi
echo "[RESULT] Final LD-pruned reference SNPs=$(count_pvar g1k.reference.final.pvar)" | tee -a "${LOG}"
echo "[RESULT] 1KG reference samples=$(count_psam g1k.reference.final.psam)" | tee -a "${LOG}"
echo "[RESULT] KoGES projected samples=$(count_psam koges.projection.final.psam)" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/pca_1kg_reference.coordinates.tsv" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/pca_1kg_reference.scaling.tsv" | tee -a "${LOG}"
echo "[OUTPUT] ${WD}/${LOG}" | tee -a "${LOG}"
echo "[PLOT] Rscript ${PLOT_RSCRIPT} all" | tee -a "${LOG}"
Rscript "${PLOT_RSCRIPT}" all \
  "${WD}/pca_1kg_reference.coordinates.tsv"
echo "[ALL DONE]" | tee -a "${LOG}"

__ALL_BODY__
