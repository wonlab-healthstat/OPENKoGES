#!/usr/bin/env bash
# ==============================================================================
# EAS COJO dual-mode driver (Step2 -> Step3a -> Step3b prepare -> reduced LD
# reference -> Step3b GCTA -> Step3c)
#
# Modes:
#   EAS_RUN_MODE=primary
#     - 1000G EAS source-specific clump-defined condition set only
#     - direct KoGES-exact rescue is audit-only
#     - rescue-only loci are NOT primary-ready
#
#   EAS_RUN_MODE=sensitivity
#     - same clump-defined set + direct KoGES-exact rescue variants
#     - rescue-only loci may be assessed as sensitivity evidence
#
# Examples:
#   EAS_RUN_MODE=primary bash run_EAS_COJO_dualmode_V5.sh
#   EAS_RUN_MODE=sensitivity bash run_EAS_COJO_dualmode_V5.sh
#   bash run_EAS_COJO_dualmode_V5.sh primary
#
# Optional resume controls:
#   START_FROM=2|3a|3b_prepare|ref|3b|3c
#   STOP_AFTER=2|3a|3b_prepare|ref|3b|3c
#
# After both modes have completed, compare them reproducibly with:
#   Rscript --vanilla step3d.compare_EAS_dualmode_V5.R
# ==============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SHARED="$(cd "$HERE/../shared" && pwd)"
cd "$HERE"

EAS_RUN_MODE="${EAS_RUN_MODE:-${1:-primary}}"
EAS_RUN_MODE="$(printf '%s' "$EAS_RUN_MODE" | tr '[:upper:]' '[:lower:]')"
case "$EAS_RUN_MODE" in
  primary)     SUFFIX="EAS_primary_clumponly" ;;
  sensitivity) SUFFIX="EAS_sensitivity_rescue" ;;
  *) echo "ERROR: EAS_RUN_MODE must be primary or sensitivity." >&2; exit 2 ;;
esac
export EAS_RUN_MODE

BASE="${BASE_DIR:-/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant}"
S1="${PRELD_DIR:-$BASE/COJO/EAS-only/step1}"
S2="${STEP2_DIR:-$BASE/COJO/EAS-only/ancestry_clumped_1000G_V5_${SUFFIX}}"
MA="${MA_DIR:-$BASE/COJO/EAS-only/GCTA_ma_V5_${SUFFIX}}"
G="${GCTA_OUT_DIR:-$BASE/COJO/EAS-only/GCTA_conditional_V5_${SUFFIX}}"
LOG="${LOG_DIR:-$BASE/COJO/EAS-only/logs_V5_${SUFFIX}}"
REF="${COJO_LD_REFERENCE:-$BASE/COJO/EAS-only/KoGES_COJO_LD_reference_V5_${SUFFIX}}"
REF_1000G_DIR="${REF_1000G_DIR:-$BASE/COJO/1000G_LD_panel_indel}"
FULL_BFILE="${KOGES_BFILE_FULL:-/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated}"

# Sized for a dedicated ~40-core / ~500G-free node per mode.
# Two nodes (primary vs sensitivity) should each export these values; they
# do not add together.  Do not run both modes on one 40-core node at 40+40.
# N_PARALLEL: 1000G clump + snplist + Step3b GCTA workers.
# N_REP_PARALLEL: representative-selection forks copy membership tables;
#   keep below 40 so one mclapply wave does not OOM.
# N_MEM_PARALLEL: KoGES 160K --r2.  Default code cap is 16; 500G free allows
#   more, but 40 concurrent full-bfile reads will stall on I/O.
# N_DT_THREADS: data.table OpenMP during Step2 prep only (forks drop to 1).
N_PARALLEL="${N_PARALLEL:-40}"
N_REP_PARALLEL="${N_REP_PARALLEL:-32}"
N_MEM_PARALLEL="${N_MEM_PARALLEL:-32}"
N_DT_THREADS="${N_DT_THREADS:-16}"
STEP2_OVERWRITE="${STEP2_OVERWRITE:-false}"
MA_OVERWRITE="${MA_OVERWRITE:-true}"
GCTA_OVERWRITE="${GCTA_OVERWRITE:-true}"
START_FROM="${START_FROM:-2}"
STOP_AFTER="${STOP_AFTER:-}"
STAGES="2 3a 3b_prepare ref 3b 3c"

STEP2_SCRIPT="step2.clump_joint_cojo_conditions_1000G_EAS_dualmode_V5.R"
STEP3A_SCRIPT="step3a.convert_bolt_sumstats_to_gcta_ma_EAS_dualmode_V5.R"
STEP3B_SCRIPT="step3b.run_joint_cojo_conditional_gcta_EAS_dualmode_V5.R"
STEP3C_SCRIPT="step3c.merge_joint_cojo_conditional_results_EAS_dualmode_V5.R"
STEP3D_SCRIPT="step3d.compare_EAS_dualmode_V5.R"
REF_SCRIPT="$SHARED/build_cojo_ld_reference.sh"

for f in "$STEP2_SCRIPT" "$STEP3A_SCRIPT" "$STEP3B_SCRIPT" "$STEP3C_SCRIPT" "$STEP3D_SCRIPT"; do
  [ -f "$f" ] || { echo "ERROR: missing script: $f" >&2; exit 1; }
done
[ -f "$REF_SCRIPT" ] || { echo "ERROR: missing script: $REF_SCRIPT" >&2; exit 1; }

mkdir -p "$LOG"

say() { echo; echo "=== [$(date '+%F %T')] $* ==="; }

stage_index() {
  local x i=0
  for x in $STAGES; do
    i=$((i+1))
    [ "$x" = "$1" ] && { echo "$i"; return 0; }
  done
  echo "ERROR: unknown stage '$1'; valid: $STAGES" >&2
  return 1
}
START_IDX="$(stage_index "$START_FROM")"
[ -z "$STOP_AFTER" ] || stage_index "$STOP_AFTER" >/dev/null
skip() { [ "$(stage_index "$1")" -lt "$START_IDX" ]; }
stop_here() { [ "$STOP_AFTER" = "$1" ] && { echo "STOP_AFTER=$1 reached."; exit 0; } || true; }

READY="$S2/joint_COJO_ready_targets.tsv"
MANIFEST="$S2/joint_COJO_condition_manifest_ancestry_clumped.tsv"

check_step2_mode() {
  [ -f "$READY" ] || { echo "ERROR: missing $READY" >&2; exit 1; }
  [ -f "$MANIFEST" ] || { echo "ERROR: missing $MANIFEST" >&2; exit 1; }

  local n_ready rescue_sum mode_values
  n_ready=$(( $(wc -l < "$READY") - 1 ))
  rescue_sum="$(awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="n_condition_snps_direct_KoGES_rescue") c=i; next }
    c && $c!="NA" && $c!="" { s += $c }
    END { print s+0 }
  ' "$READY")"
  mode_values="$(awk -F'\t' '
    NR==1 { for(i=1;i<=NF;i++) if($i=="EAS_run_mode") c=i; next }
    c && $c!="NA" && $c!="" { seen[$c]=1 }
    END { for(x in seen) print x }
  ' "$READY" | sort -u | paste -sd, -)"

  echo "  Step2 ready targets: $n_ready"
  echo "  selected direct-rescue covariates: $rescue_sum"
  echo "  Step2 EAS_run_mode: ${mode_values:-<missing>}"

  [ "$mode_values" = "$EAS_RUN_MODE" ] || {
    echo "ERROR: Step2 ready table mode '$mode_values' != requested '$EAS_RUN_MODE'." >&2
    exit 1
  }
  if [ "$EAS_RUN_MODE" = "primary" ] && [ "$rescue_sum" -ne 0 ]; then
    echo "ERROR: primary clump-only run contains selected direct KoGES rescue covariates." >&2
    exit 1
  fi
}

say "Run configuration"
echo "  EAS_RUN_MODE = $EAS_RUN_MODE"
echo "  Step1        = $S1"
echo "  Step2        = $S2"
echo "  MA           = $MA"
echo "  GCTA         = $G"
echo "  LD reference = $REF"
echo "  1000G panel  = $REF_1000G_DIR"
echo "  logs         = $LOG"
echo "  N_PARALLEL / N_REP_PARALLEL / N_MEM_PARALLEL / N_DT_THREADS = $N_PARALLEL / $N_REP_PARALLEL / $N_MEM_PARALLEL / $N_DT_THREADS"

# ------------------------------------------------------------------------- Step2
if skip 2; then
  say "Step2 skipped (START_FROM=$START_FROM); validating existing outputs"
else
  say "Step2: EAS prior-signal definition ($EAS_RUN_MODE)"
  env EAS_RUN_MODE="$EAS_RUN_MODE" \
      PRELD_DIR="$S1" \
      OUT_DIR="$S2" \
      REF_1000G_DIR="$REF_1000G_DIR" \
      N_PARALLEL="$N_PARALLEL" \
      N_REP_PARALLEL="$N_REP_PARALLEL" \
      N_MEM_PARALLEL="$N_MEM_PARALLEL" \
      N_DT_THREADS="$N_DT_THREADS" \
      OVERWRITE="$STEP2_OVERWRITE" \
      Rscript --vanilla "$STEP2_SCRIPT" \
      > "$LOG/step2.log" 2>&1
  tail -30 "$LOG/step2.log" | sed 's/^/  /'
fi
check_step2_mode
stop_here 2

mkdir -p "$MA" "$G"

# ------------------------------------------------------------------------ Step3a
if skip 3a; then
  say "Step3a skipped (START_FROM=$START_FROM)"
else
  say "Step3a: BOLT -> full GCTA .ma"
  env EAS_RUN_MODE="$EAS_RUN_MODE" \
      STEP2_DIR="$S2" \
      READY_TARGETS_FILE="$READY" \
      MA_OUT_DIR="$MA" \
      OVERWRITE="$MA_OVERWRITE" \
      Rscript --vanilla "$STEP3A_SCRIPT" \
      > "$LOG/step3a.log" 2>&1
  tail -15 "$LOG/step3a.log" | sed 's/^/  /'
fi
stop_here 3a

# --------------------------------------------------------------- Step3b prepare
if skip 3b_prepare; then
  say "Step3b prepare skipped (START_FROM=$START_FROM)"
  [ -f "$G/joint_COJO_gcta_job_manifest.EAS_primary.tsv" ] || {
    echo "ERROR: prepared Step3b manifest is missing." >&2; exit 1;
  }
else
  say "Step3b prepare: write extract/condition snplists"
  env EAS_RUN_MODE="$EAS_RUN_MODE" \
      STEP2_DIR="$S2" \
      READY_TARGETS_FILE="$READY" \
      MA_DIR="$MA" \
      OUT_DIR="$G" \
      KOGES_BFILE="$FULL_BFILE" \
      RUN_MODE=prepare RUN_GCTA=false OVERWRITE=true \
      Rscript --vanilla "$STEP3B_SCRIPT" \
      > "$LOG/step3b_prepare.log" 2>&1
  grep -E 'Prepared GCTA jobs|Assessable|Unassessable' "$LOG/step3b_prepare.log" | sed 's/^/  /' || true
fi
stop_here 3b_prepare

# ------------------------------------------------------------ reduced LD ref
if skip ref; then
  say "Reduced LD reference skipped (START_FROM=$START_FROM)"
  [ -f "${REF}.bed" ] || { echo "ERROR: missing ${REF}.bed" >&2; exit 1; }
else
  say "Build reduced KoGES LD reference"
  rm -f "${REF}.bed" "${REF}.bim" "${REF}.fam" "${REF}.log" "${REF}.nosex"
  KOGES_BFILE_FULL="$FULL_BFILE" \
    bash "$REF_SCRIPT" "$REF" "$G/work_snplists" \
    > "$LOG/build_ld_reference.log" 2>&1
  tail -20 "$LOG/build_ld_reference.log" | sed 's/^/  /'
fi
stop_here ref

# ---------------------------------------------------------------- Step3b GCTA
if skip 3b; then
  say "Step3b GCTA skipped (START_FROM=$START_FROM)"
else
  say "Step3b: GCTA --cojo-cond ($EAS_RUN_MODE)"
  env EAS_RUN_MODE="$EAS_RUN_MODE" \
      STEP2_DIR="$S2" \
      READY_TARGETS_FILE="$READY" \
      MA_DIR="$MA" \
      OUT_DIR="$G" \
      KOGES_BFILE="$REF" \
      KOGES_BFILE_FULL="$FULL_BFILE" \
      RUN_MODE=batch RUN_GCTA=true OVERWRITE="$GCTA_OVERWRITE" \
      N_PARALLEL="$N_PARALLEL" \
      Rscript --vanilla "$STEP3B_SCRIPT" \
      > "$LOG/step3b.log" 2>&1
  grep -E 'Parallel GCTA finished|Conditionally independent|Not conditionally GWS|High-LD|Conditionally unassessable' \
    "$LOG/step3b.log" | sed 's/^/  /' || true
fi
stop_here 3b

# ------------------------------------------------------------------------ Step3c
say "Step3c: integrate full denominator and final verdict"
env EAS_RUN_MODE="$EAS_RUN_MODE" \
    STEP1_MANIFEST="$S1/joint_COJO_condition_manifest_preLD.tsv" \
    STEP2_MANIFEST="$MANIFEST" \
    BATCH_RESULTS="$G/joint_COJO_conditional_results.EAS_primary.tsv" \
    READY_TARGETS_FILE="$READY" \
    REQUIRE_READY_TARGETS_QC=true \
    REQUIRE_EUR_SUPPORTING_QC=false \
    OUT_FILE="$G/joint_COJO_conditional_results_integrated.tsv" \
    OUT_SUMMARY="$G/joint_COJO_conditional_status_summary.tsv" \
    OUT_PAPER="$G/joint_COJO_conditional_results_paper.tsv" \
    Rscript --vanilla "$STEP3C_SCRIPT" \
    > "$LOG/step3c.log" 2>&1

tail -30 "$LOG/step3c.log" | sed 's/^/  /'

say "DONE"
echo "Mode       : $EAS_RUN_MODE"
echo "Integrated : $G/joint_COJO_conditional_results_integrated.tsv"
echo "Paper table: $G/joint_COJO_conditional_results_paper.tsv"
echo "Logs       : $LOG"
