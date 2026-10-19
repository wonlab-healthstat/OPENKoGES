#!/bin/bash
# ==============================================================================
# EAS-only full-run driver: Step3a -> Step3b prepare -> reduced LD reference ->
# Step3b GCTA -> Step3c.
#
# Consumes EAS-only Step2 (DIRECT_KOGES_RESCUE_IN_CORE, rescue-aware grouping,
# EAS_primary_ready without a required 1000G clump).  EUR_supporting is not
# prepared or run.  Do not point these directories at the joint V4 run.
#
# Usage:
#   bash EAS-only/run_step3_full_V4_EAS_only.sh
#   STOP_AFTER=3a bash EAS-only/run_step3_full_V4_EAS_only.sh
#   START_FROM=ref bash EAS-only/run_step3_full_V4_EAS_only.sh
#   Stages: 3a | 3b_prepare | ref | 3b | 3c
# ==============================================================================
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
SHARED="$(cd "$HERE/../shared" && pwd)"
cd "$HERE"

BASE=/data4/sjim/openKoGES_revision/DATA_cumulative_inci/forBOLT-LMM/BOLT_lead_variant
S2="$BASE/COJO/EAS-only/ancestry_clumped_1000G_V4"
S1="$BASE/COJO/EAS-only/step1"
MA="$BASE/COJO/EAS-only/GCTA_ma_V4"
G="$BASE/COJO/EAS-only/GCTA_conditional_V4"
LOG="$BASE/COJO/EAS-only/logs"
REF="$BASE/COJO/EAS-only/KoGES_COJO_LD_reference_V4"
FULL_BFILE=/data6/OpenGWAS/KBA160K/Genotype/KBAv1_KoGES_QCed_COMMON_Imputed_KRPphase2_RSupdated

N_PARALLEL="${N_PARALLEL:-32}"
STOP_AFTER="${STOP_AFTER:-}"
START_FROM="${START_FROM:-3a}"
SETS="EAS_primary"
STAGES="3a 3b_prepare ref 3b 3c"

mkdir -p "$MA" "$G" "$LOG"

say() { echo; echo "=== [$(date '+%H:%M:%S')] $* ==="; }
stop_here() { [ "$STOP_AFTER" = "$1" ] && { echo; echo "STOP_AFTER=$1 reached."; exit 0; } || true; }

stage_index() {
  local i=0 s
  for s in $STAGES; do
    i=$((i + 1))
    [ "$s" = "$1" ] && { echo "$i"; return 0; }
  done
  echo "ERROR: unknown stage '$1'. Valid stages: $STAGES" >&2
  return 1
}

START_IDX="$(stage_index "$START_FROM")"
if [ -n "$STOP_AFTER" ]; then stage_index "$STOP_AFTER" > /dev/null; fi

skip() {
  local idx
  idx="$(stage_index "$1")"
  [ "$idx" -lt "$START_IDX" ]
}

require_prepared_jobs() {
  local s manifest missing
  for s in $SETS; do
    manifest="$G/joint_COJO_gcta_job_manifest.${s}.tsv"
    [ -f "$manifest" ] || {
      echo "ERROR: START_FROM=$START_FROM skips 3b_prepare but $manifest is missing." >&2
      exit 1
    }
    missing="$(
      awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="cojo_extract_snplist_file")e=i;
                    if($i=="cojo_condition_snplist_file")c=i}; next}
                  {print $e; print $c}' "$manifest" |
      while read -r p; do [ -f "$p" ] || echo x; done | wc -l
    )"
    printf "  %-16s %5d job(s), missing snplist file(s): %s\n" \
      "$s" "$(( $(wc -l < "$manifest") - 1 ))" "$missing"
    [ "$missing" -eq 0 ] || {
      echo "ERROR: $manifest references $missing snplist file(s) that do not exist." >&2
      exit 1
    }
  done
}

# ------------------------------------------------------------------ Step2 gate
MANIFEST="$S2/joint_COJO_condition_manifest_ancestry_clumped.tsv"
GATE="$MANIFEST $S2/joint_COJO_EAS_primary_targets.tsv"

say "Step2 EAS-only output check"
waited=0
while :; do
  missing=""
  for f in $GATE; do [ -f "$f" ] || missing="$missing $f"; done
  [ -z "$missing" ] && break
  if [ "${WAIT_FOR_STEP2:-0}" = "0" ]; then
    echo "ERROR: EAS-only Step2 has not finished. Missing:" >&2
    for f in $missing; do echo "  $f" >&2; done
    echo "Re-run with WAIT_FOR_STEP2=1 to wait for it." >&2
    exit 1
  fi
  echo "  [$(date '+%H:%M:%S')] waiting for EAS-only Step2 ($(echo $missing | wc -w) file(s) missing) ..."
  waited=1
  sleep 300
done
[ "$waited" -eq 1 ] && sleep 60
for s in $SETS; do
  printf "  %-24s %s ready target(s)\n" "$s" "$(( $(wc -l < "$S2/joint_COJO_${s}_targets.tsv") - 1 ))"
done
echo "  manifest rows: $(( $(wc -l < "$MANIFEST") - 1 ))"

# --------------------------------------------------------------------- Step 3a
if skip 3a; then
  say "Step3a: skipped (START_FROM=$START_FROM)"
  echo "  existing .ma files: $(ls "$MA"/*.ma 2>/dev/null | wc -l)"
else
  say "Step3a: BOLT -> GCTA .ma (EAS-primary phenotypes)"
  env STEP2_DIR="$S2" \
      READY_TARGETS_FILE="$S2/joint_COJO_EAS_primary_targets.tsv" \
      MA_OUT_DIR="$MA" \
      OVERWRITE=FALSE \
      Rscript --vanilla step3a.convert_bolt_sumstats_to_gcta_ma_EAS_only_V4.R \
      > "$LOG/step3a_EAS_primary.log" 2>&1
  tail -3 "$LOG/step3a_EAS_primary.log" | sed 's/^/    /'
fi
stop_here 3a

# ------------------------------------------------------- Step 3b (prepare only)
if skip 3b_prepare; then
  say "Step3b prepare: skipped (START_FROM=$START_FROM), verifying existing jobs"
  require_prepared_jobs
else
  say "Step3b prepare (no GCTA): writing extract snplists"
  env READY_TARGETS_FILE="$S2/joint_COJO_EAS_primary_targets.tsv" \
      STEP2_DIR="$S2" \
      MA_DIR="$MA" \
      OUT_DIR="$G" \
      RESULTS_TAG="EAS_primary" \
      KOGES_BFILE="$FULL_BFILE" \
      RUN_MODE=prepare RUN_GCTA=false OVERWRITE=TRUE \
      Rscript --vanilla step3b.run_joint_cojo_conditional_gcta_EAS_only_V4.R \
      > "$LOG/step3b_prepare_EAS_primary.log" 2>&1
  grep -E 'Prepared GCTA jobs|Assessable|Unassessable' "$LOG/step3b_prepare_EAS_primary.log" | sed 's/^/    /'
fi
stop_here 3b_prepare

# ------------------------------------------------------- reduced LD reference
if skip ref; then
  say "Reduced LD reference: skipped (START_FROM=$START_FROM)"
  [ -f "${REF}.bed" ] || {
    echo "ERROR: START_FROM=$START_FROM skips the reference build but ${REF}.bed is missing." >&2
    exit 1
  }
  echo "  reference variants: $(wc -l < "${REF}.bim")"
else
  say "Building reduced LD reference (EAS-primary extract union)"
  KOGES_BFILE_FULL="$FULL_BFILE" \
    bash "$SHARED/build_cojo_ld_reference.sh" "$REF" "$G/work_snplists" 2>&1 | sed 's/^/  /'
fi
stop_here ref

# -------------------------------------------------------------- Step 3b (GCTA)
if skip 3b; then
  say "Step3b GCTA: skipped (START_FROM=$START_FROM)"
else
  say "Step3b GCTA (reduced reference, N_PARALLEL=$N_PARALLEL)"
  env READY_TARGETS_FILE="$S2/joint_COJO_EAS_primary_targets.tsv" \
      STEP2_DIR="$S2" \
      MA_DIR="$MA" \
      OUT_DIR="$G" \
      RESULTS_TAG="EAS_primary" \
      KOGES_BFILE="$REF" \
      KOGES_BFILE_FULL="$FULL_BFILE" \
      RUN_MODE=batch RUN_GCTA=true OVERWRITE=TRUE \
      N_PARALLEL="$N_PARALLEL" \
      Rscript --vanilla step3b.run_joint_cojo_conditional_gcta_EAS_only_V4.R \
      > "$LOG/step3b_run_EAS_primary.log" 2>&1
  grep -E 'variants (in reference|required)|all required|samples identical|position and allele|Parallel GCTA finished|Conditionally independent|Not conditionally GWS|High-LD|Conditionally unassessable' \
    "$LOG/step3b_run_EAS_primary.log" | sed 's/^/    /'
fi
stop_here 3b

# --------------------------------------------------------------------- Step 3c
say "Step3c: EAS-only integration (no EUR_supporting QC)"
env STEP1_MANIFEST="$S1/joint_COJO_condition_manifest_preLD.tsv" \
    STEP2_MANIFEST="$MANIFEST" \
    BATCH_RESULTS="$G/joint_COJO_conditional_results.EAS_primary.tsv" \
    EUR_BATCH_RESULTS="" \
    READY_TARGETS_FILE="$S2/joint_COJO_EAS_primary_targets.tsv" \
    REQUIRE_EUR_SUPPORTING_QC=FALSE \
    OUT_FILE="$G/joint_COJO_conditional_results_integrated.tsv" \
    OUT_SUMMARY="$G/joint_COJO_conditional_status_summary.tsv" \
    OUT_PAPER="$G/joint_COJO_conditional_results_paper.tsv" \
    Rscript --vanilla step3c.merge_joint_cojo_conditional_results_EAS_only_V4.R \
    > "$LOG/step3c.log" 2>&1
tail -25 "$LOG/step3c.log" | sed 's/^/  /'

say "DONE"
echo "Integrated : $G/joint_COJO_conditional_results_integrated.tsv"
echo "Paper table: $G/joint_COJO_conditional_results_paper.tsv"
echo "Logs       : $LOG"
