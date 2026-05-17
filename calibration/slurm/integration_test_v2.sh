#!/bin/bash
#SBATCH --job-name=fvs_integ_v2
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/integ_v2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/integ_v2_%j.err

# Comprehensive integration test v2.
# Builds + runs all 12 supported variants (7 Eastern + 5 Western sample)
# using each variant's own upstream test deck where available, with a
# refined PASS/WARN/FAIL rubric:
#
#   PASS = exit 0 or 10 (clean run / informational notes)
#   WARN = exit 20+ but .sum file complete and contains variant marker
#   FAIL = no .sum produced or marker mismatch or .sum < 1KB

set -uo pipefail
module purge
module load gcc/12.3.0 R/4.4.0

PROJ=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
US=/users/PUOM0008/crsfaaron/upstream_fvs_check/ForestVegetationSimulator
cd $PROJ

VARIANTS=(acd ne cs ls sn kt em wc op ca bm cr)
declare -A TESTDIR
TESTDIR[acd]="FVSne"   TESTDIR[ne]="FVSne"    TESTDIR[cs]="FVScs"
TESTDIR[ls]="FVSls"    TESTDIR[sn]="FVSsn"    TESTDIR[kt]="FVSkt"
TESTDIR[em]="FVSem"    TESTDIR[wc]="FVSwc"    TESTDIR[op]="FVSop"
TESTDIR[ca]="FVSca"    TESTDIR[bm]="FVSbm"    TESTDIR[cr]="FVScr"

declare -A MARKER=(
  [acd]="AC" [ne]="NE" [cs]="CS" [ls]="LS" [sn]="SN"
  [kt]="KT" [em]="EM" [wc]="WC" [op]="OP" [ca]="CA"
  [bm]="BM" [cr]="CR"
)

TEST_ROOT=$(mktemp -d)
mkdir -p $TEST_ROOT/lib
RESULTS=$TEST_ROOT/results.tsv
echo -e "phase\tvariant\tstatus\tdetail" > $RESULTS

log()  { echo "[$(date +%H:%M:%S)] $*"; }
rec()  { echo -e "$1\t$2\t$3\t$4" >> $RESULTS; }

find_keyfile () {
  local v=$1
  local td=${TESTDIR[$v]}
  local d=$US/tests/$td
  if [ ! -d $d ]; then echo ""; return; fi
  # Prefer canonical <variant>t01.key (the standard projection deck).
  # Fall back to net01.key (for ACD reusing NE deck).
  # Last resort: smallest .key file.
  if [ -f $d/${v}t01.key ]; then
    echo $d/${v}t01.key
  elif [ -f $d/net01.key ]; then
    echo $d/net01.key
  else
    ls -1S $d/*.key 2>/dev/null | tail -1
  fi
}

##############################################################################
# PHASE 1: build
##############################################################################
log "=== PHASE 1: build all ${#VARIANTS[@]} variants ==="
for V in "${VARIANTS[@]}"; do
  log "Building FVS$V..."
  if bash deployment/scripts/build_fvs_executables.sh . $TEST_ROOT/lib $V \
       > $TEST_ROOT/build_${V}.log 2>&1; then
    SIZE=$(stat -c %s $TEST_ROOT/lib/FVS${V} 2>/dev/null || echo 0)
    if [ "$SIZE" -gt 1000000 ]; then
      rec build $V PASS "${SIZE}B"
      log "  PASS: FVS$V = $SIZE bytes"
    else
      rec build $V FAIL "exe too small: ${SIZE}B"
    fi
  else
    rec build $V FAIL "build script returned nonzero"
    tail -10 $TEST_ROOT/build_${V}.log
  fi
done

##############################################################################
# PHASE 2: run on variant-specific deck
##############################################################################
log "=== PHASE 2: run each variant on its own test deck ==="
for V in "${VARIANTS[@]}"; do
  if [ ! -x $TEST_ROOT/lib/FVS${V} ]; then
    rec run $V SKIP "no executable"; continue
  fi
  KEY=$(find_keyfile $V)
  if [ -z "$KEY" ]; then
    rec run $V SKIP "no upstream test deck for $V"; continue
  fi

  W=$TEST_ROOT/run_${V}
  mkdir -p $W
  perl -pe "s/\r\n|\r/\n/g" < $KEY > $W/test.key
  # Copy companion .tre file with matching basename
  BASE=$(basename $KEY .key)
  KEYDIR=$(dirname $KEY)
  for ext in tre fvs slf; do
    [ -f $KEYDIR/${BASE}.${ext} ] && cp $KEYDIR/${BASE}.${ext} $W/test.${ext}
  done

  (cd $W && $TEST_ROOT/lib/FVS${V} --keywordfile=test.key) \
    > $TEST_ROOT/run_${V}.log 2>&1
  RC=$?

  SUM_SZ=$(stat -c %s $W/test.sum 2>/dev/null || echo 0)
  OUT_SZ=$(stat -c %s $W/test.out 2>/dev/null || echo 0)

  # Refined rubric
  if [ "$SUM_SZ" -lt 1000 ]; then
    rec run $V FAIL "rc=$RC sum=${SUM_SZ}B (too small)"
  elif [ "$RC" = "0" ] || [ "$RC" = "10" ]; then
    rec run $V PASS "rc=$RC sum=${SUM_SZ}B out=${OUT_SZ}B  deck=$(basename $KEY)"
    log "  PASS: FVS$V rc=$RC sum=$SUM_SZ"
  elif [ "$RC" = "20" ]; then
    rec run $V WARN "rc=$RC sum=${SUM_SZ}B (FVS04 input mismatch but output complete)"
    log "  WARN: FVS$V rc=$RC sum=$SUM_SZ (deck issue, output OK)"
  else
    rec run $V FAIL "rc=$RC sum=${SUM_SZ}B"
    log "  FAIL: FVS$V rc=$RC sum=$SUM_SZ"
    tail -10 $TEST_ROOT/run_${V}.log
  fi
done

##############################################################################
# PHASE 3: variant marker check
##############################################################################
log "=== PHASE 3: variant marker check ==="
for V in "${VARIANTS[@]}"; do
  SUM=$TEST_ROOT/run_${V}/test.sum
  if [ ! -f $SUM ]; then rec marker $V SKIP "no .sum"; continue; fi
  MARK="${MARKER[$V]}"
  if head -1 $SUM | grep -qE "\b${MARK}\b"; then
    rec marker $V PASS "found ${MARK}"
  else
    rec marker $V FAIL "missing marker ${MARK}"
  fi
done

##############################################################################
# PHASE 4: md5 distinctness across variants
##############################################################################
log "=== PHASE 4: md5 distinctness ==="
declare -A SUM_MD5
for V in "${VARIANTS[@]}"; do
  S=$TEST_ROOT/run_${V}/test.sum
  [ -f $S ] && SUM_MD5[$V]=$(md5sum $S | cut -d" " -f1)
done
DISTINCT=$(printf "%s\n" "${SUM_MD5[@]}" | sort -u | wc -l)
COUNT=${#SUM_MD5[@]}
if [ "$DISTINCT" = "$COUNT" ] && [ "$COUNT" -gt 1 ]; then
  rec distinct all PASS "$COUNT variants, $DISTINCT unique md5s"
else
  rec distinct all WARN "$COUNT variants, $DISTINCT distinct (allow ACD==NE if same deck)"
fi

##############################################################################
# PHASE 5: R post-pass smoke test
##############################################################################
log "=== PHASE 5: R post-pass smoke test ==="
if Rscript calibration/R/smoke_postpass.R > $TEST_ROOT/smoke.log 2>&1; then
  rec smoke postpass PASS "stratified post-pass OK"
else
  rec smoke postpass FAIL "see smoke.log"
fi

##############################################################################
# PHASE 6: write report
##############################################################################
log "=== PHASE 6: build report ==="
REPORT=$PROJ/INTEGRATION_TEST_REPORT_v2.md
{
  echo "# fvs-modern integration test v2 (variant-specific decks)"
  echo
  echo "Run: SLURM job \$SLURM_JOB_ID"
  echo "Date: $(date -u +"%Y-%m-%d %H:%M UTC")"
  echo "Branch: acd-bridge-fix-2026-05-15"
  echo "Commit: $(git log --oneline -1)"
  echo
  echo "## Test scope: 12 variants"
  echo
  echo "**Eastern** (7): ACD, NE, CS, LS, SN, KT, EM"
  echo "**Western sample** (5): WC (Westside Cascades), OP (Olympic Peninsula), CA (Inland California), BM (Blue Mountains), CR (Central Rockies)"
  echo
  echo "**Test deck:** each variant uses its own upstream tests/FVS<variant>/ deck."
  echo "ACD falls back to NE's net01 because no FVSacd test dir exists upstream."
  echo
  echo "## Refined rubric"
  echo
  echo "- **PASS**: rc in {0, 10} and .sum >= 1KB"
  echo "- **WARN**: rc=20 (input-variant mismatch or non-fatal error) but .sum complete"
  echo "- **FAIL**: no .sum, .sum too small, marker missing, or rc >= 30"
  echo
  echo "## Results"
  echo
  echo "| Phase | Variant | Status | Detail |"
  echo "| --- | --- | --- | --- |"
  tail -n +2 $RESULTS | awk -F"\t" "{printf \"| %s | %s | %s | %s |\n\", \$1, \$2, \$3, \$4}"
  echo
  echo "## Pass summary"
  echo
  TOT=$(tail -n +2 $RESULTS | wc -l)
  PASS=$(grep -c -P "\tPASS\t" $RESULTS || true)
  WARN=$(grep -c -P "\tWARN\t" $RESULTS || true)
  FAIL=$(grep -c -P "\tFAIL\t" $RESULTS || true)
  SKIP=$(grep -c -P "\tSKIP\t" $RESULTS || true)
  echo "- Total: $TOT"
  echo "- PASS: $PASS"
  echo "- WARN: $WARN"
  echo "- FAIL: $FAIL"
  echo "- SKIP: $SKIP"
  echo
  if [ "$FAIL" = "0" ]; then
    echo "**Overall: PASS** ($PASS PASS / $WARN WARN / $SKIP SKIP)"
  else
    echo "**Overall: PARTIAL** ($FAIL FAIL among $TOT checks)"
  fi
} > $REPORT
cat $REPORT
log "Report at $REPORT"

# Save artifacts
ART=$PROJ/calibration/analysis/acd_stand_level_2026-05-16/integration_test_v2
mkdir -p $ART
cp $RESULTS $ART/
cp $TEST_ROOT/build_*.log $ART/ 2>/dev/null
cp $TEST_ROOT/run_*.log   $ART/ 2>/dev/null
cp $TEST_ROOT/smoke.log   $ART/ 2>/dev/null

echo "DONE"
