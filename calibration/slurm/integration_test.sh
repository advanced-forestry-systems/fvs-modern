#!/bin/bash
#SBATCH --job-name=fvs_integ_test
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/integ_test_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/integ_test_%j.err

# Comprehensive integration test for the fvs-modern fork on the
# acd-bridge-fix-2026-05-15 branch.
#
# Phase 1: Build standalone FVS executables for all 7 Eastern variants
# Phase 2: Run each binary on the upstream net01 test deck
# Phase 3: Verify .sum and .out files are non-trivial and contain
#          expected variant markers
# Phase 4: Compare md5s across variants to confirm distinct outputs
# Phase 5: Run R-side smoke test (smoke_postpass.R)
# Phase 6: Build a markdown report

set -uo pipefail
module purge
module load gcc/12.3.0 R/4.4.0

PROJ=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
cd $PROJ

EASTERN_VARIANTS=(acd ne cs ls sn kt em)
NET01_KEY=/users/PUOM0008/crsfaaron/upstream_fvs_check/ForestVegetationSimulator/tests/FVSne/net01.key
NET01_TRE=/users/PUOM0008/crsfaaron/upstream_fvs_check/ForestVegetationSimulator/tests/FVSne/net01.tre

TEST_ROOT=$(mktemp -d)
RESULTS=$TEST_ROOT/results.tsv
echo -e "phase\tvariant\tstatus\tdetail" > $RESULTS

log()  { echo "[$(date +%H:%M:%S)] $*"; }
rec()  { echo -e "$1\t$2\t$3\t$4" >> $RESULTS; }

##############################################################################
# PHASE 1: build
##############################################################################
log "=== PHASE 1: build all 7 Eastern variants ==="
mkdir -p $TEST_ROOT/lib
for V in "${EASTERN_VARIANTS[@]}"; do
  log "Building FVS$V..."
  if bash deployment/scripts/build_fvs_executables.sh . $TEST_ROOT/lib $V \
       > $TEST_ROOT/build_${V}.log 2>&1; then
    SIZE=$(stat -c %s $TEST_ROOT/lib/FVS${V} 2>/dev/null || echo 0)
    if [ "$SIZE" -gt 1000000 ]; then
      rec build $V PASS "${SIZE}B"
      log "  PASS: FVS$V = $SIZE bytes"
    else
      rec build $V FAIL "exe too small or missing: ${SIZE}B"
      log "  FAIL: FVS$V too small ($SIZE bytes)"
    fi
  else
    rec build $V FAIL "build script returned nonzero"
    log "  FAIL: build script failed for $V"
    tail -10 $TEST_ROOT/build_${V}.log
  fi
done

##############################################################################
# PHASE 2: run each on net01
##############################################################################
log "=== PHASE 2: run each variant on net01 ==="
for V in "${EASTERN_VARIANTS[@]}"; do
  if [ ! -x $TEST_ROOT/lib/FVS${V} ]; then
    rec run $V SKIP "no executable"; continue
  fi
  W=$TEST_ROOT/run_${V}
  mkdir -p $W
  tr "\r" "\n" < $NET01_KEY > $W/net01.key
  cp $NET01_TRE $W/

  (cd $W && $TEST_ROOT/lib/FVS${V} --keywordfile=net01.key) \
    > $TEST_ROOT/run_${V}.log 2>&1
  RC=$?

  SUM_SZ=$(stat -c %s $W/net01.sum 2>/dev/null || echo 0)
  OUT_SZ=$(stat -c %s $W/net01.out 2>/dev/null || echo 0)

  # Pass criteria: .sum > 1KB and .out > 50KB and RC in {0,10} (10 = errgro warnings)
  if { [ "$RC" = "0" ] || [ "$RC" = "10" ]; } && \
     [ "$SUM_SZ" -gt 1000 ] && [ "$OUT_SZ" -gt 50000 ]; then
    rec run $V PASS "rc=$RC sum=${SUM_SZ}B out=${OUT_SZ}B"
    log "  PASS: FVS$V rc=$RC sum=$SUM_SZ out=$OUT_SZ"
  else
    rec run $V FAIL "rc=$RC sum=${SUM_SZ}B out=${OUT_SZ}B"
    log "  FAIL: FVS$V rc=$RC sum=$SUM_SZ out=$OUT_SZ"
    tail -10 $TEST_ROOT/run_${V}.log
  fi
done

##############################################################################
# PHASE 3: verify variant markers
##############################################################################
log "=== PHASE 3: verify each .sum carries variant marker ==="
declare -A MARKER=(
  [acd]="AC" [ne]="NE" [cs]="CS" [ls]="LS"
  [sn]="SN" [kt]="KT" [em]="EM"
)
for V in "${EASTERN_VARIANTS[@]}"; do
  SUM=$TEST_ROOT/run_${V}/net01.sum
  if [ ! -f $SUM ]; then
    rec marker $V SKIP "no .sum"; continue
  fi
  MARK="${MARKER[$V]}"
  if head -1 $SUM | grep -qE "\b${MARK}\b"; then
    rec marker $V PASS "found '${MARK}' in header"
  else
    rec marker $V FAIL "missing '${MARK}' marker"
    head -1 $SUM
  fi
done

##############################################################################
# PHASE 4: md5 diff to confirm variants produce distinct outputs
##############################################################################
log "=== PHASE 4: md5 distinctness across variants ==="
declare -A SUM_MD5
for V in "${EASTERN_VARIANTS[@]}"; do
  S=$TEST_ROOT/run_${V}/net01.sum
  [ -f $S ] && SUM_MD5[$V]=$(md5sum $S | cut -d" " -f1)
done

DISTINCT=$(printf "%s\n" "${SUM_MD5[@]}" | sort -u | wc -l)
COUNT=${#SUM_MD5[@]}
if [ "$DISTINCT" = "$COUNT" ] && [ "$COUNT" -gt 1 ]; then
  rec distinct all PASS "$COUNT variants, $DISTINCT unique .sum md5s"
  log "  PASS: $COUNT variants produced $DISTINCT distinct .sum files"
else
  rec distinct all FAIL "$COUNT variants, only $DISTINCT distinct"
  log "  FAIL: $COUNT variants produced only $DISTINCT distinct"
  for V in "${!SUM_MD5[@]}"; do echo "    $V: ${SUM_MD5[$V]}"; done
fi

##############################################################################
# PHASE 5: R post-pass smoke test
##############################################################################
log "=== PHASE 5: R post-pass smoke test ==="
if Rscript calibration/R/smoke_postpass.R > $TEST_ROOT/smoke.log 2>&1; then
  rec smoke postpass PASS "stratified post-pass logic OK"
  log "  PASS: smoke_postpass.R"
else
  rec smoke postpass FAIL "see smoke.log"
  log "  FAIL: smoke_postpass.R"
  tail -15 $TEST_ROOT/smoke.log
fi

##############################################################################
# PHASE 6: build summary report
##############################################################################
log "=== PHASE 6: build report ==="
REPORT_DST=$PROJ/INTEGRATION_TEST_REPORT.md
{
  echo "# fvs-modern integration test report"
  echo
  echo "Run: SLURM job \$SLURM_JOB_ID on Cardinal"
  echo "Date: $(date -u +"%Y-%m-%d %H:%M UTC")"
  echo "Branch: acd-bridge-fix-2026-05-15"
  echo "Commit: $(git log --oneline -1)"
  echo
  echo "## Results"
  echo
  echo "| Phase | Variant | Status | Detail |"
  echo "| --- | --- | --- | --- |"
  tail -n +2 $RESULTS | awk -F"\t" "{printf \"| %s | %s | %s | %s |\n\", \$1, \$2, \$3, \$4}"
  echo
  echo "## File sizes"
  echo
  echo "| Variant | .sum (bytes) | .out (bytes) | md5(.sum) |"
  echo "| --- | --- | --- | --- |"
  for V in "${EASTERN_VARIANTS[@]}"; do
    SUM=$TEST_ROOT/run_${V}/net01.sum
    OUT=$TEST_ROOT/run_${V}/net01.out
    if [ -f $SUM ]; then
      SS=$(stat -c %s $SUM)
      OS=$(stat -c %s $OUT 2>/dev/null || echo 0)
      MD=$(md5sum $SUM | cut -d" " -f1)
      echo "| $V | $SS | $OS | $MD |"
    else
      echo "| $V | - | - | - |"
    fi
  done
  echo
  echo "## Pass summary"
  echo
  TOT=$(tail -n +2 $RESULTS | wc -l)
  PASS=$(grep -c -P "\tPASS\t" $RESULTS || true)
  FAIL=$(grep -c -P "\tFAIL\t" $RESULTS || true)
  SKIP=$(grep -c -P "\tSKIP\t" $RESULTS || true)
  echo "- Total checks: $TOT"
  echo "- PASS: $PASS"
  echo "- FAIL: $FAIL"
  echo "- SKIP: $SKIP"
  echo
  if [ "$FAIL" = "0" ]; then
    echo "**Overall: PASS** — all $TOT checks passed."
  else
    echo "**Overall: PARTIAL** — $FAIL of $TOT checks failed."
  fi
} > $REPORT_DST
cat $REPORT_DST
log "Report written to $REPORT_DST"

# Copy test artifacts back for later debugging
ARTIFACTS=$PROJ/calibration/analysis/acd_stand_level_2026-05-16/integration_test
mkdir -p $ARTIFACTS
cp $RESULTS $ARTIFACTS/results.tsv
cp $TEST_ROOT/build_*.log $ARTIFACTS/ 2>/dev/null
cp $TEST_ROOT/run_*.log   $ARTIFACTS/ 2>/dev/null
cp $TEST_ROOT/smoke.log   $ARTIFACTS/ 2>/dev/null

# Capture stand S248112 year-2090 line from each variant
echo "Stand S248112 UNTHINNED CONTROL year 2090, all variants:" > $ARTIFACTS/year2090_per_variant.txt
for V in "${EASTERN_VARIANTS[@]}"; do
  S=$TEST_ROOT/run_${V}/net01.sum
  if [ -f $S ]; then
    YEAR=$(grep "^2090 160" $S | head -1)
    echo "$V: $YEAR" >> $ARTIFACTS/year2090_per_variant.txt
  fi
done
cat $ARTIFACTS/year2090_per_variant.txt

echo "DONE"
