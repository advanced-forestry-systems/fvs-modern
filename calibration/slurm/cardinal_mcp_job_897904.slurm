#!/bin/bash
#SBATCH --job-name=mortcap_build_run
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_build_run_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_build_run_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0
set -euo pipefail

echo "=== MORTCAP BUILD+RUN START $(date) ==="
echo "Node: $(hostname)"

# ── Step 1: build FVSne with patched gompmort ──────────────────────────────
BINDIR=/users/PUOM0008/crsfaaron/fvs-modern/bin-mortcap
mkdir -p "$BINDIR"
cd /users/PUOM0008/crsfaaron/fvs-modern
bash deployment/scripts/build_fvs_executables.sh . "$BINDIR" ne 2>&1 | tail -40
echo "=== BUILD DONE $(date) ==="
ls -lh "$BINDIR/FVSne"
nm "$BINDIR/FVSne" | grep -i "gompsurv\|gmortmult\|gsfloor" | head -10

# ── Step 2: Bakuzis long-horizon test with mortcap allon config ───────────
BIN="$BINDIR/FVSne"
CONFDIR=/fs/scratch/PUOM0008/crsfaaron/wt-rescore/config
SRC=/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/fvs_input_harness/stands
DG_COEF=$CONFDIR/greg_dg_coefficients_decel.csv
HT_COEF=$CONFDIR/conus_htdbh_coefficients.csv
GOMP_COEF=/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv
BASE=/fs/scratch/PUOM0008/crsfaaron/rescore_mortcap
mkdir -p "$BASE/logs" "$BASE/work"

echo ""
echo "=== BAKUZIS MORTCAP RUN START $(date) ==="
echo "Binary: $BIN"
echo "FVS_GOMP_MORTMULT=0.5  FVS_GOMP_SFLOOR=0.95"

for CFG in native allon_mortcap; do
  WORK="$BASE/work/$CFG"
  rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
  for k in $SRC/*.key; do
    b=$(basename "${k%.key}")
    cp "$SRC/$b.tre" "$WORK/"
    sed -e "s/^NUMCYCLE.*/NUMCYCLE          30/" \
        -e "s/^TIMEINT.*/TIMEINT           1        10/" \
        -e "/^MORTMULT/d" "$k" > "$WORK/$b.key"
  done
  unset FVS_GREGDG FVS_GREGDG_COEF FVS_GREGHTDBH FVS_GREGHTDBH_COEF \
        FVS_GOMPIT FVS_GOMPIT_COEF FVS_GOMP_MORTMULT FVS_GOMP_SFLOOR \
        FVS_GOMP_DBHMIN FVS_GOMP_SFLOOR 2>/dev/null || true
  case "$CFG" in
    native)
      ;;
    allon_mortcap)
      export FVS_GREGDG=1 FVS_GREGDG_COEF="$DG_COEF"
      export FVS_GREGHTDBH=1 FVS_GREGHTDBH_COEF="$HT_COEF"
      export FVS_GOMPIT=1 FVS_GOMPIT_COEF="$GOMP_COEF"
      export FVS_GOMP_MORTMULT=0.5
      export FVS_GOMP_SFLOOR=0.95
      ;;
  esac
  ok=0; fail=0
  for k in *.key; do
    b="${k%.key}"
    "$BIN" --keywordfile="$k" > /dev/null 2>"$b.err"
    if [ -s "$b.out" ]; then ok=$((ok+1)); else fail=$((fail+1)); fi
  done
  echo "CFG=$CFG  OK=$ok  FAIL=$fail"
  # Quick BA check: last-cycle BA from .out files
  echo "  BA at final cycle (sample 5 stands):"
  ls "$WORK"/*.out 2>/dev/null | head -5 | while read f; do
    tail -20 "$f" | grep -E "Basal|BASAL|BA " | tail -1 || true
  done
done

echo "=== ALL DONE $(date) ==="
echo "Results in: $BASE/work/"

