#!/bin/bash
#SBATCH --job-name=mortcap_bakuzis_v2
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_bakuzis_v2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_bakuzis_v2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0
set -uo pipefail   # NOTE: no -e; FVS exit 10 (ICCODE 1, informational) is benign

echo "=== BAKUZIS MORTCAP v2 START $(date) ==="
echo "Node: $(hostname)"

BIN=/users/PUOM0008/crsfaaron/fvs-modern/bin-mortcap/FVSne
CONFDIR=/fs/scratch/PUOM0008/crsfaaron/wt-rescore/config
SRC=/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/fvs_input_harness/stands
DG_COEF=$CONFDIR/greg_dg_coefficients_decel.csv
HT_COEF=$CONFDIR/conus_htdbh_coefficients.csv
GOMP_COEF=/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv
BASE=/fs/scratch/PUOM0008/crsfaaron/rescore_mortcap
mkdir -p "$BASE/logs" "$BASE/work"

echo "Binary: $BIN"; ls -lh "$BIN"
echo "Coef check:"; for c in "$DG_COEF" "$HT_COEF" "$GOMP_COEF"; do [ -f "$c" ] && echo "  OK  $c" || echo "  MISSING $c"; done

for CFG in native allon_mortcap; do
  WORK="$BASE/work/$CFG"
  rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
  for k in "$SRC"/*.key; do
    b=$(basename "${k%.key}")
    cp "$SRC/$b.tre" "$WORK/" 2>/dev/null || true
    sed -e "s/^NUMCYCLE.*/NUMCYCLE          30/" \
        -e "s/^TIMEINT.*/TIMEINT           1        10/" \
        -e "/^MORTMULT/d" "$k" > "$WORK/$b.key"
  done
  unset FVS_GREGDG FVS_GREGDG_COEF FVS_GREGHTDBH FVS_GREGHTDBH_COEF \
        FVS_GOMPIT FVS_GOMPIT_COEF FVS_GOMP_MORTMULT FVS_GOMP_SFLOOR \
        FVS_GOMP_DBHMIN 2>/dev/null || true
  case "$CFG" in
    native) ;;
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
    "$BIN" --keywordfile="$k" > /dev/null 2>"$b.err" || true
    if [ -s "$b.out" ]; then ok=$((ok+1)); else fail=$((fail+1)); fi
  done
  echo "CFG=$CFG  OK=$ok  FAIL=$fail"
done

echo "=== SCORING $(date) ==="
module load R/4.4.0
Rscript /users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/score_mortcap_v2.R "$BASE"
echo "=== ALL DONE $(date) ==="
