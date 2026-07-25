#!/bin/bash
#SBATCH --job-name=mortcap_sweep
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_sweep_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/mortcap_sweep_%j.err

module load gcc/12.3.0
set -uo pipefail
echo "=== MORTCAP SWEEP START $(date) ==="

BIN=/users/PUOM0008/crsfaaron/fvs-modern/bin-mortcap/FVSne
CONFDIR=/fs/scratch/PUOM0008/crsfaaron/wt-rescore/config
SRC=/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/fvs_input_harness/stands
DG_COEF=$CONFDIR/greg_dg_coefficients_decel.csv
HT_COEF=$CONFDIR/conus_htdbh_coefficients.csv
GOMP_COEF=/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv
BASE=/fs/scratch/PUOM0008/crsfaaron/rescore_mortcap/sweep
mkdir -p "$BASE"
SUMMARY="$BASE/sweep_median_BA.csv"
echo "sfloor,gmortmult,year,median_BA_m2ha" > "$SUMMARY"

FT2AC_TO_M2HA=0.2295684

# sweep grid: SFLOOR values at GMORTMULT 0.5, plus one tighter GMORTMULT combo
for COMBO in "0.97 0.5" "0.98 0.5" "0.99 0.5" "0.995 0.5" "0.99 0.3"; do
  SFLOOR=$(echo $COMBO | cut -d' ' -f1); GMM=$(echo $COMBO | cut -d' ' -f2)
  TAG="sf${SFLOOR}_gm${GMM}"
  WORK="$BASE/$TAG"; rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK"
  for k in "$SRC"/*.key; do
    b=$(basename "${k%.key}")
    cp "$SRC/$b.tre" "$WORK/" 2>/dev/null || true
    sed -e "s/^NUMCYCLE.*/NUMCYCLE          30/" \
        -e "s/^TIMEINT.*/TIMEINT           1        10/" \
        -e "/^MORTMULT/d" "$k" > "$WORK/$b.key"
  done
  export FVS_GREGDG=1 FVS_GREGDG_COEF="$DG_COEF"
  export FVS_GREGHTDBH=1 FVS_GREGHTDBH_COEF="$HT_COEF"
  export FVS_GOMPIT=1 FVS_GOMPIT_COEF="$GOMP_COEF"
  export FVS_GOMP_MORTMULT=$GMM FVS_GOMP_SFLOOR=$SFLOOR
  for k in *.key; do "$BIN" --keywordfile="$k" > /dev/null 2>/dev/null || true; done
  # score inline with awk: for each .out, pull BA (field4) at AGE 100/200/300
  for YR in 100 200 300; do
    vals=$(for f in *.out; do
      awk -v yr=$YR 'f&&/^[0-9]{4}[ \t]+[0-9]+/{if($2==yr)print $4} /SUMMARY STATISTICS/{f=1}' "$f"
    done | sort -n)
    n=$(echo "$vals" | grep -c .); 
    med=$(echo "$vals" | awk '{a[NR]=$1} END{if(NR==0){print "NA"}else if(NR%2){print a[(NR+1)/2]}else{print (a[NR/2]+a[NR/2+1])/2}}')
    medm=$(awk -v b="$med" -v c=$FT2AC_TO_M2HA 'BEGIN{if(b=="NA")print"NA";else printf "%.2f", b*c}')
    echo "$SFLOOR,$GMM,$YR,$medm" >> "$SUMMARY"
  done
  echo "done $TAG"
done
echo "=== SWEEP RESULTS ==="
cat "$SUMMARY"
echo "=== DONE $(date) ==="
