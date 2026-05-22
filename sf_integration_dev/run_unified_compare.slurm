#!/bin/bash
#SBATCH --job-name=unified_compare
#SBATCH --time=03:00:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/unified_compare_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/unified_compare_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
# Unified vs species-free LOO + RE-scale decomposition. Run AFTER unified
# production fits land (hcb 10302998, htdbh 10302999).
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
BASE=/users/PUOM0008/crsfaaron/fvs-conus/output/conus
OUT=$BASE/sf_integration/unified_compare
mkdir -p $OUT; cd $OUT
run_one () {
  local comp=$1 unidir=$2 sfdir=$3
  local uloo=$(ls -t $unidir/*_loo.rds 2>/dev/null | head -1)
  local usum=$(ls -t $unidir/*_summary.csv 2>/dev/null | head -1)
  local sfit=$(ls -tL $sfdir/*_fit.rds 2>/dev/null | head -1)
  local ssum=$(ls -t $sfdir/*_summary.csv 2>/dev/null | head -1)
  echo "### $comp"; echo "  uni_loo=$uloo"; echo "  sf_fit=$sfit"
  if [ -z "$uloo" ] || [ -z "$sfit" ]; then echo "  SKIP ($comp): missing inputs"; return; fi
  Rscript --vanilla /users/PUOM0008/crsfaaron/fvs-modern/sf_integration_dev/compare_unified_vs_sf.R \
    --comp=$comp --unified_loo=$uloo --unified_summary=$usum \
    --sf_fit=$sfit --sf_summary=$ssum --out=$OUT/compare_$comp 2>&1
}
run_one hcb   $BASE/hcb/unified_v2split   $BASE/hcb/sf_v2split
run_one htdbh $BASE/htdbh/unified_v2split $BASE/ht_dbh/sf_v2split
echo "DONE unified_compare"
