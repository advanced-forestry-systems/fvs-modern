#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J sc_bakuzis
#SBATCH -t 2:00:00
#SBATCH --mem=96G
#SBATCH -c 8
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/sc_bakuzis_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
OUT=bakuzis_singlecohort
mkdir -p $OUT
BUILDS="engine_default engine_calibrated conus_b2 conus_b1 greg gompit sdicon_b2 sdicon_b1"
echo "=================== BUILD TRAJECTORIES (single-cohort) ==================="
for b in $BUILDS; do
  if [ -f "$OUT/${b}_sc_all.csv" ]; then echo "SKIP traj $b (exists)"; continue; fi
  echo "########## TRAJ $b ##########"
  Rscript build_singlecohort_traj_v2.R $b 10 40 || echo "TRAJ FAILED: $b"
done
echo
echo "=================== BAKUZIS CHECK + MATRIX per build ==================="
for b in $BUILDS; do
  for sfx in _sc_all _sc_ne _sc_sn _sc_pn; do
    csv=$OUT/${b}${sfx}.csv
    if [ -f "$csv" ]; then
      echo "########## CHECK ${b}${sfx} ##########"
      Rscript bakuzis_check.R "$csv" 0.40 > $OUT/${b}${sfx}_check.txt 2>&1 || echo "check failed $csv"
      Rscript bakuzis_matrix.R "$csv" $OUT/${b}${sfx}_matrix.png 0.40 >/dev/null 2>&1 || echo "matrix failed $csv"
      echo "wrote $OUT/${b}${sfx}_check.txt"
    fi
  done
done
echo "ALL DONE"
