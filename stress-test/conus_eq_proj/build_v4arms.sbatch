#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J build_v4arms
#SBATCH -t 1:30:00
#SBATCH -N 1 -n 8 --mem=96G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/build_v4arms_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
for b in v4_b2 engine_default_ht engine_calibrated_ht; do
  echo "########## TRAJ $b ##########"
  Rscript build_singlecohort_traj_v3.R $b 10 40 || echo "TRAJ FAILED: $b"
done
echo BUILD_ALL_DONE
