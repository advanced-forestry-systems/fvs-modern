#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J build_v4only
#SBATCH -t 0:40:00
#SBATCH -N 1 -n 8 --mem=96G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/build_v4only_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript build_singlecohort_traj_v3.R v4_b2 10 40
echo V4_BUILD_DONE
