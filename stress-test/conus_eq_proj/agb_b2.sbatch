#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:20:00
#SBATCH -N 1 -n 4
#SBATCH --mem=16G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/agb_b2_%j.out
module load gcc/12.3.0 python/3.12
export PYTHONNOUSERSITE=1
export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
python conus_eq_agb.py --metrics conus_eq_ne_conus_b2_metrics.csv --treelists conus_eq_ne_conus_b2_treelists.csv
