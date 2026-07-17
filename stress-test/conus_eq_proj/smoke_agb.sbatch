#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:30:00
#SBATCH -N 1 -n 4
#SBATCH --mem=24G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/smoke_agb_%j.out
module load gcc/12.3.0 python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
python conus_eq_agb.py --metrics smoketest/conus_eq_ne_conus_b2_metrics.csv --treelists smoketest/conus_eq_ne_conus_b2_treelists.csv
