#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 4
#SBATCH --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/reval_%j.out
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
cd $D
module load gcc/12.3.0 python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python conus_eq_agb.py --metrics $D/val_v2/conus_eq_ne_conus_b2_metrics.csv --treelists $D/val_v2/conus_eq_ne_conus_b2_treelists.csv
python conus_eq_agb.py --metrics $D/val_v2/conus_eq_ne_conus_b1_metrics.csv --treelists $D/val_v2/conus_eq_ne_conus_b1_treelists.csv
module load R/4.4.0
Rscript validate_ne_v2.R val_v2
echo REVAL_DONE
