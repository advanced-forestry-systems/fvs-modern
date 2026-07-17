#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -p hugemem
#SBATCH -N 1 -n 46
#SBATCH --mem=900G
#SBATCH -t 20:00:00
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/greg_agb_ne_%j.out
module load gcc/12.3.0 python/3.12
export PYTHONNOUSERSITE=1
export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_agb.py --metrics /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq_greg/conus_eq_ne_greg_metrics.csv --treelists /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_conus_eq_greg/conus_eq_ne_greg_treelists.csv
echo GREG_AGB_DONE_ne
