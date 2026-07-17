#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -p hugemem
#SBATCH -N 1 -n 46
#SBATCH --mem=900G
#SBATCH -t 20:00:00
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/agb_patch_sn_%j.out
module load gcc/12.3.0 python/3.12
export PYTHONNOUSERSITE=1
export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
cd $D
python conus_eq_agb.py --metrics $D/out_conus_eq_b2_gompit/conus_eq_sn_conus_b2_gompit_metrics.csv --treelists $D/out_conus_eq_b2_gompit/conus_eq_sn_conus_b2_gompit_treelists.csv
echo AGB_PATCH_DONE_SN
