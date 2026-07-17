#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/val_b2_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript conus_eq_projector.R --variant=NE --mode=dependent --nstands=300 --seed=7
module load python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python conus_eq_agb.py --metrics conus_eq_ne_conus_b2_metrics.csv --treelists conus_eq_ne_conus_b2_treelists.csv
