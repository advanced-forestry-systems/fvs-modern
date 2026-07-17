#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/val_v2_%j.out
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
cd $D
module load gcc/12.3.0 R/4.4.0
Rscript conus_eq_projector_v2.R --variant=NE --mode=dependent --nstands=1500 --seed=7 --outdir=$D/val_v2
Rscript conus_eq_projector_v2.R --variant=NE --mode=free      --nstands=1500 --seed=7 --outdir=$D/val_v2
module load python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python conus_eq_agb.py --metrics $D/val_v2/conus_eq_ne_conus_b2_metrics.csv --treelists $D/val_v2/conus_eq_ne_conus_b2_treelists.csv
python conus_eq_agb.py --metrics $D/val_v2/conus_eq_ne_conus_b1_metrics.csv --treelists $D/val_v2/conus_eq_ne_conus_b1_treelists.csv
module load R/4.4.0
Rscript validate_ne_v2.R val_v2
echo ALL_DONE
