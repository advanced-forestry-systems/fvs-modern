#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/smoke_htbio_%j.out
set -e
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
cd $D
module load gcc/12.3.0 R/4.4.0
Rscript conus_eq_projector_v2.R --variant=NE --mode=dependent --nstands=300 --seed=7 --outdir=$D/smoke_htbio
module load python/3.12; export PYTHONNOUSERSITE=1; export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
python conus_eq_agb.py --metrics $D/smoke_htbio/conus_eq_ne_conus_b2_metrics.csv --treelists $D/smoke_htbio/conus_eq_ne_conus_b2_treelists.csv
echo "=== HT_BIO_M present? ==="
head -1 $D/smoke_htbio/conus_eq_ne_conus_b2_treelists.csv
echo SMOKE_DONE
