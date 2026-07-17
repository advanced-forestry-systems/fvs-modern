#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 4
#SBATCH --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_smoke_ne_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript --vanilla $D/conus_eq_projector_greg.R --variant=NE --nstands=200 --outdir=$D/greg_smoke
echo SMOKE_DONE
