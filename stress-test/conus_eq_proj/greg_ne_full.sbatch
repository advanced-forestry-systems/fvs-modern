#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 03:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=248G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_ne_full_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
D=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript --vanilla $D/conus_eq_projector_greg.R --variant=NE --nstands=0 --outdir=$D/out_conus_eq_greg
echo NE_FULL_DONE
Rscript --vanilla $D/validate_ne_greg.R
echo VALIDATE_DONE
