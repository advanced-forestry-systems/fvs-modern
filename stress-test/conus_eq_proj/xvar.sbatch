#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J xvar_check
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 2 --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/xvar_%j.log
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/xvar.R
echo XVAR_DONE
