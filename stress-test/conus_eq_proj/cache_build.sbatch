#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v3_cache_build
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 8 --mem=160G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest/cache_%j.log
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v3.R --variant=AK --mode=dependent --nstands=50 --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest/cachebuild 2>&1 | grep -E 'cached covariate|stands projected|Realized'
echo CACHE_DONE
