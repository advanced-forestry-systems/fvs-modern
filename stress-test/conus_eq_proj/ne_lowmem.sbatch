#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v3_ne_lowmem
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 8 --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest/ne_lowmem_%j.log
module load gcc/12.3.0 R/4.4.0
/usr/bin/time -v Rscript --vanilla /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=0 --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest/nelowmem 2>&1 | grep -E 'loaded cached|stands projected|Realized self-thinning|Maximum resident set'
echo NE_LOWMEM_DONE
