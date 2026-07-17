#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v3_memtest_ne
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest/ne_%j.log
module load gcc/12.3.0 R/4.4.0
P=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
O=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/v3_memtest; mkdir -p $O
/usr/bin/time -v Rscript --vanilla $P/conus_eq_projector_v3.R --variant=NE --mode=dependent --nstands=0 --cyclelen=1 --outdir=$O 2>&1 | grep -E "stands projected|Realized self-thinning|Maximum resident|Wrote:|Error"
echo MEMTEST_DONE
