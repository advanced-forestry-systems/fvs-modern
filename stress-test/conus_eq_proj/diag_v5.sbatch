#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J diagv5
#SBATCH -t 1:00:00
#SBATCH -N 1 -n 4 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/diagv5_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
export DIAG_THIN=815194084290487
# v4-equivalent: st_bias=0, st_rate=0.34 st_exp=2.5; run full NE seed so the trace stand is hit
Rscript --vanilla conus_eq_projector_v5.R --variant=NE --mode=dependent --nstands=0 --cyclelen=1 \
  --st_rate=0.34 --st_exp=2.5 --st_bias=0.0 \
  --standinit_dir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4 \
  --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/diag_v5 2>&1 | grep -E 'DIAGTHIN|self-thinning slope|year 0' 
echo DIAGDONE
