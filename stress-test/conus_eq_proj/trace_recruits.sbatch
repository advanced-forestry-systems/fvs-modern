#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v4trace
#SBATCH -t 0:30:00
#SBATCH -N 1 -n 4 --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/v4trace_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
export TRACE_RECRUITS=72745554010538
Rscript --vanilla conus_eq_projector_v4.R --variant=NE --mode=dependent --nstands=40 --seed=7 --cyclelen=1 \
  --standinit_dir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4 \
  --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_v4_trace
echo TRACE_DONE
