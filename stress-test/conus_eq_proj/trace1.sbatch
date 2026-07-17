#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v4trace1
#SBATCH -t 0:15:00
#SBATCH -N 1 -n 2 --mem=24G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/v4trace1_%j.out
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
export TRACE_RECRUITS=815194084290487
Rscript --vanilla conus_eq_projector_v4.R --variant=NE --mode=dependent --nstands=0 --cyclelen=1 \
  --standinit_dir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_trace1 \
  --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_v4_trace1
echo TRACE1_DONE
