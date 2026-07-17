#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J v5seed
#SBATCH -t 4:00:00
#SBATCH -N 1 -n 8 --mem=96G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/v5seed_%A_%a.out
#SBATCH --array=1-3
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
VARS=(NONE NE SN PN)
V=${VARS[$SLURM_ARRAY_TASK_ID]}
Rscript --vanilla conus_eq_projector_v5.R --variant=$V --mode=dependent --nstands=0 --cyclelen=1 \
  --st_rate=0.20 --st_exp=2.20 --st_bias=0.0 \
  --standinit_dir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4 \
  --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/out_v5_seed
echo DONE $V
