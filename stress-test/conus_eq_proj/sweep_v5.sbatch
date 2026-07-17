#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J sweepv5
#SBATCH -t 2:00:00
#SBATCH -N 1 -n 4 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/sweepv5_%A_%a.out
#SBATCH --array=1-5
set -e
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
RATE=(x 0.28 0.24 0.20 0.28 0.24)
EXP=(x 2.40 2.30 2.20 2.40 2.30)
BIAS=(x 0.0 0.0 0.0 0.5 0.5)
TAG=(x A B C D E)
i=$SLURM_ARRAY_TASK_ID
OD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/sweep_v5/set_${TAG[$i]}
mkdir -p $OD
Rscript --vanilla conus_eq_projector_v5.R --variant=NE --mode=dependent --nstands=0 --cyclelen=1 \
  --st_rate=${RATE[$i]} --st_exp=${EXP[$i]} --st_bias=${BIAS[$i]} \
  --standinit_dir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_seed_v4 \
  --outdir=$OD > $OD/run.log 2>&1
Rscript --vanilla score_v5.R $OD/conus_eq_ne_conus_b2_metrics.csv NE "${TAG[$i]}_rate${RATE[$i]}_exp${EXP[$i]}_bias${BIAS[$i]}"
echo SWEEPDONE ${TAG[$i]}
