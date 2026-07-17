#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:35:00
#SBATCH -N 1 -n 4 --mem=32G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/tune_%A_%a.out
set -e; cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla conus_eq_projector_v4.R --variant=NE --mode=dependent --nstands=300 --seed=7 --cyclelen=1   --st_rate=$1 --st_exp=$2 --ig_rdcap=$3 --ig_p=$4 --outdir=$5
