#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -J greg_smoke
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 4 --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/greg_smoke_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript conus_eq_projector_greg.R --variant=NE --nstands=40 --seed=7 --outdir=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/greg_smoke
