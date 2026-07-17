#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:25:00
#SBATCH -N 1 -n 8
#SBATCH --mem=48G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/smoke_dep_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript conus_eq_projector.R --variant=NE --mode=dependent --nstands=30 --seed=7
