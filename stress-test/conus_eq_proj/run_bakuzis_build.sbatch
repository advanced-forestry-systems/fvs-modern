#!/bin/bash
#SBATCH -J bak_build
#SBATCH -A PUOM0008
#SBATCH -p cpu
#SBATCH -t 02:00:00
#SBATCH -N 1 -n 1 -c 8
#SBATCH --mem=110G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bak_build2_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript build_bakuzis_traj.R
