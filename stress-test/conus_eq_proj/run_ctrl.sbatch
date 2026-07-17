#!/bin/bash
#SBATCH -J bak_ctrl
#SBATCH -A PUOM0008
#SBATCH -p cpu
#SBATCH -t 00:40:00
#SBATCH -N 1 -n 1 -c 4
#SBATCH --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bak_ctrl_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj && Rscript build_bakuzis_controlled.R
