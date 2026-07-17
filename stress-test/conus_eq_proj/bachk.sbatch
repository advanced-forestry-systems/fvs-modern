#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:08:00
#SBATCH -N 1 -n 4
#SBATCH --mem=24G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bachk_%j.out
module load gcc/12.3.0 R/4.4.0
Rscript /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/bachk.R
