#!/bin/bash
#SBATCH -A PUOM0008
#SBATCH -t 00:30:00
#SBATCH -N 1 -n 8
#SBATCH --mem=80G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/extract_b1_%j.out
module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
Rscript extract_b1_re.R
