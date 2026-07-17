#!/bin/bash
#SBATCH -A PUOM0008 -t 01:00:00 -N 1 -n 4 --mem=64G
#SBATCH -o /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/logs/matrix_%j.out
module load gcc/12.3.0 R/4.4.0
Rscript ~/growth_mort_matrix.R
