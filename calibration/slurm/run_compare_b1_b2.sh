#!/bin/bash
#SBATCH --job-name=compare_b1_b2
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/MEMORY/compare_b1_b2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/MEMORY/compare_b1_b2_%j.err

module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
Rscript --vanilla calibration/R/compare_b1_b2_hg.R
