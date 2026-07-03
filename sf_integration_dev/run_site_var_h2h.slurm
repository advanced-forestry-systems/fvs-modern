#!/bin/bash
#SBATCH --job-name=site_h2h
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/site_h2h_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/site_h2h_%j.err
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla /users/PUOM0008/crsfaaron/site_var_headtohead.R
