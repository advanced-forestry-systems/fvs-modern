#!/bin/bash
#SBATCH --job-name=cspi4_compare
#SBATCH --time=00:20:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/cspi4_compare_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/cspi4_compare_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
Rscript --vanilla /users/PUOM0008/crsfaaron/compare_cspi4_vs_climsi.R
