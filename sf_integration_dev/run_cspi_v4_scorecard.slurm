#!/bin/bash
#SBATCH --job-name=cspi_v4_score
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/cspi_v4_score_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/cspi_v4_score_%j.err
set -uo pipefail
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
Rscript --vanilla /users/PUOM0008/crsfaaron/cspi_v4_scorecard.R
