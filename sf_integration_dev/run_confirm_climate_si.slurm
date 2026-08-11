#!/bin/bash
#SBATCH --job-name=confirm_csi
#SBATCH --time=00:30:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/confirm_csi_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/confirm_csi_%j.err
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
mkdir -p /users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration
Rscript --vanilla /users/PUOM0008/crsfaaron/fvs-modern/sf_integration_dev/confirm_climate_si.R
