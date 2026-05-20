#!/bin/bash
#SBATCH --job-name=sanity_traits
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

## Decision-quality summary for species-free model design.
## Reads each banked base-model fit and extracts sigma_sp, sigma_eco,
## sigma, and gamma_trait coefficients. Outputs CSV under MEMORY/.

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

Rscript calibration/R/sanity_check_traits.R
