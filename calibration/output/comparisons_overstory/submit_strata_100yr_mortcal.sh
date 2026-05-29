#!/bin/bash
#SBATCH --job-name=silc_strata_100yr_mc
#SBATCH --account=PUOM0008
#SBATCH --time=03:00:00
#SBATCH --mem=16G
#SBATCH --output=silc_strata_100yr_mc_%j.out
#SBATCH --error=silc_strata_100yr_mc_%j.err
set -euo pipefail
module load gcc/12.3.0 R/4.5.2
cd ~/silc_strata
Rscript run_silc_strata_100yr_mortcal.R
