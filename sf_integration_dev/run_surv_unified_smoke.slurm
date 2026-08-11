#!/bin/bash
#SBATCH --job-name=surv_unified_smoke
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/surv_unified_smoke_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/surv_unified_smoke_%j.err
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
O=calibration/output/conus/sf_integration/unified_smoke
mkdir -p \$O
Rscript --vanilla calibration/R/34_fit_survival_speciesfree.R \
  --stan_file=calibration/stan/survival_unified.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --subsample=40000 --smoke --outdir=\$O --outname=surv_unified_smoke 2>&1 | tail -30
