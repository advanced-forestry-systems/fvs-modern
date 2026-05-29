#!/bin/bash
#SBATCH --job-name=hg_v6_relht_smoke
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_v6_relht_smoke_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_v6_relht_smoke_%j.err
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
mkdir -p /users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration/unified_smoke
cd /users/PUOM0008/crsfaaron/fvs-modern
Rscript --vanilla calibration/R/32_fit_hg_speciesfree_v5.R \
  --stan_file=calibration/stan/hg_v6_relht.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --subsample=40000 --smoke --compact \
  --outdir=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration/unified_smoke \
  --outname=hg_v6_relht_smoke
echo DONE
