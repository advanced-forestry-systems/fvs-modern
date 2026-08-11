#!/bin/bash
#SBATCH --job-name=hg_v6_relht_100k
#SBATCH --time=36:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_v6_relht_100k_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_v6_relht_100k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
mkdir -p /users/PUOM0008/crsfaaron/fvs-conus/output/conus/hg/v6_relht_v2split
cd /users/PUOM0008/crsfaaron/fvs-modern
Rscript --vanilla calibration/R/32_fit_hg_speciesfree_v5.R \
  --stan_file=calibration/stan/hg_v6_relht.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --outdir=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/hg/v6_relht_v2split \
  --outname=hg_v6_relht_100k \
  --subsample=100000 --compact
echo DONE
