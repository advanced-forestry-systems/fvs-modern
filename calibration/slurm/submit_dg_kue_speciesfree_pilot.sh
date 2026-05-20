#!/bin/bash
#SBATCH --job-name=fvs_dg_kue_sf_b1_pilot
#SBATCH --time=18:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/logs/fvs_dg_kue_sf_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/logs/fvs_dg_kue_sf_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## DG_Kuehne species-free B1 pilot.
##
## Subsample of 200K observations (close to v2's 160K). z_sp removed from
## the linear predictor; trait_effect = W * gamma is the only species signal.
##
## Decision criterion when fit lands:
##   - v2 banked sigma = 1.913 [1.907, 1.918]
##   - B1 sigma within +/- 5 percent of 1.913 validates species-free on DG_Kue
##   - Use compare_b1_b2_dg_kue.R (to be written) for full posterior comparison.
##
## Expected wall time: 7-9 hours (matches v2 7.2h benchmark).

set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern

mkdir -p logs calibration/output/conus/dg/speciesfree_pilot

Rscript --vanilla calibration/R/32_fit_dg_kuehne_speciesfree.R \
  --stan_file=calibration/stan/dg_kuehne2022_speciesfree.stan \
  --outdir=calibration/output/conus/dg/speciesfree_pilot \
  --outname=dg_kuehne_cspi_traits1_b1 \
  --subsample=100000
