#!/bin/bash
#SBATCH --job-name=fvs_prod_dg_org_ln_v3
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=2-12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## DG ORGANON lognormal v3 (2026-05-08). v2 (9217968) completed cleanly
## (0 divergent, 0 treedepth, 0 BFMI) but rhat 1.86-3.79 across 16 fixed
## effects with ESS 4.4-5.8. Diagnosed as architectural identifiability:
## K1*b1, K2*b2, and the b3 terms defined trade-off ridges along which
## chains found different modes. Three changes from v2:
##   1. K1, K2 fixed to 1.0 / 1.0 in the Stan model (was: estimated).
##      Eliminates the worst trade-off ridges.
##   2. Tightened priors on mu_b0, b2 (rescaled for K2=1), and sigma_eco.
##   3. Warmup 1000 -> 2000, adapt_delta 0.9 -> 0.95.

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

Rscript calibration/R/31_fit_dg_organon.R \
  --data           calibration/data/conus_remeasurement_pairs_metric_cond.rds \
  --stan_file      calibration/stan/dg_organon_lognormal_v3.stan \
  --site           cspi \
  --out            calibration/output/conus/dg \
  --n              500000 \
  --chains         4 \
  --warmup         2000 \
  --sampling       1000 \
  --max_treedepth  10 \
  --adapt_delta    0.95 \
  --seed           42 \
  --save_draws
