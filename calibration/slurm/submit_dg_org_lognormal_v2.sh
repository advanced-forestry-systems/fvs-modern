#!/bin/bash
#SBATCH --job-name=fvs_prod_dg_org_ln_v2
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=192G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## DG ORGANON lognormal v2 (2026-05-05, recovered 2026-05-06). R driver
## patched to honor --site flag; data file is metric_cond.rds (post-COND
## join) to match other production fits. Memory 192G.

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

Rscript calibration/R/31_fit_dg_organon.R \
  --data           calibration/data/conus_remeasurement_pairs_metric_cond.rds \
  --stan_file      calibration/stan/dg_organon_lognormal.stan \
  --site           cspi \
  --out            calibration/output/conus/dg \
  --n              500000 \
  --chains         4 \
  --warmup         1000 \
  --sampling       1000 \
  --max_treedepth  8 \
  --adapt_delta    0.9 \
  --seed           42 \
  --save_draws
