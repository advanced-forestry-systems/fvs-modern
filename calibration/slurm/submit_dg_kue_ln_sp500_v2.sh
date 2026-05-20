#!/bin/bash
#SBATCH --job-name=fvs_prod_dg_kue_ln_sp500_v2
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --time=2-00:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## DG Kuehne refit v2 (2026-05-05, recovered 2026-05-06). Stan likelihood
## scale guarded via fmin(fmax(sigma / sqrt_years[i], 1e-4), 50.0).
## Memory 256G to accommodate min_sp_n=500 species count (~50+).

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

Rscript calibration/R/31_fit_dg_multimodel.R \
  --form           kuehne \
  --stan_file      calibration/stan/dg_kuehne2022_lognormal.stan \
  --data           calibration/data/conus_remeasurement_pairs_metric.rds \
  --stan_dir       calibration/stan \
  --traits_rds     calibration/traits/species_traits.rds \
  --traits         1 \
  --site           cspi \
  --out            calibration/output/conus/dg \
  --n              500000 \
  --chains         4 \
  --warmup         1000 \
  --sampling       1000 \
  --max_treedepth  8 \
  --adapt_delta    0.9 \
  --seed           42 \
  --min_sp_n       500 \
  --save_draws
