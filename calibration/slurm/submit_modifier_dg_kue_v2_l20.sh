#!/bin/bash
#SBATCH --job-name=fvs_mod_dg_kue_l20
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

Rscript calibration/R/51_fit_modifier.R \
  --residuals calibration/output/conus/dg/dg_kuehne_residuals.rds \
  --stan_file calibration/stan/modifier_common.stan \
  --lambda    20 \
  --chains    2 \
  --warmup    500 \
  --sampling  500 \
  --max_treedepth 10 \
  --adapt_delta   0.95 \
  --seed      42 \
  --out       calibration/output/conus/dg/modifier
