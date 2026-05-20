#!/bin/bash
#SBATCH --job-name=fvs_mod_dg_kue_v2
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --time=04:00:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

## Modifier extraction for the May 7 banked DG Kuehne v2 fit
## (output/conus/dg/dg_kuehne_cspi_traits1_fit.rds, 13.6 GB).
## Re-extracts residuals (existing ones are stale May 3) then fits the
## modifier_common.stan at lambda = 5, 10, 20 so 60b_pivot_modifier_table.R
## can include Kuehne in the cross-model wide table.
##
## Output goes to output/conus/dg/modifier/ (NOT modifier_kuehne/) so the
## 60b pivot script picks it up automatically — fixes the May 5 anti-pattern
## where modifier_kuehne/ was excluded from the pivot.

set -uo pipefail
source /etc/profile.d/lmod.sh
module load gcc/12.3.0
module load R/4.4.0
cd $HOME/fvs-modern

mkdir -p calibration/output/conus/dg/modifier

#--- Stage A: extract fresh residuals from the May 7 banked Kuehne v2 fit ---
echo "=== Extracting residuals from Kuehne v2 fit ==="
Rscript calibration/R/50_extract_base_residuals.R \
  --fit      calibration/output/conus/dg/dg_kuehne_cspi_traits1_fit.rds \
  --pairs    calibration/data/conus_remeasurement_pairs_metric_cond.rds \
  --family   log \
  --response dg_obs_a \
  --model    dg_kuehne \
  --out      calibration/output/conus/dg/dg_kuehne_residuals.rds

#--- Stage B: fit modifier coefficients at three lambdas ---
for LAMBDA in 5 10 20; do
  echo "=== Fitting modifier at lambda=$LAMBDA ==="
  Rscript calibration/R/51_fit_modifier.R \
    --residuals calibration/output/conus/dg/dg_kuehne_residuals.rds \
    --stan_file calibration/stan/modifier_common.stan \
    --lambda    $LAMBDA \
    --chains    2 \
    --warmup    500 \
    --sampling  500 \
    --max_treedepth 10 \
    --adapt_delta   0.95 \
    --seed      42 \
    --out       calibration/output/conus/dg/modifier
done

echo "=== Done. Outputs in output/conus/dg/modifier/ ==="
ls -lah calibration/output/conus/dg/modifier/
