#!/bin/bash
#SBATCH --job-name=acd_dg_refit
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_dg_refit_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_dg_refit_%j.err

set -euo pipefail
module purge
module load gcc/12.3.0 R/4.4.0

export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
cd $FVS_PROJECT_ROOT

# Longer-warmup HMC config aimed at rhat < 1.05 for ACD diameter growth.
# Prior run with warmup=500 sampling=500 adapt_delta=0.90 produced rhat=1.95.
# Quadruple warmup, triple sampling, tighten step size, deepen tree.
export FVS_HMC_WARMUP=2000
export FVS_HMC_SAMPLING=1500
export FVS_HMC_ADAPT_DELTA=0.99
export FVS_HMC_TREEDEPTH=12
export FVS_HMC_CHAINS=4
export FVS_MAX_OBS=10000

# Backup existing posterior before overwrite
SNAP=$FVS_PROJECT_ROOT/calibration/output/variants/acd/diameter_growth_posterior.csv.refit_pre_$(date +%Y%m%d_%H%M%S)
cp $FVS_PROJECT_ROOT/calibration/output/variants/acd/diameter_growth_posterior.csv $SNAP

Rscript calibration/R/02c_fit_dg_hmc_small.R --variant acd 2>&1 | tee calibration/logs/02c_acd_${SLURM_JOB_ID}.console.log

# Report rhat status
Rscript -e '
post <- read.csv("calibration/output/variants/acd/diameter_growth_posterior.csv")
key <- post[post$variable %in% c("mu_b0","sigma_b0","sigma","mu_b1","mu_b2"), ]
print(key)
cat("\n=== max rhat across all variables ===\n")
cat(max(post$rhat, na.rm=TRUE), "\n")
cat("=== converged fraction ===\n")
cat(mean(post$converged, na.rm=TRUE), "\n")
'

echo "DONE"
