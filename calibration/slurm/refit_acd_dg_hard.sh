#!/bin/bash
#SBATCH --job-name=acd_dg_hard
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_dg_hard_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_dg_hard_%j.err

# Hardened HMC re-fit. Use ONLY if 9812192 finishes with rhat > 1.05.
# Doubles warmup, tightens adapt_delta to 0.995, deepens tree to 13.

set -euo pipefail
module purge
module load gcc/12.3.0 R/4.4.0

export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
cd $FVS_PROJECT_ROOT

export FVS_HMC_WARMUP=2500
export FVS_HMC_SAMPLING=1500
export FVS_HMC_ADAPT_DELTA=0.995
export FVS_HMC_TREEDEPTH=13
export FVS_HMC_CHAINS=4
export FVS_MAX_OBS=10000

SNAP=$FVS_PROJECT_ROOT/calibration/output/variants/acd/diameter_growth_posterior.csv.hard_pre_$(date +%Y%m%d_%H%M%S)
cp $FVS_PROJECT_ROOT/calibration/output/variants/acd/diameter_growth_posterior.csv $SNAP

Rscript calibration/R/02c_fit_dg_hmc_small.R --variant acd 2>&1 \
  | tee calibration/logs/02c_acd_hard_${SLURM_JOB_ID}.console.log

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
