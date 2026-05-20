#!/bin/bash
#SBATCH --job-name=acd_ab_chain
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_ab_chain_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/logs/acd_ab_chain_%j.err

# Runs three back-to-back benchmark passes after the longer-warmup ACD HMC
# re-fit completes. Each pass writes a tagged CSV under
# calibration/output/comparisons/manuscript_tables/
#
# Pass 1: refit_only         — converged ACD posterior, no post-pass
# Pass 2: refit_postpass_pop — converged ACD posterior + population multipliers
# Pass 3: refit_postpass_strat + NY-Adirondacks — full stack
#
# Submit with --dependency=afterok:<hmc-jobid> so it waits for the HMC job.

set -euo pipefail
module purge
module load gcc/12.3.0 R/4.4.0

export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
cd $FVS_PROJECT_ROOT

# Pre-check the new ACD posterior converged before doing anything
Rscript -e '
post <- read.csv("calibration/output/variants/acd/diameter_growth_posterior.csv")
mx <- max(post$rhat, na.rm=TRUE)
conv <- mean(post$converged, na.rm=TRUE)
cat(sprintf("ACD posterior rhat_max=%.3f converged_frac=%.3f\n", mx, conv))
if (is.na(mx) || mx >= 1.10) {
  cat("WARNING: rhat still high; benchmarks will still run but quality is suspect\n")
}
' 2>&1 | tee calibration/logs/acd_posterior_check_${SLURM_JOB_ID}.log

run_pass () {
  local tag="$1"; shift
  echo "==== PASS: $tag ===="
  # Wrap Rscript in subshell with || true so figure-step crashes do not abort
  # the whole chain. Step 5 CSVs are saved before step 6 figures run.
  ( env "$@" Rscript calibration/R/19_fia_benchmark_engine.R 2>&1 \
      | tee calibration/logs/ab_${tag}_${SLURM_JOB_ID}.log ) || \
    echo "(Rscript exited non-zero; checking for tagged CSV anyway)"

  # Save tagged CSVs. The engine writes fia_benchmark_results.csv reliably;
  # fia_benchmark_pctrmse.csv may not be produced in newer engine versions.
  for kind in pctrmse results; do
    SRC=calibration/output/comparisons/manuscript_tables/fia_benchmark_${kind}.csv
    DST=calibration/output/comparisons/manuscript_tables/fia_benchmark_${kind}_${tag}.csv
    if [ -f "$SRC" ]; then
      mv "$SRC" "$DST"
      echo "  saved: $DST"
    else
      echo "  note: $SRC not produced this pass (continuing)"
    fi
  done
}

# Pass 1: refit-only A/B
run_pass refit_only \
  FVS_ACD_RELABEL=TRUE \
  FVS_ACD_POSTPASS=FALSE

# Pass 2: refit + population post-pass
run_pass refit_postpass_pop \
  FVS_ACD_RELABEL=TRUE \
  FVS_ACD_POSTPASS=TRUE \
  FVS_ACD_POSTPASS_MODE=population

# Pass 3: refit + stratified post-pass + NY Adirondacks
run_pass refit_postpass_strat_ny \
  FVS_ACD_RELABEL=TRUE \
  FVS_ACD_POSTPASS=TRUE \
  FVS_ACD_POSTPASS_MODE=stratified \
  FVS_ACD_NY_COUNTIES="19,31,33,35,41,49,89,115"

echo "==== all three A/Bs complete ===="
ls -la calibration/output/comparisons/manuscript_tables/fia_benchmark_*_refit_*.csv
