#!/bin/bash
# Run after SLURM job 9912896 (or any equivalent A/B chain) completes.
# Snapshots the tagged CSVs to a versioned artifact directory, runs the
# comparison reporter, and commits everything to the branch.

set -euo pipefail
PROJ=/users/PUOM0008/crsfaaron/fvs-modern-acdbridge
cd $PROJ

# 1. Check the chain completed
LATEST_JOB=${1:-9912896}
STATE=$(sacct -j $LATEST_JOB --format=State -P -n 2>/dev/null | head -1)
echo "Job $LATEST_JOB state: $STATE"

# 2. Locate tagged CSVs
TAGS=(refit_only refit_postpass_pop refit_postpass_strat_ny)
TABLES=$PROJ/calibration/output/comparisons/manuscript_tables

echo "===tagged CSVs==="
for t in "${TAGS[@]}"; do
  for kind in pctrmse results; do
    f=$TABLES/fia_benchmark_${kind}_${t}.csv
    if [ -f $f ]; then
      ls -la $f
    else
      echo "MISSING: $f"
    fi
  done
done

# 3. Run the comparison reporter
echo ""
echo "===running compare_post_refit_ab.R==="
module load gcc/12.3.0 R/4.4.0
Rscript $PROJ/calibration/R/compare_post_refit_ab.R 2>&1 | tee $PROJ/calibration/logs/compare_ab.log

# 4. Snapshot to v2 artifact directory
ART=$PROJ/calibration/analysis/acd_stand_level_2026-05-16/calibrated_ne_vs_acd_v2
mkdir -p $ART
for t in "${TAGS[@]}"; do
  for kind in pctrmse results; do
    f=$TABLES/fia_benchmark_${kind}_${t}.csv
    [ -f $f ] && cp $f $ART/
  done
done
cp $PROJ/calibration/analysis/acd_stand_level_2026-05-16/post_refit_comparison/comparison.md $ART/ 2>/dev/null

# 5. Write a README
cat > $ART/README.md <<MDEOF
# Calibrated NE vs ACD A/B v2 (post-HMC re-fit) — 2026-05-17

## Source

SLURM job $LATEST_JOB, A/B chain after HMC re-fit 9812192. With the
HMC sigma_b0 dropping 3.5x but the z_b0 silent-NA bug not yet
diagnosed, this run uses ACD's NE-fallback path (the canonical
round-5/6 path).

## Files

- fia_benchmark_pctrmse_refit_only.csv : ACD via NE-fallback, no post-pass
- fia_benchmark_pctrmse_refit_postpass_pop.csv : + population multipliers
- fia_benchmark_pctrmse_refit_postpass_strat_ny.csv : + stratified + NY counties
- fia_benchmark_results_*.csv : full-metric versions
- comparison.md : side-by-side compare across all three configurations
                  + baseline (28.52% ACD vs 23.19% NE pre-refit)

## Interpretation

The current ACD posterior is the round-5/6 baseline. The HMC refit
posterior (round 10) exists in
calibration/output/variants/acd/diameter_growth_samples.zb0_refit.rds
but is not yet usable due to the z_b0 silent-NA bug (task #104).

Next round: diagnose z_b0 NA, switch ACD to use the refit posterior,
and rerun this chain for the converged-A/B numbers.
MDEOF

# 6. Commit and push
cd $PROJ
git add $ART/
git commit -m "Round 13 close: calibrated NE vs ACD A/B v2 results

A/B chain $LATEST_JOB results harvested:
- 3 tagged CSVs (refit_only, postpass_pop, postpass_strat_ny)
- comparison.md side-by-side
- Snapshotted to calibrated_ne_vs_acd_v2/

ACD uses NE-fallback path (the round-5/6 baseline) due to z_b0
silent-NA bug in the new HMC refit posterior (tracked as task #104).
Refit posterior preserved as samples.zb0_refit.rds for future use."

git push origin acd-bridge-fix-2026-05-15

echo "DONE. Branch is now ready for PR review."
