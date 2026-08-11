#!/bin/bash
# run_fvs_perseus_pipeline.sh -- end-to-end FVS -> PERSEUS pipeline, one place.
#
# Folds the post-hoc fixes (treeinit TPA expansion + height completion) into the
# front of the pipeline so a campaign always starts from correct tree lists, and
# chains through to the dashboard merge + stress test. Each STAGE is idempotent
# and can be run alone (comment the others). Heavy stages submit SLURM arrays;
# wait for them before the dependent stage.
#
# Env: FVS_PROJECT_ROOT=$HOME/fvs-modern, FVS_LIB_DIR, FVS_GOMPIT lib for gompit.
set -euo pipefail
SD=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
FIA=/fs/scratch/PUOM0008/crsfaaron/FIA
FRESH=/fs/scratch/PUOM0008/crsfaaron/FIA_fresh
PW=/fs/scratch/PUOM0008/crsfaaron/perseus_wire
module load gcc/12.3.0 python/3.12 2>/dev/null || true

echo "== STAGE 1: treeinit fix (TPA_UNADJ expansion + height completion) =="
# The DataMart FVS_TREEINIT_PLOT under-expands overstory ~6.5x and is ~30% short
# on heights; this restores both from the raw FIA TREE table. ALWAYS run first.
python3 $SD/treeinit_fix_v2.py --in-dir $FRESH/treeinit --tree-dir $FIA \
    --out-dir $FRESH/treeinit_h

echo "== STAGE 2: strata inputs (plantation flag, harvest rates) =="
python3 $SD/build_plantation_flag.py                       # plt_plantation.csv
# state_harvest_rates.csv needs gdal/R (conus_hcs raster sampling):
# module load gdal/3.7.3 R/4.4.0; Rscript $SD/build_state_harvest_rates.R

echo "== STAGE 3: campaign (SLURM arrays vs treeinit_h) =="
# sbatch --array=0-380%40 $SD/submit_conus_fvs_v3.slurm     # default+calibrated
# sbatch --array=0-380%40 $SD/submit_conus_gompit_v3.slurm  # gompit
echo "   (submit the two arrays, wait for completion, then continue)"

echo "== STAGE 4: aggregate per-state density series =="
for c in default calibrated; do
  python3 $SD/fvs_perseus_aggregate.py --in-dir $SD/out_fvs_v3 \
      --out-dir $SD/perseus_series_${c}_v3 --config $c --start 2025 --engine fvs_$c
done
python3 $SD/fvs_perseus_aggregate.py --in-dir $SD/out_gompit_v3 \
    --out-dir $SD/perseus_series_gompit_v3 --config gompit --start 2025 --engine fvs_gompit

echo "== STAGE 5: plantation-aware managed scenarios =="
for c in default calibrated; do
  python3 $SD/fvs_managed_v2.py --campaign $SD/out_fvs_v3 --config $c \
      --plantation $SD/plt_plantation.csv --rates $SD/state_harvest_rates.csv \
      --start 2030 --k 1.9 --window 20 --out $SD/managed3_$c
done
python3 $SD/fvs_managed_v2.py --campaign $SD/out_gompit_v3 --config gompit \
    --plantation $SD/plt_plantation.csv --rates $SD/state_harvest_rates.csv \
    --start 2030 --k 1.9 --window 20 --out $SD/managed3_gompit

echo "== STAGE 6: posterior parameter CI (anchored states, SLURM array) =="
# sbatch --array=0-179%45 $SD/submit_post.slurm  (manifest_post.tsv) -> aggregate
echo "   (optional; produces posterior_ci_all.csv)"

echo "== STAGE 7: merge onto a CLEAN main + ribbon + stress test =="
# Run on a fresh checkout of origin/main so unreleased feature work never leaks:
#   git checkout main && git reset --hard origin/main
#   FVS_MANAGED_ROOT=<root> python3 fvs_perseus_merge.py $PW <series_root>
#   python3 fvs_posterior_ribbon.py $PW posterior_ci_all.csv
#   python3 fvs_dashboard_stress.py $PW/public/api   # MUST report 0 violations
#   git add public/api && git commit && git push origin main
echo "== DONE =="
