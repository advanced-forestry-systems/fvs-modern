#!/bin/bash
#SBATCH --job-name=cat_final7 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/finalize_v6_%j.out
# Catalog v6 (fifth catalog version, October 30, 2026): runs_v5 default rows + runs_v6 calibrated rows (fixed
# overlay, PR 135), then score and gate exactly as finalize_catalog_v6.slurm did for v5.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog; cd $C
echo "== concatenate $(date)"
Rscript concat_catalog_v6.R
echo "== score $(date)"
mkdir -p score_catalog_v6
Rscript $J/bakuzis_stress_v2.R --input=catalog_long_v6.csv --source=long --out=score_catalog_v6 --sdimax=catalog_sdimax.csv --obs_ba_file=catalog_obs_ba.csv --label=catalog_20261030v6 > score_catalog_v6/run_report.txt 2>&1
grep -E '^Reineke axis|Error' score_catalog_v6/run_report.txt
echo "== gates $(date)"
P=$J/proj_rescore2
Rscript catalog_gates.R --long=catalog_long_v6.csv --traj=score_catalog_v6/catalog_20261030v6_scorecard_by_trajectory.csv \
  --long=$P/b1_native/proj_b1_native_long_metric.csv --traj=$P/b1_native/proj_b1_native_scorecard_by_trajectory.csv \
  --long=$P/b2_native/proj_b2_native_long_metric.csv --traj=$P/b2_native/proj_b2_native_scorecard_by_trajectory.csv \
  --long=$P/b1_gompit/proj_b1_gompit_long_metric.csv --traj=$P/b1_gompit/proj_b1_gompit_scorecard_by_trajectory.csv \
  --long=$P/b2_gompit/proj_b2_gompit_long_metric.csv --traj=$P/b2_gompit/proj_b2_gompit_scorecard_by_trajectory.csv \
  --thresholds=fia_variant_thresholds_plot_acd.csv --seedgate=catalog_seedgate_v6.csv --out=gates_v6 > gates_v6_report.txt 2>&1
tail -120 gates_v6_report.txt
echo CATALOG_FINAL_DONE
