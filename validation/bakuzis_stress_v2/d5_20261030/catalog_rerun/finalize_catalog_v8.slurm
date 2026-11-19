#!/bin/bash
#SBATCH --job-name=finalize8 --account=PUOM0008 --time=03:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/finalize8_%j.out
# Catalog v8: concatenate, score and gate exactly as finalize_catalog_v7.slurm did for v6, then compare with v6.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog; F=$J/fix_b_20261030; cd $F
Rscript concat_catalog_v8.R
mkdir -p score_catalog_v8
Rscript $J/bakuzis_stress_v2.R --input=catalog_long_v8.csv --source=long --out=score_catalog_v8 --sdimax=$C/catalog_sdimax.csv --obs_ba_file=$C/catalog_obs_ba.csv --label=catalog_20261030v8 > score_catalog_v8/run_report.txt 2>&1
grep -E '^Reineke axis|Error' score_catalog_v8/run_report.txt
P=$J/proj_rescore2
Rscript $C/catalog_gates.R --long=catalog_long_v8.csv --traj=score_catalog_v8/catalog_20261030v8_scorecard_by_trajectory.csv \
  --long=$P/b1_native/proj_b1_native_long_metric.csv --traj=$P/b1_native/proj_b1_native_scorecard_by_trajectory.csv \
  --long=$P/b2_native/proj_b2_native_long_metric.csv --traj=$P/b2_native/proj_b2_native_scorecard_by_trajectory.csv \
  --long=$P/b1_gompit/proj_b1_gompit_long_metric.csv --traj=$P/b1_gompit/proj_b1_gompit_scorecard_by_trajectory.csv \
  --long=$P/b2_gompit/proj_b2_gompit_long_metric.csv --traj=$P/b2_gompit/proj_b2_gompit_scorecard_by_trajectory.csv \
  --thresholds=$C/fia_variant_thresholds_plot_acd.csv --seedgate=catalog_seedgate_v8.csv --out=gates_v8 > gates_v8_report.txt 2>&1
tail -40 gates_v8_report.txt
Rscript compare_gates_v6_v8.R
echo CATALOG_V8_DONE
