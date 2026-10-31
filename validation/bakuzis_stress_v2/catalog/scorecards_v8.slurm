#!/bin/bash
#SBATCH --job-name=bak_score8 --account=PUOM0008 --time=03:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=64G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/scorecards8_%j.out
# Catalog v8 per-variant Bakuzis scorecards, same scripts and settings as scorecards_v6.slurm.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog; P=$J/proj_rescore2; F=$J/fix_b_20261030; cd $F
echo "== prep $(date)"
Rscript $C/prep_scorecard_long.R --long=catalog_long_v8.csv --seedgate=catalog_seedgate_v8.csv --out=catalog_long_scorecard_v8.csv \
  --proj=$P/b1_native/proj_b1_native_long_metric.csv,$P/b1_gompit/proj_b1_gompit_long_metric.csv,$P/b2_native/proj_b2_native_long_metric.csv,$P/b2_gompit/proj_b2_gompit_long_metric.csv
echo "== score $(date)"
mkdir -p score_scorecard_v8
Rscript $J/bakuzis_stress_v2.R --input=catalog_long_scorecard_v8.csv --source=long --out=score_scorecard_v8 --sdimax=$C/catalog_sdimax.csv --obs_ba_file=$C/catalog_obs_ba.csv --label=scorecard_20261030v8 > score_scorecard_v8/run_report.txt 2>&1
grep -E '^Reineke axis|Error|WARNING' score_scorecard_v8/run_report.txt
echo "== scorecards $(date)"
Rscript $C/bakuzis_scorecards.R --traj=score_scorecard_v8/scorecard_20261030v8_scorecard_by_trajectory.csv --gates=gates_v8/catalog_gates_by_arm.csv --thresholds=$C/fia_variant_thresholds_plot_acd.csv --out=scorecards_v8 --step=5 --qmd_rate=0.10 --label="catalog v8, sixth version (D5 plot-size fix, all arms rerun), October 30, 2026"
ls scorecards_v8 | head; echo SCORECARDS_DONE
