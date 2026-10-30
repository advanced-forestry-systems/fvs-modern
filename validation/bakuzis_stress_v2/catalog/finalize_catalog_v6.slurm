#!/bin/bash
#SBATCH --job-name=cat_final5 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/finalize_v4_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog; cd $C
echo "== concatenate $(date)"
Rscript -e 'suppressMessages(library(data.table)); fs <- list.files("runs_v5", "^catalog_[a-z]+_(native|gompit)\\.csv$", full.names = TRUE); cat("from runs_v5:", length(fs), "\n"); d <- rbindlist(lapply(fs, fread), fill = TRUE); d[, traj := paste(group, traj, sep = ":")]; cat("files", length(fs), "rows", nrow(d), "arms", paste(unique(d$arm), collapse=","), "groups", uniqueN(d$group), "\n"); fwrite(d, "catalog_long_v5.csv"); th <- fread("fia_variant_thresholds_plot.csv"); if (!"ACD" %in% th$variant) { acd <- copy(th[variant == "NE"]); acd[, variant := "ACD"]; th <- rbind(th, acd) }; fwrite(th, "fia_variant_thresholds_plot_acd.csv"); fwrite(th[, .(group = variant, ba_p99, ba_p999)], "catalog_obs_ba.csv"); fwrite(th[, .(group = variant, sdimax = round(sdi_p999_ha))], "catalog_sdimax.csv"); sgf <- list.files("runs_v5", "_seedgate\\.csv$", full.names = TRUE); sg <- rbindlist(lapply(sgf, fread), fill = TRUE); fwrite(sg, "catalog_seedgate_v5.csv"); print(sg[, .N, by = .(status)]); print(sg[status %in% c("OK","SEED_FAIL"), .(ok = mean(status == "OK"), n = .N), by = .(variant, arm, config)][order(variant, arm, config)])'
echo "== score $(date)"
mkdir -p score_catalog_v5
Rscript $J/bakuzis_stress_v2.R --input=catalog_long_v5.csv --source=long --out=score_catalog_v5 --sdimax=catalog_sdimax.csv --obs_ba_file=catalog_obs_ba.csv --label=catalog_20261029v5 > score_catalog_v5/run_report.txt 2>&1
grep -E '^Reineke axis|Error' score_catalog_v5/run_report.txt
echo "== gates $(date)"
P=$J/proj_rescore2
Rscript catalog_gates.R --long=catalog_long_v5.csv --traj=score_catalog_v5/catalog_20261029v5_scorecard_by_trajectory.csv \
  --long=$P/b1_native/proj_b1_native_long_metric.csv --traj=$P/b1_native/proj_b1_native_scorecard_by_trajectory.csv \
  --long=$P/b2_native/proj_b2_native_long_metric.csv --traj=$P/b2_native/proj_b2_native_scorecard_by_trajectory.csv \
  --long=$P/b1_gompit/proj_b1_gompit_long_metric.csv --traj=$P/b1_gompit/proj_b1_gompit_scorecard_by_trajectory.csv \
  --long=$P/b2_gompit/proj_b2_gompit_long_metric.csv --traj=$P/b2_gompit/proj_b2_gompit_scorecard_by_trajectory.csv \
  --thresholds=fia_variant_thresholds_plot_acd.csv --seedgate=catalog_seedgate_v5.csv --out=gates_v5 > gates_v5_report.txt 2>&1
tail -120 gates_v5_report.txt
echo CATALOG_FINAL_DONE
