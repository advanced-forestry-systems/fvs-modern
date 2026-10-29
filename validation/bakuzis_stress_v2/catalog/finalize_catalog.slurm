#!/bin/bash
#SBATCH --job-name=cat_final --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/finalize_v4_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog; cd $C
echo "== concatenate $(date)"
Rscript -e 'suppressMessages(library(data.table)); fs <- list.files("runs", "^catalog_[a-z]+_(native|gompit)\\.csv$", full.names = TRUE); fs <- sapply(fs, function(f) if (file.exists(file.path("runs_v4", basename(f)))) file.path("runs_v4", basename(f)) else if (file.exists(file.path("runs_v3", basename(f)))) file.path("runs_v3", basename(f)) else f); cat("from runs_v4:", sum(grepl("runs_v4", fs)), "\n"); d <- rbindlist(lapply(fs, fread), fill = TRUE); d[, traj := paste(group, traj, sep = ":")]; cat("files", length(fs), "rows", nrow(d), "arms", paste(unique(d$arm), collapse=","), "groups", uniqueN(d$group), "\n"); fwrite(d, "catalog_long_v4.csv"); th <- fread("fia_variant_thresholds_plot.csv"); if (!"ACD" %in% th$variant) { acd <- copy(th[variant == "NE"]); acd[, variant := "ACD"]; th <- rbind(th, acd) }; fwrite(th, "fia_variant_thresholds_plot_acd.csv"); fwrite(th[, .(group = variant, ba_p99, ba_p999)], "catalog_obs_ba.csv"); fwrite(th[, .(group = variant, sdimax = round(sdi_p999_ha))], "catalog_sdimax.csv"); sgf <- list.files("runs", "_seedgate\\.csv$", full.names = TRUE); sgf <- sapply(sgf, function(f) if (file.exists(file.path("runs_v4", basename(f)))) file.path("runs_v4", basename(f)) else if (file.exists(file.path("runs_v3", basename(f)))) file.path("runs_v3", basename(f)) else f); sg <- rbindlist(lapply(sgf, fread), fill = TRUE); fwrite(sg, "catalog_seedgate_v4.csv"); print(sg[, .N, by = .(status)]); print(sg[status %in% c("OK","SEED_FAIL"), .(ok = mean(status == "OK"), n = .N), by = .(variant, arm, config)][order(variant, arm, config)])'
echo "== score $(date)"
mkdir -p score_catalog_v4
Rscript $J/bakuzis_stress_v2.R --input=catalog_long_v4.csv --source=long --out=score_catalog_v4 --sdimax=catalog_sdimax.csv --obs_ba_file=catalog_obs_ba.csv --label=catalog_20261028 > score_catalog_v4/run_report.txt 2>&1
grep -E '^Reineke axis|Error' score_catalog_v4/run_report.txt
echo "== gates $(date)"
P=$J/proj_rescore2
Rscript catalog_gates.R --long=catalog_long_v4.csv --traj=score_catalog_v4/catalog_20261028_scorecard_by_trajectory.csv \
  --long=$P/b1_native/proj_b1_native_long_metric.csv --traj=$P/b1_native/proj_b1_native_scorecard_by_trajectory.csv \
  --long=$P/b2_native/proj_b2_native_long_metric.csv --traj=$P/b2_native/proj_b2_native_scorecard_by_trajectory.csv \
  --long=$P/b1_gompit/proj_b1_gompit_long_metric.csv --traj=$P/b1_gompit/proj_b1_gompit_scorecard_by_trajectory.csv \
  --long=$P/b2_gompit/proj_b2_gompit_long_metric.csv --traj=$P/b2_gompit/proj_b2_gompit_scorecard_by_trajectory.csv \
  --thresholds=fia_variant_thresholds_plot_acd.csv --seedgate=catalog_seedgate_v4.csv --out=gates_v4 > gates_v4_report.txt 2>&1
tail -120 gates_v4_report.txt
echo CATALOG_FINAL_DONE
