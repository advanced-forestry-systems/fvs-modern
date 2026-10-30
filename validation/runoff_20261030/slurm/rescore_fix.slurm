#!/bin/bash
#SBATCH --job-name=ro_rescore_fix --account=PUOM0008 --time=04:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=96G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/logs/rescore_fix_%j.out
# Amendment A3 rescoring: primary = A calibrated near from engines_fix + A long from catalog_long_v6; sensitivity adds _inproc arms.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog
echo "== engines_fix status"; for f in $F/engines_fix/a_near/*_gate.csv; do echo "$(basename $f): $(tail -n +2 $f | awk -F, '{print $6}' | sort | uniq -c | tr '\n' ' ')"; done
Rscript $F/prep_fix.R
NEAROTHER="$(ls $F/c_near/runoff_near_*.csv | grep -v _gate) $(ls $F/b_out/v*/near/runoff_near_*.csv)"
for S in score_fix score_fix_sens; do
  O=$F/$S; mkdir -p $O/near $O/eng $O/gates $O/decision
  EXTRA=""; [ $S = score_fix_sens ] && EXTRA=$O/inputs/a_near_calibrated_inproc.csv
  Rscript $RO/scripts/runoff_score_near.R $RO/obs/nearterm_obs_frame.csv $O/near $F/score_fix/inputs/a_near_default.csv $F/score_fix/inputs/a_near_calibrated_fix.csv $EXTRA $NEAROTHER > $O/near/report.txt 2>&1; tail -2 $O/near/report.txt
  Rscript $J/bakuzis_stress_v2.R --input=$O/eng_long.csv --source=long --out=$O/eng --sdimax=$C/catalog_sdimax.csv --obs_ba_file=$C/catalog_obs_ba.csv --label=runoff_eng > $O/eng/run_report.txt 2>&1; grep -E 'Reineke axis|Error' $O/eng/run_report.txt
  Rscript $C/catalog_gates.R --long=$O/eng/runoff_eng_long_metric.csv --traj=$O/eng/runoff_eng_scorecard_by_trajectory.csv \
    --long=$F/score/b/runoff_b_long_metric.csv --traj=$F/score/b/runoff_b_scorecard_by_trajectory.csv \
    --thresholds=$C/fia_variant_thresholds_plot_acd.csv --out=$O/gates > $O/gates/report.txt 2>&1; grep -c . $O/gates/catalog_gates_by_arm.csv
  [ $S = score_fix ] && Rscript $RO/scripts/runoff_decide.R $O/gates/catalog_gates_by_arm.csv $O/near/near_joined.csv $O/decision > $O/decision/report.txt 2>&1
done
echo "== fix effect on A calibrated arms"
Rscript -e '
suppressMessages(library(data.table)); F <- "'$F'"
n <- fread(file.path(F, "score_fix_sens/near/near_scores.csv"))[set == "common" & metric == "BA" & grepl("^calibrated", arm)]
n[, base := sub("_inproc$", "", arm)][, path := ifelse(grepl("_inproc$", arm), "inproc", "fixed")]
w <- dcast(n, variant + base ~ path, value.var = c("pbias", "pbias_lo", "pbias_hi", "n"))
g <- fread(file.path(F, "score_fix_sens/gates/catalog_gates_by_arm.csv"))[grepl("^calibrated", arm)]
g[, base := sub("_inproc$", "", arm)][, path := ifelse(grepl("_inproc$", arm), "inproc", "fixed")]
gw <- dcast(g, group + base ~ path, value.var = c("share_ba999_50", "share_sdi_50", "reineke_flag_rate", "med_peak_ba100", "admissible"))
setnames(gw, "group", "variant"); out <- merge(w, gw, by = c("variant", "base"), all = TRUE)
fwrite(out, file.path(F, "score_fix/fix_effect_calibrated.csv")); print(out[, .(variant, base, pb_fixed = round(pbias_fixed, 1), pb_inproc = round(pbias_inproc, 1),
  ba999_50_fixed = round(share_ba999_50_fixed, 3), ba999_50_inproc = round(share_ba999_50_inproc, 3), adm_fixed = admissible_fixed, adm_inproc = admissible_inproc)], nrows = 100)'
md5sum $F/score_fix/decision/runoff_decision_by_variant.csv $F/score_fix/gates/catalog_gates_by_arm.csv $F/score_fix/near/near_scores.csv $F/score_fix/fix_effect_calibrated.csv $F/score_fix_sens/near/near_scores.csv $F/score_fix_sens/gates/catalog_gates_by_arm.csv $F/score/decision/runoff_decision_by_variant.csv
echo RESCORE_FIX_DONE
