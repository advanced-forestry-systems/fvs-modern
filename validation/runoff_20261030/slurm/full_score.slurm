#!/bin/bash
#SBATCH --job-name=ro_full_score --account=PUOM0008 --time=04:00:00 --nodes=1 --ntasks=1 --cpus-per-task=4 --mem=96G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/logs/score_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; C=$J/catalog
S=$F/score; mkdir -p $S/eng $S/b $S/gates $S/near $S/decision
G=ACD,BM,CA,CI,CR,CS,EC,EM,IE,LS,NE,PN,SN,SO,TT,UT,WC,WS
echo "== run status"; for f in $F/a_near/*_gate.csv $F/c_long/*_gate.csv $F/c_near/*_gate.csv; do echo "$(basename $f): $(tail -n +2 $f | awk -F, '{print $6}' | sort | uniq -c | tr '\n' ' ')"; done
Rscript $RO/scripts/runoff_score_near.R $RO/obs/nearterm_obs_frame.csv $S/near $F/a_near/runoff_near_*.csv $F/c_near/runoff_near_*.csv $F/b_out/v*/near/runoff_near_*.csv
Rscript $RO/scripts/runoff_long_prep.R $S/eng_long.csv $G 0 $F/c_long
Rscript $J/bakuzis_stress_v2.R --input=$S/eng_long.csv --source=long --out=$S/eng --sdimax=$C/catalog_sdimax.csv --obs_ba_file=$C/catalog_obs_ba.csv --label=runoff_eng > $S/eng/run_report.txt 2>&1
BF=$(ls $F/b_out/v*/long/runoff_long_b_*_metrics.csv | tr '\n' ',' | sed 's/,$//')
Rscript $J/bakuzis_stress_v2.R --input=$BF --source=projector --out=$S/b --sdimax=$C/catalog_sdimax.csv --obs_ba_file=$C/catalog_obs_ba.csv --label=runoff_b > $S/b/run_report.txt 2>&1
Rscript $C/catalog_gates.R --long=$S/eng/runoff_eng_long_metric.csv --traj=$S/eng/runoff_eng_scorecard_by_trajectory.csv \
  --long=$S/b/runoff_b_long_metric.csv --traj=$S/b/runoff_b_scorecard_by_trajectory.csv \
  --thresholds=$C/fia_variant_thresholds_plot_acd.csv --out=$S/gates > $S/gates/report.txt 2>&1
Rscript $RO/scripts/runoff_decide.R $S/gates/catalog_gates_by_arm.csv $S/near/near_joined.csv $S/decision
echo FULL_SCORE_DONE
