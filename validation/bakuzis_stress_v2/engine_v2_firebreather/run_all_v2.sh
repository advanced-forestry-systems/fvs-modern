#!/bin/bash
J=~/jobs/bakuzis_fvsconus_20261027; cd $J
rm -rf engine_v2; bash build/run_engine_v2.sh > build/run_engine_v2.log 2>&1
for t in gompit native; do
  rm -rf score_engine_v2_$t; mkdir -p score_engine_v2_$t
  Rscript bakuzis_stress_v2.R --input=engine_v2/engine_bakuzis_100yr_$t.csv --source=engine --out=score_engine_v2_$t --sdimax=sdimax_ne_engine_tph.csv --obs_ba_p99=70.57 --label=engine_v2_$t > score_engine_v2_$t/run_report.txt 2>&1
done
echo ALL_DONE
