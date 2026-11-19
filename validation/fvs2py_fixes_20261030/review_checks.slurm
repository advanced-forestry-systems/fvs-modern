#!/bin/bash
#SBATCH --job-name=rev143 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/eng_20261030/rev_%j.out
module load gcc/12.3.0
export PATH=$HOME/.local/bin:$PATH
E=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/eng_20261030
R=$E/rev; mkdir -p $R
export FIA_DATA_DIR=$HOME/fia_data
OLD=/fs/scratch/PUOM0008/crsfaaron/wt-eng-20261030
NEW=/fs/scratch/PUOM0008/crsfaaron/wt-rev-20261030
setenv() { W=$1; export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python; }
python3 -c "import platform;print('glibc', platform.libc_ver())"
for T in old new; do
  if [ $T = old ]; then W=$OLD; else W=$NEW; fi
  setenv $W
  echo "== $T smoke_summary (buc.run_scenario) $(date)"; (cd $R && python3 $E/smoke_summary.py $R/smoke_$T.csv 2>&1 | tail -n 2)
  echo "== $T smoke_both (buc and b100 run_scenario) $(date)"; (cd $R && python3 $R/smoke_both.py $R/smokeboth_$T.csv 2>&1 | tail -n 2)
  echo "== $T engine_bakuzis_v2 $(date)"; (cd $R && python3 $R/engine_bakuzis_v2.py --wt $W --lib $W/lib/FVSne.so --tag ne --out $R/eb_$T 2>&1 | tail -n 3)
  echo "== $T confirm_real d12 $(date)"; sed -e "s#^F = .*#F = \"$W\"#" -e "s#^X = .*#X = \"$R/cr_$T\"#" $R/confirm_real.py > $R/confirm_real_$T.py
  (cd $R && python3 $R/confirm_real_$T.py d12 2>&1 | tail -n 4)
done
setenv $NEW; cd $NEW
echo "== new tests $(date)"; python3 -m pytest deployment/fvs2py/fvs2py/tests/test_engine_fixes.py -v -p no:cacheprovider 2>&1 | grep -E "PASSED|FAILED|ERROR|passed|failed"
echo "== new tests, 5 repeats $(date)"; for i in 1 2 3 4 5; do python3 -m pytest deployment/fvs2py/fvs2py/tests/test_engine_fixes.py -q -p no:cacheprovider 2>&1 | tail -n 1; done
echo "== config tests $(date)"; python3 -m pytest config/tests/test_draws_overlay_standgen.py -q -p no:cacheprovider 2>&1 | tail -n 2
cd $NEW/deployment/fvs2py/fvs2py/tests && sed "s#/usr/local/lib/FVSso.so#$NEW/lib/FVSso.so#" test_base.py > tmp_lib_test_base.py && cd $NEW
echo "== existing suite (TEST_DLL=lib/FVSso.so) $(date)"; PYTHONPATH=$NEW:$NEW/deployment/fvs2py timeout 900 python3 -m pytest deployment/fvs2py/fvs2py/tests/tmp_lib_test_base.py deployment/fvs2py/fvs2py/tests/test_core.py -q -p no:cacheprovider 2>&1 | grep -E "^FAILED" | sed 's/ - .*//' | sort > $R/fails_rev.txt; wc -l < $R/fails_rev.txt
PYTHONPATH=$NEW:$NEW/deployment/fvs2py timeout 900 python3 -m pytest deployment/fvs2py/fvs2py/tests/tmp_lib_test_base.py deployment/fvs2py/fvs2py/tests/test_core.py -q -p no:cacheprovider 2>&1 | tail -n 1
rm -f deployment/fvs2py/fvs2py/tests/tmp_lib_test_base.py
echo "fails identical to PR 143 head run:"; diff <(sort $E/fails_eng.txt) $R/fails_rev.txt && echo SAME
echo "== comparisons $(date)"
python3 $R/compare_pairs.py $R/rev_comparisons.csv \
  smoke_pr143head_vs_review $R/smoke_old.csv $R/smoke_new.csv \
  smoke_asrecorded_after_vs_review $E/smoke_after.csv $R/smoke_new.csv \
  smokeboth_pr143head_vs_review $R/smokeboth_old.csv $R/smokeboth_new.csv \
  engine_bakuzis_v2_pr143head_vs_review $R/eb_old/engine_bakuzis_100yr_ne.csv $R/eb_new/engine_bakuzis_100yr_ne.csv \
  confirm_real_d12_runs $R/cr_old/confirm/d12_runs.csv $R/cr_new/confirm/d12_runs.csv \
  confirm_real_d1_echo $R/cr_old/confirm/d1_echo_decode.csv $R/cr_new/confirm/d1_echo_decode.csv 2>&1
python3 $R/compare_pairs.py $R/smoke_comparison.csv \
  before_fix_vs_pr143_head $E/smoke_before.csv $E/smoke_after.csv \
  pr143_head_vs_review_fix $E/smoke_after.csv $R/smoke_new.csv 2>&1
echo REV_DONE
