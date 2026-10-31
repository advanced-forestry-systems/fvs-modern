#!/bin/bash
#SBATCH --job-name=rev143b --account=PUOM0008 --time=00:45:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/eng_20261030/rev2_%j.out
# Second review round of PR 143 (October 30, 2026), run in the worktree holding the PR head.
# Writes into validation/fvs2py_fixes_20261030/: test_engine_fixes_run.txt, concurrency_20_repeats.txt
# and smoke_comparison.csv (regenerated from the CSVs in that folder).
module load gcc/12.3.0
export PATH=$HOME/.local/bin:$PATH
W=${W:-/fs/scratch/PUOM0008/crsfaaron/wt-rev-20261030}
V=$W/validation/fvs2py_fixes_20261030
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python
cd $W
T=deployment/fvs2py/fvs2py/tests/test_engine_fixes.py
{ echo "# $T, $(date '+%B %d, %Y %H:%M'), $(python3 -c 'import platform;print(*platform.libc_ver())')"
  python3 -m pytest $T -v -p no:cacheprovider 2>&1 | grep -E "PASSED|FAILED|ERROR|passed|failed"; } > $V/test_engine_fixes_run.txt
{ echo "# test_concurrent_construction_registers_one_owner, 20 separate pytest runs, $(date '+%B %d, %Y %H:%M')"
  for i in $(seq 1 20); do echo "run $i: $(python3 -m pytest $T -q -p no:cacheprovider -k test_concurrent_construction_registers_one_owner 2>&1 | tail -n 1)"; done; } > $V/concurrency_20_repeats.txt
cd $V && python3 compare_pairs.py smoke_comparison.csv \
  before_fix_vs_pr143_first_commit smoke_before.csv smoke_after.csv \
  first_commit_vs_review_fix smoke_after.csv smoke_after_review.csv
cat $V/test_engine_fixes_run.txt $V/concurrency_20_repeats.txt
echo REV2_DONE
