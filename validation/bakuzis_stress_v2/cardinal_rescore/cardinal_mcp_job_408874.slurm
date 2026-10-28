#!/bin/bash
#SBATCH --job-name=bak_build_smoke
#SBATCH --account=PUOM0008
#SBATCH --time=01:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/build_smoke_%j.out
# Build FVSne.so and FVSacd.so from origin/main (895f5aa) in the worktree, then rerun the ensemble smoke
# test against that library (the ~/fvs-modern/lib build from the modified Aug 20 tree segfaulted on the
# fixed-generator keyfile).
module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main
cd $W
bash deployment/scripts/build_fvs_libraries.sh src-converted $W/lib ne acd 2>&1 | grep -v 'COMPILE FAIL' | tail -15
ls -la $W/lib/*.so
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export BAKUZIS_OUTPUT_DIR=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/ens_smoke
export PYTHONPATH=$W:$W/deployment/fvs2py:${PYTHONPATH:-}
mkdir -p $BAKUZIS_OUTPUT_DIR
python3 -X faulthandler calibration/python/bakuzis_uncertainty_comparison.py --variant ne --n-draws 0 --seed 42 --batch-id 0 --batch-size 6 --use-fia --fia-n-plots 5 --output-dir $BAKUZIS_OUTPUT_DIR
echo "rc=$?"; ls -la $BAKUZIS_OUTPUT_DIR

