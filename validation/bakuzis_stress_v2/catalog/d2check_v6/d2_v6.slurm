#!/bin/bash
#SBATCH --job-name=d2v6 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=0-11
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/d2check_v6/logs/d2_%A_%a.out
# D2 recheck on the fixed code: calibrated_native, first 20 catalog stands, the 6 variants where the
# zero-fill was material (CA CR CS EM LS WC). inproc vs kwonly should now be identical.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030; C=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog; D=$C/d2check_v6
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
VARS=(ca cr cs em ls wc); MODES=(inproc kwonly)
V=${VARS[$((SLURM_ARRAY_TASK_ID / 2))]}; M=${MODES[$((SLURM_ARRAY_TASK_ID % 2))]}
unset FVS_GOMPIT FVS_NOAUTOES; case $V in em|ie) export FVS_NOAUTOES=1;; esac
export D2_MODE=$M D2_SPSDI_LOG=$D/$M/spsdi_${V}.csv; mkdir -p $D/$M; rm -f $D2_SPSDI_LOG
cd $C
echo ">>> D2 $V $M $(date)"
python3 $D/d2_wrapper_v6.py --variant $V --arm native --assign $C/fia_variant_assignment.csv --n-stands 300 --seed 20261028 --configs calibrated --first 20 --out $D/$M 2>&1 | grep -E 'wrote|Traceback|Error|CATALOG_TASK_DONE|spsdi read failed'
echo "rc=$? $(date)"
