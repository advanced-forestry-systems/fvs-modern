#!/bin/bash
#SBATCH --job-name=ro_d2 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=0-23
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/d2check/logs/d2_%A_%a.out
# D2 audit: calibrated_native, 100 yr, first 20 catalog stands, 12 variants where the in-process fvs2py
# overlay sets spsdi (besides NE). Mode "inproc" = as the catalog ran it (FVS config_version calibrated:
# keyword overlay + apply_calibrated_config). Mode "kwonly" = same keyword block, FVS config_version default.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main; RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; D=$RO/d2check
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
VARS=(ca ci cr cs ec em ls sn tt ut wc ws); MODES=(inproc kwonly)
V=${VARS[$((SLURM_ARRAY_TASK_ID / 2))]}; M=${MODES[$((SLURM_ARRAY_TASK_ID % 2))]}
unset FVS_GOMPIT FVS_NOAUTOES; case $V in em|ie) export FVS_NOAUTOES=1;; esac
export D2_MODE=$M D2_SPSDI_LOG=$D/$M/spsdi_${V}.csv; mkdir -p $D/$M; rm -f $D2_SPSDI_LOG
echo ">>> D2 $V $M $(date)"
python3 $D/d2_wrapper.py --variant $V --arm native --mode long --configs calibrated --first 20 --out $D/$M 2>&1 | grep -E 'status counts|Traceback|Error|RUNOFF_TASK_DONE|spsdi read failed'
echo "rc=$? $(date)"

