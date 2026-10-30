#!/bin/bash
#SBATCH --job-name=defv6 --account=PUOM0008 --time=01:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=0-1
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/d2check_v6/logs/def_%A_%a.out
# Default-arm reuse check: rerun default_native with the fixed code on the first 10 catalog stands of CA and WC
# and compare with the runs_v5 default rows.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030; C=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog; D=$C/d2check_v6
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
VARS=(ca wc); V=${VARS[$SLURM_ARRAY_TASK_ID]}
unset FVS_GOMPIT FVS_NOAUTOES
cd $C
echo ">>> default check $V $(date)"
python3 $C/conus_catalog_runs.py --variant $V --arm native --assign $C/fia_variant_assignment.csv --n-stands 300 --seed 20261028 --configs default --first 10 --out $D/default 2>&1 | grep -E 'wrote|Traceback|Error|CATALOG_TASK_DONE'
echo "rc=$? $(date)"
