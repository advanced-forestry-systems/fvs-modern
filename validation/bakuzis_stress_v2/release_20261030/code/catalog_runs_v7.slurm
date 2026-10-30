#!/bin/bash
#SBATCH --job-name=cat_runs7 --account=PUOM0008 --time=10:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G
#SBATCH --array=0-1,2-3,4-5,6-7,8-9,10-11,12-13,14-15,16-17,20-21,24-25,30-31,32-33,34-35,36-37,38-39,40-41,42-43
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/cat7_%A_%a.out
# Catalog v6 (fifth catalog version), October 30, 2026: calibrated arms only, rerun with the fixed code in
# wt-fix-20261030 (fvs-modern main d08b8b96, PR 135: the in-process MAP overlay no longer writes SDImax 0 for
# species without a calibrated value). Same 18 variants x 2 survival arms, same sample and seed as runs_v5;
# default arms are unchanged by the fix (no overlay on default) and are reused from runs_v5.
# task = variant index * 2 + arm index
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
VARS=(acd bm ca ci cr cs ec em ie kt ls nc ne oc op pn so sn tt ut wc ws)
V=${VARS[$((SLURM_ARRAY_TASK_ID / 2))]}; A=$((SLURM_ARRAY_TASK_ID % 2))
cd $J
if [ $A -eq 1 ]; then
  export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$W/config/greg_mortality_coefficients.csv; ARM=gompit
else
  unset FVS_GOMPIT; ARM=native
fi
unset FVS_NOAUTOES; case $V in em|ie) export FVS_NOAUTOES=1;; esac
echo ">>> $V $ARM NOAUTOES=${FVS_NOAUTOES:-0} $(date)"
python3 $J/conus_catalog_runs.py --variant $V --arm $ARM --assign $J/fia_variant_assignment.csv --n-stands 300 --seed 20261028 --configs calibrated --out $J/runs_v6
echo "rc=$? $(date)"
