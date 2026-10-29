#!/bin/bash
#SBATCH --job-name=cat_runs --account=PUOM0008 --time=10:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G
#SBATCH --array=0,1,14,15,16,17
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/logs/cat_%A_%a.out
# 22 variants x 2 survival arms; task = variant index * 2 + arm index
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog
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
case $V in em|ie) export FVS_NOAUTOES=1;; esac
echo ">>> $V $ARM NOAUTOES=${FVS_NOAUTOES:-0} $(date)"
python3 $J/conus_catalog_runs.py --variant $V --arm $ARM --assign $J/fia_variant_assignment.csv --n-stands 300 --seed 20261028 --out $J/runs_v4
echo "rc=$? $(date)"
