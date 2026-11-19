#!/bin/bash
#SBATCH --job-name=cat_runs8 --account=PUOM0008 --time=10:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=24G
#SBATCH --array=0-35
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/cat8_%A_%a.out
# Catalog v8 (October 30, 2026): all four arms (default and calibrated x native and GOMPIT) on the catalog sample
# (--n-stands 300 --seed 20261028), from wt-fix2-20261030 with D5: each tree expanded by the condition share of
# its own plot size (MICRPROP/SUBPPROP/MACRPROP_UNADJ) instead of CONDPROP_UNADJ (266 of 4,748 stands move > 2%).
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix2-20261030; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; F=$J/fix_b_20261030
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
VARS=(acd bm ca ci cr cs ec em ie ls ne pn sn so tt ut wc ws)
V=${VARS[$((SLURM_ARRAY_TASK_ID / 2))]}; A=$((SLURM_ARRAY_TASK_ID % 2))
if [ $A -eq 1 ]; then export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$W/config/greg_mortality_coefficients.csv; ARM=gompit; else unset FVS_GOMPIT; ARM=native; fi
case $V in em|ie) export FVS_NOAUTOES=1;; *) unset FVS_NOAUTOES;; esac
cd $F
echo ">>> $V $ARM NOAUTOES=${FVS_NOAUTOES:-0} $(date)"
python3 $F/conus_catalog_runs.py --variant $V --arm $ARM --assign $J/catalog/fia_variant_assignment.csv --n-stands 300 --seed 20261028 --out $F/runs_v8
echo "rc=$? $(date)"
