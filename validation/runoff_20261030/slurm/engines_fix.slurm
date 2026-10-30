#!/bin/bash
#SBATCH --job-name=ro_fix_eng --account=PUOM0008 --time=06:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=0-35
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/engines_fix/logs/eng_%A_%a.out
# Amendment A3: A calibrated near-term arms (calibrated_native, calibrated_gompit) rerun with the fixed fvs-modern
# worktree wt-fix-20261030 (fvs-modern main d08b8b96: D2 zero fill and three other defects fixed).
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030; RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full/engines_fix
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python
VARS=(acd bm ca ci cr cs ec em ie ls ne pn sn so tt ut wc ws); ARMS=(native gompit)
V=${VARS[$((SLURM_ARRAY_TASK_ID / 2))]}; ARM=${ARMS[$((SLURM_ARRAY_TASK_ID % 2))]}
unset FVS_GOMPIT FVS_GOMPIT_COEF FVS_NOAUTOES
[ $ARM = gompit ] && export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$W/config/greg_mortality_coefficients.csv
case $V in em|ie) export FVS_NOAUTOES=1;; esac
echo ">>> fix $V $ARM lib=$FVS_LIB_DIR $(date)"
python3 -c "import bakuzis_uncertainty_comparison as b; print('buc from', b.__file__, 'lib', b.FVS_LIB_DIR)"
python3 $RO/scripts/runoff_engine_runs.py --variant $V --arm $ARM --mode near --configs calibrated --out $F/a_near 2>&1 | grep -E 'status counts|Traceback|Error|RUNOFF_TASK_DONE'
echo "rc=$? $(date)"

