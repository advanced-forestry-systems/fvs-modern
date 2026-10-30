#!/bin/bash
#SBATCH --job-name=ro_full_eng --account=PUOM0008 --time=06:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --array=0-71
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/logs/eng_%A_%a.out
# Run-off full run, engines. 72 tasks (full/tasks_engines.tsv): A = fvs-modern near-term (18 variants x
# native/gompit, default+calibrated configs per task); C = unmodified USDA FVS b183c10f long-term and
# near-term (18 variants, default parameters, native mortality; ACD stands run on the upstream NE library).
# Every catalog stand (--n-stands 300 --seed 20261028, no --first). NOT SUBMITTED by the pilot agent.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main; RO=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff; F=$RO/full
export FVS_PROJECT_ROOT=$W FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export PYTHONPATH=$W:$W/deployment/fvs2py:$W/calibration/python:${PYTHONPATH:-}
read I V C M ARM EV <<< "$(awk -v t=$SLURM_ARRAY_TASK_ID '$1==t' $F/tasks_engines.tsv)"
unset FVS_GOMPIT FVS_GOMPIT_COEF FVS_NOAUTOES
if [ "$C" = A ]; then export FVS_LIB_DIR=$W/lib; CFG=default,calibrated; OUT=$F/a_$M
else export FVS_LIB_DIR=$RO/upstream_build/lib; CFG=usfs; OUT=$F/c_$M; fi
[ "$ARM" = gompit ] && export FVS_GOMPIT=1 FVS_GOMPIT_COEF=$W/config/greg_mortality_coefficients.csv
case $V in em|ie) export FVS_NOAUTOES=1;; esac
echo ">>> task $I $V contestant=$C mode=$M arm=$ARM engine_variant=$EV lib=$FVS_LIB_DIR configs=$CFG $(date)"
python3 $RO/scripts/runoff_engine_runs.py --variant $V --engine-variant $EV --arm $ARM --mode $M --configs $CFG --out $OUT
echo "rc=$? $(date)"
