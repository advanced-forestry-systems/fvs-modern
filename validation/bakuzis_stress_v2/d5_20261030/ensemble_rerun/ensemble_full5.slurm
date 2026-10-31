#!/bin/bash
#SBATCH --job-name=bak_ens_full5
#SBATCH --account=PUOM0008
#SBATCH --time=08:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=0-5
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/ens5_%A_%a.out
# Paper 1 Bakuzis ensemble, fifth regeneration (October 30, 2026): ens_full4 design, rerun from wt-fix2-20261030 with D5 (trees expanded by the condition share of their own plot size: MICRPROP/SUBPPROP/MACRPROP_UNADJ)
# (NE and ACD, 36 cells x 5 real FIA stands x (default + calibrated MAP + 100 posterior draws), seed 42,
# six batches of six cells), run from the fix tree wt-fix-20261030 (branch fix/draws-overlay-standgen-20261030):
#   D1 draw SDIMAX/BAIMULT/MORTMULT written with the MAP block's fixed-field formatter
#   D2 no MAP in-process overlay on draw runs; the MAP overlay skips species with no calibrated SDImax
#   D3 stand trees filtered to the sampled condition and expanded by 1/CONDPROP_UNADJ
#   D4 draw multipliers on the MAP serializer basis (calibration_multipliers.py), MORTMULT included
module load gcc/12.3.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix2-20261030
export FVS_PROJECT_ROOT=$W FVS_LIB_DIR=$W/lib FVS_CONFIG_DIR=$W/config FIA_DATA_DIR=$HOME/fia_data
export BAKUZIS_OUTPUT_DIR=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/ens_full5
export PYTHONPATH=$W:$W/deployment/fvs2py:${PYTHONPATH:-}
mkdir -p $BAKUZIS_OUTPUT_DIR; cd $W
for v in ne acd; do
  echo ">>> $v batch $SLURM_ARRAY_TASK_ID $(date)"
  python3 calibration/python/bakuzis_uncertainty_comparison.py --variant $v --n-draws 100 --seed 42 --batch-id $SLURM_ARRAY_TASK_ID --batch-size 6 --use-fia --fia-n-plots 5 --output-dir $BAKUZIS_OUTPUT_DIR
  echo "rc=$? $(date)"
done
echo ENS_FULL_TASK_DONE
