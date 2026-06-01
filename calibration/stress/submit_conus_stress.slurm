#!/bin/bash
#SBATCH --job-name=fvs_conus_stress
#SBATCH --time=01:30:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/logs/%x_%A_%a.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/logs/%x_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
#
# Full-CONUS FVS stress test, one array task per (variant, batch).
# STAGED, NOT AUTO-SUBMITTED. Submit only after prep + manifest are built:
#
#   STRESS=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
#   python $STRESS/prep_split_standinit.py --standinit $FIA/ENTIRE_FVS_STANDINIT_PLOT.csv \
#          --out-dir $STRESS/standinit_by_variant --batch-size 5000      # one-time (run as small sbatch)
#   python $STRESS/build_manifest.py --counts $STRESS/standinit_by_variant/counts.tsv \
#          --manifest $STRESS/manifest.tsv --batch-size 5000             # prints the --array range
#   sbatch --array=0-<N>%16 $STRESS/submit_conus_stress.slurm            # %16 backfills around other jobs
#
# Throttle (%16) keeps this polite while the asym_agb jobs are running.

set -euo pipefail
module load gcc/12.3.0 >/dev/null 2>&1
module load python/3.12 >/dev/null 2>&1

STRESS=/fs/scratch/PUOM0008/crsfaaron/fvs_stress
export FVS_PROJECT_ROOT=$HOME/fvs-modern
export FVS_LIB_DIR=$HOME/fvs-modern/lib
export FVS_FIA_DATA_DIR=/fs/scratch/PUOM0008/crsfaaron/FIA
mkdir -p "$STRESS/logs" "$STRESS/out"

echo "task ${SLURM_ARRAY_TASK_ID} on $(hostname) $(date)"
python3 "$STRESS/run_stress_task.py" \
    --manifest      "$STRESS/manifest.tsv" \
    --standinit-dir "$STRESS/standinit_by_variant" \
    --output-dir    "$STRESS/out" \
    --task-id       "${SLURM_ARRAY_TASK_ID}"
echo "task ${SLURM_ARRAY_TASK_ID} complete $(date)"
