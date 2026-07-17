#!/bin/bash
#SBATCH --job-name=v4_ne_redo --account=PUOM0008 --time=24:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8 --mem=64G
#SBATCH --array=288-291
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_v4_logs/v4redo_%A_%a.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_v4_logs/v4redo_%A_%a.err
module load gcc/12.3.0 python/3.12 >/dev/null 2>&1
export PYTHONNOUSERSITE=1
export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
export FVS_LIB_DIR=/users/PUOM0008/crsfaaron/fvs-modern/lib
export NSBE_ROOT=/users/PUOM0008/crsfaaron/fvs-modern/data/NSBE
export FIA_DATA_DIR=/fs/scratch/PUOM0008/crsfaaron/FIA
cd /fs/scratch/PUOM0008/crsfaaron/fvs_stress
python3 run_conus_task_wo1_v4.py --manifest manifest.tsv --standinit-dir standinit_by_variant \
  --treeinit-dir /fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit_h \
  --output-dir /fs/scratch/PUOM0008/crsfaaron/fvs_stress/out_conus_engine_v4 --task-id ${SLURM_ARRAY_TASK_ID}
