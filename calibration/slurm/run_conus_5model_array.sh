#!/bin/bash
#SBATCH --job-name=conus_5model
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --array=0-4
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_5model_%a_%A.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_5model_%a_%A.err

# CONUS 5-model scorecard array — one task per compiled FVS variant.
# Compiled variants as of 2026-07-10: ne sn wc em pn
# Expand VARIANTS array and update --array range when additional variants are compiled.
#
# States covered:
#   ne  -> CT DE MA MD ME NH NJ NY PA RI VT WV
#   sn  -> AL GA MS SC
#   wc  -> OR WA
#   em  -> MT ND SD
#   pn  -> OR WA
#
# States NOT yet covered (no compiled binary as of 2026-07-10):
#   AK, AZ, AR, CA, CO, FL, IA, ID, IL, IN, KS, KY, LA, MI, MN, MO,
#   NC, NE, NM, NV, OH, OK, TN, TX, UT, VA, WI, WY
# Pending variants: ls (MI MN WI), cs (IL IN MO), ca (CA), so (FL LA NC TX VA),
#   cr (CO WY), ie (ID MT WA), kt (MT ID), ut (UT NV)

module load gcc/12.3.0
module load python/3.11

VARIANTS=(ne sn wc em pn)
VARIANT=${VARIANTS[$SLURM_ARRAY_TASK_ID]}

FIA_DIR=/fs/scratch/PUOM0008/crsfaaron/FIA
OUT_DIR=/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus
SCRIPT=/users/PUOM0008/crsfaaron/fvs-modern/calibration/run_5model_scorecard.py

echo "[conus_array] task=${SLURM_ARRAY_TASK_ID} variant=${VARIANT}"
echo "[conus_array] node=$(hostname)  cpus=${SLURM_CPUS_PER_TASK}"

cd /users/PUOM0008/crsfaaron/fvs-modern/calibration

python3 "$SCRIPT" \
    --variant   "$VARIANT" \
    --max-pairs 5000 \
    --fia-dir   "$FIA_DIR" \
    --outdir    "$OUT_DIR"

EXIT=$?
echo "[conus_array] variant=${VARIANT} exit=${EXIT}"
exit $EXIT
