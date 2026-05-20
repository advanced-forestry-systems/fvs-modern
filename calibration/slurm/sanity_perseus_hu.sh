#!/bin/bash
# =============================================================================
# Sanity test for perseus_harvest_uncertainty.py
# =============================================================================
# 1 batch (50 plots), 5 draws, 1 scenario (RCP 4.5 harvest), both variants.
# Verifies the new harvest+climate+uncertainty driver runs end-to-end before
# launching the full 500-draw factorial.
# =============================================================================
set -euo pipefail
FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/users/PUOM0008/crsfaaron/fvs-modern}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-PUOM0008}"
FVS_FIA_DATA_DIR="${FVS_FIA_DATA_DIR:-/users/PUOM0008/crsfaaron/fia_data}"
SI_RATIO_CSV="${SI_RATIO_CSV:-/users/PUOM0008/crsfaaron/SiteIndex/future_SI/perseus_si_ratios.csv}"

PROJECT_DIR="${FVS_PROJECT_ROOT}"
SCRIPT="${PROJECT_DIR}/calibration/python/perseus_harvest_uncertainty.py"
PERSEUS_CSV="${PROJECT_DIR}/calibration/data/perseus_plots.csv"
BASE_OUTPUT_DIR="${PROJECT_DIR}/calibration/output/perseus_hu_sanity"
LOG_DIR="${PROJECT_DIR}/calibration/logs/perseus_hu_sanity"
mkdir -p "${BASE_OUTPUT_DIR}" "${LOG_DIR}"

if [ ! -f "$SCRIPT" ]; then
  echo "ERROR: $SCRIPT not found"; exit 1
fi
if [ ! -f "$SI_RATIO_CSV" ]; then
  echo "ERROR: SI ratio CSV not found: $SI_RATIO_CSV"; exit 1
fi

JOB_ID=$(sbatch --parsable <<SBATCH
#!/bin/bash
#SBATCH --job-name=hu_sanity
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${LOG_DIR}/sanity.out
#SBATCH --error=${LOG_DIR}/sanity.err

module purge
module load python/3.12

export FVS_PROJECT_ROOT="${PROJECT_DIR}"
export FVS_FIA_DATA_DIR="${FVS_FIA_DATA_DIR}"
export FVS_LIB_DIR="${PROJECT_DIR}/lib"
export LD_LIBRARY_PATH="${PROJECT_DIR}/lib:\${LD_LIBRARY_PATH:-}"
export NSBE_ROOT="${PROJECT_DIR}/data/NSBE"
export FVS_CONFIG_DIR="${PROJECT_DIR}/config"
export PERSEUS_OUTPUT_DIR="${BASE_OUTPUT_DIR}"
export PERSEUS_CSV="${PERSEUS_CSV}"
export SI_RATIO_CSV="${SI_RATIO_CSV}"

echo "=== HU SANITY ==="
echo "Start: \$(date)"
echo "Node:  \$(hostname)"

python3 "${SCRIPT}" \
    --batch-id 1 --batch-size 50 \
    --n-draws 5 \
    --variants ne acd \
    --rcp 4.5 --harvest \
    --output-dir "${BASE_OUTPUT_DIR}" \
    --perseus-csv "${PERSEUS_CSV}" \
    --si-ratio-csv "${SI_RATIO_CSV}" \
    --seed 42

echo ""
echo "Files written:"
ls -la "${BASE_OUTPUT_DIR}/harvest_rcp45/" 2>&1 || true
echo "End: \$(date)"
SBATCH
)
echo "Submitted sanity job: ${JOB_ID}"
echo "Tail logs: tail -f ${LOG_DIR}/sanity.{out,err}"
