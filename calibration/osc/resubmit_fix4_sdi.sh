#!/bin/bash
#SBATCH --job-name=fvs-fix4-sdi
#SBATCH --account=PUOM0008
#SBATCH --array=0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=2:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix4_sdi_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix4_sdi_%A_%a.err

# =============================================================================
# FVS Fix Round 4: Stand Density Only (v4 elev patch) + JSON re-export
# Patch: 08_fetch_stand_data.R v4 - fixed elev = first(c(na.omit(.elev_col), NA_real_))
# Runs AFTER fix3 has completed H-D for all 25 variants
# =============================================================================

VARIANTS=( acd ak bc bm ca ci cr cs ec em ie kt ls nc ne oc on op pn sn so tt ut wc ws )
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Fix4 SDI: ${VARIANT}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname), Cores: ${SLURM_CPUS_PER_TASK}, Mem: 128G"
echo "Start: $(date)"
echo "==========================================="

module purge
module load gcc/12.3.0
module load R/4.4.0
module load gdal/3.7.3
module load proj/9.2.1
module load geos/3.12.0

export R_LIBS_USER="$(Rscript -e 'cat(.libPaths()[1])' 2>/dev/null)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export CMDSTAN_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export FVS_FIA_DATA_DIR="/path/to/user/path"
export FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}"
export FVS_MAX_OBS=30000

SCRIPTS_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/R"
OUTPUT_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/${VARIANT}"
DATA_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/processed/${VARIANT}"

mkdir -p "${OUTPUT_DIR}"

COMPLETED=0
FAILED=0

# ---- Skip if stand density already done ----
if [ -f "${OUTPUT_DIR}/stand_density_summary.csv" ]; then
    echo "[$(date)] Stand density already complete for ${VARIANT}, re-exporting JSON only"
else
    # ---- Step 1: Force re-extract stand data with v4 elev patch ----
    echo "[$(date)] Removing old stand_density.csv if exists..."
    rm -f "${DATA_DIR}/stand_density.csv" "${DATA_DIR}/species_sdi.csv"

    echo "[$(date)] STEP 08: Stand data extraction (v4 elev patched)..."
    if Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R" \
        --variant "$VARIANT" --fia-dir "/path/to/user/path" 2>&1; then
        echo "[$(date)] STEP 08: SUCCESS"
        COMPLETED=$((COMPLETED+1))
    else
        echo "[$(date)] STEP 08: FAILED (non fatal)"
        FAILED=$((FAILED+1))
    fi

    # ---- Step 2: Fit stand density / SDIMAX ----
    if [ -f "${DATA_DIR}/stand_density.csv" ] && [ "$(wc -l < ${DATA_DIR}/stand_density.csv)" -gt 1 ]; then
        echo "[$(date)] STEP 09: SDIMAX calibration..."
        if Rscript "${SCRIPTS_DIR}/09_fit_stand_density.R" \
            --variant "$VARIANT" 2>&1; then
            echo "[$(date)] STEP 09: SUCCESS"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 09: FAILED (non fatal)"
            FAILED=$((FAILED+1))
        fi
    else
        echo "[$(date)] STEP 09: SKIPPED (no stand data or empty)"
    fi
fi

# ---- Step 3: JSON Export (always rerun to pick up new sdimax) ----
echo ""
echo "[$(date)] STEP 06: Posterior to JSON..."
if Rscript "${SCRIPTS_DIR}/06_posterior_to_json.R" \
    --variant "$VARIANT" 2>&1; then
    echo "[$(date)] STEP 06: SUCCESS"
    COMPLETED=$((COMPLETED+1))
else
    echo "[$(date)] STEP 06: FAILED"
    FAILED=$((FAILED+1))
fi

echo ""
echo "==========================================="
echo "DONE: ${VARIANT} | ok=${COMPLETED} fail=${FAILED}"
echo "End: $(date)"
ls -la "${OUTPUT_DIR}/" 2>/dev/null
echo "==========================================="
