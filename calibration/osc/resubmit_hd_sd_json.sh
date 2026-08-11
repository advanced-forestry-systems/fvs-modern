#!/bin/bash
#SBATCH --job-name=fvs-fix2
#SBATCH --account=PUOM0008
#SBATCH --array=0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix2_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix2_%A_%a.err

# =============================================================================
# FVS Fix Round 2: Height-Diameter + Stand Density + JSON Export
# Patches applied:
#   - 03_fit_height_diameter.R: MAP-based init for HMC (fixes NaN at initialization)
#   - 08_fetch_stand_data.R: any_of() for COND columns (fixes missing ELEV)
# =============================================================================

VARIANTS=( acd ak bc bm ca ci cr cs ec em ie kt ls nc ne oc on op pn sn so tt ut wc ws )
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Fix2: ${VARIANT}"
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

# ---- Step 1: Height-Diameter (patched: MAP-based init) ----
if [ ! -f "${OUTPUT_DIR}/height_diameter_map.csv" ] && \
   [ ! -f "${OUTPUT_DIR}/height_diameter_posterior.csv" ]; then
    echo ""
    echo "[$(date)] STEP 03: Height-diameter (patched MAP init)..."
    if [ -f "${DATA_DIR}/height_diameter.csv" ]; then
        if Rscript "${SCRIPTS_DIR}/03_fit_height_diameter.R" \
            --variant "$VARIANT" 2>&1; then
            echo "[$(date)] STEP 03: SUCCESS"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 03: FAILED (non fatal)"
            FAILED=$((FAILED+1))
        fi
    else
        echo "[$(date)] STEP 03: SKIPPED (no data file)"
    fi
else
    echo "[$(date)] STEP 03: SKIPPED (output already exists)"
fi

# ---- Step 2: Stand Data Extraction (patched: any_of for ELEV) ----
# Force re-extraction since previous run produced empty files
if [ -f "${DATA_DIR}/stand_density.csv" ]; then
    NLINES=$(wc -l < "${DATA_DIR}/stand_density.csv")
    if [ "$NLINES" -le 1 ]; then
        echo "[$(date)] Removing empty stand_density.csv (header only)"
        rm -f "${DATA_DIR}/stand_density.csv"
        rm -f "${DATA_DIR}/species_sdi.csv"
    fi
fi

if [ ! -f "${DATA_DIR}/stand_density.csv" ] || [ "$(wc -l < ${DATA_DIR}/stand_density.csv)" -le 1 ]; then
    echo ""
    echo "[$(date)] STEP 08: Stand data extraction (patched ELEV)..."
    if Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R" \
        --variant "$VARIANT" --fia-dir "/path/to/user/path" 2>&1; then
        echo "[$(date)] STEP 08: SUCCESS"
        COMPLETED=$((COMPLETED+1))
    else
        echo "[$(date)] STEP 08: FAILED (non fatal)"
        FAILED=$((FAILED+1))
    fi
else
    echo "[$(date)] STEP 08: SKIPPED (stand data exists with data)"
fi

# ---- Step 3: Stand Density / SDIMAX Calibration ----
if [ -f "${DATA_DIR}/stand_density.csv" ] && [ "$(wc -l < ${DATA_DIR}/stand_density.csv)" -gt 1 ]; then
    if [ ! -f "${OUTPUT_DIR}/stand_density_summary.csv" ]; then
        echo ""
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
        echo "[$(date)] STEP 09: SKIPPED (output exists)"
    fi
else
    echo "[$(date)] STEP 09: SKIPPED (no stand data)"
fi

# ---- Step 4: JSON Export ----
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
