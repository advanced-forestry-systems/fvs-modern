#!/bin/bash
#SBATCH --job-name=fvs-resub
#SBATCH --account=PUOM0008
#SBATCH --array=0-2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/resub_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/resub_%A_%a.err

VARIANTS=( oc tt ws )
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Phase 1 RESUBMISSION: ${VARIANT}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname)"
echo "Memory: 64G (increased from 32G)"
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
LOG="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/resub_${VARIANT}.log"

mkdir -p "${OUTPUT_DIR}"

echo "[$(date)] RESUBMISSION for ${VARIANT} with 64G memory" | tee "$LOG"
echo "[$(date)] Scripts have been patched for posterior::summarise_draws" | tee -a "$LOG"

# Step 03: H-D (will likely fail again for these variants, but try)
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 03: Height diameter model..." | tee -a "$LOG"
if [ -f "${DATA_DIR}/height_diameter.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/03_fit_height_diameter.R" \
        --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 03: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 03: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 03: SKIPPED (no data)" | tee -a "$LOG"
fi

# Step 03b: Height increment
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 03b: Height increment model..." | tee -a "$LOG"
if [ -f "${DATA_DIR}/height_growth.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/03b_fit_height_increment.R" \
        --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 03b: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 03b: SKIPPED or FAILED (non fatal)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 03b: SKIPPED (no data)" | tee -a "$LOG"
fi

# Step 04: Mortality (FIXED: posterior extraction)
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 04: Mortality model (patched)..." | tee -a "$LOG"
if [ -f "${DATA_DIR}/mortality.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/04_fit_mortality.R" \
        --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 04: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 04: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 04: SKIPPED (no data)" | tee -a "$LOG"
fi

# Step 05: Crown ratio (FIXED: posterior extraction)
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 05: Crown ratio model (patched)..." | tee -a "$LOG"
if [ -f "${DATA_DIR}/crown_ratio_change.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/05_fit_crown_ratio.R" \
        --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 05: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 05: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 05: SKIPPED (no data)" | tee -a "$LOG"
fi

# Step 08 + 09: Stand density (with 64G should not OOM)
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 08: Stand level data extraction (64G)..." | tee -a "$LOG"
if Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R" \
    --variant "$VARIANT" --fia-dir "/path/to/user/path" >> "$LOG" 2>&1; then
    echo "[$(date)] STEP 08: SUCCESS" | tee -a "$LOG"

    echo "[$(date)] STEP 09: Stand density calibration..." | tee -a "$LOG"
    if Rscript "${SCRIPTS_DIR}/09_fit_stand_density.R" \
        --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 09: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 09: FAILED (non fatal)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 08: FAILED (non fatal, skipping step 09)" | tee -a "$LOG"
fi

# Step 06: JSON conversion
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 06: Posterior to JSON..." | tee -a "$LOG"
if Rscript "${SCRIPTS_DIR}/06_posterior_to_json.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "[$(date)] STEP 06: SUCCESS" | tee -a "$LOG"
else
    echo "[$(date)] STEP 06: FAILED" | tee -a "$LOG"
fi

# Summary
echo "" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
echo "RESUBMISSION COMPLETE: ${VARIANT}" | tee -a "$LOG"
echo "End: $(date)" | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
ls -la "${OUTPUT_DIR}/" 2>/dev/null | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
