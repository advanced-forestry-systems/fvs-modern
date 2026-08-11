#!/bin/bash
#
# Resubmit ALL variants with fixes for:
#   1. Step 03: H-D Chapman-Richards NaN at init (added lb bounds + safe init)
#   2. Steps 04/05: posterior extraction p95 column naming (unname + gsub)
#   3. Recovery: extract posteriors from existing RDS files
#   4. Step 08: 64G memory for rFIA
#
# This runs extraction first (fast, uses existing RDS), then re-runs H-D,
# then 04/05 with fixes, then 08/09 with more memory, then 06 (JSON).
#
# Usage: bash calibration/osc/resubmit_all_fixes.sh

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BATCH_SCRIPT="${SCRIPT_DIR}/_resubmit_all_generated.sh"

cat > "$BATCH_SCRIPT" << 'BATCH_EOF'
#!/bin/bash
#SBATCH --job-name=fvs-fix
#SBATCH --account=PUOM0008
#SBATCH --array=0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix_%A_%a.err

VARIANTS=( acd ak bc bm ca ci cr cs ec em ie kt ls nc ne oc on op pn sn so tt ut wc ws )
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Fix Resubmission: ${VARIANT}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID}"
echo "Node: $(hostname), Cores: ${SLURM_CPUS_PER_TASK}, Mem: 64G"
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
LOG="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/fix_${VARIANT}.log"

mkdir -p "${OUTPUT_DIR}"

echo "[$(date)] FIX RESUBMISSION for ${VARIANT}" | tee "$LOG"
COMPLETED=0
FAILED=0

# ---- Phase 1: Extract posteriors from existing RDS files ----
echo "" | tee -a "$LOG"
echo "[$(date)] EXTRACT: Recovering posteriors from RDS files..." | tee -a "$LOG"
if Rscript "${SCRIPTS_DIR}/extract_posteriors.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "[$(date)] EXTRACT: SUCCESS" | tee -a "$LOG"
    COMPLETED=$((COMPLETED+1))
else
    echo "[$(date)] EXTRACT: FAILED (no RDS files yet, will fit fresh)" | tee -a "$LOG"
    FAILED=$((FAILED+1))
fi

# ---- Phase 2: Re-run H-D with bounded params (only if no H-D posterior yet) ----
if [ ! -f "${OUTPUT_DIR}/height_diameter_posterior.csv" ] && \
   [ ! -f "${OUTPUT_DIR}/height_diameter_map.csv" ]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEP 03: Height-diameter (fixed init + bounds)..." | tee -a "$LOG"
    if [ -f "${DATA_DIR}/height_diameter.csv" ]; then
        if Rscript "${SCRIPTS_DIR}/03_fit_height_diameter.R" \
            --variant "$VARIANT" >> "$LOG" 2>&1; then
            echo "[$(date)] STEP 03: SUCCESS" | tee -a "$LOG"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 03: FAILED (non fatal)" | tee -a "$LOG"
            FAILED=$((FAILED+1))
        fi
    else
        echo "[$(date)] STEP 03: SKIPPED (no data)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 03: SKIPPED (output already exists)" | tee -a "$LOG"
fi

# ---- Phase 3: Re-run mortality if no posterior CSV yet ----
if [ ! -f "${OUTPUT_DIR}/mortality_posterior.csv" ] || \
   [ ! -f "${OUTPUT_DIR}/mortality_map.csv" ]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEP 04: Mortality (fixed extraction)..." | tee -a "$LOG"
    if [ -f "${DATA_DIR}/mortality.csv" ]; then
        if Rscript "${SCRIPTS_DIR}/04_fit_mortality.R" \
            --variant "$VARIANT" >> "$LOG" 2>&1; then
            echo "[$(date)] STEP 04: SUCCESS" | tee -a "$LOG"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 04: FAILED (non fatal)" | tee -a "$LOG"
            FAILED=$((FAILED+1))
        fi
    fi
else
    echo "[$(date)] STEP 04: SKIPPED (output already exists)" | tee -a "$LOG"
fi

# ---- Phase 4: Re-run crown ratio if no posterior CSV yet ----
if [ ! -f "${OUTPUT_DIR}/crown_ratio_posterior.csv" ] || \
   [ ! -f "${OUTPUT_DIR}/crown_ratio_map.csv" ]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEP 05: Crown ratio (fixed extraction)..." | tee -a "$LOG"
    if [ -f "${DATA_DIR}/crown_ratio_change.csv" ]; then
        if Rscript "${SCRIPTS_DIR}/05_fit_crown_ratio.R" \
            --variant "$VARIANT" >> "$LOG" 2>&1; then
            echo "[$(date)] STEP 05: SUCCESS" | tee -a "$LOG"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 05: FAILED (non fatal)" | tee -a "$LOG"
            FAILED=$((FAILED+1))
        fi
    fi
else
    echo "[$(date)] STEP 05: SKIPPED (output already exists)" | tee -a "$LOG"
fi

# ---- Phase 5: Stand density (64G memory) ----
if [ ! -f "${OUTPUT_DIR}/stand_density.csv" ]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEP 08: Stand data extraction (64G)..." | tee -a "$LOG"
    if Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R" \
        --variant "$VARIANT" --fia-dir "/path/to/user/path" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 08: SUCCESS" | tee -a "$LOG"
        COMPLETED=$((COMPLETED+1))

        echo "[$(date)] STEP 09: SDIMAX calibration..." | tee -a "$LOG"
        if Rscript "${SCRIPTS_DIR}/09_fit_stand_density.R" \
            --variant "$VARIANT" >> "$LOG" 2>&1; then
            echo "[$(date)] STEP 09: SUCCESS" | tee -a "$LOG"
            COMPLETED=$((COMPLETED+1))
        else
            echo "[$(date)] STEP 09: FAILED (non fatal)" | tee -a "$LOG"
            FAILED=$((FAILED+1))
        fi
    else
        echo "[$(date)] STEP 08: FAILED (non fatal)" | tee -a "$LOG"
        FAILED=$((FAILED+1))
    fi
else
    echo "[$(date)] STEPS 08+09: SKIPPED (data exists)" | tee -a "$LOG"
fi

# ---- Phase 6: JSON conversion ----
echo "" | tee -a "$LOG"
echo "[$(date)] STEP 06: Posterior to JSON..." | tee -a "$LOG"
if Rscript "${SCRIPTS_DIR}/06_posterior_to_json.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "[$(date)] STEP 06: SUCCESS" | tee -a "$LOG"
    COMPLETED=$((COMPLETED+1))
else
    echo "[$(date)] STEP 06: FAILED" | tee -a "$LOG"
    FAILED=$((FAILED+1))
fi

echo "" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
echo "DONE: ${VARIANT} | ok=${COMPLETED} fail=${FAILED}" | tee -a "$LOG"
echo "End: $(date)" | tee -a "$LOG"
ls -la "${OUTPUT_DIR}/" 2>/dev/null | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
BATCH_EOF

echo "Generated batch script: $BATCH_SCRIPT"
echo ""
echo "NOTE: Wait for current job 7429482 to finish (all variants on Step 08)"
echo "      before submitting to avoid file conflicts."
echo ""
echo "To submit now:  sbatch $BATCH_SCRIPT"
echo "To submit later: sbatch $BATCH_SCRIPT"
