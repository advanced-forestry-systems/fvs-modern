#!/bin/bash
#SBATCH --job-name=fvs-oom-fix
#SBATCH --account=PUOM0008
#SBATCH --array=0-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/oomfix_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/oomfix_%A_%a.err

VARIANTS=( ls sn )
VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS OOM Resubmission: ${VARIANT} (128G)"
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

export R_LIBS_USER=$(Rscript -e 'cat(.libPaths()[1])' 2>/dev/null)
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export CMDSTAN_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export FVS_FIA_DATA_DIR="/path/to/user/path"
export FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}"
export FVS_MAX_OBS=30000

SCRIPTS_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/R"
OUTPUT_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/${VARIANT}"
DATA_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/processed/${VARIANT}"
LOG="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/oomfix_${VARIANT}.log"

echo "[$(date)] OOM RESUBMISSION for ${VARIANT} with 128G" | tee "$LOG"
COMPLETED=0
FAILED=0

# Core model files already exist; just need remaining steps

# ---- Extract posteriors (if needed) ----
echo "[$(date)] EXTRACT: Recovering posteriors..." | tee -a "$LOG"
if Rscript "${SCRIPTS_DIR}/extract_posteriors.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "[$(date)] EXTRACT: SUCCESS" | tee -a "$LOG"
    COMPLETED=$((COMPLETED+1))
else
    echo "[$(date)] EXTRACT: FAILED" | tee -a "$LOG"
    FAILED=$((FAILED+1))
fi

# ---- Stand data extraction + SDIMAX (the OOM step) ----
if [ ! -f "${OUTPUT_DIR}/stand_density.csv" ]; then
    echo "[$(date)] STEP 08: Stand data extraction (128G)..." | tee -a "$LOG"
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
        echo "[$(date)] STEP 08: FAILED" | tee -a "$LOG"
        FAILED=$((FAILED+1))
    fi
else
    echo "[$(date)] STEPS 08+09: SKIPPED (data exists)" | tee -a "$LOG"
fi

# ---- JSON conversion ----
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
