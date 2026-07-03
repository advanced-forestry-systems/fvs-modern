#!/bin/bash
#SBATCH --job-name=fvs-cal
#SBATCH --account=PUOM0008
#SBATCH --array=0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/slurm_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/slurm_%A_%a.err

# ============================================================================
# Variant Array Mapping
# ============================================================================

VARIANTS=(
  acd ak bc bm ca ci cr cs ec em
  ie kt ls nc ne oc on op pn sn
  so tt ut wc ws
)

VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Calibration: Variant ${VARIANT}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID} of ${#VARIANTS[@]}"
echo "Node: $(hostname)"
echo "Start: $(date)"
echo "==========================================="

# ============================================================================
# Load Modules
# ============================================================================

module purge
module load gcc/12.3.0
module load R/4.4.0
module load gdal/3.7.3
module load proj/9.2.1
module load geos/3.12.0

# ============================================================================
# Set Environment
# ============================================================================

# OSC uses per cluster R library paths: ~/R/<cluster>/<R_version>
# Detect dynamically from R itself to handle any cluster/version combo
export R_LIBS_USER="$(Rscript -e 'cat(.libPaths()[1])' 2>/dev/null)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export CMDSTAN_NUM_THREADS=${SLURM_CPUS_PER_TASK}

# Tell the R scripts where the FIA data lives
export FVS_FIA_DATA_DIR="/path/to/user/path"
export FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}"

# ============================================================================
# Create Directories
# ============================================================================

mkdir -p "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs"
mkdir -p "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/${VARIANT}"
mkdir -p "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/processed/${VARIANT}"

# ============================================================================
# Run Pipeline
# ============================================================================

SCRIPTS_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/R"
LOG="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/variant_${VARIANT}.log"

echo "Starting pipeline for variant ${VARIANT}..." | tee "$LOG"
echo "FIA data: /path/to/user/path" | tee -a "$LOG"
echo "Cores: ${SLURM_CPUS_PER_TASK}" | tee -a "$LOG"

# Step 01: Fetch and prepare FIA data from local files
echo "[$(date)] STEP 01: Preparing FIA data..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/01_fetch_fia_data.R" \
    --variant "$VARIANT" --fia-dir "/path/to/user/path" >> "$LOG" 2>&1
echo "[$(date)] STEP 01: Done" | tee -a "$LOG"

# Step 02: Diameter growth (REQUIRED)
echo "[$(date)] STEP 02: Diameter growth model..." | tee -a "$LOG"
if ! Rscript "${SCRIPTS_DIR}/02_fit_diameter_growth.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "FATAL: Diameter growth failed for ${VARIANT}" | tee -a "$LOG"
    exit 1
fi
echo "[$(date)] STEP 02: Done" | tee -a "$LOG"

# Step 03: Height diameter
echo "[$(date)] STEP 03: Height diameter model..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/03_fit_height_diameter.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 03: Done" | tee -a "$LOG"

# Step 03b: Height increment (auto skips if variant lacks HG params)
echo "[$(date)] STEP 03b: Height increment model..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/03b_fit_height_increment.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 03b: Done" | tee -a "$LOG"

# Step 04: Mortality
echo "[$(date)] STEP 04: Mortality model..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/04_fit_mortality.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 04: Done" | tee -a "$LOG"

# Step 05: Crown ratio
echo "[$(date)] STEP 05: Crown ratio model..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/05_fit_crown_ratio.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 05: Done" | tee -a "$LOG"

# Step 08: Stand level data extraction
echo "[$(date)] STEP 08: Stand level data..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R" \
    --variant "$VARIANT" --fia-dir "/path/to/user/path" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 08: Done" | tee -a "$LOG"

# Step 09: SDIMAX / BAMAX / self thinning calibration
echo "[$(date)] STEP 09: Stand density calibration..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/09_fit_stand_density.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 09: Done" | tee -a "$LOG"

# Step 06: Convert all posteriors to calibrated JSON (REQUIRED)
echo "[$(date)] STEP 06: Posterior to JSON..." | tee -a "$LOG"
if ! Rscript "${SCRIPTS_DIR}/06_posterior_to_json.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "FATAL: Posterior to JSON failed for ${VARIANT}" | tee -a "$LOG"
    exit 1
fi
echo "[$(date)] STEP 06: Done" | tee -a "$LOG"

# Step 07: Diagnostics
echo "[$(date)] STEP 07: Diagnostics..." | tee -a "$LOG"
Rscript "${SCRIPTS_DIR}/07_diagnostics.R" \
    --variant "$VARIANT" >> "$LOG" 2>&1 || true
echo "[$(date)] STEP 07: Done" | tee -a "$LOG"

# ============================================================================
# Summary
# ============================================================================

echo "" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
echo "Pipeline complete for variant: ${VARIANT}" | tee -a "$LOG"
echo "End time: $(date)" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Outputs:" | tee -a "$LOG"
ls -lh "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/${VARIANT}/" 2>/dev/null | tee -a "$LOG"

if [ -f "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/config/calibrated/${VARIANT}.json" ]; then
    echo "Calibrated config: config/calibrated/${VARIANT}.json" | tee -a "$LOG"
else
    echo "WARNING: No calibrated config produced" | tee -a "$LOG"
fi
