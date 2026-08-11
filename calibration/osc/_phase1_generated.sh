#!/bin/bash
#SBATCH --job-name=fvs-p1
#SBATCH --account=PUOM0008
#SBATCH --array=0-24
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=36:00:00
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/phase1_%A_%a.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/phase1_%A_%a.err

# ============================================================================
# Variant Mapping
# ============================================================================

VARIANTS=( acd ak bc bm ca ci cr cs ec em ie kt ls nc ne oc on op pn sn so tt ut wc ws)

VARIANT="${VARIANTS[$SLURM_ARRAY_TASK_ID]}"

if [ -z "$VARIANT" ]; then
    echo "ERROR: Invalid array task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "==========================================="
echo "FVS Phase 1: ${VARIANT}"
echo "Array Task: ${SLURM_ARRAY_TASK_ID} of ${#VARIANTS[@]}"
echo "Node: $(hostname)"
echo "Cores: ${SLURM_CPUS_PER_TASK}"
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
# Environment
# ============================================================================

export R_LIBS_USER="$(Rscript -e 'cat(.libPaths()[1])' 2>/dev/null)"
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export CMDSTAN_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export FVS_FIA_DATA_DIR="/path/to/user/path"
export FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}"

# Max observations for subsampling (controls runtime)
export FVS_MAX_OBS=30000

# ============================================================================
# Directories
# ============================================================================

SCRIPTS_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/R"
OUTPUT_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/output/variants/${VARIANT}"
DATA_DIR="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/data/processed/${VARIANT}"
LOG="${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/phase1_${VARIANT}.log"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs"

# ============================================================================
# Preflight Checks
# ============================================================================

echo "[$(date)] Preflight checks for ${VARIANT}..." | tee "$LOG"

# Verify FIA data exists (produced by step 01 in prior run)
if [ ! -f "${DATA_DIR}/diameter_growth.csv" ]; then
    echo "FATAL: No diameter_growth.csv for ${VARIANT}. Run 01_fetch_fia_data.R first." | tee -a "$LOG"
    exit 1
fi

# Verify diameter growth posteriors exist (produced by step 02 in prior run)
if [ ! -f "${OUTPUT_DIR}/diameter_growth_map.csv" ] &&    [ ! -f "${OUTPUT_DIR}/diameter_growth_posterior.csv" ]; then
    echo "WARNING: No diameter growth output for ${VARIANT}. Step 06 may be incomplete." | tee -a "$LOG"
fi

# Count available data files
echo "Available data files:" | tee -a "$LOG"
ls -1 "${DATA_DIR}/"*.csv 2>/dev/null | while read f; do
    echo "  $(basename $f): $(wc -l < $f) lines" | tee -a "$LOG"
done

# ============================================================================
# Step 03: Height Diameter Model (Chapman Richards via brms)
# ============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] STEP 03: Height diameter model..." | tee -a "$LOG"

if [ -f "${DATA_DIR}/height_diameter.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/03_fit_height_diameter.R"         --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 03: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 03: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 03: SKIPPED (no height_diameter.csv)" | tee -a "$LOG"
fi

# ============================================================================
# Step 03b: Height Increment Model (CmdStanR)
# Automatically skips if variant lacks HG parameters in config
# ============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] STEP 03b: Height increment model..." | tee -a "$LOG"

if [ -f "${DATA_DIR}/height_growth.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/03b_fit_height_increment.R"         --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 03b: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 03b: SKIPPED or FAILED (non fatal)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 03b: SKIPPED (no height_growth.csv)" | tee -a "$LOG"
fi

# ============================================================================
# Step 04: Mortality Model (Logistic via brms)
# ============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] STEP 04: Mortality model..." | tee -a "$LOG"

if [ -f "${DATA_DIR}/mortality.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/04_fit_mortality.R"         --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 04: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 04: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 04: SKIPPED (no mortality.csv)" | tee -a "$LOG"
fi

# ============================================================================
# Step 05: Crown Ratio Change Model (Linear via brms)
# ============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] STEP 05: Crown ratio model..." | tee -a "$LOG"

if [ -f "${DATA_DIR}/crown_ratio_change.csv" ]; then
    if Rscript "${SCRIPTS_DIR}/05_fit_crown_ratio.R"         --variant "$VARIANT" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 05: SUCCESS" | tee -a "$LOG"
    else
        echo "[$(date)] STEP 05: FAILED (non fatal, continuing)" | tee -a "$LOG"
    fi
else
    echo "[$(date)] STEP 05: SKIPPED (no crown_ratio_change.csv)" | tee -a "$LOG"
fi

# ============================================================================
# Step 08 + 09: Stand Density (conditional)
# ============================================================================

SKIP_SDI=0

if [ "$SKIP_SDI" -eq 0 ]; then
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEP 08: Stand level data extraction..." | tee -a "$LOG"

    if Rscript "${SCRIPTS_DIR}/08_fetch_stand_data.R"         --variant "$VARIANT" --fia-dir "/path/to/user/path" >> "$LOG" 2>&1; then
        echo "[$(date)] STEP 08: SUCCESS" | tee -a "$LOG"

        echo "[$(date)] STEP 09: Stand density calibration..." | tee -a "$LOG"
        if Rscript "${SCRIPTS_DIR}/09_fit_stand_density.R"             --variant "$VARIANT" >> "$LOG" 2>&1; then
            echo "[$(date)] STEP 09: SUCCESS" | tee -a "$LOG"
        else
            echo "[$(date)] STEP 09: FAILED (non fatal)" | tee -a "$LOG"
        fi
    else
        echo "[$(date)] STEP 08: FAILED (non fatal, skipping step 09)" | tee -a "$LOG"
    fi
else
    echo "" | tee -a "$LOG"
    echo "[$(date)] STEPS 08+09: SKIPPED (--skip-stand-density)" | tee -a "$LOG"
fi

# ============================================================================
# Step 06: Posterior to JSON (combines ALL components)
# ============================================================================

echo "" | tee -a "$LOG"
echo "[$(date)] STEP 06: Posterior to JSON..." | tee -a "$LOG"

if Rscript "${SCRIPTS_DIR}/06_posterior_to_json.R"     --variant "$VARIANT" >> "$LOG" 2>&1; then
    echo "[$(date)] STEP 06: SUCCESS" | tee -a "$LOG"
else
    echo "[$(date)] STEP 06: FAILED" | tee -a "$LOG"
    echo "WARNING: Calibrated JSON not produced for ${VARIANT}" | tee -a "$LOG"
fi

# ============================================================================
# Summary
# ============================================================================

echo "" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"
echo "Phase 1 complete for variant: ${VARIANT}" | tee -a "$LOG"
echo "End time: $(date)" | tee -a "$LOG"
echo "==========================================" | tee -a "$LOG"

echo "" | tee -a "$LOG"
echo "Output files:" | tee -a "$LOG"
ls -lh "${OUTPUT_DIR}/" 2>/dev/null | tee -a "$LOG"

# Count successful components
N_DONE=0
for f in height_diameter_posterior.csv height_increment_posterior.csv          mortality_posterior.csv crown_ratio_posterior.csv          stand_density_posterior.csv; do
    if [ -f "${OUTPUT_DIR}/${f}" ]; then
        N_DONE=$((N_DONE + 1))
    fi
done

echo "" | tee -a "$LOG"
echo "Components completed: ${N_DONE} of 5 possible" | tee -a "$LOG"

if [ -f "${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/config/calibrated/${VARIANT}.json" ]; then
    echo "Calibrated config: YES" | tee -a "$LOG"
else
    echo "Calibrated config: NO" | tee -a "$LOG"
fi
