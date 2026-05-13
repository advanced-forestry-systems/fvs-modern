#!/bin/bash
# =============================================================================
# Full PERSEUS Harvest x Climate x Posterior Uncertainty factorial submitter
# =============================================================================
# Submits the complete 6-scenario uncertainty factorial for FVS NE and ACD.
#
# Decomposition: SLURM array task = (scenario x batch).
#   6 scenarios x N_BATCHES = total array tasks.
#   Each task runs BATCH_SIZE plots through one (rcp x harvest) corner with
#   N_DRAWS posterior draws plus default and calibrated MAP point estimates.
#
# Aggregation: chained per-scenario aggregator that runs after the array
#   completes, calling perseus_uncertainty_aggregate.py on each scenario subdir.
#
# Sizing notes (from the sanity test 9207452 once it lands):
#   * Per-plot wall time = (sanity wall) / 50 plots / 5 draws
#   * Per-task wall time = BATCH_SIZE x N_DRAWS x per_plot_per_draw_sec
#   * Pick BATCH_SIZE so per-task wall stays under 24 hr (the cpu partition cap)
#
# Defaults below assume ~2 sec per FVS run scaled from prior baseline pipelines.
# Override via env vars before running (BATCH_SIZE=25 N_DRAWS=500 bash submit_...).
#
# Usage:
#   bash submit_perseus_hu_factorial.sh                         # all 6 scenarios
#   SCENARIOS="harvest_rcp45 noharvest_rcp45" bash submit_...   # subset
#   BATCH_SIZE=25 N_DRAWS=500 bash submit_...                   # override sizing
#   DRY_RUN=1 bash submit_...                                   # print only
# =============================================================================
set -euo pipefail

# ---------------------------- Tunable parameters -----------------------------
BATCH_SIZE="${BATCH_SIZE:-25}"          # plots per task
N_DRAWS="${N_DRAWS:-500}"               # posterior draws per plot per scenario
N_PLOTS="${N_PLOTS:-3586}"              # PERSEUS plot count after year filter
TASK_WALL="${TASK_WALL:-24:00:00}"      # per-task wall (max 24 hr on cpu)
AGG_WALL="${AGG_WALL:-04:00:00}"
TASK_MEM="${TASK_MEM:-16G}"
AGG_MEM="${AGG_MEM:-32G}"
TASK_CPUS="${TASK_CPUS:-4}"
SEED="${SEED:-42}"
START_YEAR="${START_YEAR:-1999}"
END_YEAR="${END_YEAR:-2004}"
VARIANTS="${VARIANTS:-ne acd}"

# Six default scenarios. Comment out to skip the no-harvest no-climate cell
# if you want to reuse the existing baseline uncertainty_1999_2004 outputs.
SCENARIOS_DEFAULT="noharvest_noclimate noharvest_rcp45 noharvest_rcp85 harvest_noclimate harvest_rcp45 harvest_rcp85"
SCENARIOS="${SCENARIOS:-$SCENARIOS_DEFAULT}"

DRY_RUN="${DRY_RUN:-0}"

# ---------------------------- HPC config -------------------------------------
FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-/users/PUOM0008/crsfaaron/fvs-modern}"
SLURM_ACCOUNT="${SLURM_ACCOUNT:-PUOM0008}"
FVS_FIA_DATA_DIR="${FVS_FIA_DATA_DIR:-/users/PUOM0008/crsfaaron/fia_data}"
SI_RATIO_CSV="${SI_RATIO_CSV:-/users/PUOM0008/crsfaaron/SiteIndex/future_SI/perseus_si_ratios.csv}"

PROJECT_DIR="${FVS_PROJECT_ROOT}"
SCRIPT="${PROJECT_DIR}/calibration/python/perseus_harvest_uncertainty.py"
PERSEUS_CSV="${PROJECT_DIR}/calibration/data/perseus_plots.csv"
BASE_OUTPUT_DIR="${PROJECT_DIR}/calibration/output/perseus_hu_factorial"
LOG_DIR="${PROJECT_DIR}/calibration/logs/perseus_hu_factorial"
mkdir -p "${BASE_OUTPUT_DIR}" "${LOG_DIR}"

# ---------------------------- Validation -------------------------------------
[ -f "${SCRIPT}" ] || { echo "ERROR: script missing: ${SCRIPT}"; exit 1; }
[ -f "${SI_RATIO_CSV}" ] || { echo "ERROR: SI ratio CSV missing"; exit 1; }
[ -f "${PERSEUS_CSV}" ] || { echo "ERROR: PERSEUS CSV missing"; exit 1; }

N_BATCHES=$(( (N_PLOTS + BATCH_SIZE - 1) / BATCH_SIZE ))

cat <<EOF
=== PERSEUS HU Factorial submission ===
Script:        ${SCRIPT}
Scenarios:     ${SCENARIOS}
Variants:      ${VARIANTS}
Plots:         ${N_PLOTS}
Batch size:    ${BATCH_SIZE} plots / task
Batches each:  ${N_BATCHES}
Draws / plot:  ${N_DRAWS}
Task wall:     ${TASK_WALL}
Output base:   ${BASE_OUTPUT_DIR}
DRY_RUN:       ${DRY_RUN}
EOF

# ---------------------------- Per-scenario submit ----------------------------
ARRAY_IDS=()
for scen in ${SCENARIOS}; do
    case "${scen}" in
        harvest_noclimate)   RCP="none"; HARV="--harvest" ;;
        harvest_rcp45)       RCP="4.5";  HARV="--harvest" ;;
        harvest_rcp85)       RCP="8.5";  HARV="--harvest" ;;
        noharvest_noclimate) RCP="none"; HARV="" ;;
        noharvest_rcp45)     RCP="4.5";  HARV="" ;;
        noharvest_rcp85)     RCP="8.5";  HARV="" ;;
        *) echo "ERROR: unknown scenario ${scen}"; exit 1 ;;
    esac

    JOB_NAME="hu_${scen}"
    SBATCH_BODY="$(cat <<SBATCH
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --array=1-${N_BATCHES}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${TASK_CPUS}
#SBATCH --mem=${TASK_MEM}
#SBATCH --time=${TASK_WALL}
#SBATCH --output=${LOG_DIR}/${scen}_%a.out
#SBATCH --error=${LOG_DIR}/${scen}_%a.err

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

echo "=== HU ${scen} batch \${SLURM_ARRAY_TASK_ID} ==="
echo "Start: \$(date)"
echo "Node:  \$(hostname)"

# Build args as a bash array so empty HARV does not become a positional arg
ARGS=(--batch-id "\${SLURM_ARRAY_TASK_ID}" \
      --batch-size ${BATCH_SIZE} \
      --n-draws ${N_DRAWS} \
      --variants ${VARIANTS} \
      --rcp ${RCP} \
      --output-dir "${BASE_OUTPUT_DIR}" \
      --perseus-csv "${PERSEUS_CSV}" \
      --si-ratio-csv "${SI_RATIO_CSV}" \
      --seed ${SEED} \
      --start-year ${START_YEAR} \
      --end-year ${END_YEAR})
if [ -n "${HARV}" ]; then ARGS+=(--harvest); fi
python3 "${SCRIPT}" "\${ARGS[@]}"
PY_EXIT=\$?
echo "python3 exited \${PY_EXIT}"
[ "\${PY_EXIT}" -eq 0 ] || exit \${PY_EXIT}

echo "End: \$(date)"
SBATCH
)"

    if [ "${DRY_RUN}" = "1" ]; then
        echo
        echo "----- DRY RUN sbatch for ${scen} -----"
        echo "${SBATCH_BODY}" | head -20
        echo "..."
        continue
    fi

    JID=$(echo "${SBATCH_BODY}" | sbatch --parsable)
    echo "${scen}: array job ${JID} (${N_BATCHES} tasks)"
    ARRAY_IDS+=("${JID}")
done

if [ "${DRY_RUN}" = "1" ]; then
    echo "DRY_RUN: skipping aggregator submit."
    exit 0
fi

# ---------------------------- Per-scenario aggregator ------------------------
DEP_LIST=$(IFS=:; echo "${ARRAY_IDS[*]}")
AGG_ID=$(sbatch --parsable --dependency=afterok:${DEP_LIST} <<SBATCH
#!/bin/bash
#SBATCH --job-name=hu_agg
#SBATCH --account=${SLURM_ACCOUNT}
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=${AGG_MEM}
#SBATCH --time=${AGG_WALL}
#SBATCH --output=${LOG_DIR}/aggregate.out
#SBATCH --error=${LOG_DIR}/aggregate.err

module purge
module load python/3.12

export FVS_PROJECT_ROOT="${PROJECT_DIR}"

echo "=== HU factorial aggregation ==="
echo "Start: \$(date)"

cd "${PROJECT_DIR}"

for scen in ${SCENARIOS}; do
    SCEN_DIR="${BASE_OUTPUT_DIR}/\${scen}"
    [ -d "\${SCEN_DIR}" ] || { echo "skip \${scen} (missing dir)"; continue; }
    echo
    echo "--- aggregating \${scen} ---"
    POINT_FILES="\$(ls \${SCEN_DIR}/perseus_\${scen}_unc_point_batch*.csv 2>/dev/null | wc -l)"
    DRAW_FILES="\$(ls \${SCEN_DIR}/perseus_\${scen}_unc_draws_batch*.csv 2>/dev/null | wc -l)"
    echo "  point batches: \${POINT_FILES}; draw batches: \${DRAW_FILES}"
    python3 "${PROJECT_DIR}/calibration/python/perseus_uncertainty_aggregate.py" \\
        --output-dir "\${SCEN_DIR}" \\
        --perseus-csv "${PERSEUS_CSV}" || true
done

echo
echo "End: \$(date)"
SBATCH
)
echo "Aggregator: ${AGG_ID} (depends on ${DEP_LIST})"
echo
echo "Monitor: squeue -u \$USER --name=hu_*"
