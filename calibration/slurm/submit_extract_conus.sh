#!/bin/bash
#SBATCH --job-name=extract_conus
#SBATCH --account=PUOM0008
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=00:30:00
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

# =============================================================================
# Run 61_extract_conus_summaries.R for one component (or all).
#
# Usage:
#   sbatch --export=COMPONENT=hg calibration/slurm/submit_extract_conus.sh
#   sbatch --export=COMPONENT=all calibration/slurm/submit_extract_conus.sh
#   sbatch --export=COMPONENT=dg,MODEL=dg_kuehne_cspi_traits1 ...
#
# Why a SLURM wrapper:
#   `module load` doesn't initialize correctly under `ssh ... 'cmd'` so
#   R isn't on PATH there. SLURM batch jobs get a proper login environment
#   where module loading works.
#
# Memory note: 32G should cover the largest fit (mortality at 14G) plus
# the as_draws_df() expansion. Bump --mem=64G if extraction OOMs.
# =============================================================================

set -uo pipefail
source /etc/profile.d/lmod.sh
# R/4.4.0 requires gcc/12.3.0 to be loaded first (per `module spider R/4.4.0`)
module load gcc/12.3.0
module load R/4.4.0

COMPONENT="${COMPONENT:-hg}"
MODEL="${MODEL:-}"

# FVS_PROJECT_ROOT must point AT the fvs-modern dir, not its parent.
# 61_extract_conus_summaries.R uses Sys.getenv("FVS_PROJECT_ROOT", "..").
# The default ".." is wrong when launched from inside fvs-modern, so we
# export it explicitly here.
export FVS_PROJECT_ROOT="${FVS_PROJECT_ROOT:-$HOME/fvs-modern}"
cd "$FVS_PROJECT_ROOT"
echo "=== Extraction job ==="
echo "host:     $(hostname)"
echo "time:     $(date)"
echo "component: $COMPONENT"
echo "model:    ${MODEL:-<default>}"
echo "Rscript:  $(which Rscript)"
echo "======================="

ARGS=("--component" "$COMPONENT")
if [ -n "$MODEL" ]; then
    ARGS+=("--model" "$MODEL")
fi

Rscript calibration/R/61_extract_conus_summaries.R "${ARGS[@]}"
RC=$?

echo "=== Done at $(date) (rc=$RC) ==="

if [ $RC -eq 0 ]; then
    echo "--- output CSVs (top of calibration/output/conus/) ---"
    ls -lh calibration/output/conus/*.csv 2>/dev/null | tail -20
fi

exit $RC
