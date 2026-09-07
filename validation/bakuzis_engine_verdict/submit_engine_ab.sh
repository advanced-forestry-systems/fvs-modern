#!/bin/bash
#SBATCH --job-name=eng_bakab
#SBATCH --account=PUOM0008
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/wt-engine/logs/engbak_%x_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/wt-engine/logs/engbak_%x_%j.err

# Engine A/B + Bakuzis: run the 36-scenario 100yr Bakuzis harness with a
# specified gompmort lib (GOMPIT or EXP), gompit mortality enabled.
# Usage: sbatch --job-name=<gompit|exp> submit_engine_ab.sh <libdir> <tag>
set -x
W=/fs/scratch/PUOM0008/crsfaaron/wt-engine
LIBVARIANT_DIR="$1"     # e.g. $W/lib_gompit
TAG="$2"                # gompit | exp
export FVS_PROJECT_ROOT="$HOME/fvs-modern"      # fvs2py + config live here (stable)
export FVS_LIB_DIR="$LIBVARIANT_DIR"            # <-- the A/B toggle
export FVS_CONFIG_DIR="$W/config"
export FVS_GOMPIT=1
export FVS_GOMPIT_COEF="$W/config/greg_mortality_coefficients.csv"
export BAKUZIS_OUTPUT_DIR="$W/ab_engine/bakuzis_${TAG}"
export PYTHONPATH="$W:$FVS_PROJECT_ROOT:$FVS_PROJECT_ROOT/deployment/fvs2py:$FVS_PROJECT_ROOT/deployment/microfvs:${PYTHONPATH:-}"
export PYTHONNOUSERSITE=1
module load gcc/12.3.0
mkdir -p "$BAKUZIS_OUTPUT_DIR" "$W/logs"
echo "=== ENGINE BAKUZIS/AB  tag=$TAG lib=$FVS_LIB_DIR  $(date) ==="
ls -la "$FVS_LIB_DIR/FVSne.so"
"$W/fvsvenv/bin/python" "$W/calibration/python/bakuzis_100yr_comparison.py" \
    --variant ne --seed 42 --output-dir "$BAKUZIS_OUTPUT_DIR"
echo "=== DONE $(date) rc=$? ==="
ls -la "$BAKUZIS_OUTPUT_DIR"
