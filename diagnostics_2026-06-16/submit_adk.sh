#!/bin/bash
#SBATCH --job-name=adk_build_calib
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_build_calib_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_build_calib_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "=== BUILD FVSadk $(date) ==="
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib adk 2>&1 | tail -20
echo "=== adk artifacts ==="; ls -la lib/FVSadk lib/FVSadk.so 2>/dev/null
echo "=== CALIBRATE adk on NY FIA $(date) ==="
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
VAR=adk STATES=NY NSAMP=150 SEED=5 OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_baimult_calib.json python3 calib_ne.py 2>&1 | tail -30
echo "=== DONE $(date) ==="
