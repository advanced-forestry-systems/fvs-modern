#!/bin/bash
#SBATCH --job-name=adk_calib
#SBATCH --time=01:00:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_calib_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_calib_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
# ADK = Acadian engine recalibrated to NY (Adirondack) FIA growth; uses the working FVSacd executable
VAR=acd STATES=NY NSAMP=150 SEED=5 OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_calibration.json python3 calib_ne.py 2>&1 | tail -30
echo "=== DONE $(date) ==="
