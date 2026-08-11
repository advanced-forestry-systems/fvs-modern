#!/bin/bash
#SBATCH --job-name=adk_build2
#SBATCH --time=00:40:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_build2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/adk_build2_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib adk 2>&1 | tail -15
echo "=== artifacts ==="; ls -la lib/FVSadk.so 2>/dev/null && echo SO_OK || echo SO_MISSING
echo "=== load test ==="; python3 -c "import ctypes; ctypes.CDLL('lib/FVSadk.so', mode=1); print('ADK_SO_LOADS')" 2>&1 | tail -3
