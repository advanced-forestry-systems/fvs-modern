#!/bin/bash
#SBATCH --job-name=prtexm_fix
#SBATCH --time=00:25:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/prtexm_fix_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/prtexm_fix_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib-test ne 2>&1 | tail -5
echo "=== RERUN addtrees with patched prtexm ==="
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
LIB=/users/PUOM0008/crsfaaron/fvs-modern/lib-test/FVSne.so python3 test_addtrees.py 2>&1 | tail -15
