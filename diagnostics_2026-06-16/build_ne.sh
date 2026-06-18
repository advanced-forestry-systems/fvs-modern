#!/bin/bash
#SBATCH --job-name=build_ne_so
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/build_ne_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/build_ne_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "=== BUILD START $(date) ==="
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib-test ne 2>&1 | tail -25
echo "=== BUILT SO ==="; ls -la /users/PUOM0008/crsfaaron/fvs-modern/lib-test/FVSne.so 2>/dev/null
echo "=== IN-PROCESS TEST against fresh .so ==="
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
LIB=/users/PUOM0008/crsfaaron/fvs-modern/lib-test/FVSne.so python3 test_treeattr_injection.py 2>&1 | tail -20
echo "=== JOB DONE $(date) ==="
