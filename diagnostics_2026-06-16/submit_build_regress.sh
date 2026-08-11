#!/bin/bash
#SBATCH --job-name=fvs_build_regress
#SBATCH --time=05:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/build_regress_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/build_regress_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "=== FULL BUILD all variants START $(date) ==="
env FC=gfortran CC=gcc bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib 2>&1 | tail -50
echo "=== BUILT .so count ==="; ls lib/*.so 2>/dev/null | wc -l
echo "=== BUILT executables ==="; ls lib/FVS* 2>/dev/null | grep -vE '\.so|\.bak|\.mod' | wc -l
echo "=== REGRESSION SUITE START $(date) ==="
bash deployment/scripts/run_regression_tests.sh ./lib src-converted/tests 2>&1 | tail -80
echo "=== ALL DONE $(date) ==="
