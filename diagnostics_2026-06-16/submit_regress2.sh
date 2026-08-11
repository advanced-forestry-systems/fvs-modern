#!/bin/bash
#SBATCH --job-name=fvs_regress2
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/regress2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/regress2_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "=== REGRESSION (executables present) $(date) ==="
bash deployment/scripts/run_regression_tests.sh ./lib src-converted/tests 2>&1 | tail -90
echo "=== DONE $(date) ==="
