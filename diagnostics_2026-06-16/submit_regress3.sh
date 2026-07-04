#!/bin/bash
#SBATCH --job-name=fvs_regress3
#SBATCH --time=01:30:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/regress3_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/regress3_%j.err
module load gcc/12.3.0
cd /users/PUOM0008/crsfaaron/fvs-modern
echo "=== REGRESSION with ABSOLUTE paths $(date) ==="
bash deployment/scripts/run_regression_tests.sh /users/PUOM0008/crsfaaron/fvs-modern/lib /users/PUOM0008/crsfaaron/fvs-modern/src-converted/tests 2>&1 | tail -90
echo "=== DONE $(date) ==="
