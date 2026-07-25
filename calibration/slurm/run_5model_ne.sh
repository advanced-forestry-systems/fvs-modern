#!/usr/bin/env bash
#SBATCH --job-name=5model_ne
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/5model_ne_%j.log

# Systems 3 (greg_dg) and 5 (greg_full) are NE-only; skipped automatically for ne.
# Systems 1, 2, 4 run on all variants.

module load gcc/12.3.0

echo "[5model_ne] starting at $(date)"
echo "[5model_ne] host: $(hostname)"

/usr/bin/python3 /users/PUOM0008/crsfaaron/fvs-modern/calibration/run_5model_scorecard.py \
    --variant ne \
    --max-pairs 200

echo "[5model_ne] finished at $(date)"
