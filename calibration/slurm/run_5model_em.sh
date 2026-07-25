#!/usr/bin/env bash
#SBATCH --job-name=5model_em
#SBATCH --account=PUOM0008
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/5model_em_%j.log

# Systems 3 (greg_dg) and 5 (greg_full) are NE-only; skipped automatically for em.
# Systems 1, 2, 4 run on all variants.

module load gcc/12.3.0

echo "[5model_em] starting at $(date)"
echo "[5model_em] host: $(hostname)"

/usr/bin/python3 /users/PUOM0008/crsfaaron/fvs-modern/calibration/run_5model_scorecard.py \
    --variant em \
    --max-pairs 200

echo "[5model_em] finished at $(date)"
