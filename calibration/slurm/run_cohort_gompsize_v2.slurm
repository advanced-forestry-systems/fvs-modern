#!/bin/bash
#SBATCH --job-name=score_gs_v2
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/score_gs_v2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/score_gs_v2_%j.err
module load gcc/12.3.0
module load python/3.11 2>/dev/null || true
cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
echo "=== gompsize v2 start $(date) ==="
python3 run_cohort_gompsize_v2.py --variant ne --max-pairs 5000 \
  --outdir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_cohort_gompsize_v2 \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA
echo "=== done $(date) ==="
