#!/bin/bash
#SBATCH --job-name=conus_ne_pcfix
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_pcfix_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_pcfix_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0
module load python/3.11

echo "=== CONUS NE PC-fix scorecard started: $(date) ==="

cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python run_5model_scorecard.py \
  --variant ne \
  --max-pairs 5000 \
  --output-dir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_pcfix \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA

echo "=== CONUS NE PC-fix scorecard complete: $(date) ==="
