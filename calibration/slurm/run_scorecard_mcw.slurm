#!/bin/bash
#SBATCH --job-name=scorecard_mcw
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=3:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/scorecard_mcw_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/scorecard_mcw_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
module load gcc/12.3.0
module load python/3.11
echo "=== MCW-active bias-fixed scorecard start $(date) ==="
cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 run_5model_scorecard_biasfix.py --variant ne --max-pairs 5000 \
  --outdir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_mcw \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA
echo "=== done $(date) ==="
