#!/bin/bash
#SBATCH --account=puom0008
#SBATCH --partition=condo-ccapp-backfill
#SBATCH --job-name=biascheck
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_biascheck/biascheck_%j.log
set -e
cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 run_5model_scorecard_biasfix.py \
  --variant ne --max-pairs 500 --seed 42 \
  --outdir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_biascheck \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA
echo "SCORECARD_DONE"
