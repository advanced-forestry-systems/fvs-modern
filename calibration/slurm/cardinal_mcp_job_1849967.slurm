#!/bin/bash
#SBATCH --job-name=conus_5model_ne
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_%j.err

module load gcc/12.3.0
module load python/3.11

cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 run_5model_scorecard.py \
  --variant ne \
  --max-pairs 5000 \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA \
  --outdir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus

EXIT=$?
echo "[conus_ne] exit=${EXIT}"
exit $EXIT

