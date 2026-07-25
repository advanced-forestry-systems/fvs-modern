#!/bin/bash
#SBATCH --job-name=conus_ne_pcfix_v2
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=6:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_pcfix_v2_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/conus_ne_pcfix_v2_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

# Module setup copied verbatim from the known-good scorecard job 12306461
module load gcc/12.3.0
module load python/3.11

echo "=== CONUS NE PC-fix scorecard v2 started: $(date) ==="

cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 run_5model_scorecard.py \
  --variant ne \
  --max-pairs 5000 \
  --outdir /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_pcfix \
  --fia-dir /fs/scratch/PUOM0008/crsfaaron/FIA

EXIT=$?
echo "[conus_ne_pcfix] exit=${EXIT}"
echo "=== CONUS NE PC-fix scorecard v2 complete: $(date) ==="
exit $EXIT
