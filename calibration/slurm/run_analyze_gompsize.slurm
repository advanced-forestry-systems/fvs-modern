#!/bin/bash
#SBATCH --job-name=analyze_gompsize
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=00:15:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/analyze_gompsize_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/analyze_gompsize_%j.err
module load python/3.11
cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 analyze_gompsize.py > /fs/scratch/PUOM0008/crsfaaron/fvs-modern/scorecard_conus_cohort_gompsize/ANALYSIS.txt 2>&1
echo "ANALYSIS_DONE"
