#!/bin/bash
#SBATCH --job-name=ba_decomp_ne
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=1:30:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/ba_decomp_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/ba_decomp_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

module load gcc/12.3.0
module load python/3.11

cd /users/PUOM0008/crsfaaron/fvs-modern/calibration
python3 run_ba_decomp.py --variant ne --max-pairs 200 \
  --outdir /users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard

