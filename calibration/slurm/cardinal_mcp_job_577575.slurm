#!/bin/bash
#SBATCH --job-name=stress_test
#SBATCH --account=PUOM0008
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=0:30:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/stress_test_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/logs/stress_test_%j.err

module load gcc/12.3.0
module load R/4.4.0

Rscript /users/PUOM0008/crsfaaron/fvs-modern/calibration/stress_test_components.R

