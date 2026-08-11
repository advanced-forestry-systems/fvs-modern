#!/bin/bash
#SBATCH --job-name=inspect
#SBATCH --account=PUOM0008
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect3_%j.out
#SBATCH --error=${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/logs/inspect3_%j.err

module load gcc/12.3.0
module load R/4.4.0
export R_LIBS_USER="/path/to/user/path"

Rscript ${FVS_PROJECT_ROOT:-/path/to/fvs-modern}/calibration/slurm/inspect_script.R
