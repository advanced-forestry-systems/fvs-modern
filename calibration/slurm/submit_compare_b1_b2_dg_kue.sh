#!/bin/bash
#SBATCH --job-name=fvs_cmp_b1_b2_dg_kue
#SBATCH --time=00:30:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

set -e
module load gcc/12.3.0 R/4.4.0
cd $HOME/fvs-modern
mkdir -p logs MEMORY
echo "Host: $(hostname)"
echo "Started: $(date)"
Rscript --vanilla calibration/R/compare_b1_b2_dg_kue.R
echo "Finished: $(date)"
