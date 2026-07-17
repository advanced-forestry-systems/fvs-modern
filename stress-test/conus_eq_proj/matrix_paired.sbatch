#!/bin/bash
#SBATCH --job-name=matrix_paired
#SBATCH --account=PUOM0008
#SBATCH --time=02:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=%x_%j.out
module load gcc/12.3.0 gdal/3.7.3 geos/3.12.0 proj/9.2.1 R/4.4.0
cd $SLURM_SUBMIT_DIR
Rscript matrix_paired.R
