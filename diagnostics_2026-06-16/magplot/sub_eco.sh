#!/bin/bash
#SBATCH --account=PUOM0008
#SBATCH --job-name=ak_eco2
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/akwork/ak_eco2_%j.out
cd /fs/scratch/PUOM0008/crsfaaron/akwork
export TMPDIR=/tmp
python3 ak_eco_par2.py
