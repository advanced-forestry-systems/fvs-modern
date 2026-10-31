#!/bin/bash
#SBATCH --job-name=diagB --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/diagB/logs/b_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
export FIA_DATA_DIR=$HOME/fia_data
python3 /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/diagB/b_expansion.py
