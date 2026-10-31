#!/bin/bash
#SBATCH --job-name=diagB2 --account=PUOM0008 --time=02:00:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/r1/logs/b_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0
export FIA_DATA_DIR=$HOME/fia_data
python3 /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/r1/b_expansion.py
