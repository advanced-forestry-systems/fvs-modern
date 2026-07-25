#!/bin/bash
#SBATCH --job-name=fullsys_ne
#SBATCH --account=PUOM0008
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/fullsys_ne_%j.log
module load gcc/12.3.0
/usr/bin/python3 /users/PUOM0008/crsfaaron/fvs-modern/calibration/run_fullsystem_fia_scorecard.py --variant ne --max-pairs 200
