#!/bin/bash
#SBATCH --job-name=7cfg_em_fix
#SBATCH --account=PUOM0008
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/calibration/slurm/7cfg_em_fix_%j.log

module load gcc/12.3.0
/usr/bin/python3 /users/PUOM0008/crsfaaron/fvs-modern/calibration/run_7config_fia_scorecard.py --variant em --max-pairs 200 --out-suffix _fix
