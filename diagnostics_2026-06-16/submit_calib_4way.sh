#!/bin/bash
#SBATCH --job-name=calib_4way
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/calib_4way_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/calib_4way_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export NSAMP=250 SEED=5 OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/calib_4way_20260618.csv
python3 calib_4way.py
