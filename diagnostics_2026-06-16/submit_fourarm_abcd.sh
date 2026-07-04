#!/bin/bash
#SBATCH --job-name=fourarm_abcd
#SBATCH --time=10:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_abcd_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_abcd_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export NSAMP=600 SEED=5 OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_abcd_20260618.csv
python3 fourarm_abcd.py
