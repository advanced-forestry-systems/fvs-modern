#!/bin/bash
#SBATCH --job-name=heldout_smoke
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=2
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/heldout_smoke_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/heldout_smoke_%j.err

cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export VARS=ne
export NSAMP=120
export SEED=5
export OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/held_out_dd_smoke.csv
python3 held_out_validation.py
