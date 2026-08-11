#!/bin/bash
#SBATCH --job-name=heldout_dd
#SBATCH --time=05:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/heldout_dd_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/heldout_dd_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export VARS=ne,sn,kt,pn
export NSAMP=400
export SEED=5
export OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/held_out_density_dependent_20260618.csv
python3 held_out_validation.py
