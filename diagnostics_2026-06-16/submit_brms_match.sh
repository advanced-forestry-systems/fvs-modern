#!/bin/bash
#SBATCH --job-name=brms_match
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/brms_match_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/brms_match_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/brms_match_rate_20260618.csv
python3 brms_match_rate.py
