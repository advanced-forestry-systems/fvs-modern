#!/bin/bash
#SBATCH --job-name=fourarm_proj
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_proj_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_proj_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_projector_NE_20260618.csv
python3 fourarm_projector.py
