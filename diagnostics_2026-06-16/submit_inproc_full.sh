#!/bin/bash
#SBATCH --job-name=inproc_full
#SBATCH --time=00:15:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/inproc_full_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/inproc_full_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
LIB=/users/PUOM0008/crsfaaron/fvs-modern/lib-test/FVSne.so python3 test_inproc_full.py
