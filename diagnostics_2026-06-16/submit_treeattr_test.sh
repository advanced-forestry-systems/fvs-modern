#!/bin/bash
#SBATCH --job-name=treeattr_test
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/treeattr_test_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/treeattr_test_%j.err
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
python3 test_treeattr_injection.py
