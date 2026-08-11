#!/bin/bash
#SBATCH --job-name=fourarm_engine
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_engine_%j.out
#SBATCH --error=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_engine_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
cd /users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16
export VARS=ne,sn,kt,pn,cr,ut,nc,ec,wc; export NSAMP=400; export SEED=5
export OUT=/users/PUOM0008/crsfaaron/fvs-modern/diagnostics_2026-06-16/fourarm_engine_20260618.csv
python3 fourarm_engine.py
