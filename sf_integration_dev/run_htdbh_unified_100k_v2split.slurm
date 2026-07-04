#!/bin/bash
#SBATCH --job-name=htdbh_unified_100k
#SBATCH --time=36:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/htdbh_unified_100k_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/htdbh_unified_100k_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
# UNIFIED HT-DBH production fit (species-free B1 hierarchy + species RE z_sp).
# Mirrors run_hcb_sf_100k_v2split.slurm (same data, seed 42, 100k subsample ->
# identical rows as the species-free production fit, so LOO is directly
# comparable). --compact saves the LOO object + summary (incl sigma_sp) and
# skips the multi-GB save_object. DEV path only; no production config writes.
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
SCRATCH=/fs/scratch/PUOM0008/crsfaaron/fvs-conus
OUTDIR=$SCRATCH/output/conus/htdbh/unified_v2split
mkdir -p "$OUTDIR"
cd /users/PUOM0008/crsfaaron/fvs-modern

Rscript --vanilla calibration/R/36_fit_htdbh_speciesfree.R \
  --stan_file=calibration/stan/htdbh_unified.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --outdir=$OUTDIR \
  --outname=htdbh_unified_100k_v2split_${SLURM_JOB_ID} \
  --subsample=100000 --compact

HOME_OUT=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/htdbh/unified_v2split
mkdir -p "$HOME_OUT"
cp -v $OUTDIR/htdbh_unified_100k_v2split_${SLURM_JOB_ID}_summary.csv "$HOME_OUT/" 2>/dev/null
cp -v $OUTDIR/htdbh_unified_100k_v2split_${SLURM_JOB_ID}_meta.rds    "$HOME_OUT/" 2>/dev/null
cp -v $OUTDIR/htdbh_unified_100k_v2split_${SLURM_JOB_ID}_loo.rds     "$HOME_OUT/" 2>/dev/null
echo "DONE htdbh_unified"
