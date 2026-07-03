#!/bin/bash
#SBATCH --job-name=htdbh_cspi4_100k
#SBATCH --time=36:00:00
#SBATCH --mem=128G
#SBATCH --cpus-per-task=8
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/htdbh_cspi4_100k_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/htdbh_cspi4_100k_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
# HT-DBH unified fit with CSPI v4 swapped in for climate_si as the productivity
# covariate. Same data/seed/subsample/--compact as the climate_si unified fit
# (10302999) so LOO is directly comparable. DEV path only.
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
SCRATCH=/fs/scratch/PUOM0008/crsfaaron/fvs-conus
OUTDIR=$SCRATCH/output/conus/htdbh/unified_cspi4_v2split
mkdir -p "$OUTDIR"
cd /users/PUOM0008/crsfaaron/fvs-modern
Rscript --vanilla calibration/R/36_fit_htdbh_speciesfree.R \
  --stan_file=calibration/stan/htdbh_unified.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --cspi_file=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/sf_integration/cspi_v4_at_calib_plots.csv \
  --outdir=$OUTDIR --outname=htdbh_unified_cspi4_100k_${SLURM_JOB_ID} \
  --subsample=100000 --compact
HOME_OUT=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/htdbh/unified_cspi4_v2split
mkdir -p "$HOME_OUT"
cp -v $OUTDIR/htdbh_unified_cspi4_100k_${SLURM_JOB_ID}_summary.csv "$HOME_OUT/" 2>/dev/null
cp -v $OUTDIR/htdbh_unified_cspi4_100k_${SLURM_JOB_ID}_meta.rds    "$HOME_OUT/" 2>/dev/null
cp -v $OUTDIR/htdbh_unified_cspi4_100k_${SLURM_JOB_ID}_loo.rds     "$HOME_OUT/" 2>/dev/null
echo "DONE htdbh_cspi4"
