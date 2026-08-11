#!/bin/bash
#SBATCH --job-name=unified_smoke
#SBATCH --time=01:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/unified_smoke_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/unified_smoke_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
# Smoke test of the unified models (species-free hierarchy + species random
# effect). Reports sigma_sp (species signal beyond traits?) alongside the
# L1/L2/L3 ecoregion and forest-type RE scales. Runs with --compact so the
# drivers save a small LOO + summary CSV instead of the multi-GB save_object
# (that save was the cause of the earlier smoke hang).
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
O=calibration/output/conus/sf_integration/unified_smoke
V2=calibration/traits/species_traits_v2.rds
mkdir -p $O

echo "########## HCB unified smoke ##########"
Rscript --vanilla calibration/R/33_fit_hcb_speciesfree.R \
  --stan_file=calibration/stan/hcb_unified.stan --traits=$V2 \
  --subsample=40000 --smoke --compact --outdir=$O --outname=hcb_unified_smoke 2>&1 | tail -30
echo "----- HCB sigma_sp / RE scales / phi (from summary csv) -----"
awk -F, 'NR==1 || $1 ~ /^sigma|^phi/' $O/hcb_unified_smoke_summary.csv 2>&1

echo
echo "########## HT-DBH unified smoke ##########"
Rscript --vanilla calibration/R/36_fit_htdbh_speciesfree.R \
  --stan_file=calibration/stan/htdbh_unified.stan --traits=$V2 \
  --subsample=40000 --smoke --compact --outdir=$O --outname=htdbh_unified_smoke 2>&1 | tail -30
echo "----- HT-DBH sigma_sp / RE scales (from summary csv) -----"
awk -F, 'NR==1 || $1 ~ /^sigma/' $O/htdbh_unified_smoke_summary.csv 2>&1
echo "DONE"
