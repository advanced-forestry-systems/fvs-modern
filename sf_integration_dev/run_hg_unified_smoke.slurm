#!/bin/bash
#SBATCH --job-name=hg_unified_smoke
#SBATCH --time=02:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_unified_smoke_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/hg_unified_smoke_%j.err
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
O=calibration/output/conus/sf_integration/unified_smoke
mkdir -p \$O
Rscript --vanilla calibration/R/32_fit_hg_speciesfree_v5.R \
  --stan_file=calibration/stan/hg_unified_v5.stan \
  --traits=calibration/traits/species_traits_v2.rds \
  --subsample=40000 --smoke --compact --outdir=\$O --outname=hg_unified_smoke 2>&1 | tail -30
echo \"----- HG sigma_sp / RE scales -----\"
awk -F, 'NR==1 || \$1 ~ /^sigma/' \$O/hg_unified_smoke_summary.csv 2>&1
echo DONE
