#!/bin/bash
#SBATCH --job-name=bench_sf_all
#SBATCH --time=01:40:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/bench_sf_all_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/bench_sf_all_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
B=calibration/dev_sf_integration/benchmark_sf_vs_legA.R
D=calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds
OUT=calibration/output/conus/sf_integration/benchmark
V1=calibration/traits/species_traits.rds
V2=calibration/traits/species_traits_v2.rds

run () {  # comp prefix traits split vtag region vjson
  echo "##### $1 / $5 #####"
  Rscript --vanilla $B --component=$1 --bundle_prefix=$2 --data=$D \
    --traits=$3 --variety_split=$4 --variant_tag=$5 --region=$6 \
    --variant_json=config/calibrated/$7 --test_n=200000 --outdir=$OUT
}

# Northeast variant (region = eastern ecoregions)
run hg    hg_v5_prod    $V1 FALSE ne east ne.json
run hcb   hcb_v2split   $V2 TRUE  ne east ne.json
run htdbh htdbh_v2split $V2 TRUE  ne east ne.json
# Pacific Northwest variant (region = western ecoregions)
run hg    hg_v5_prod    $V1 FALSE pn west pn.json
run hcb   hcb_v2split   $V2 TRUE  pn west pn.json
run htdbh htdbh_v2split $V2 TRUE  pn west pn.json
