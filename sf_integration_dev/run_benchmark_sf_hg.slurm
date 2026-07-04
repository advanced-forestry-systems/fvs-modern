#!/bin/bash
#SBATCH --job-name=bench_sf_hg
#SBATCH --time=00:50:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/bench_sf_hg_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/bench_sf_hg_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu

set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern

# HG benchmark against the hg_v5_prod bundle (v1 traits, no variety split, to
# match the fit that produced the bundle). DEV/TEST only.
Rscript --vanilla calibration/dev_sf_integration/benchmark_sf_vs_legA.R \
  --component=hg \
  --bundle_dir=calibration/output/conus/sf_integration \
  --bundle_prefix=hg_v5_prod \
  --data=calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds \
  --traits=calibration/traits/species_traits.rds \
  --variant_json=config/calibrated/ne.json \
  --variety_split=FALSE \
  --test_n=200000 \
  --outdir=calibration/output/conus/sf_integration/benchmark
