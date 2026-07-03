#!/bin/bash
#SBATCH --job-name=sf_extract_bench
#SBATCH --time=01:30:00
#SBATCH --mem=96G
#SBATCH --cpus-per-task=4
#SBATCH --account=PUOM0008
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/sf_extract_bench_%j.out
#SBATCH --error=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/logs/sf_extract_bench_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=aaron.weiskittel@maine.edu
# Args: $1=component(hcb|htdbh|cr) $2=fit_jobid $3=scratch_subdir(hcb|ht_dbh|cr)
#       $4=fit_prefix(e.g. hcb_sf_100k_v2split) $5=bundle_prefix(e.g. hcb_v2split)
# Recovers fit from scratch, extracts the SF bundle (61b), checks convergence,
# then runs the three-way benchmark. DEV/TEST only.
set -uo pipefail
module load gcc/12.3.0 R/4.4.0
cd /users/PUOM0008/crsfaaron/fvs-modern
COMP=$1; JOB=$2; SUB=$3; FITPRE=$4; BPRE=$5
S=/fs/scratch/PUOM0008/crsfaaron/fvs-conus/output/conus/$SUB/sf_v2split
H=/users/PUOM0008/crsfaaron/fvs-conus/output/conus/$SUB/sf_v2split
FIT=$S/${FITPRE}_${JOB}_fit.rds; META=$S/${FITPRE}_${JOB}_meta.rds
echo "Component=$COMP fit=$FIT"
mkdir -p "$H"
cp -v "$META" "$H/" 2>/dev/null || true
ln -sf "$FIT" "$H/${FITPRE}_${JOB}_fit.rds" 2>/dev/null || true

echo "=== convergence check ==="
Rscript -e 'm<-readRDS(commandArgs(TRUE)[1]); r<-m$summary$rhat; cat("MAXRHAT",round(max(r,na.rm=TRUE),3)," FRAC_GT_1.05",round(mean(r>1.05,na.rm=TRUE),2),"\n"); if(max(r,na.rm=TRUE)>1.05) cat("WARNING: NOT CONVERGED — benchmark numbers for this component are unreliable.\n")' "$META" 2>&1

echo "=== 61b extraction ==="
Rscript --vanilla calibration/R/61b_extract_speciesfree_summaries.R \
  --component=$COMP --fit="$FIT" --meta="$META" \
  --traits=calibration/traits/species_traits_v2.rds \
  --outdir=calibration/output/conus/sf_integration --outname=$BPRE 2>&1

echo "=== benchmark ==="
Rscript --vanilla calibration/dev_sf_integration/benchmark_sf_vs_legA.R \
  --component=$COMP \
  --bundle_dir=calibration/output/conus/sf_integration --bundle_prefix=$BPRE \
  --data=calibration/data/conus_remeasurement_pairs_metric_cond_v2.rds \
  --traits=calibration/traits/species_traits_v2.rds \
  --variant_json=config/calibrated/ne.json --variety_split=TRUE --test_n=200000 \
  --outdir=calibration/output/conus/sf_integration/benchmark 2>&1
