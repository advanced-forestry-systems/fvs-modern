#!/bin/bash
#SBATCH --job-name=bak_proj_rescore2
#SBATCH --account=PUOM0008
#SBATCH --time=03:00:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=8
#SBATCH --mem=120G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/proj_rescore2_%j.out
# Rescore the CONUS projector arms (the Paper 2 Table 6 inputs) with bakuzis_stress_v2.R, with the FIA
# observed BA and SDI percentiles by variant as the magnitude test. 1,500 stands per variant per arm
# (seed 20261027). The metrics files carry no dominant height, so height and volume relations report
# NOT EVALUATED; what is live is the shock, QMD monotonicity, Reineke (RD-gated), frontier, SDI and BA
# magnitude, plateau and BA inversion tests.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027
B=/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj
PAIRS=/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds
cd $J; mkdir -p proj_rescore proj_rescore2
Rscript -e 'suppressMessages(library(data.table)); p<-fread("proj_rescore/fia_ba_sdi_percentiles.csv"); fwrite(p[,.(group=variant,ba_p99,ba_p999)],"proj_rescore/obs_ba_by_variant.csv")'
echo "== observed percentiles $(date)"
echo "percentiles already computed (per hectare SDI columns added by the seedcheck job)"
# SDImax lookup by variant for the scorer (group = VARIANT), from the observed p99.9 of SDI
P999=$(Rscript -e 'suppressMessages(library(data.table)); f<-"proj_rescore/fia_ba_sdi_percentiles.csv"; if(file.exists(f)) cat(round(max(fread(f)$ba_p999),1)) else cat("NA")')
P99=$(Rscript -e 'suppressMessages(library(data.table)); f<-"proj_rescore/fia_ba_sdi_percentiles.csv"; if(file.exists(f)) cat(round(max(fread(f)$ba_p99),1)) else cat("NA")')
echo "BA p99=$P99 p999=$P999 (CONUS-wide max over variants; per-variant values in the csv)"
declare -A DIRS=( [b1_native]="$B/out_conus_eq|_conus_b1_metrics.csv" [b2_native]="$B/out_conus_eq|_conus_b2_metrics.csv" [b1_gompit]="$B/out_conus_eq_b1_gompit|_conus_b1_gompit_metrics.csv" [b2_gompit]="$B/out_conus_eq_b2_gompit|_conus_b2_gompit_metrics.csv" )
for arm in b1_native b2_native b1_gompit b2_gompit; do
  IFS='|' read -r dir pat <<< "${DIRS[$arm]}"
  files=$(ls $dir/*$pat 2>/dev/null | tr '\n' ',' | sed 's/,$//')
  n=$(echo "$files" | tr ',' '\n' | grep -c .)
  echo "== $arm: $n variant files $(date)"
  [ "$n" -gt 0 ] || continue
  Rscript bakuzis_stress_v2.R --input="$files" --source=projector --out=proj_rescore2/$arm --sdimax=proj_rescore/sdimax_by_variant.csv --obs_ba_file=proj_rescore/obs_ba_by_variant.csv --max_traj_per_group=1500 --label=proj_$arm > proj_rescore2/${arm}_run_report.txt 2>&1
  grep -E '^Reineke axis|subsampled|Error|Wrote' proj_rescore2/${arm}_run_report.txt
done
echo "PROJ_RESCORE_DONE $(date)"

