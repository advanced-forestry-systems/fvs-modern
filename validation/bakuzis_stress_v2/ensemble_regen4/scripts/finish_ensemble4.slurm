#!/bin/bash
#SBATCH --job-name=bak_finish4
#SBATCH --account=PUOM0008
#SBATCH --time=01:30:00
#SBATCH --nodes=1 --ntasks=1 --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/finish4_%j.out
# Seed gate, aggregation (pooled legacy band + per-stand within-stand bands), scoring (variant and
# species-group thresholds), Paper 1 figures and Section 4.6 numbers for ens_full4 (October 30, 2026).
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027
export PYTHONPATH=$W:$W/deployment/fvs2py:${PYTHONPATH:-}; export MPLCONFIGDIR=/tmp/mpl_$$
cd $J
ls ens_full4/bakuzis_uncertainty_*_n100_batch0*.csv | wc -l
echo "== seed gate $(date)"; Rscript seedgate_full4.R
echo "== aggregate $(date)"; mkdir -p ens_full4_agg
python3 $W/calibration/python/bakuzis_uncertainty_aggregate.py --input-dir $J/ens_full4 --output-dir $J/ens_full4_agg 2>&1 | tail -12
ls -la ens_full4_agg
echo "== score (variant thresholds) $(date)"; mkdir -p score_ens_full4
Rscript bakuzis_stress_v2.R --input=ens_full4_agg/bakuzis_uncertainty_summary_long.csv --source=uncert --out=score_ens_full4 --sdimax=ens_sdimax.csv --obs_ba_file=ens_obs_ba.csv --label=ens_full4_20261030 > score_ens_full4/run_report.txt 2>&1
grep -E '^Reineke axis|Error|Wrote' score_ens_full4/run_report.txt
echo "== rescore with species-group thresholds $(date)"; mkdir -p score_ens_full4b
Rscript bakuzis_stress_v2.R --input=ens_full4_agg/bakuzis_uncertainty_summary_long.csv --source=uncert --out=score_ens_full4b --sdimax=ens_sdimax_spgrp.csv --obs_ba_file=ens_obs_ba_spgrp.csv --label=ens_full4_20261030b > score_ens_full4b/run_report.txt 2>&1
grep -E '^Reineke axis|Error|Wrote' score_ens_full4b/run_report.txt
echo "== paper 1 figures $(date)"; rm -rf paper1_figs_20261030b; mkdir -p paper1_figs_20261030b
python3 $W/calibration/python/bakuzis_uncertainty_figures.py --summary ens_full4_agg/bakuzis_uncertainty_summary_long.csv --bench ens_full4_agg/bakuzis_benchmark_wide.csv --laws ens_full4_agg/bakuzis_laws_compliance.csv --cell-bands ens_full4_agg/bakuzis_uncertainty_cell_bands.csv --outdir paper1_figs_20261030b 2>&1 | tail -12
ls paper1_figs_20261030b
echo "== Section 4.6 numbers $(date)"
Rscript ens4_numbers.R ens_full3 ens_full4 > ens4_numbers.out 2>&1; tail -3 ens4_numbers.out
cp ens4_numbers.R ens4_numbers.out ens4_numbers_old_vs_new.csv $W/calibration/python/bakuzis_uncertainty_figures.py paper1_figs_20261030b/
echo FINISH4_DONE
