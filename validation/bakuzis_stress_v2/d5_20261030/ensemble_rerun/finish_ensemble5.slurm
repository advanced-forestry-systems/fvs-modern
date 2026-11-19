#!/bin/bash
#SBATCH --job-name=finish5 --account=PUOM0008 --time=03:00:00 --nodes=1 --ntasks=1 --cpus-per-task=2 --mem=48G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/finish5_%j.out
# Seed gate, aggregation, scoring, Paper 1 figures and Section 4.6 numbers for ens_full5 (D5), mirroring finish_ensemble4.slurm.
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
W=/fs/scratch/PUOM0008/crsfaaron/wt-fix2-20261030; J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027
export PYTHONPATH=$W:$W/deployment/fvs2py:${PYTHONPATH:-}; export MPLCONFIGDIR=/tmp/mpl_$$
cd $J
ls ens_full5/bakuzis_uncertainty_*_n100_batch0*.csv | wc -l
echo "== seed gate $(date)"; Rscript seedgate_full5.R
echo "== aggregate $(date)"; mkdir -p ens_full5_agg
python3 $W/calibration/python/bakuzis_uncertainty_aggregate.py --input-dir $J/ens_full5 --output-dir $J/ens_full5_agg 2>&1 | tail -12
echo "== score (variant thresholds) $(date)"; mkdir -p score_ens_full5
Rscript bakuzis_stress_v2.R --input=ens_full5_agg/bakuzis_uncertainty_summary_long.csv --source=uncert --out=score_ens_full5 --sdimax=ens_sdimax.csv --obs_ba_file=ens_obs_ba.csv --label=ens_full5_20261030 > score_ens_full5/run_report.txt 2>&1
grep -E '^Reineke axis|Error|Wrote' score_ens_full5/run_report.txt
echo "== rescore with species-group thresholds $(date)"; mkdir -p score_ens_full5b
Rscript bakuzis_stress_v2.R --input=ens_full5_agg/bakuzis_uncertainty_summary_long.csv --source=uncert --out=score_ens_full5b --sdimax=ens_sdimax_spgrp.csv --obs_ba_file=ens_obs_ba_spgrp.csv --label=ens_full5_20261030b > score_ens_full5b/run_report.txt 2>&1
grep -E '^Reineke axis|Error|Wrote' score_ens_full5b/run_report.txt
echo "== paper 1 figures $(date)"; mkdir -p paper1_figs_20261030c
python3 $W/calibration/python/bakuzis_uncertainty_figures.py --summary ens_full5_agg/bakuzis_uncertainty_summary_long.csv --bench ens_full5_agg/bakuzis_benchmark_wide.csv --laws ens_full5_agg/bakuzis_laws_compliance.csv --cell-bands ens_full5_agg/bakuzis_uncertainty_cell_bands.csv --outdir paper1_figs_20261030c 2>&1 | tail -12
echo "== Section 4.6 numbers $(date)"
Rscript ens5_numbers.R ens_full4 ens_full5 > ens5_numbers.out 2>&1; tail -3 ens5_numbers.out
echo "== within-stand ACD band check (part A) and year-0 BA vs BALIVE (part B)"
Rscript -e 'suppressMessages(library(data.table)); for (tg in c("ens_full4","ens_full5")) { d <- rbindlist(lapply(Sys.glob(sprintf("%s/bakuzis_uncertainty_*_n100_batch*.csv", tg)), fread), fill=TRUE); y <- d[year == max(year) & config == "posterior", .(u = uniqueN(atba)), by = .(variant, scenario, replicate)]; cat(tg, ": stands with identical year-100 BA across draws:", paste(y[, .(s = paste0(variant, " ", sum(u == 1), "/", .N)), by = variant]$s, collapse = "; "), "\n"); y0 <- d[config == "default", .SD[which.min(year)], by = .(variant, scenario, replicate)]; fwrite(y0[, .(variant, scenario, replicate, ba0 = atba, tpa0 = tpa)], sprintf("%s_year0_default.csv", tg)) }'
Rscript -e 'suppressMessages(library(data.table)); a <- fread("ens_full4_year0_default.csv"); b <- fread("ens_full5_year0_default.csv"); m <- merge(a, b, by = c("variant","scenario","replicate"), suffixes = c("_4","_5")); m[, ch := 100 * (ba0_5 / ba0_4 - 1)]; cat("year-0 BA change ens_full4 -> ens_full5: stands", nrow(m), "changed", m[ch != 0, .N], "moved >2%", m[abs(ch) > 2, .N], "max |%|", round(max(abs(m$ch)), 2), "\n"); fwrite(m, "ens5_year0_change.csv")'
echo FINISH5_DONE
