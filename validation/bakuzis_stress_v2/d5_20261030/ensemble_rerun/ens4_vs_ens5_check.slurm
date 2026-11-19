#!/bin/bash
#SBATCH --job-name=ens45cmp --account=PUOM0008 --time=00:20:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/fix_b_20261030/logs/ens45cmp_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
cd /fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027
Rscript -e '
suppressMessages(library(data.table)); options(width = 220)
for (b in c("", "b")) {
  a <- fread(sprintf("score_ens_full4%s/ens_full4_20261030%s_scorecard_summary.csv", b, b)); c <- fread(sprintf("score_ens_full5%s/ens_full5_20261030%s_scorecard_summary.csv", b, b))
  m <- merge(a, c, by = c("arm", "relation"), suffixes = c("_4", "_5"))
  cat("== scorer", b, ": rows", nrow(m), " differing rows:\n"); print(m[evaluated_4 != evaluated_5 | flagged_4 != flagged_5 | verdict_4 != verdict_5, .(arm, relation, evaluated_4, evaluated_5, flagged_4, flagged_5, verdict_4, verdict_5)])
}
for (tg in c("ens_full4", "ens_full5")) {
  d <- rbindlist(lapply(Sys.glob(sprintf("%s/bakuzis_uncertainty_*_n100_batch*.csv", tg)), fread), fill = TRUE)
  y <- d[year == max(year) & config %in% c("default", "calibrated_map"), .(ba = median(atba)), by = .(variant, scenario, species_group, site_class, density_class, config)]
  w <- dcast(y, variant + scenario + species_group + site_class + density_class ~ config, value.var = "ba"); w[, div := 100 * (calibrated_map - default) / default]
  cat(tg, ": NE cells MAP below default", w[variant == "NE", sum(calibrated_map < default)], "of", w[variant == "NE", .N], "; NE |div| median", round(w[variant == "NE", median(abs(div))], 1), "max", round(w[variant == "NE", max(abs(div))], 1),
      "at", w[variant == "NE"][which.max(abs(div)), paste(species_group, site_class, density_class, calibrated_map, "vs", default)], "; ACD identical cells", w[variant == "ACD", sum(calibrated_map == default)], "max |div|", round(w[variant == "ACD", max(abs(div))], 2), "\n")
  s <- fread(sprintf("%s_seedgate.csv", tg))
  print(s[, .(med_ba0 = median(atba), inbin = sum(inbin), n = .N), by = .(variant, density_class)][order(density_class, variant)])
}'

