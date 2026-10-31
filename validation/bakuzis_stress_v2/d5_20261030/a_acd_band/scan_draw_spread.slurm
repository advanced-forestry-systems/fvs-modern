#!/bin/bash
#SBATCH --job-name=d4a_scan --account=PUOM0008 --time=00:20:00 --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=16G
#SBATCH --output=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/diagA/logs/scan_%j.out
source /etc/profile.d/lmod.sh; module load gcc/12.3.0 R/4.4.0
J=/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027; cd $J/diagA
Rscript -e '
suppressMessages(library(data.table))
fs <- Sys.glob("'$J'/ens_full4/bakuzis_uncertainty_*_n100_batch*.csv")
d <- rbindlist(lapply(fs, fread), fill = TRUE)
cat("rows", nrow(d), "configs", paste(unique(d$config), collapse=","), "\n")
y <- d[year == max(year), .(variant, scenario, replicate, config, draw_id, atba, tpa)]
s <- y[config == "posterior", .(n = .N, uBA = uniqueN(atba), uTPA = uniqueN(tpa), sdBA = sd(atba)), by = .(variant, scenario, replicate)]
m <- y[config == "calibrated_map", .(variant, scenario, replicate, mapBA = atba)]; df <- y[config == "default", .(variant, scenario, replicate, defBA = atba)]
s <- merge(merge(s, m, by = c("variant","scenario","replicate")), df, by = c("variant","scenario","replicate"))
print(s[, .(stands = .N, ident_BA = sum(uBA == 1), ident_BA_TPA = sum(uBA == 1 & uTPA == 1), map_eq_default = sum(mapBA == defBA)), by = variant])
fwrite(s, "yr100_draw_spread_by_stand.csv")
cat("cross-tab ACD: identical draws vs map==default\n"); print(s[variant == "ACD", table(ident = uBA == 1, map_eq_def = mapBA == defBA)])
# when do draws diverge? per stand first year with >1 unique BA
f <- d[config == "posterior", .(u = uniqueN(atba)), by = .(variant, scenario, replicate, year)][u > 1, .(first_div = min(year)), by = .(variant, scenario, replicate)]
print(merge(s[, .(variant, scenario, replicate, uBA)], f, all.x = TRUE)[, .(n = .N, never_div = sum(is.na(first_div))), by = .(variant)])
print(d[variant == "ACD" & config == "posterior", .N, by = .(species_group)])
z <- s[variant == "ACD"][, ident := uBA == 1]; print(merge(z, unique(d[variant=="ACD", .(scenario, species_group, site_class, density_class)]), by = "scenario")[, .(n=.N, ident=sum(ident)), by = .(species_group)])
'
echo SCAN_DONE

