suppressMessages(library(data.table))
fs <- list.files("ens_full5", pattern = "^bakuzis_uncertainty_(ne|acd)_n100_batch0[0-5]\\.csv$", full.names = TRUE)
d <- rbindlist(lapply(fs, fread), fill = TRUE)
cat("raw rows", nrow(d), "variants", paste(unique(d$variant), collapse=","), "configs", paste(unique(d$config), collapse=","), "draws", uniqueN(d$draw_id), "\n")
bins <- data.table(density_class = c("Low","Medium","High"), lo = c(50,90,150), hi = c(90,150,240))
s <- d[config == "default", .SD[which.min(year)], by = .(variant, scenario, replicate)][, .(variant, scenario, species_group, site_class, density_class, replicate, year, tpa, atba, attopht)]
s <- merge(s, bins, by = "density_class"); s[, inbin := atba >= lo & atba <= hi]
cat("\nYEAR-0 SEED GATE (default config, each replicate vs its BALIVE selection bin)\n")
print(s[, .(replicates = .N, in_bin = sum(inbin), zero_stems = sum(tpa == 0), ba_min = min(atba), ba_med = median(atba), ba_max = max(atba)), by = .(variant, density_class)][order(variant, density_class)])
print(s[, .(cells = uniqueN(scenario), cells_populated = uniqueN(scenario[tpa > 0])), by = variant])
fwrite(s, "ens_full5_seedgate.csv")
cat("\nF5 bound: cell-median year-100 BA (ft2/ac) by variant and config:\n")
e <- d[, .SD[which.max(year)], by = .(variant, scenario, replicate, config, draw_id)][, .(ba100 = median(atba)), by = .(variant, scenario, config)]
print(e[, .(max_cell_median_ft2ac = round(max(ba100)), m2ha = round(max(ba100)*0.2295684,1)), by = .(variant, config)][order(variant, config)])
