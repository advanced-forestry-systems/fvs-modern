## ens3_numbers.R (2026-10-30): every number quoted in Paper 1 Section 4.6 rebuild, from the null-TPA-fixed ensemble.
suppressMessages(library(data.table))
for (tag in c("ens_full", "ens_full3")) {
  cat("\n######", tag, "\n")
  sg <- fread(sprintf("%s_seedgate.csv", tag)); W <- fread(sprintf("%s_agg/bakuzis_benchmark_wide.csv", tag)); L <- fread(sprintf("%s_agg/bakuzis_uncertainty_summary_long.csv", tag))
  cat("seedgate rows", nrow(sg), "inbin", sum(sg$inbin), sprintf("(%.1f%%)", 100*mean(sg$inbin)), "min tpa", min(sg$tpa), "\n")
  print(sg[, .(med_atba = median(atba)), by = .(variant, density_class)][order(variant, density_class)])
  print(sg[inbin == FALSE, .N, by = .(variant, density_class, below = atba < lo)][order(variant, density_class)])
  ba <- W[variable == "atba"]
  cat("populated cells:", ba[, .N, by = variant][, paste(variant, N)], "\n")
  ## default vs calibrated MAP divergence at year 100 (diff_pct_y100 is (cal - def) / def x 100)
  print(ba[, .(med_abs_diff = median(abs(diff_pct_y100), na.rm = TRUE), max_abs_diff = max(abs(diff_pct_y100), na.rm = TRUE), n_within1 = sum(abs(diff_pct_y100) <= 1, na.rm = TRUE), n = .N, identical = sum(def_median_y100 == cal_median_y100, na.rm = TRUE)), by = variant])
  print(ba[order(diff_pct_y100)][1:3, .(variant, species_group, site_class, density_class, def = def_median_y100, cal = cal_median_y100, diff_pct_y100)])
  ## posterior band at year 100 as percent of posterior median
  print(ba[, .(band_pct_med = median(100 * post_band_y100 / post_median_y100, na.rm = TRUE), band_pct_q25 = quantile(100 * post_band_y100 / post_median_y100, .25, na.rm = TRUE), band_pct_q75 = quantile(100 * post_band_y100 / post_median_y100, .75, na.rm = TRUE)), by = variant])
  print(ba[, .(band25 = median(100 * post_band_y025 / post_median_y025, na.rm = TRUE), band50 = median(100 * post_band_y050 / post_median_y050, na.rm = TRUE), band100 = median(100 * post_band_y100 / post_median_y100, na.rm = TRUE)), by = variant])
  ## stems retained at year 100 under default
  tp <- L[variable == "tpa" & config == "default", .(y0 = median[year == min(year)], y100 = median[year == max(year)]), by = .(variant, scenario)]
  print(tp[, .(med_ret100 = median(y100 / y0)), by = variant])
  ## cell-median BA maximum over 100 years (ft2/ac and m2/ha)
  bm <- L[variable == "atba" & config %in% c("default", "calibrated_map", "posterior_band")]
  print(bm[, .(max_ft2ac = max(median), max_m2ha = round(max(median) * 0.229568, 1)), by = .(variant, config)])
}
