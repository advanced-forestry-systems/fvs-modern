## ens4_numbers.R (2026-10-30): every number quoted in Paper 1 Section 4.6, old (ens_full3, as run on
## October 29) next to new (ens_full4, the D1 to D4 fixed code). Updated copy of ens3_numbers.R.
## The posterior band is reported two ways and never pooled:
##   within-stand  = parameter only: 95 percent band of a stand's 100 draws, as percent of that stand's
##                   posterior median, then the median across the cell's stands, then across cells
##   between-stand = spread of the 5 stands' calibrated MAP values in a cell (range and CV, percent of median)
## The legacy pooled band (draws x stands) is printed once more, labelled as such, for comparison only.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
tags <- if (length(args)) args else c("ens_full3", "ens_full4")
score_dir <- c(ens_full3 = "score_ens_full3b", ens_full4 = "score_ens_full4b")
score_lab <- c(ens_full3 = "ens_full3_20261029b", ens_full4 = "ens_full4_20261030b")
out_rows <- list()
put <- function(tag, metric, variant, value) out_rows[[length(out_rows) + 1]] <<- data.table(tag = tag, metric = metric, variant = variant, value = value)

for (tag in tags) {
  cat("\n######", tag, "\n")
  sg <- fread(sprintf("%s_seedgate.csv", tag)); W <- fread(sprintf("%s_agg/bakuzis_benchmark_wide.csv", tag)); L <- fread(sprintf("%s_agg/bakuzis_uncertainty_summary_long.csv", tag))
  cat("seedgate rows", nrow(sg), "inbin", sum(sg$inbin), sprintf("(%.1f%%)", 100 * mean(sg$inbin)), "min tpa", min(sg$tpa), "\n")
  put(tag, "seedgate_inbin", "ALL", sum(sg$inbin))
  print(sg[, .(med_atba = median(atba), in_bin = sum(inbin), n = .N), by = .(variant, density_class)][order(variant, density_class)])
  ba <- W[variable == "atba"]
  cat("populated cells:", ba[, .N, by = variant][, paste(variant, N)], "\n")

  cat("\n-- default vs calibrated MAP divergence at year 100 (cell medians; diff_pct_y100 = (cal - def) / def x 100)\n")
  dv <- ba[, .(med_abs_diff = median(abs(diff_pct_y100), na.rm = TRUE), max_abs_diff = max(abs(diff_pct_y100), na.rm = TRUE),
               n_within1 = sum(abs(diff_pct_y100) <= 1, na.rm = TRUE), n = .N, identical = sum(def_median_y100 == cal_median_y100, na.rm = TRUE)), by = variant]
  print(dv); for (i in seq_len(nrow(dv))) { put(tag, "div_med_abs_pct", dv$variant[i], dv$med_abs_diff[i]); put(tag, "div_max_abs_pct", dv$variant[i], dv$max_abs_diff[i]); put(tag, "div_n_within1", dv$variant[i], dv$n_within1[i]) }
  print(ba[order(diff_pct_y100)][1:3, .(variant, species_group, site_class, density_class, def = def_median_y100, cal = cal_median_y100, diff_pct_y100)])

  ## ------------------------------------------------------------------ bands from the raw per-draw runs
  fs <- list.files(tag, pattern = "^bakuzis_uncertainty_(ne|acd)_n100_batch0[0-5]\\.csv$", full.names = TRUE)
  d <- rbindlist(lapply(fs, fread, select = c("year", "atba", "tpa", "variant", "config", "draw_id", "scenario", "replicate", "species_group", "site_class", "density_class")), fill = TRUE)
  d[, h := year - min(year), by = .(variant, scenario, replicate, config, draw_id)]
  d[, h := 5L * as.integer(round(h / 5))]
  for (hz in c(25L, 50L, 100L)) {
    x <- d[h == hz]
    p <- x[config == "posterior", .(q025 = quantile(atba, .025), q500 = median(atba), q975 = quantile(atba, .975), nuniq = uniqueN(atba), n = .N),
           by = .(variant, scenario, replicate)]
    p[, band_pct := 100 * (q975 - q025) / q500]
    m <- x[config == "calibrated_map", .(map = median(atba)), by = .(variant, scenario, replicate)]
    p <- merge(p, m, by = c("variant", "scenario", "replicate"), all.x = TRUE)
    cellw <- p[, .(within = median(band_pct, na.rm = TRUE)), by = .(variant, scenario)]
    btw <- m[, .(range_pct = 100 * (max(map) - min(map)) / median(map), cv_pct = 100 * sd(map) / mean(map), n_st = .N), by = .(variant, scenario)]
    r <- merge(cellw[, .(within_cellmed_med = median(within), within_q25 = quantile(within, .25), within_q75 = quantile(within, .75)), by = variant],
               btw[, .(between_range_med = median(range_pct), between_cv_med = median(cv_pct)), by = variant], by = "variant")
    r <- merge(r, p[, .(stands = .N, within_stand_med = median(band_pct, na.rm = TRUE), stands_identical_draws = sum(nuniq == 1),
                        post_med_over_map_med = median(q500 / map, na.rm = TRUE), post_med_over_map_q05 = quantile(q500 / map, .05, na.rm = TRUE),
                        post_med_over_map_q95 = quantile(q500 / map, .95, na.rm = TRUE)), by = variant], by = "variant")
    pooled <- x[config == "posterior", .(pb = 100 * (quantile(atba, .975) - quantile(atba, .025)) / median(atba)), by = .(variant, scenario)][, .(pooled_legacy_med = median(pb)), by = variant]
    r <- merge(r, pooled, by = "variant")
    cat(sprintf("\n-- BA bands at year %d (percent of median). within = parameter only; between = 5 stands' MAP values; pooled_legacy mixes both (not used)\n", hz))
    print(r, digits = 4)
    for (i in seq_len(nrow(r))) for (k in setdiff(names(r), "variant")) put(tag, sprintf("y%03d_%s", hz, k), r$variant[i], r[[k]][i])
  }

  cat("\n-- stems retained at year 100 under default (cell medians)\n")
  tp <- L[variable == "tpa" & config == "default", .(y0 = median[year == min(year)], y100 = median[year == max(year)]), by = .(variant, scenario)]
  tr <- tp[, .(med_ret100 = median(y100 / y0)), by = variant]; print(tr); for (i in seq_len(nrow(tr))) put(tag, "stems_retained_y100", tr$variant[i], tr$med_ret100[i])
  cat("\n-- cell-median BA maximum over 100 years\n")
  bm <- L[variable == "atba" & config %in% c("default", "calibrated_map")][, .(max_ft2ac = max(median), max_m2ha = round(max(median) * 0.229568, 1)), by = .(variant, config)]
  print(bm); for (i in seq_len(nrow(bm))) put(tag, paste0("max_cell_med_ba_ft2ac_", bm$config[i]), bm$variant[i], bm$max_ft2ac[i])

  sd_ <- score_dir[[tag]]; f <- file.path(sd_, sprintf("%s_scorecard_summary.csv", score_lab[[tag]]))
  if (file.exists(f)) {
    sc <- fread(f)
    cat("\n-- magnitude (T4) and Reineke (T3, T3b) counts from", f, "\n")
    s4 <- sc[grepl("^T4|^T3", relation), .(arm, relation, evaluated, not_evaluated, not_applicable, flagged, verdict)]
    print(s4)
    for (i in seq_len(nrow(s4))) { put(tag, paste0("score_flagged|", s4$relation[i]), s4$arm[i], s4$flagged[i]); put(tag, paste0("score_evaluated|", s4$relation[i]), s4$arm[i], s4$evaluated[i]) }
  } else cat("\n(no scorecard at", f, ")\n")
  lf <- sprintf("%s_agg/bakuzis_laws_compliance.csv", tag)
  if (file.exists(lf)) {
    la <- fread(lf); cat("\n-- Bakuzis law fractions (posterior_band row is computed on the pooled band medians)\n"); print(la)
    for (i in seq_len(nrow(la))) for (k in grep("fraction$", names(la), value = TRUE)) put(tag, paste0(k, "|", la$config[i]), la$variant[i], la[[k]][i])
  }
}
o <- rbindlist(out_rows)
if (length(tags) == 2) {
  w <- dcast(o, metric + variant ~ tag, value.var = "value", fun.aggregate = function(z) z[1])
  fwrite(w, "ens4_numbers_old_vs_new.csv"); cat("\nWrote ens4_numbers_old_vs_new.csv (", nrow(w), "rows )\n")
} else fwrite(o, sprintf("ens4_numbers_%s.csv", tags[1]))
