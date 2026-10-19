# Title:   fvs-modern manuscript Figures 2, 4, and 5 regenerated from the 19-variant ablation_full benchmark
# Input:   calibration/output/comparisons/manuscript_tables/fia_benchmark_results_ablation_full.csv
# Output:  manuscript/figures_20261019/fig2_benchmark_rmse.(png|pdf), fig4_composite_pctrmse_dumbbell.(png|pdf),
#          fig5_r2_heatmap.(png|pdf), and fig_data_20261019.csv (every plotted value)
# Run on:  OSC Cardinal, R 4.5.2 (Aaron named Cardinal, October 19, 2026)

suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- "/users/PUOM0008/crsfaaron/fvs-modern"
f_in <- file.path(root, "calibration/output/comparisons/manuscript_tables/fia_benchmark_results_ablation_full.csv")
out  <- file.path(root, "manuscript/figures_20261019")
dir.create(out, showWarnings = FALSE)
d <- fread(f_in)
v <- d[VARIANT != "OVERALL"]
stopifnot(nrow(v) == 19)

names_v <- c(AK = "Alaska", BM = "Blue Mountains", CA = "Inland California and Southern Cascades",
  CI = "Central Idaho", CR = "Central Rockies", CS = "Central States", EC = "East Cascades",
  EM = "Eastern Montana", IE = "Inland Empire", LS = "Lake States", NC = "Klamath Mountains",
  NE = "Northeast", PN = "Pacific Northwest Coast", SN = "Southern", SO = "South Central Oregon and Northeast California",
  TT = "Tetons", UT = "Utah", WC = "Westside Cascades", WS = "Western Sierra Nevada")
stopifnot(all(v$VARIANT %in% names(names_v)))

theme_j <- theme_classic(base_size = 9, base_family = "sans") +
  theme(panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_text(face = "bold", size = 9, hjust = 0),
        legend.position = "bottom", legend.title = element_blank(), axis.text = element_text(colour = "black"),
        plot.margin = margin(4, 8, 4, 4))
col_def <- "grey55"; col_cal <- "#134a31"

# ---------- Figure 2: BA and gross volume RMSE, default vs calibrated, 95% bootstrap interval on calibrated ----------
mk <- function(p, lab) v[, .(VARIANT, n = n_conditions, default = get(paste0(p, "_RMSE_default")),
                              calib = get(paste0(p, "_RMSE_calib")), lo = get(paste0(p, "_RMSE_lo_calib")),
                              hi = get(paste0(p, "_RMSE_hi_calib")), panel = lab)]
f2 <- rbind(mk("BA", "a"), mk("VOL_CFGRS", "b"))
f2[, red := 100 * (default - calib) / default]
ord <- f2[panel == "a"][order(red), VARIANT]
f2[, VARIANT := factor(VARIANT, levels = ord)]
f2[, panel := factor(panel, levels = c("a", "b"),
     labels = c("(a)~Basal~area~RMSE~(ft^2~ac^-1)", "(b)~Gross~cubic~volume~RMSE~(ft^3~ac^-1)"))]
p2 <- ggplot(f2, aes(y = VARIANT)) +
  geom_segment(aes(x = default, xend = calib, yend = VARIANT), colour = "grey80", linewidth = 0.5) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, colour = col_cal, linewidth = 0.5) +
  geom_point(aes(x = default, shape = "Default"), colour = col_def, size = 1.8, stroke = 0.6) +
  geom_point(aes(x = calib, shape = "Calibrated"), colour = col_cal, size = 1.8) +
  scale_shape_manual(values = c(Default = 1, Calibrated = 16),
                     labels = c(Calibrated = "Calibrated (95% bootstrap interval)", Default = "Default")) +
  facet_wrap(~panel, scales = "free_x", labeller = label_parsed) +
  labs(x = "RMSE", y = "Variant (ordered by basal area RMSE reduction)") + theme_j
ggsave(file.path(out, "fig2_benchmark_rmse.png"), p2, width = 6.5, height = 4.6, dpi = 600)
ggsave(file.path(out, "fig2_benchmark_rmse.pdf"), p2, width = 6.5, height = 4.6)

# ---------- Figure 4: composite percent RMSE dumbbell ----------
attrs <- c("BA", "SDI", "QMD", "VOL_CFGRS")
v[, comp_default := rowMeans(.SD), .SDcols = paste0(attrs, "_RMSE_pct_default")]
v[, comp_calib   := rowMeans(.SD), .SDcols = paste0(attrs, "_RMSE_pct_calib")]
f4 <- v[, .(VARIANT, comp_default, comp_calib, worse = comp_calib > comp_default)]
f4[, VARIANT := factor(VARIANT, levels = f4[order(comp_calib - comp_default, decreasing = TRUE), VARIANT])]
p4 <- ggplot(f4, aes(y = VARIANT)) +
  geom_segment(aes(x = comp_default, xend = comp_calib, yend = VARIANT, colour = worse), linewidth = 0.6) +
  geom_point(aes(x = comp_default, shape = "Default"), colour = col_def, size = 1.9, stroke = 0.6) +
  geom_point(aes(x = comp_calib, shape = "Calibrated"), colour = col_cal, size = 1.9) +
  scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#b2182b"), guide = "none") +
  scale_shape_manual(values = c(Default = 1, Calibrated = 16)) +
  labs(x = "Composite RMSE (% of observed mean; equal weights on BA, SDI, QMD, gross volume)", y = "Variant") + theme_j
ggsave(file.path(out, "fig4_composite_pctrmse_dumbbell.png"), p4, width = 6.5, height = 4.4, dpi = 600)
ggsave(file.path(out, "fig4_composite_pctrmse_dumbbell.pdf"), p4, width = 6.5, height = 4.4)

# ---------- Figure 5: R2 by variant and metric, default and calibrated ----------
mets <- c(BA = "BA", TPA = "TPA", QMD = "QMD", SDI = "SDI", VOL_CFGRS = "Gross cubic", VOL_CFNET = "Net cubic",
          VOL_BFNET = "Net board", HT_top = "Top height")
f5 <- rbindlist(lapply(names(mets), function(m) rbind(
  v[, .(VARIANT, metric = mets[[m]], approach = "Default", R2 = get(paste0(m, "_r2_default")))],
  v[, .(VARIANT, metric = mets[[m]], approach = "Calibrated", R2 = get(paste0(m, "_r2_calib")))])))
f5[, metric := factor(metric, levels = unname(mets))]
f5[, approach := factor(approach, levels = c("Default", "Calibrated"))]
f5[, VARIANT := factor(VARIANT, levels = v[order(BA_r2_calib), VARIANT])]
p5 <- ggplot(f5, aes(metric, VARIANT, fill = R2)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = ifelse(is.na(R2), "", sprintf("%.2f", R2)),
                colour = !is.na(R2) & R2 > 0.6), size = 2.1) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
  scale_fill_viridis_c(option = "D", limits = c(0, 1), na.value = "grey90", direction = -1,
                       name = expression(italic(R)^2)) +
  facet_wrap(~approach) + labs(x = NULL, y = "Variant") + theme_j +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "right", legend.title = element_text(size = 9),
        axis.line = element_blank(), panel.grid.major.x = element_blank())
ggsave(file.path(out, "fig5_r2_heatmap.png"), p5, width = 6.5, height = 5.0, dpi = 600)
ggsave(file.path(out, "fig5_r2_heatmap.pdf"), p5, width = 6.5, height = 5.0)


# ---------- Figure 3: percent BA RMSE reduction with 95% interval, grouped by region ----------
reg <- list(Pacific = c("AK", "PN", "WC", "NC", "CA", "SO", "WS"),
            `Interior West` = c("BM", "EC", "IE", "CI", "EM", "UT", "CR", "TT"), East = c("SN", "CS", "LS", "NE"))
stopifnot(setequal(unlist(reg), v$VARIANT))
f3 <- v[, .(VARIANT, red = 100 * (BA_RMSE_default - BA_RMSE_calib) / BA_RMSE_default,
             lo = 100 * (BA_RMSE_default - BA_RMSE_hi_calib) / BA_RMSE_default,
             hi = 100 * (BA_RMSE_default - BA_RMSE_lo_calib) / BA_RMSE_default)]
f3[, region := factor(names(reg)[sapply(VARIANT, function(x) which(sapply(reg, function(r) x %in% r)))], levels = names(reg))]
f3[, inc0 := lo <= 0]
f3[, VARIANT := factor(VARIANT, levels = f3[order(red), VARIANT])]
p3 <- ggplot(f3, aes(y = VARIANT, colour = inc0)) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.5) +
  geom_point(aes(x = red), size = 1.8) +
  scale_colour_manual(values = c(`FALSE` = col_cal, `TRUE` = col_def), guide = "none") +
  facet_grid(region ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Reduction in basal area RMSE (%), 95% bootstrap interval", y = NULL) + theme_j +
  theme(strip.text.y = element_text(angle = 0, face = "bold", hjust = 0))
ggsave(file.path(out, "fig3_ba_rmse_reduction_by_region.png"), p3, width = 6.5, height = 5.6, dpi = 600)
ggsave(file.path(out, "fig3_ba_rmse_reduction_by_region.pdf"), p3, width = 6.5, height = 5.6)

# ---------- Plotted values ----------
fwrite(rbind(
  f2[, .(figure = "2", VARIANT, panel = as.character(panel), default, calib, lo, hi, n)],
  f4[, .(figure = "4", VARIANT, panel = "composite_pct", default = comp_default, calib = comp_calib, lo = NA, hi = NA, n = NA)],
  f5[, .(figure = "5", VARIANT, panel = paste(metric, approach), default = NA, calib = R2, lo = NA, hi = NA, n = NA)],
  f3[, .(figure = "3", VARIANT, panel = as.character(region), default = NA, calib = red, lo, hi, n = NA)],
  fill = TRUE), file.path(out, "fig_data_20261019.csv"))
cat("worse composite:", paste(f4[worse == TRUE, as.character(VARIANT)], collapse = ","), "\n")
print(f4[order(comp_calib - comp_default)])
cat("done\n")
