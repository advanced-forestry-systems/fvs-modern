## runoff_figures.R (2026-10-30): figures for the CONUS model run-off technical report.
## Inputs: full/score_fix/near/near_scores.csv (common set), full/score_fix/gates/catalog_gates_by_arm.csv,
## full/score_fix/decision/runoff_decision_by_variant.csv. Headless (ragg), 300 dpi.
suppressMessages({ library(data.table); library(ggplot2); library(ragg) })
F <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/full/score_fix"; O <- "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/report/figures"
lab <- c("usfs_native" = "USDA FVS", "default_native" = "FVS-modern default, native", "default_gompit" = "FVS-modern default, GOMPIT",
  "calibrated_native" = "FVS-modern calibrated, native", "calibrated_gompit" = "FVS-modern calibrated, GOMPIT",
  "conus_v2_b1" = "fvs-conus v2 b1, logit", "conus_v2_b1_gompit" = "fvs-conus v2 b1, GOMPIT", "conus_v2_b2" = "fvs-conus v2 b2, logit", "conus_v2_b2_gompit" = "fvs-conus v2 b2, GOMPIT",
  "conus_v6_b1" = "fvs-conus v6 b1, logit (diagnostic)", "conus_v6_b1_gompit" = "fvs-conus v6 b1, GOMPIT (diagnostic)",
  "conus_v6_b2" = "fvs-conus v6 b2, logit (diagnostic)", "conus_v6_b2_gompit" = "fvs-conus v6 b2, GOMPIT (diagnostic)")
fam <- function(a) fifelse(a == "usfs_native", "USDA FVS", fifelse(grepl("^default", a), "FVS-modern default", fifelse(grepl("^calibrated", a), "FVS-modern calibrated",
  fifelse(grepl("conus_v2", a), "fvs-conus v2", "fvs-conus v6 (diagnostic)"))))
famcol <- c("USDA FVS" = "#4D4D4D", "FVS-modern default" = "#2C7FB8", "FVS-modern calibrated" = "#08519C", "fvs-conus v2" = "#D95F0E", "fvs-conus v6 (diagnostic)" = "#BDBDBD")
norm <- function(d, col = "arm") { d[get(col) == "usfs_native(ne_lib)", (col) := "usfs_native"]; d[, arm_lab := factor(lab[get(col)], levels = rev(lab))]; d[, family := factor(fam(get(col)), levels = names(famcol))]; d }
th <- theme_classic(base_size = 9) + theme(panel.background = element_rect(fill = "white", colour = NA), strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 9), panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
  legend.position = "bottom", legend.title = element_blank(), axis.title.y = element_blank())
## (a) near-term BA percent bias
N <- norm(fread(file.path(F, "near/near_scores.csv"))[set == "common" & metric == "BA"])
N[, vlab := paste0(variant, " (n = ", n, ")")]
pa <- ggplot(N, aes(pbias, arm_lab, colour = family)) + geom_vline(xintercept = 0, linewidth = 0.4) +
  geom_errorbarh(aes(xmin = pbias_lo, xmax = pbias_hi), height = 0, linewidth = 0.5) + geom_point(size = 1.4) +
  facet_wrap(~vlab, ncol = 6) + scale_colour_manual(values = famcol) +
  labs(x = "Near-term basal area bias (% of observed survivor basal area, m² ha⁻¹ basis), 95% stand-bootstrap interval") + th
ggsave(file.path(O, "fig_a_nearterm_ba_bias.png"), pa, width = 13, height = 8.5, dpi = 300, device = ragg::agg_png)
## (b) share above FIA BA p99.9 by year 50, Wilson 95% interval
G <- norm(fread(file.path(F, "gates/catalog_gates_by_arm.csv")))
wil <- function(x, n, z = 1.959964) { p <- x / n; c0 <- (p + z^2 / (2 * n)) / (1 + z^2 / n); h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n); list(pmax(c0 - h, 0), pmin(c0 + h, 1)) }
G[, x := round(share_ba999_50 * n_traj)][, c("lo", "hi") := wil(x, n_traj)][, vlab := paste0(group, " (n = ", n_traj, ")")]
pb <- ggplot(G, aes(100 * share_ba999_50, arm_lab, colour = family)) + geom_vline(xintercept = 5, linetype = "dashed", linewidth = 0.3) + geom_vline(xintercept = 10, linetype = "dotted", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = 100 * lo, xmax = 100 * hi), height = 0, linewidth = 0.5) + geom_point(size = 1.4) +
  facet_wrap(~vlab, ncol = 6) + scale_colour_manual(values = famcol) +
  labs(x = "Stands whose maximum basal area (m² ha⁻¹) exceeds the FIA plot p99.9 within 50 years (%), Wilson 95% interval; dashed 5%, dotted 10%") + th
ggsave(file.path(O, "fig_b_longterm_ba_exceedance.png"), pb, width = 13, height = 8.5, dpi = 300, device = ragg::agg_png)
## (c) decision map (v6 excluded)
D <- fread(file.path(F, "decision/runoff_decision_by_variant.csv"))
arms <- names(lab)[!grepl("v6", names(lab))]
M <- CJ(variant = D$variant, arm = arms)
M <- merge(M, D[, .(variant, admissible, chosen, noadm = admissible == "" | is.na(admissible))], by = "variant")
M[, adm := mapply(function(a, s) a %in% strsplit(sub("usfs_native\\(ne_lib\\)", "usfs_native", s), ";")[[1]], arm, admissible)]
M[, win := !is.na(chosen) & sub("usfs_native\\(ne_lib\\)", "usfs_native", chosen) == arm]
M[, status := fifelse(win, "Winner", fifelse(adm, "Admissible", fifelse(noadm, "No admissible arm (not recommended)", "Not admissible")))]
M[, status := factor(status, levels = c("Winner", "Admissible", "Not admissible", "No admissible arm (not recommended)"))]
M[, arm_lab := factor(lab[arm], levels = lab[arms])][, variant := factor(variant, levels = rev(sort(unique(variant))))]
pc <- ggplot(M, aes(arm_lab, variant, fill = status)) + geom_tile(colour = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("Winner" = "#08519C", "Admissible" = "#9ECAE1", "Not admissible" = "#F0F0F0", "No admissible arm (not recommended)" = "#FDD0A2")) +
  geom_text(data = M[win == TRUE], aes(label = "W"), colour = "white", size = 3, fontface = "bold") +
  labs(x = NULL, y = "FVS variant") + theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "bottom", legend.title = element_blank(), axis.line = element_blank(), axis.ticks = element_blank())
ggsave(file.path(O, "fig_c_decision_map.png"), pc, width = 8, height = 7, dpi = 300, device = ragg::agg_png)
cat("figures written\n")
