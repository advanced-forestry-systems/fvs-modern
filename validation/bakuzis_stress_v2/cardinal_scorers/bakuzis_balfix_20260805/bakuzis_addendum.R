#!/usr/bin/env Rscript
# bakuzis_addendum.R
#  (a) Eichhorn on the RAW plot-level trajectories (the age-gridded aggregate
#      has no common height window between site class 1 and site class 4).
#  (b) 55 %-of-SDImax gate diagnostic -> vls/morts.f90 RIP=RI branch.
#  (c) BAMAX ceiling diagnostic       -> vls/morts.f90 lines 750-786.

suppressPackageStartupMessages(library(ggplot2))
has_ggsci <- requireNamespace("ggsci", quietly = TRUE)
setwd(Sys.getenv("BAKDIR", "/fs/scratch/PUOM0008/crsfaaron/bakuzis_20260804"))

K_TPH <- 2.4710538; K_BAPH <- 0.2295684; K_QMD <- 2.54
K_HT  <- 0.3048;    K_VOL  <- 0.06997
PMSDIU <- 0.85; KSDI <- 0.5454154

raw <- read.csv("out/bakuzis_traj_raw.csv", stringsAsFactors = FALSE)
raw <- raw[is.finite(raw$Tpa) & raw$Tpa > 0 & is.finite(raw$QMD) & raw$QMD > 0, ]
raw$HTm  <- raw$TopHt * K_HT
raw$VOLm <- raw$TCuFt * K_VOL
raw$BAPH <- raw$BA * K_BAPH
raw$RD   <- ifelse(raw$SDIMax > 0, raw$SDI / raw$SDIMax, NA)
raw$BAMAX_ft2 <- raw$SDIMax * KSDI * PMSDIU     # sdical.f90 line 202
raw$BA_frac   <- raw$BA / raw$BAMAX_ft2

ARMS <- c("fvs_base", "fvs_regional", "organon")
ARMS <- ARMS[ARMS %in% unique(raw$arm)]
out <- list(); add <- function(...) out[[length(out) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

sink("bakuzis_addendum.txt")
cat(strrep("=", 74), "\nADDENDUM: Eichhorn (raw), SDI gate, BAMAX ceiling\n",
    strrep("=", 74), "\n", sep = "")

## ---- (a) Eichhorn on raw data, height-binned ----------------------------
cat("\n\n[8 EICHHORN, raw plot-level, 2 m height bins, min 20 plot-cycles per",
    "site class per bin]\n")
eich_rows <- list()
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$HTm) & is.finite(raw$VOLm) & raw$VOLm > 0, ]
  d$hbin <- floor(d$HTm / 2) * 2 + 1
  tb <- aggregate(list(VOL = d$VOLm), by = list(hbin = d$hbin, site = d$site_cls),
                  FUN = mean)
  nb <- aggregate(list(n = d$VOLm), by = list(hbin = d$hbin, site = d$site_cls),
                  FUN = length)
  tb <- merge(tb, nb); tb <- tb[tb$n >= 20, ]
  ok <- names(which(table(tb$hbin) == length(unique(d$site_cls))))
  tb <- tb[as.character(tb$hbin) %in% ok, ]
  cv <- tapply(tb$VOL, tb$hbin, function(x) sd(x) / mean(x))
  mcv <- mean(cv, na.rm = TRUE)
  cat(sprintf("\n  %-13s bins with all 4 site classes: %d (HT %.0f-%.0f m)\n",
              a, length(cv), min(tb$hbin), max(tb$hbin)))
  cat(sprintf("  %-13s per-bin CV: %s\n", "",
              paste(sprintf("%.0fm=%.1f%%", as.numeric(names(cv)), 100 * cv),
                    collapse = " ")))
  cat(sprintf("  %-13s MEAN Eichhorn CV = %.1f%%  -> %s\n", "", 100 * mcv,
              ifelse(mcv < 0.15, "PASS (site-independent)",
                     "FLAG (site-dependent V-H)")))
  add(arm = a, relation = "8 Eichhorn (volume-height)",
      statistic = sprintf("mean CV of VOL at matched HT = %.1f%% over %d common 2 m bins (threshold 15%%)",
                          100 * mcv, length(cv)),
      decision = ifelse(mcv < 0.15, "PASS", "FLAG"))
  tb$arm <- a; eich_rows[[a]] <- tb
}

## ---- (b) the 55 % SDImax gate -------------------------------------------
cat("\n\n[SDI GATE] vls/morts.f90: density mortality engages only above",
    "PMSDIL = 55% of SDImax\n")
cat("           (line ~530: IF(T .LE. TEM .OR. RN .LE. 0.0) RIP=RI)\n")
cat("           Below the gate ONLY the halved background rate applies",
    "(line 512: RI = 0.5 * RI)\n\n")
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$RD), ]
  f <- mean(d$RD < 0.55, na.rm = TRUE)
  fl <- mean(d$RD[d$dens_cls == "low"] < 0.55, na.rm = TRUE)
  fh <- mean(d$RD[d$dens_cls == "high"] < 0.55, na.rm = TRUE)
  bysite <- tapply(d$RD, d$site_cls, function(x) mean(x < 0.55, na.rm = TRUE))
  cat(sprintf("  %-13s %% of plot-cycles BELOW the 55%% gate: all=%.1f%%  low-density=%.1f%%  high-density=%.1f%%\n",
              a, 100 * f, 100 * fl, 100 * fh))
  cat(sprintf("  %-13s by site class: %s\n", "",
              paste(sprintf("S%s=%.0f%%", names(bysite), 100 * bysite), collapse = "  ")))
  add(arm = a, relation = "DIAG: SDI gate (55% SDImax)",
      statistic = sprintf("%.1f%% of plot-cycles below the gate (low-density %.1f%%, high-density %.1f%%): density mortality inactive",
                          100 * f, 100 * fl, 100 * fh),
      decision = "DIAGNOSTIC")
}

## ---- (c) the BAMAX ceiling ----------------------------------------------
cat("\n\n[BAMAX CEILING] vls/morts.f90 lines 750-786 inflate mortality until\n")
cat("                residual BA <= BAMAX + 1 ft2/ac, where BAMAX =\n")
cat("                SDImax * 0.5454154 * 0.85 (vbase/sdical.f90 line 202)\n\n")
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$BA_frac), ]
  pinned <- mean(abs(d$BA - d$BAMAX_ft2) <= 1.0, na.rm = TRUE)
  over   <- mean(d$BA_frac > 0.99, na.rm = TRUE)
  cat(sprintf("  %-13s max BA/BAMAX = %.3f ; %% of plot-cycles pinned within 1 ft2/ac of BAMAX = %.1f%% ; %% above 0.99*BAMAX = %.1f%%\n",
              a, max(d$BA_frac, na.rm = TRUE), 100 * pinned, 100 * over))
  add(arm = a, relation = "DIAG: BAMAX ceiling",
      statistic = sprintf("%.1f%% of plot-cycles pinned within 1 ft2 ac-1 of BAMAX (max BA/BAMAX %.3f)",
                          100 * pinned, max(d$BA_frac, na.rm = TRUE)),
      decision = "DIAGNOSTIC")
}
RES <- do.call(rbind, out)
cat("\n\n", strrep("=", 74), "\n", sep = "")
print(RES, right = FALSE)
sink()
write.csv(RES, "bakuzis_addendum_table.csv", row.names = FALSE)

## ---- figure: Eichhorn on raw data ---------------------------------------
E <- do.call(rbind, eich_rows)
E$site_f <- factor(E$site, levels = 1:4,
                   labels = c("S1 (SI 47)", "S2 (SI 61)", "S3 (SI 70)", "S4 (SI 83)"))
E$arm_f <- factor(E$arm, levels = ARMS,
                  labels = c(fvs_base = "fvs_base", fvs_regional = "fvs_regional",
                             organon = "organon (Greg)")[ARMS])
pal_site <- if (has_ggsci) ggsci::pal_npg()(4) else
  c("#3B4992", "#008B45", "#EE7733", "#BB0021")
names(pal_site) <- levels(E$site_f)
p <- ggplot(E, aes(hbin, VOL, colour = site_f)) +
  geom_line(linewidth = 0.9) + geom_point(size = 1.1) +
  facet_wrap(~arm_f) +
  scale_colour_manual(values = pal_site, name = "Site class") +
  labs(x = expression(Dominant~height~(m)),
       y = expression(Gross~cubic~volume~(m^3~ha^-1))) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))
png("fig/fig5_eichhorn.png", width = 9, height = 4, units = "in", res = 300)
print(p); dev.off()

## ---- figure: BA relative to the BAMAX ceiling ---------------------------
B <- raw[is.finite(raw$BA_frac), ]
B$arm_f <- factor(B$arm, levels = ARMS,
                  labels = c(fvs_base = "fvs_base", fvs_regional = "fvs_regional",
                             organon = "organon (Greg)")[ARMS])
B$site_f <- factor(B$site_cls, levels = 1:4,
                   labels = c("S1 (SI 47)", "S2 (SI 61)", "S3 (SI 70)", "S4 (SI 83)"))
agg <- aggregate(list(BA_frac = B$BA_frac),
                 by = list(Age = B$Age, arm_f = B$arm_f, site_f = B$site_f), FUN = mean)
agg <- agg[agg$Age <= 200, ]
p7 <- ggplot(agg, aes(Age, BA_frac, colour = site_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  geom_hline(yintercept = 1, colour = "grey25", linetype = "dashed") +
  scale_colour_manual(values = pal_site, name = "Site class") +
  labs(x = "Stand age (yr)",
       y = expression(Basal~area~relative~to~BAMAX~ceiling)) +
  theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"), legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))
png("fig/fig7_bamax_ceiling.png", width = 9, height = 4, units = "in", res = 300)
print(p7); dev.off()
cat("DONE\n")
