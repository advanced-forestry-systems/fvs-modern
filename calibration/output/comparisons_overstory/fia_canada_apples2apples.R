#!/usr/bin/env Rscript
# fia_canada_apples2apples.R   (base R)
# =====================================================================
# Apples-to-apples cross-region scorecard for basal area, using the
# overstory FIA recompute (SLURM 10591333; DIA>=5", FVS_ACD_RELABEL=TRUE)
# plus the Canadian MAGPlot NB benchmark from earlier sessions.
#
# Models compared on BA bias (negative = under-projects):
#                       FIA (overstory)    Canada NB (262)
#   FVS-ACD calibrated      -0.06%             -0.04%
#   OSM-ACD                 -2.18%             -3.00%
#   FVS-NE calibrated       -3.95%             (no ACD overlay on NB)
#   FVS-NE default          -8.20%             +9.10%
#   FVS-ACD default         -6.07%              n/a
#
# Output: fia_canada_apples2apples.png + fia_canada_apples2apples.csv
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

# Pull live FIA biases from the validation file so the figure stays in sync
d <- read.csv(file.path(od, "validation_data_overstory.csv"))
m <- function(x) mean(x, na.rm = TRUE)
bias <- function(p, o) 100*(m(p)/m(o) - 1)
fmt <- function(b) sprintf("%+.2f%%", b)

acd <- d[d$VARIANT == "ACD" & is.finite(d$BA_t2) & d$BA_t2 > 0, ]
ne  <- d[d$VARIANT == "NE"  & is.finite(d$BA_t2) & d$BA_t2 > 0, ]

fia_acd_cal <- bias(acd$BA_pred_calib,   acd$BA_t2)
fia_acd_def <- bias(acd$BA_pred_default, acd$BA_t2)
fia_ne_cal  <- bias(ne$BA_pred_calib,    ne$BA_t2)
fia_ne_def  <- bias(ne$BA_pred_default,  ne$BA_t2)
fia_osm     <- -2.18  # OSM-ACD prior benchmark (12,029 plots, R2 0.96)

# Canada NB (from prior MAGPlot work)
ca_acd_cal  <- -0.04   # AcadianGY R on 262 NB pairs
ca_osm      <- -3.00   # OSM-ACD on 262 NB pairs
ca_ne_def   <-  9.10   # FVS-NE Fortran default on 262 NB pairs (this session's earlier work)

models <- c("FVS-ACD\ncalibrated", "OSM-ACD", "FVS-NE\ncalibrated", "FVS-NE\ndefault", "FVS-ACD\ndefault")
fia    <- c(fia_acd_cal, fia_osm, fia_ne_cal, fia_ne_def, fia_acd_def)
canada <- c(ca_acd_cal,  ca_osm,  NA,         ca_ne_def,  NA)
n_fia  <- c(nrow(acd), 12029, nrow(ne), nrow(ne), nrow(acd))
n_ca   <- c(262, 262, NA, 262, NA)
cols   <- c("#2ca02c", "#7DB5D5", "#C5A55A", "#5A5A5A", "#a07050")

# CSV
summ <- data.frame(model = c("FVS-ACD calibrated","OSM-ACD","FVS-NE calibrated","FVS-NE default","FVS-ACD default"),
                   FIA_bias_pct = round(fia, 2), n_FIA = n_fia,
                   Canada_bias_pct = round(canada, 2), n_Canada = n_ca,
                   row.names = NULL)
write.csv(summ, file.path(od, "fia_canada_apples2apples.csv"), row.names = FALSE)
cat("=== Apples-to-apples cross-region scorecard (BA bias) ===\n"); print(summ)

png(file.path(od, "fia_canada_apples2apples.png"), width = 2300, height = 980, res = 170)
par(mfrow = c(1, 2), mar = c(5.5, 4.8, 4.0, 1.0), mgp = c(2.8, 0.7, 0))

## Panel A: grouped bars FIA vs Canada
M <- rbind(FIA = fia, Canada = canada)
colnames(M) <- models
bp <- barplot(M, beside = TRUE, col = c("#5A5A5A", "#2ca02c"), border = NA,
              ylim = c(-11, 12), ylab = "Basal-area bias vs observed (%)",
              cex.names = 0.78, las = 1)
abline(h = 0, lwd = 1.2)
v <- as.vector(M)
for (i in seq_along(v)) if (!is.na(v[i]))
  text(as.vector(bp)[i], v[i] + ifelse(v[i] >= 0, 0.7, -0.9), sprintf("%+.1f", v[i]), cex = 0.78, font = 2)
nax <- which(is.na(v))
if (length(nax)) text(as.vector(bp)[nax], 0.6, "n/a", cex = 0.7, col = "#888")
legend("topright", legend = c("FIA (DIA >= 5\")", "Canada NB (262)"), fill = c("#5A5A5A", "#2ca02c"),
       border = NA, bty = "n", cex = 0.9)
title(main = "Apples-to-apples FIA + Canada (overstory)\nFVS-ACD calibrated essentially unbiased BOTH regions",
      cex.main = 1.1, font.main = 2)

## Panel B: scatter FIA vs Canada for models with both regions
ok <- !is.na(canada)
fx <- fia[ok]; cy <- canada[ok]; nm <- models[ok]; co <- cols[ok]
lim <- c(-10, 11)
plot(NA, xlim = lim, ylim = lim, xlab = "FIA Maine/NH/VT BA bias (%)",
     ylab = "Canada NB BA bias (%)", las = 1)
abline(h = 0, v = 0, col = "#bbb"); abline(0, 1, lty = 3, col = "#888")
rect(-5, -5, 5, 5, border = "#2ca02c", lty = 2)
points(fx, cy, pch = 19, col = co, cex = 1.8)
text(fx, cy, nm, pos = c(4, 4, 1), cex = 0.85, font = 2, col = co, xpd = NA)
text(0, -8, "dashed box = accurate both regions (|bias| <= 5%)", cex = 0.78, col = "#2ca02c")
text(9, 9.5, "1:1 = consistent", cex = 0.75, col = "#888", srt = 25)
title(main = "Cross-region consistency:\nFVS-ACD calibrated and OSM-ACD both inside the accuracy box",
      cex.main = 1.05, font.main = 2)
dev.off()
cat("wrote", file.path(od, "fia_canada_apples2apples.png"), "\n")
