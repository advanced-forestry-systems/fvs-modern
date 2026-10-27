#!/usr/bin/env Rscript
# bakuzis_matrix.R -- draw the modified Bakuzis matrix (Leary triangular form)
# Usage: Rscript bakuzis_matrix.R trajectories.csv out.png [form_factor=0.40]
# Input columns: site, age, HT, TPH, BAPH, QMD (origin, VOL optional). Base R graphics.

args <- commandArgs(trailingOnly = TRUE)
csv <- args[1]; out <- if (length(args) >= 2) args[2] else "bakuzis_matrix.png"
ff  <- if (length(args) >= 3) as.numeric(args[3]) else 0.40
d <- read.csv(csv, stringsAsFactors = FALSE)
if (!"origin" %in% names(d)) d$origin <- "all"
if (!"VOL" %in% names(d)) d$VOL <- d$BAPH * d$HT * ff
d <- d[order(d$origin, d$site, d$age), ]
sites <- sort(unique(d$site)); origins <- unique(d$origin)
pal <- setNames(grDevices::hcl.colors(length(sites), "Zissou 1"), sites)
lty <- setNames(seq_along(origins), origins)

png(out, width = 2400, height = 2200, res = 220)
par(mfrow = c(3, 3), mar = c(4, 4, 2.5, 1), cex.axis = 0.8, cex.lab = 0.9)
panel <- function(xk, yk, title, log = "") {
  plot(d[[xk]], d[[yk]], type = "n", xlab = xk, ylab = yk, main = title, log = log)
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]
    lines(g[[xk]], g[[yk]], col = pal[as.character(s)], lty = lty[[o]], lwd = 2)
  }
  grid()
}
panel("age","HT","1. Site curves (HT-age)")
panel("QMD","HT","2. Height-diameter")
plot.new(); legend("center", title = "Bakuzis matrix",
  legend = c(paste("site", sites), origins),
  col = c(pal, rep(1, length(origins))), lty = c(rep(1, length(sites)), lty), lwd = 2, bty = "n", cex = 0.9)
panel("age","TPH","3. Sukachev (stems-age)")
panel("QMD","TPH","4. Reineke (log-log)", log = "xy")
panel("HT","TPH","5. Spacing (stems-height)")
panel("TPH","BAPH","6. Stocking guide (BA-stems)")
panel("HT","VOL","8. Eichhorn (volume-height)")
panel("TPH","VOL","9. Yield-density (volume-stems)")
dev.off()
cat("wrote", out, "\n")
