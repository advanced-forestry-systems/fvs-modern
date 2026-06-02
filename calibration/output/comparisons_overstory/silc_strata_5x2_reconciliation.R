#!/usr/bin/env Rscript
# silc_strata_5x2_reconciliation.R
# =====================================================================
# 11-byStrata-stand to 10-cell reconciliation. Shows how the prior
# SILC byStrata 11-stand grouping (used in v25..v28 decks) maps onto
# the new 5-type x 2-density break. One-page reference for SILC.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory"

# The 11 byStrata stands present in GrownDB and their year-0 metrics
# (from the AGM cross-check earlier this session). Hard-coded since
# GrownDB stand set is fixed.
stands <- data.frame(
  byStrata_id = c("C4C-N","H3B-N","H3C-N","H3D-N","H4C-N","HS3C-N",
                  "S2A-N","S3B-N","S3C-N","S3D-N","S4C-N"),
  BA_0  = c(161.6, 97.2, 75.9, 68.7, 79.2, 97.7, 88.4, 164.3, 100.6, 89.7, 138.1),
  TPA_0 = c(359.8, 234.3, 225.7, 229.1, 167.7, 267.0, 347.0, 394.3, 263.4, 328.5, 283.3),
  QMD_0 = c(9.07, 8.72, 7.85, 7.42, 9.30, 8.19, 6.83, 8.74, 8.37, 7.07, 9.45),
  stringsAsFactors = FALSE
)
# parse the same way the regen script does
parse_id <- function(s) {
  pre <- regmatches(s, regexpr("^[A-Z]+", s))
  rest <- substring(s, nchar(pre)+1)
  dens_letter <- regmatches(rest, regexpr("[A-D]", rest))
  ft <- switch(pre,
    "C"  = "Cedar", "H" = "Hardwood", "HS"="Mixedwood",
    "SH" = "Mixedwood", "S"="Commercial Softwood",
    "OS" = "Other Softwood", NA)
  dc <- if (dens_letter %in% c("A","B")) "A+B (high)" else "C+D (low)"
  c(forest_type = ft, density_class = dc,
    prefix = pre, density_letter = dens_letter)
}
m <- t(sapply(stands$byStrata_id, parse_id))
stands <- cbind(stands, as.data.frame(m, stringsAsFactors=FALSE))

# Cell rollup
cell_summary <- aggregate(BA_0 ~ forest_type + density_class,
                          data = stands, FUN = length)
names(cell_summary)[3] <- "n_stands"
ba_mean <- aggregate(BA_0 ~ forest_type + density_class,
                     data = stands, FUN = mean)
cell_summary$BA_0_mean <- round(ba_mean$BA_0, 1)
stand_ids <- aggregate(byStrata_id ~ forest_type + density_class,
                       data = stands,
                       FUN = function(x) paste(x, collapse=", "))
cell_summary$byStrata_ids <- stand_ids$byStrata_id
write.csv(cell_summary,
          file.path(od, "silc_strata_5x2_reconciliation_table.csv"),
          row.names = FALSE)

cat("=== 11 byStrata stands -> 10-cell reconciliation ===\n")
print(cell_summary, row.names = FALSE)

# === Figure: Sankey-like flow (left: 11 byStrata stands; right: 10 cells)
png(file.path(od, "silc_strata_5x2_reconciliation.png"),
    width = 2300, height = 1200, res = 165)
par(mar = c(2.0, 1.0, 3.2, 1.0))
plot(NA, xlim = c(0, 10), ylim = c(0, 12), axes = FALSE,
     xlab = "", ylab = "")

# Left column: 11 stands (sorted by forest_type to match cell order)
type_order <- c("Cedar","Hardwood","Mixedwood",
                 "Commercial Softwood","Other Softwood")
dens_order <- c("A+B (high)","C+D (low)")
stands$forest_type <- factor(stands$forest_type, levels = type_order)
stands$density_class <- factor(stands$density_class, levels = dens_order)
stands <- stands[order(stands$forest_type, stands$density_class), ]
left_y <- seq(11, 1, length.out = nrow(stands))
ftype_col <- c("Cedar" = "#7b3294", "Hardwood" = "#c2a5cf",
                "Mixedwood" = "#a6dba0",
                "Commercial Softwood" = "#008837",
                "Other Softwood" = "#5e4fa2")

# Right column: 10 cells in canonical order
cells <- expand.grid(forest_type = type_order,
                      density_class = dens_order,
                      stringsAsFactors = FALSE)
cells$forest_type <- factor(cells$forest_type, levels = type_order)
cells$density_class <- factor(cells$density_class, levels = dens_order)
cells <- cells[order(cells$forest_type, cells$density_class), ]
right_y <- seq(11, 1, length.out = nrow(cells))
cells$n <- mapply(function(ft, dc) {
  sum(stands$forest_type == ft & stands$density_class == dc)
}, cells$forest_type, cells$density_class)

# Draw left labels
text(0.05, 11.7, "byStrata stand (n=11)", cex = 1.05, font = 2, adj = 0)
for (i in seq_len(nrow(stands))) {
  ft <- as.character(stands$forest_type[i])
  rect(0.2, left_y[i] - 0.18, 2.0, left_y[i] + 0.18,
       col = ftype_col[ft], border = NA)
  text(1.1, left_y[i],
       sprintf("%s   BA=%.0f", stands$byStrata_id[i], stands$BA_0[i]),
       cex = 0.78, col = "white", font = 2)
}

# Draw right labels
text(7.5, 11.7, "5x2 cell (n=10, four empty)", cex = 1.05, font = 2, adj = 0)
for (j in seq_len(nrow(cells))) {
  ft <- as.character(cells$forest_type[j])
  empty <- cells$n[j] == 0
  rect(7.7, right_y[j] - 0.20, 9.95, right_y[j] + 0.20,
       col = if (empty) "#dddddd" else ftype_col[ft],
       border = if (empty) "#aa3333" else NA,
       lty = if (empty) 2 else 1)
  text(8.8, right_y[j],
       sprintf("%s / %s%s",
               cells$forest_type[j], cells$density_class[j],
               if (empty) " (no AGM data)"
               else sprintf(" (n=%d)", cells$n[j])),
       cex = 0.78,
       col = if (empty) "#aa3333" else "white",
       font = 2)
}

# Connector lines: each byStrata stand to its cell
for (i in seq_len(nrow(stands))) {
  ft <- as.character(stands$forest_type[i])
  dc <- as.character(stands$density_class[i])
  jrow <- which(as.character(cells$forest_type) == ft &
                as.character(cells$density_class) == dc)
  lines(c(2.05, 7.65), c(left_y[i], right_y[jrow]),
        col = ftype_col[ft], lwd = 1.4)
}

mtext("11 byStrata stands -> 5-type x 2-density 10-cell rollup  (AGM source)",
      side = 3, line = 1.4, cex = 1.15, font = 2)
mtext("Colour codes forest type; dashed grey-red blocks are cells with no AGM byStrata coverage",
      side = 3, line = 0.0, cex = 0.83, col = "#444")
dev.off()

cat("\nwrote silc_strata_5x2_reconciliation_table.csv\n")
cat("wrote silc_strata_5x2_reconciliation.png\n")
