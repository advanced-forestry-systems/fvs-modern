##=============================================================================
## eval_gamma_heatmap.R
##
## Build a multi-panel γ heatmap from modifier_traitmed fits.
## Rows = 8 traits, columns = 8 modifier types, one panel per component.
## Color = γ_mean; significance overlay = filled point if 90% CI excludes zero.
##
## Inputs: gamma_summary CSVs from one or more modifier_traitmed fits.
## CLI:
##   --csv_dir=PATH (directory with *_gamma_summary.csv files)
##   --components=DG,HG,HCB,CR (comma-separated, optional)
##   --out=FILE.png (output figure)
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

CSV_DIR <- get_arg("csv_dir", "calibration/output/conus")
OUT_FILE <- get_arg("out", "modifier_traitmed_gamma_heatmap.png")

# Hard-coded paths to traitmed gamma summaries
SUMMARIES <- list(
  list(component = "DG_Kuehne v8",
       path = file.path(CSV_DIR, "dg/modifier_traitmed_v8/dg_kuehne_v8_traitmed_lambda10_gamma_summary.csv")),
  list(component = "HG_Organon v5",
       path = file.path(CSV_DIR, "hg/modifier_traitmed_v5_for_loo/hg_organon_v5_traitmed_lambda10_gamma_summary.csv")),
  list(component = "HCB",
       path = file.path(CSV_DIR, "hcb/modifier_traitmed/hcb_speciesfree_traitmed_lambda10_gamma_summary.csv")),
  list(component = "CR",
       path = file.path(CSV_DIR, "cr/modifier_traitmed/cr_speciesfree_traitmed_lambda10_gamma_summary.csv")),
  list(component = "Mortality",
       path = file.path(CSV_DIR, "mort/modifier_traitmed_plotlevel/mort_speciesfree_traitmed_plotlevel_lambda10_gamma_summary.csv")),
  list(component = "Ingrowth",
       path = file.path(CSV_DIR, "ingrowth/modifier_traitmed_plotlevel/ingrowth_v4_traitmed_plotlevel_lambda10_gamma_summary.csv"))
)

cat("Building gamma heatmap...\n")
rows <- list()
for (s in SUMMARIES) {
  if (!file.exists(s$path)) {
    cat("  MISSING:", s$path, "\n")
    next
  }
  d <- fread(s$path)
  d[, component := s$component]
  d[, modifier := sub("gamma_alpha_(.+)\\[.*", "\\1", variable)]
  d[, trait_idx := as.integer(sub(".*\\[(\\d+)\\]", "\\1", variable))]
  # Significant if 90% CI excludes zero
  d[, sig := q5 > 0 | q95 < 0]
  rows[[length(rows) + 1]] <- d
  cat("  loaded:", s$component, "n =", nrow(d), "\n")
}

dt <- rbindlist(rows, fill = TRUE)

# Order traits consistently
trait_order <- c("wood_specific_gravity", "shade_tolerance_num", "softwood",
                 "leaf_longevity_months", "max_ht_m", "max_dbh_cm",
                 "vulnerability_score", "sensitivity")
if ("trait" %in% names(dt)) {
  dt[, trait_label := trait]
} else {
  dt[, trait_label := trait_order[trait_idx]]
}

mod_order <- c("plant", "fire", "insect", "disease", "wind",
               "harvest", "cutting", "siteprep")
dt[, modifier := factor(modifier, levels = mod_order)]
dt[, trait_label := factor(trait_label, levels = rev(trait_order))]
dt[, component := factor(component,
                          levels = c("DG_Kuehne v8", "HG_Organon v5",
                                     "HCB", "CR", "Mortality", "Ingrowth"))]

# Build heatmap
p <- ggplot(dt, aes(x = modifier, y = trait_label, fill = mean)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_point(data = dt[sig == TRUE], aes(x = modifier, y = trait_label),
             color = "black", size = 1.6, shape = 16, inherit.aes = FALSE) +
  facet_wrap(~ component, ncol = 3, scales = "fixed") +
  scale_fill_gradient2(low = "#3B7DC4", mid = "white", high = "#D45050",
                        midpoint = 0, name = expression(gamma)) +
  labs(x = "Modifier type",
       y = "Species trait",
       title = "Trait-mediated species heterogeneity in modifier coefficients",
       subtitle = "Dots indicate 90% CI excludes zero") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold", size = 11),
        legend.position = "right")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
ggsave(OUT_FILE, p, width = 11, height = 8, dpi = 150)
cat("\nSaved:", OUT_FILE, "\n")

# Also save the table
tbl_out <- sub("\\.png$", "_table.csv", OUT_FILE)
fwrite(dt[, .(component, modifier, trait_label, mean, q5, q95, sig)],
       tbl_out)
cat("Saved table:", tbl_out, "\n")

# Summary
cat("\n=== Summary: significant γ pairs per component ===\n")
print(dt[, .(n_sig = sum(sig)), by = component])
