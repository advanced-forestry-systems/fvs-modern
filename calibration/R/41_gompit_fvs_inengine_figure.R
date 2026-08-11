#!/usr/bin/env Rscript
# In-engine gompit vs native FVS mortality: AGB trajectories, NE / CS / LS.
#
# Reads the per-variant, per-mode validation CSVs produced by
# gompit_fvs_validate.py (native vs gompit, same FVS binary, env-toggled) and
# draws mean above-ground biomass over a 100-yr no-harvest projection. The
# point is that gompit mortality, substituted INSIDE the FVS growth loop, gives
# bounded, realistic biomass (no runaway), trimming late-rotation stocking as
# crown closure rises.
#
# Usage:
#   Rscript 41_gompit_fvs_inengine_figure.R --indir <dir> --out gompit_fvs_inengine.png

suppressWarnings(suppressMessages(library(tidyverse)))

args <- commandArgs(trailingOnly = TRUE)
ga <- function(f, d = NULL) {
  i <- grep(paste0("^--", f, "="), args, value = TRUE)
  if (!length(i)) return(d)
  sub(paste0("^--", f, "="), "", i[1])
}
INDIR <- ga("indir", ".")
OUT   <- ga("out", "gompit_fvs_inengine.png")

files <- list.files(INDIR, pattern = "^val_(ne|cs|ls|sn)_(native|gompit)\\.csv$",
                     full.names = TRUE)
dat <- map_dfr(files, read_csv, show_col_types = FALSE)

vlab <- c(ne = "NE  (Northeast)", cs = "CS  (Central States)",
          ls = "LS  (Lake States)", sn = "SN  (Southern, n=10*)")
summ <- dat %>%
  mutate(variant = factor(tolower(VARIANT), levels = c("ne","cs","ls","sn")),
         mode = factor(MODE, levels = c("native","gompit"),
                       labels = c("Native FVS","Gompit-in-FVS"))) %>%
  group_by(variant, mode, PROJ_YEAR) %>%
  summarise(AGB = mean(AGB), .groups = "drop")

col_nat <- "#999999"; col_gom <- "#0072B2"

p <- ggplot(summ, aes(PROJ_YEAR, AGB, colour = mode)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.6) +
  facet_wrap(~variant, nrow = 1, scales = "free_y",
             labeller = labeller(variant = vlab)) +
  scale_colour_manual(values = c("Native FVS" = col_nat,
                                 "Gompit-in-FVS" = col_gom), name = NULL) +
  labs(x = "Projection year", y = "Mean AGB (short tons / ac)",
       title = "Gompit mortality substituted inside the FVS growth loop",
       subtitle = "100-yr no-harvest projection, same binary, env-toggled; bounded and realistic, no runaway") +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top",
        panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))

ggsave(OUT, p, width = 11, height = 4.2, dpi = 300, bg = "white")
cat("wrote", OUT, "\n")
# terminal-year table
summ %>% filter(PROJ_YEAR == max(PROJ_YEAR)) %>%
  pivot_wider(names_from = mode, values_from = AGB) %>% print()
