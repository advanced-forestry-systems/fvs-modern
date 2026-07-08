## Analyze the default-vs-localized max-SDI side-by-side projection output.
## Reproduces the density-error-by-stand-density table and the two-panel figure
## from the paired FVS runs (short.csv = short-interval vs observed; traj.csv = 100-yr).
##
## Usage: Rscript analyze_sidebyside.R <demo_dir> [variant_label]
##   <demo_dir> holds short.csv and traj.csv from var_sdimax_sidebyside.py

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})

args <- commandArgs(trailingOnly = TRUE)
dir  <- if (length(args) >= 1) args[1] else "."
vlab <- if (length(args) >= 2) args[2] else "variant"

S <- read.csv(file.path(dir, "short.csv"))
T <- read.csv(file.path(dir, "traj.csv"))

## percent RMSE / bias vs observed
pstat <- function(pred, obs) {
  ok <- is.finite(pred) & is.finite(obs); e <- pred[ok] - obs[ok]; mo <- mean(obs[ok])
  c(rmse = 100 * sqrt(mean(e^2)) / mo, bias = 100 * mean(e) / mo, n = sum(ok))
}

strata <- list(
  "All stands"         = S,
  "Binding (RD>0.45)"  = dplyr::filter(S, RD_brms > 0.45),
  "Dense (RD>0.6)"     = dplyr::filter(S, RD_brms > 0.6)
)

cat(sprintf("\n=== %s: density (TPH) error vs observed ===\n", vlab))
tab <- lapply(names(strata), function(nm) {
  s <- strata[[nm]]
  d <- pstat(s$default_TPH,   s$oTPH)
  l <- pstat(s$localized_TPH, s$oTPH)
  data.frame(stratum = nm, n = d["n"],
             default_rmse = round(d["rmse"], 1), localized_rmse = round(l["rmse"], 1),
             default_bias = round(d["bias"], 1), localized_bias = round(l["bias"], 1))
}) %>% bind_rows()
print(tab, row.names = FALSE)

## year-100 per-stand density difference (localized - default)
end <- T %>% group_by(sid, mode) %>% arrange(year) %>%
  summarise(TPH_end = dplyr::last(TPH), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = mode, values_from = TPH_end) %>%
  mutate(dTPH100 = localized - default)
S2 <- dplyr::left_join(S, dplyr::select(end, sid, dTPH100), by = "sid")

## two-panel figure
c_def <- "#b0682f"; c_loc <- "#2f6db0"
pA_df <- tab %>%
  tidyr::pivot_longer(c(default_rmse, localized_rmse),
                      names_to = "source", values_to = "rmse") %>%
  mutate(source = ifelse(grepl("default", source), "Default (species-weighted)",
                         "Localized (FIA-derived)"),
         stratum = factor(stratum, levels = names(strata)))
pA <- ggplot(pA_df, aes(stratum, rmse, fill = source)) +
  geom_col(position = position_dodge(0.7), width = 0.65) +
  geom_text(aes(label = round(rmse)), position = position_dodge(0.7),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = c("Default (species-weighted)" = c_def,
                               "Localized (FIA-derived)" = c_loc)) +
  labs(x = NULL, y = "Density error vs observed (% RMSE)",
       title = "A. Localized max SDI cuts density error, most where the limit binds",
       fill = NULL) +
  theme_minimal(base_size = 11) + theme(legend.position = "top")

pB <- ggplot(S2, aes(RD_brms, dTPH100, color = RD_brms)) +
  geom_hline(yintercept = 0, color = "grey60") +
  geom_point(size = 2) +
  geom_smooth(method = "loess", se = FALSE, color = "#c0392b") +
  scale_color_viridis_c(guide = "none") +
  labs(x = "Relative density (SDI / localized max SDI)",
       y = "Year-100 density difference\nlocalized - default (trees/ha)",
       title = "B. The 100-year effect is real, signed, and grows with density") +
  theme_minimal(base_size = 11)

if (requireNamespace("patchwork", quietly = TRUE)) {
  fig <- patchwork::wrap_plots(pA, pB, ncol = 2)
} else { fig <- pA }
ggsave(file.path(dir, paste0(vlab, "_maxSDI_demo_R.png")), fig,
       width = 12, height = 5, dpi = 160)
cat(sprintf("\nFigure written to %s\n", file.path(dir, paste0(vlab, "_maxSDI_demo_R.png"))))
