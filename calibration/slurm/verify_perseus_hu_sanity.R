#!/usr/bin/env Rscript
# =============================================================================
# verify_perseus_hu_sanity.R
# -----------------------------------------------------------------------------
# Inspects the sanity batch outputs from perseus_harvest_uncertainty.py and
# reports whether the schema, row counts, AGB ranges, and per draw spread
# look right before scaling to the full 500 draw factorial.
#
# Usage:
#   Rscript verify_perseus_hu_sanity.R \
#     [--scen-dir /users/PUOM0008/crsfaaron/fvs-modern/calibration/output/perseus_hu_sanity/harvest_rcp45]
#     [--n-plots 50] [--n-draws 5] [--n-variants 2]
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

opt_list <- list(
  make_option("--scen-dir", type = "character",
              default = "/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/perseus_hu_sanity/harvest_rcp45",
              help = "Scenario subdir"),
  make_option("--n-plots", type = "integer", default = 50,
              help = "Plots in the sanity batch"),
  make_option("--n-draws", type = "integer", default = 5,
              help = "Draws per plot per variant"),
  make_option("--n-variants", type = "integer", default = 2,
              help = "Variants exercised (NE + ACD)"),
  make_option("--cycles", type = "integer", default = 21,
              help = "Cycles emitted (year 0 plus 20 FVS cycles of 5 yr)")
)
opt <- parse_args(OptionParser(option_list = opt_list))

scen_dir <- opt[["scen-dir"]]
scen_name <- basename(scen_dir)

cat("=== verify perseus_hu sanity ===\n")
cat("scenario dir: ", scen_dir, "\n")
cat("scenario:     ", scen_name, "\n\n")

# Locate batch CSVs
point_path <- file.path(scen_dir, sprintf("perseus_%s_unc_point_batch1.csv", scen_name))
draws_path <- file.path(scen_dir, sprintf("perseus_%s_unc_draws_batch1.csv", scen_name))
meta_path  <- file.path(scen_dir, sprintf("perseus_%s_unc_meta_batch1.json",  scen_name))

ok <- TRUE
for (p in c(point_path, draws_path)) {
  if (!file.exists(p)) {
    cat(sprintf("MISSING: %s\n", p))
    ok <- FALSE
  }
}
if (!ok) {
  cat("\nFAIL: required output files missing. Check sanity.err.\n")
  quit(status = 1)
}

point <- read_csv(point_path, show_col_types = FALSE)
draws <- read_csv(draws_path, show_col_types = FALSE)

# ----- Schema checks ---------------------------------------------------------
expected_point <- c("PLOT", "FIRST_PLTCN", "YEAR", "PROJ_YEAR",
                    "VARIANT", "CONFIG", "AGB_TONS_AC",
                    "SCENARIO", "RCP", "HARVEST", "SI_DELTA_FT")
expected_draws <- c("PLOT", "FIRST_PLTCN", "YEAR", "PROJ_YEAR",
                    "VARIANT", "DRAW", "AGB_TONS_AC",
                    "SCENARIO", "RCP", "HARVEST", "SI_DELTA_FT")

check_cols <- function(df, expected, label) {
  miss <- setdiff(expected, names(df))
  extra <- setdiff(names(df), expected)
  cat(sprintf("[%s] columns: %d", label, ncol(df)))
  if (length(miss)) cat(sprintf("  MISSING: %s", paste(miss, collapse = ", ")))
  if (length(extra)) cat(sprintf("  EXTRA: %s", paste(extra, collapse = ", ")))
  cat("\n")
  length(miss) == 0
}
ok_p <- check_cols(point, expected_point, "point")
ok_d <- check_cols(draws, expected_draws, "draws")

# ----- Row counts ------------------------------------------------------------
n_plots <- opt[["n-plots"]]
n_draws <- opt[["n-draws"]]
n_var   <- opt[["n-variants"]]
n_cyc   <- opt$cycles

expected_point_rows <- n_plots * n_var * 2 * n_cyc        # default + calibrated
expected_draws_rows <- n_plots * n_var * n_draws * n_cyc

cat(sprintf("\n[point] rows: %d (expected ~%d)\n", nrow(point), expected_point_rows))
cat(sprintf("[draws] rows: %d (expected ~%d)\n", nrow(draws), expected_draws_rows))

# Allow up to ~10% loss to FVS run failures or empty plots
tolerable <- function(actual, expected) actual >= 0.85 * expected
ok_pn <- tolerable(nrow(point), expected_point_rows)
ok_dn <- tolerable(nrow(draws), expected_draws_rows)
cat(sprintf("  point rows >= 85%% expected? %s\n", ok_pn))
cat(sprintf("  draws rows >= 85%% expected? %s\n", ok_dn))

# ----- Scenario tagging ------------------------------------------------------
scen_levels <- unique(point$SCENARIO)
cat(sprintf("\n[point] SCENARIO values: %s\n", paste(scen_levels, collapse = ", ")))
cat(sprintf("[point] RCP values:      %s\n", paste(unique(point$RCP), collapse = ", ")))
cat(sprintf("[point] HARVEST values:  %s\n", paste(unique(point$HARVEST), collapse = ", ")))

# ----- AGB ranges (sanity) ---------------------------------------------------
cat("\n[point] AGB_TONS_AC summary by VARIANT and CONFIG:\n")
point %>%
  group_by(VARIANT, CONFIG) %>%
  summarise(n = n(), mean = mean(AGB_TONS_AC, na.rm = TRUE),
            min = min(AGB_TONS_AC, na.rm = TRUE),
            max = max(AGB_TONS_AC, na.rm = TRUE),
            .groups = "drop") %>%
  print(n = Inf)

# ----- Draw spread: mean across draws should be near calibrated point --------
if (nrow(draws) > 0) {
  cat("\n[draws] mean AGB across draws by VARIANT and PROJ_YEAR (compared to point calibrated):\n")
  draws_m <- draws %>%
    group_by(VARIANT, PROJ_YEAR) %>%
    summarise(draw_mean = mean(AGB_TONS_AC, na.rm = TRUE),
              draw_sd = sd(AGB_TONS_AC, na.rm = TRUE),
              n_draws_obs = dplyr::n_distinct(DRAW),
              .groups = "drop")
  point_cal <- point %>%
    filter(CONFIG == "calibrated") %>%
    group_by(VARIANT, PROJ_YEAR) %>%
    summarise(point_cal_mean = mean(AGB_TONS_AC, na.rm = TRUE), .groups = "drop")
  comp <- left_join(draws_m, point_cal, by = c("VARIANT", "PROJ_YEAR")) %>%
    mutate(rel_diff_pct = round(100 * (draw_mean - point_cal_mean) / point_cal_mean, 1))
  print(comp, n = Inf)
}

# ----- ThinBBA harvest signal ------------------------------------------------
if (any(point$HARVEST)) {
  cat("\n[point] AGB drop signal by PROJ_YEAR (calibrated, harvest scenario):\n")
  hsig <- point %>%
    filter(CONFIG == "calibrated", HARVEST == TRUE) %>%
    group_by(VARIANT, PROJ_YEAR) %>%
    summarise(mean_agb = round(mean(AGB_TONS_AC), 2), .groups = "drop") %>%
    arrange(VARIANT, PROJ_YEAR)
  print(hsig, n = Inf)
}

# ----- Wall time (from meta) -------------------------------------------------
if (file.exists(meta_path)) {
  meta <- jsonlite::fromJSON(meta_path)
  cat(sprintf("\n[meta] wall_clock_s: %.1f (~%.1f min)\n",
              meta$wall_clock_s, meta$wall_clock_s / 60))
  per_plot <- meta$wall_clock_s / max(meta$n_plots_in_batch, 1)
  per_draw_per_plot <- per_plot / max(meta$n_draws_requested, 1)
  cat(sprintf("[meta] per-plot wall: %.1f s\n", per_plot))
  cat(sprintf("[meta] per-draw-per-plot wall: %.2f s\n\n", per_draw_per_plot))

  # Project full factorial cost
  cat("=== Projected full factorial cost ===\n")
  full_n_plots <- 3586
  full_n_draws <- 500
  n_scenarios <- 6
  n_variants_full <- 2
  per_plot_per_scenario_full_s <- (1 + 1 + full_n_draws) * n_variants_full * (per_draw_per_plot * 1.0)
  total_core_s <- per_plot_per_scenario_full_s * full_n_plots * n_scenarios
  cat(sprintf("  per plot per scenario (500 draws, 2 variants): %.0f s\n",
              per_plot_per_scenario_full_s))
  cat(sprintf("  total task time at 25 plots/batch and 24h cap: ~%.1f hr per batch\n",
              25 * per_plot_per_scenario_full_s / 3600))
  cat(sprintf("  total core hours across 6 scenarios: %.0f core-hr\n",
              total_core_s / 3600))
  cat(sprintf("  wall clock at 50 concurrent tasks: ~%.1f days\n",
              total_core_s / 3600 / 50 / 24))
}

cat("\n=== verdict ===\n")
verdict_ok <- all(c(ok_p, ok_d, ok_pn, ok_dn))
cat(if (verdict_ok) "PASS\n" else "INVESTIGATE\n")
quit(status = if (verdict_ok) 0 else 2)
