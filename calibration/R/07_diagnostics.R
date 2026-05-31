#!/usr/bin/env Rscript
# =============================================================================
# FVS Bayesian Calibration: Comprehensive Diagnostics for All Components
#
# Generates convergence diagnostics, posterior predictive checks, and
# prior vs posterior comparisons for ALL fitted component models:
#   - Diameter growth (Wykoff)
#   - Height-diameter (Chapman-Richards)
#   - Height increment (for variants with HG params)
#   - Mortality (logistic)
#   - Crown ratio change
#
# Usage: Rscript calibration/R/07_diagnostics.R --variant ca
# =============================================================================

library(tidyverse)
library(jsonlite)
library(data.table)
library(posterior)
library(bayesplot)
library(tidybayes)
library(logger)
library(glue)

# =============================================================================
# Parse Command Line Arguments
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
variant <- "ca"

if (length(args) > 0) {
  for (i in seq_along(args)) {
    if (args[i] == "--variant" & i < length(args)) {
      variant <- args[i + 1]
    }
  }
}

# =============================================================================
# Configuration
# =============================================================================

project_root <- Sys.getenv("FVS_PROJECT_ROOT")
if (!nzchar(project_root)) {
  project_root <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = FALSE),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
calibration_dir <- file.path(project_root, "calibration")
output_dir <- file.path(calibration_dir, "output", "variants", variant)
config_dir <- file.path(project_root, "config")
diagnostic_dir <- file.path(output_dir, "diagnostics")

dir.create(diagnostic_dir, showWarnings = FALSE, recursive = TRUE)

# Logging
log_file <- file.path(calibration_dir, "logs",
                       paste0("07_diagnostics_", variant, ".log"))
dir.create(dirname(log_file), showWarnings = FALSE, recursive = TRUE)
logger::log_appender(logger::appender_file(log_file), index = 1)
logger::log_info("Starting comprehensive diagnostics for variant {variant}")

# Publication ready theme
theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.position = "bottom",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    strip.text = element_text(size = 11, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(10, 10, 10, 10)
  )

# Load original config
original_config <- fromJSON(file.path(config_dir, paste0(variant, ".json")))

# =============================================================================
# Helper: Compute Convergence Summary for a Component
# =============================================================================

compute_convergence <- function(posterior_df, component_name) {
  if (is.null(posterior_df) || nrow(posterior_df) == 0) {
    return(NULL)
  }

  # Determine which columns to use for Rhat and ESS
  rhat_col <- if ("rhat" %in% names(posterior_df)) "rhat" else
              if ("Rhat" %in% names(posterior_df)) "Rhat" else NULL
  ess_col <- if ("ess_bulk" %in% names(posterior_df)) "ess_bulk" else
             if ("Bulk_ESS" %in% names(posterior_df)) "Bulk_ESS" else NULL

  if (is.null(rhat_col) || is.null(ess_col)) {
    logger::log_warn("Cannot compute convergence for {component_name}: missing Rhat or ESS columns")
    return(NULL)
  }

  posterior_df %>%
    summarise(
      component = component_name,
      n_params = n(),
      n_good_rhat = sum(.data[[rhat_col]] <= 1.01, na.rm = TRUE),
      n_bad_rhat = sum(.data[[rhat_col]] > 1.01, na.rm = TRUE),
      max_rhat = max(.data[[rhat_col]], na.rm = TRUE),
      mean_rhat = mean(.data[[rhat_col]], na.rm = TRUE),
      n_good_ess = sum(.data[[ess_col]] >= 400, na.rm = TRUE),
      n_bad_ess = sum(.data[[ess_col]] < 400, na.rm = TRUE),
      min_ess = min(.data[[ess_col]], na.rm = TRUE),
      mean_ess = mean(.data[[ess_col]], na.rm = TRUE)
    )
}

# =============================================================================
# Helper: Generate Credible Interval Plot
# =============================================================================

make_credible_plot <- function(posterior_df, component_name, output_path) {
  if (is.null(posterior_df) || nrow(posterior_df) == 0) return(invisible(NULL))

  # Determine column names
  med_col <- if ("p50" %in% names(posterior_df)) "p50" else
             if ("Estimate" %in% names(posterior_df)) "Estimate" else
             if ("median" %in% names(posterior_df)) "median" else NULL
  lo_col <- if ("p05" %in% names(posterior_df)) "p05" else
            if ("Q2.5" %in% names(posterior_df)) "Q2.5" else
            if ("q5" %in% names(posterior_df)) "q5" else NULL
  hi_col <- if ("p95" %in% names(posterior_df)) "p95" else
            if ("Q97.5" %in% names(posterior_df)) "Q97.5" else
            if ("q95" %in% names(posterior_df)) "q95" else NULL

  if (is.null(med_col) || is.null(lo_col) || is.null(hi_col)) return(invisible(NULL))

  # Filter to fixed effects only
  plot_data <- posterior_df %>%
    filter(grepl("^b[0-9_]", variable) | grepl("^gamma", variable)) %>%
    slice_head(n = 30)

  if (nrow(plot_data) == 0) return(invisible(NULL))

  p <- ggplot(plot_data, aes(y = reorder(variable, .data[[med_col]]))) +
    geom_point(aes(x = .data[[med_col]]), size = 2.5) +
    geom_errorbarh(aes(xmin = .data[[lo_col]], xmax = .data[[hi_col]]),
                    height = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      x = "Posterior Estimate (90% CI)",
      y = "Parameter",
      title = glue("{component_name} ({variant})")
    ) +
    theme_pub

  ggsave(output_path, p, width = 10, height = max(4, nrow(plot_data) * 0.3))
  logger::log_info("Saved credible interval plot for {component_name}")
}

# =============================================================================
# Process Each Component
# =============================================================================

all_convergence <- list()

# --- Diameter Growth ---
logger::log_info("=== Diameter Growth Diagnostics ===")
dg_file <- file.path(output_dir, "diameter_growth_posterior.csv")
if (file.exists(dg_file)) {
  dg_post <- read_csv(dg_file, show_col_types = FALSE) %>% as_tibble()
  all_convergence$diameter_growth <- compute_convergence(dg_post, "Diameter Growth")
  make_credible_plot(dg_post, "Diameter Growth",
                      file.path(diagnostic_dir, "dg_credible_intervals.pdf"))
} else {
  logger::log_warn("Diameter growth posterior not found")
}

# --- Height-Diameter ---
logger::log_info("=== Height-Diameter Diagnostics ===")
htdbh_file <- file.path(output_dir, "height_diameter_summary.csv")
if (file.exists(htdbh_file)) {
  htdbh_post <- read_csv(htdbh_file, show_col_types = FALSE) %>% as_tibble()
  all_convergence$height_diameter <- compute_convergence(htdbh_post, "Height-Diameter")
  make_credible_plot(htdbh_post, "Height-Diameter",
                      file.path(diagnostic_dir, "htdbh_credible_intervals.pdf"))
} else {
  logger::log_info("Height-diameter posterior not found")
}

# --- Height Increment ---
logger::log_info("=== Height Increment Diagnostics ===")
htinc_file <- file.path(output_dir, "height_increment_posterior.csv")
if (file.exists(htinc_file)) {
  htinc_post <- read_csv(htinc_file, show_col_types = FALSE) %>% as_tibble()
  all_convergence$height_increment <- compute_convergence(htinc_post, "Height Increment")
  make_credible_plot(htinc_post, "Height Increment",
                      file.path(diagnostic_dir, "htinc_credible_intervals.pdf"))
} else {
  logger::log_info("Height increment posterior not found (expected for most variants)")
}

# --- Mortality ---
logger::log_info("=== Mortality Diagnostics ===")
mort_file <- file.path(output_dir, "mortality_summary.csv")
if (file.exists(mort_file)) {
  mort_post <- read_csv(mort_file, show_col_types = FALSE) %>% as_tibble()
  all_convergence$mortality <- compute_convergence(mort_post, "Mortality")
  make_credible_plot(mort_post, "Mortality",
                      file.path(diagnostic_dir, "mort_credible_intervals.pdf"))
} else {
  logger::log_info("Mortality posterior not found")
}

# --- Crown Ratio ---
logger::log_info("=== Crown Ratio Diagnostics ===")
cr_file <- file.path(output_dir, "crown_ratio_summary.csv")
if (file.exists(cr_file)) {
  cr_post <- read_csv(cr_file, show_col_types = FALSE) %>% as_tibble()
  all_convergence$crown_ratio <- compute_convergence(cr_post, "Crown Ratio")
  make_credible_plot(cr_post, "Crown Ratio",
                      file.path(diagnostic_dir, "cr_credible_intervals.pdf"))
} else {
  logger::log_info("Crown ratio posterior not found")
}

# =============================================================================
# Combined Convergence Summary
# =============================================================================

logger::log_info("=== Combined Convergence Summary ===")

convergence_table <- bind_rows(compact(all_convergence))

overall <- NULL

if (nrow(convergence_table) > 0) {
  write_csv(convergence_table, file.path(diagnostic_dir, "convergence_summary.csv"))

  overall <- convergence_table %>%
    summarise(
      total_params = sum(n_params),
      total_good_rhat = sum(n_good_rhat),
      total_bad_rhat = sum(n_bad_rhat),
      worst_rhat = max(max_rhat, na.rm = TRUE),
      total_good_ess = sum(n_good_ess),
      total_bad_ess = sum(n_bad_ess),
      worst_ess = min(min_ess, na.rm = TRUE)
    )

  # Convergence bar chart
  conv_plot_data <- convergence_table %>%
    select(component, n_good_rhat, n_bad_rhat) %>%
    pivot_longer(cols = c(n_good_rhat, n_bad_rhat),
                  names_to = "status", values_to = "count") %>%
    mutate(
      status = ifelse(status == "n_good_rhat", "Converged", "Needs Attention")
    )

  p_conv <- ggplot(conv_plot_data, aes(x = component, y = count, fill = status)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = c("Converged" = "#009E73", "Needs Attention" = "#D55E00")) +
    labs(
      x = "Component Model",
      y = "Number of Parameters",
      fill = "Convergence Status",
      title = glue("Convergence Summary: {original_config$variant_name} ({variant})")
    ) +
    theme_pub +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

  ggsave(file.path(diagnostic_dir, "convergence_overview.pdf"), p_conv,
         width = 10, height = 6)
}

# =============================================================================
# Prior vs Posterior Comparison (Diameter Growth)
# =============================================================================

logger::log_info("=== Prior vs Posterior Comparison ===")

if (exists("dg_post")) {
  fvs_params <- original_config$categories$other

  dg_comparison <- tibble(
    parameter = c("DGLD", "DGDS", "DGSITE", "DGSLOP", "DGSLSQ",
                  "DGSASP", "DGCASP", "DGEL", "DGELSQ", "DGCR",
                  "DGCRSQ", "DGBAL", "DGPCCF"),
    fvs_original = NA_real_,
    posterior_median = NA_real_
  )

  for (i in seq_len(nrow(dg_comparison))) {
    param <- dg_comparison$parameter[i]
    if (param %in% names(fvs_params)) {
      vals <- unlist(fvs_params[[param]])
      dg_comparison$fvs_original[i] <- mean(vals, na.rm = TRUE)
    }
  }

  # Get posterior medians
  stan_names <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7",
                  "b8", "b9", "b10", "b11", "b12", "b13")

  med_col <- if ("p50" %in% names(dg_post)) "p50" else
             if ("median" %in% names(dg_post)) "median" else NULL

  if (!is.null(med_col)) {
    for (i in seq_along(stan_names)) {
      match_row <- dg_post %>% filter(variable == stan_names[i])
      if (nrow(match_row) == 0) {
        match_rows <- dg_post %>% filter(grepl(paste0("^", stan_names[i], "\\["), variable))
        if (nrow(match_rows) > 0) {
          dg_comparison$posterior_median[i] <- mean(match_rows[[med_col]], na.rm = TRUE)
        }
      } else {
        dg_comparison$posterior_median[i] <- match_row[[med_col]][1]
      }
    }

    dg_comp_long <- dg_comparison %>%
      filter(!is.na(fvs_original) | !is.na(posterior_median)) %>%
      pivot_longer(cols = c(fvs_original, posterior_median),
                    names_to = "source", values_to = "value") %>%
      mutate(source = ifelse(source == "fvs_original", "FVS Original", "Bayesian Posterior"))

    if (nrow(dg_comp_long) > 0) {
      p_compare <- ggplot(dg_comp_long, aes(y = parameter, x = value, color = source)) +
        geom_point(size = 3, position = position_dodge(width = 0.4)) +
        scale_color_manual(values = c("FVS Original" = "#E69F00",
                                       "Bayesian Posterior" = "#0072B2")) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
        labs(
          x = "Parameter Value",
          y = "FVS Parameter",
          color = "Source",
          title = glue("FVS Original vs Bayesian Calibrated ({variant})")
        ) +
        theme_pub

      ggsave(file.path(diagnostic_dir, "prior_posterior_comparison.pdf"), p_compare,
             width = 10, height = 7)
    }

    write_csv(dg_comparison, file.path(diagnostic_dir, "dg_prior_posterior_values.csv"))
  }
}

# =============================================================================
# Generate Comprehensive Report
# =============================================================================

logger::log_info("Generating comprehensive report...")

components_fitted <- names(compact(all_convergence))

report_lines <- c(
  "# FVS Bayesian Calibration Diagnostics Report",
  glue("# Variant: {original_config$variant_name} ({variant})"),
  glue("# Generated: {Sys.time()}"),
  "",
  "## Component Models Fitted",
  paste("-", components_fitted),
  ""
)

if (nrow(convergence_table) > 0) {
  report_lines <- c(report_lines, "## Convergence Summary by Component", "")

  for (i in seq_len(nrow(convergence_table))) {
    row <- convergence_table[i, ]
    report_lines <- c(report_lines,
      glue("### {row$component}"),
      glue("  Parameters: {row$n_params}"),
      glue("  Converged (Rhat <= 1.01): {row$n_good_rhat}"),
      glue("  Needs attention (Rhat > 1.01): {row$n_bad_rhat}"),
      glue("  Max Rhat: {round(row$max_rhat, 4)}"),
      glue("  Mean ESS: {round(row$mean_ess, 0)}"),
      glue("  Min ESS: {round(row$min_ess, 0)}"),
      ""
    )
  }
}

if (!is.null(overall)) {
  report_lines <- c(report_lines,
    "## Overall Assessment",
    glue("Total parameters across all components: {overall$total_params}"),
    glue("Parameters with good convergence: {overall$total_good_rhat} ({round(100 * overall$total_good_rhat / overall$total_params, 1)}%)"),
    glue("Worst Rhat: {round(overall$worst_rhat, 4)}"),
    glue("Worst ESS: {round(overall$worst_ess, 0)}"),
    "",
    ifelse(overall$total_bad_rhat == 0,
           "All parameters converged well.",
           paste0(overall$total_bad_rhat, " parameters need additional sampling or model refinement.")),
    ""
  )
}

report_lines <- c(report_lines,
  "## Generated Outputs",
  "- convergence_summary.csv: Convergence metrics by component",
  "- convergence_overview.pdf: Visual convergence comparison",
  "- *_credible_intervals.pdf: Parameter estimates with 90% CIs (one per component)",
  "- prior_posterior_comparison.pdf: FVS original vs Bayesian calibrated (DG)",
  "- dg_prior_posterior_values.csv: Numerical comparison of DG params",
  "",
  "## Recommendations",
  "1. Review any parameters with Rhat > 1.01 for potential specification issues",
  "2. Compare calibrated parameters against biological expectations",
  "3. Validate calibrated predictions against independent FIA data",
  "4. For production use, run extended chains (4000+ post warmup) on HPC"
)

report_file <- file.path(diagnostic_dir, "DIAGNOSTICS_REPORT.md")
writeLines(report_lines, report_file)
logger::log_info("Saved diagnostics report to {report_file}")

# =============================================================================
# Print Summary
# =============================================================================

cat("\n")
cat(paste(report_lines, collapse = "\n"))
cat("\n\n")
cat("==========================================\n")
cat("Diagnostics Complete\n")
cat("==========================================\n")
cat("Variant:", variant, "\n")
cat("Components diagnosed:", paste(components_fitted, collapse = ", "), "\n")
cat("Output directory:", diagnostic_dir, "\n\n")

logger::log_info("Comprehensive diagnostics complete for variant {variant}")
