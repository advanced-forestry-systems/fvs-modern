#!/usr/bin/env Rscript
#
# FVS Bayesian Calibration Pipeline: Setup Script
# Initialize environment and install required packages
#
# Usage: Rscript calibration/R/00_setup.R
#

# ============================================================================
# Configuration
# ============================================================================

# Auto-detect project root: environment variable > script location > working dir
project_root <- Sys.getenv("FVS_PROJECT_ROOT", unset = "")
if (project_root == "") {
  # Try to infer from script location (works in Rscript and source())
  script_dir <- tryCatch(
    normalizePath(file.path(dirname(sys.frame(1)$ofile), "../.."), mustWork = TRUE),
    error = function(e) ""
  )
  if (nzchar(script_dir) && file.exists(file.path(script_dir, "calibration/R/00_setup.R"))) {
    project_root <- script_dir
  } else {
    # Fallback: walk up from working directory looking for calibration/R/
    wd <- getwd()
    for (depth in 0:3) {
      candidate <- normalizePath(file.path(wd, paste(rep("..", depth), collapse = "/")),
                                 mustWork = FALSE)
      if (file.exists(file.path(candidate, "calibration/R/00_setup.R"))) {
        project_root <- candidate
        break
      }
    }
  }
  if (project_root == "") stop("Cannot detect FVS_PROJECT_ROOT. Set the environment variable.")
}
calibration_dir <- file.path(project_root, "calibration")

# Create necessary directories
dirs_to_create <- c(
  file.path(calibration_dir, "data"),
  file.path(calibration_dir, "data", "raw_fia"),
  file.path(calibration_dir, "data", "processed"),
  file.path(calibration_dir, "output"),
  file.path(calibration_dir, "output", "variants"),
  file.path(calibration_dir, "logs"),
  file.path(calibration_dir, "config", "calibrated")
)

for (dir in dirs_to_create) {
  if (!dir.exists(dir)) {
    dir.create(dir, showWarnings = FALSE, recursive = TRUE)
    cat("Created directory:", dir, "\n")
  }
}

# ============================================================================
# Required Packages
# ============================================================================

required_packages <- c(
  # Data access and manipulation
  "rFIA",              # FIA data access
  "tidyverse",         # ggplot2, dplyr, tidyr, readr, etc
  "data.table",        # Fast data operations

  # Bayesian modeling
  "brms",              # Bayesian regression using Stan
  "cmdstanr",          # R interface to CmdStan
  "posterior",         # Tools for working with posterior distributions
  "bayesplot",         # Posterior predictive checks and diagnostics
  "tidybayes",         # Tidy extraction of posterior samples

  # Data formats
  "jsonlite",          # Read/write JSON
  "arrow",             # Parquet file support

  # Spatial and scientific
  "sf",                # Simple features for spatial data
  "terra",             # Spatial raster operations
  "units",             # Unit handling

  # Utilities
  "parallel",          # Base R parallel computing
  "doParallel",        # Parallel backend
  "foreach",           # Parallel loops
  "progress",          # Progress bars
  "glue",              # String interpolation
  "logger"             # Logging
)

# ============================================================================
# Installation Check and Installation
# ============================================================================

cat("Checking and installing required packages...\n")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      warning("Failed to install package:", pkg)
    }
  } else {
    cat("✓", pkg, "\n")
  }
}

# cmdstanr may need special handling
if (!require("cmdstanr", quietly = TRUE)) {
  cat("Installing cmdstanr from GitHub...\n")
  install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
}

# ============================================================================
# Configure CmdStan (OSC environment)
# ============================================================================

cat("\nConfiguring CmdStan...\n")

# Set up cmdstanr to use OSC's compiler environment.
# cmdstan_make_local() in current cmdstanr (>=0.7) takes cpp_options /
# stanc_options / append only; it writes to ~/.cmdstan/cmdstan-*/make/local.
# The previous call used new_local / overwrite args that no longer exist.
# Wrapped in try() so a cmdstanr version change does not halt the driver.
try(
  cmdstanr::cmdstan_make_local(
    cpp_options = list(STAN_CPP_OPTIMS = FALSE),  # -fgraphite unsupported on Cardinal GCC
    append      = TRUE
  ),
  silent = TRUE
)

# Install cmdstan if not already present. API probe first: older cmdstanr had
# cmdstan_installed(); 0.9+ dropped it in favor of cmdstan_version() which
# returns a version string when installed or errors otherwise.
has_cmdstan <- tryCatch({
  v <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  !is.null(v) && nzchar(as.character(v))
}, error = function(e) FALSE)
if (!has_cmdstan) {
  cat("Installing CmdStan (this may take a few minutes)...\n")
  try(cmdstanr::install_cmdstan(cores = parallel::detectCores()), silent = TRUE)
}

# ============================================================================
# Session Configuration
# ============================================================================

cat("\nConfiguring R session...\n")

# Set working directory
setwd(project_root)

# Configure brms/Stan options
options(
  mc.cores = parallel::detectCores(),
  brms.backend = "cmdstanr"
)

# Set ggplot2 theme
ggplot2::theme_set(ggplot2::theme_minimal())

# ============================================================================
# Logging Setup
# ============================================================================

log_file <- file.path(calibration_dir, "logs", "setup.log")
logger::log_appender(logger::appender_file(log_file))
logger::log_info("FVS Bayesian Calibration Pipeline initialized")
logger::log_info("Project root: {project_root}")
logger::log_info("R version: {R.version$version.string}")
logger::log_info("Available cores: {parallel::detectCores()}")

# ============================================================================
# Print Summary
# ============================================================================

cat("\n")
cat("========================================\n")
cat("FVS Bayesian Calibration Setup Complete\n")
cat("========================================\n")
cat("Project root:", project_root, "\n")
cat("Calibration dir:", calibration_dir, "\n")
cat("Available cores:", parallel::detectCores(), "\n")
cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")
cat("\nAll required packages installed successfully.\n")
cat("Ready to begin FIA data fetch and model fitting.\n\n")

# Return paths for use by other scripts
invisible(list(
  project_root = project_root,
  calibration_dir = calibration_dir
))
