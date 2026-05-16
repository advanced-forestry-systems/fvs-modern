#!/usr/bin/env Rscript
################################################################################
# 27_acd_post_pass.R
# ------------------------------------------------------------------------------
# Run FVS-ACD calibrated and default projections on NE-tagged plots from the
# existing observed_changes.csv intermediate file. FVS-ACD is a subvariant of
# FVS-NE: it shares the NE Fortran binary and falls back to NE calibrated
# parameters where ACD-specific posteriors are unavailable.
#
# Reads:
#   intermediate/observed_changes.csv           (from a prior engine run)
#   calibration/output/variants/acd/*.rds       (ACD posteriors, if present)
#   calibration/output/variants/ne/*.rds        (NE fallback posteriors)
#
# Writes:
#   intermediate/validation_data_acd_post.csv   (per-plot residuals)
#   manuscript_tables/fia_benchmark_results_acd_row.csv  (single ACD row)
#
# Optionally appends to manuscript_tables/fia_benchmark_results.csv.
################################################################################

suppressPackageStartupMessages({
  library(data.table)
})

project_root <- Sys.getenv("FVS_PROJECT_ROOT",
                            unset = "/users/PUOM0008/crsfaaron/fvs-modern")
fia_root <- Sys.getenv("FVS_FIA_DATA_DIR",
                        unset = file.path(dirname(project_root), "FIA"))
output_root <- file.path(project_root, "calibration/output/comparisons")
inter_dir <- file.path(output_root, "intermediate")

cat("FVS-ACD post-pass on NE plots\n")
cat(strrep("=", 80), "\n\n")
cat("Project root:", project_root, "\n")
cat("FIA root:", fia_root, "\n\n")

# ------------------------------------------------------------------------------
# Engine sourcing: extract the function-definition prefix + the validation
# metrics block. Robust to engine refactors: each marker is looked up by name
# and verified for ordering before slicing.
# ------------------------------------------------------------------------------

engine_path <- file.path(project_root, "calibration/R/19_fia_benchmark_engine.R")
cat("Sourcing engine for projection + metric functions...\n")
cat("  Engine path:", engine_path, "\n")
engine_text <- readLines(engine_path)

# Find marker line; stop with clear message if absent.
find_marker <- function(pattern, label) {
  hit <- grep(pattern, engine_text)
  if (length(hit) == 0)
    stop("Engine missing marker '", label, "' (pattern: ", pattern,
         "). Engine may have been refactored.", call. = FALSE)
  hit[1]
}

# Function-definitions prefix: helper functions, STEP 1 (variant params load),
# and STEP 2 (projection functions). Terminating just before STEP 2: Data
# Assembly keeps it cheap to source (no FIA reads, no raster lookup).
fn_start  <- find_marker("^# Helper Functions",       "Helper Functions")
fn_end    <- find_marker("^# STEP 2: Data Assembly",  "STEP 2: Data Assembly") - 1L

# Validation-metrics block: the calc_* / compute_* helpers used downstream.
# Driver loop starts at '# Compute by variant', so slice ends just before it.
stats_start <- find_marker("^# --- Validation metrics ---", "Validation metrics")
stats_end   <- find_marker("^# Compute by variant",         "Compute by variant") - 1L

# Sanity: each slice must be a forward-going, non-empty range.
stopifnot(fn_start  < fn_end)
stopifnot(stats_start < stats_end)

funcs_part <- paste(engine_text[fn_start:fn_end], collapse = "\n")
stats_part <- paste(engine_text[stats_start:stats_end], collapse = "\n")

cat("  Slicing prefix    : lines", fn_start,  "to", fn_end,
    "(", fn_end  - fn_start  + 1L, "lines )\n")
cat("  Slicing stats     : lines", stats_start, "to", stats_end,
    "(", stats_end - stats_start + 1L, "lines )\n")

# Hoist symbols the engine prefix expects but that live OUTSIDE the slice
# (defined before line `load_variant_params <-`). The prefix references
# `calib_root`, `data_root`, mortality/CI toggles, etc., and runs the eager
# `for (v in ALL_VARIANTS)` loop. Restricting ALL_VARIANTS to {acd, ne} keeps
# that loop fast (don't touch every variant's posterior directory).
ALL_VARIANTS            <- c("acd", "ne")
MORT_STRATEGY           <- "no_si"
# INGROWTH: prefer a cached lookup written by the engine's STEP 2b (added in
# the 2026-05-15 patch). When it's present we get full-fidelity ingrowth; when
# absent (older engine runs) we silently disable and document below.
ingrowth_cache_path     <- file.path(project_root,
                                     "calibration/output/comparisons",
                                     "intermediate/ingrowth_lookup.rds")
if (file.exists(ingrowth_cache_path)) {
  ingrowth_lookup    <- readRDS(ingrowth_cache_path)
  INGROWTH_ENABLED   <- TRUE
  cat("  Ingrowth lookup loaded from cache (", length(ingrowth_lookup), "entries)\n")
} else {
  ingrowth_lookup    <- NULL
  INGROWTH_ENABLED   <- FALSE
  cat("  Ingrowth lookup NOT cached; running without ingrowth contribution.\n")
  cat("  (To enable, re-run 19_fia_benchmark_engine.R so STEP 2b writes the cache.)\n")
}
CLIMATE_SI_ENABLED      <- TRUE
SDIMAX_RASTER_ENABLED   <- TRUE
DG_BACKTRANSFORM        <- "median"
CI_BRACKETS_ENABLED     <- TRUE
BOOTSTRAP_CI_ENABLED    <- TRUE
BOOTSTRAP_N             <- 1000
DG_VARIANT_MULT_ENABLED <- FALSE
DG_VARIANT_MULT         <- c()
OUTPUT_TAG              <- ""
calib_root        <- file.path(project_root, "calibration/output/variants")
data_root         <- file.path(project_root, "calibration/data/processed")
RASTER_LOOKUP_CSV <- file.path(project_root, "calibration/data/plot_raster_lookup.csv")

cat("  Evaluating function definitions...\n")
eval(parse(text = paste(funcs_part, stats_part, sep = "\n\n")))

required_fns <- c("load_variant_params", "project_condition_calibrated",
                  "project_condition_default", "compute_metrics",
                  "bootstrap_stat_ci", "compute_boot_ci",
                  "compute_variant_stats")
missing_fns <- required_fns[!vapply(required_fns, exists, logical(1))]
if (length(missing_fns) > 0)
  stop("Engine sourcing did not produce required functions: ",
       paste(missing_fns, collapse = ", "), call. = FALSE)
cat("  Functions available:", paste(required_fns, collapse = ", "), "\n\n")

# ------------------------------------------------------------------------------
# Toggles above are set before eval so the engine prefix can reference them.
# If you need to tweak any, do so above this block, not here.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Bridge: load ACD params; fall back to NE params for any component missing
# from the ACD posterior set. This is the "subvariant of NE" contract.
# ------------------------------------------------------------------------------

cat("Loading ACD calibrated parameters (with NE fallback)...\n")
acd_params <- load_variant_params("acd")
ne_params  <- load_variant_params("ne")

component_status <- function(p) {
  c(dg = !is.null(p$dg), mort = !is.null(p$mort),
    hd = !is.null(p$hd), htg = !is.null(p$htg))
}
acd_have <- component_status(acd_params)
ne_have  <- component_status(ne_params)

fallback_components <- character(0)
for (k in names(acd_have)) {
  if (!acd_have[[k]] && ne_have[[k]]) {
    acd_params[[k]]              <- ne_params[[k]]
    acd_params[[paste0(k, "_lo")]] <- ne_params[[paste0(k, "_lo")]]
    acd_params[[paste0(k, "_hi")]] <- ne_params[[paste0(k, "_hi")]]
    fallback_components <- c(fallback_components, k)
  }
}
if (length(fallback_components) > 0) {
  warning("ACD posterior missing for ", paste(fallback_components, collapse = ", "),
          "; falling back to NE.", call. = FALSE)
  cat("  Bridge fallback (ACD <- NE) for:",
      paste(fallback_components, collapse = ", "), "\n")
}

calibration_source <- if (length(fallback_components) == length(acd_have)) "NE_full" else
                       if (length(fallback_components) > 0)            "ACD_partial" else
                                                                         "ACD_native"
final_have <- component_status(acd_params)
cat("  Final components available:",
    paste(names(final_have)[final_have], collapse = "+"),
    "[", calibration_source, "]\n\n")

# Guard rail: refuse to run with completely empty params.
if (!any(final_have))
  stop("Neither ACD nor NE posteriors are available; cannot run post-pass. ",
       "Restore calibration/output/variants/{acd,ne}/ from archive first.",
       call. = FALSE)

# ------------------------------------------------------------------------------
# Read intermediate observed_changes.csv (built by the engine) and subset to
# the rows we want to re-project under ACD. After a relabel run, NE rows are
# already retagged ACD; before that, ACD rows are NE-tagged. Handle both.
# ------------------------------------------------------------------------------

obs_path <- file.path(inter_dir, "observed_changes.csv")
cat("Reading observed_changes.csv ...\n")
matched <- fread(obs_path)
cat("  Total matched conditions:", nrow(matched), "\n")
present <- sort(unique(matched$VARIANT))
cat("  Variants present:", paste(present, collapse = ", "), "\n")

if ("ACD" %in% present) {
  acd_subset <- matched[VARIANT == "ACD"]
  cat("  Using ACD-tagged conditions (engine ran with FVS_ACD_RELABEL):",
      nrow(acd_subset), "\n\n")
} else if ("NE" %in% present) {
  acd_subset <- matched[VARIANT == "NE"]
  cat("  Treating NE-tagged conditions as ACD candidates:",
      nrow(acd_subset), "\n\n")
} else {
  stop("observed_changes.csv has neither ACD nor NE rows; nothing to project.",
       call. = FALSE)
}

# ------------------------------------------------------------------------------
# Load FIA tree data once, then project each condition.
# ------------------------------------------------------------------------------

cat("Reading ENTIRE_TREE.csv for tree-level lookup...\n")
tree_cols <- c("CN", "PLT_CN", "INVYR", "STATECD", "CONDID", "STATUSCD",
                "SPCD", "DIA", "HT", "CR", "TPA_UNADJ",
                "VOLCFGRS", "VOLCFNET", "VOLBFNET",
                if (INGROWTH_ENABLED) "PREV_TRE_CN" else NULL)
fia_trees <- fread(file.path(fia_root, "ENTIRE_TREE.csv"),
                    select = tree_cols, showProgress = FALSE)
fia_trees <- fia_trees[STATUSCD == 1]
setkey(fia_trees, PLT_CN, CONDID)
cat("  Live trees loaded:", nrow(fia_trees), "\n\n")

cat("Projecting", nrow(acd_subset), "conditions with ACD parameters...\n")
projection_results <- vector("list", nrow(acd_subset))
t_start <- proc.time()
for (i in seq_len(nrow(acd_subset))) {
  row <- acd_subset[i]
  trees <- fia_trees[PLT_CN == row$PLT_CN_t1 & CONDID == row$CONDID]
  pr_calib <- project_condition_calibrated(
    trees = trees,
    params = acd_params,
    interval_years = row$interval_years,
    si = row$FVS_SITE_INDEX,
    sdimax_raster_val = row$SDIMAX_imperial
  )
  pr_default <- project_condition_default(
    trees = trees,
    interval_years = row$interval_years,
    variant_code = "ACD"
  )
  projection_results[[i]] <- list(
    PLT_CN_t1 = row$PLT_CN_t1,
    CONDID = row$CONDID,
    VARIANT = "ACD",
    BA_pred_calib = pr_calib$BA_pred,
    TPA_pred_calib = pr_calib$TPA_pred,
    QMD_pred_calib = pr_calib$QMD_pred,
    SDI_pred_calib = pr_calib$SDI_pred,
    VOL_CFGRS_pred_calib = pr_calib$VOL_CFGRS_pred,
    VOL_CFNET_pred_calib = pr_calib$VOL_CFNET_pred,
    VOL_BFNET_pred_calib = pr_calib$VOL_BFNET_pred,
    HT_top_calib = pr_calib$HT_top,
    BA_pred_default = pr_default$BA_pred,
    TPA_pred_default = pr_default$TPA_pred,
    QMD_pred_default = pr_default$QMD_pred,
    SDI_pred_default = pr_default$SDI_pred,
    VOL_CFGRS_pred_default = pr_default$VOL_CFGRS_pred,
    VOL_CFNET_pred_default = pr_default$VOL_CFNET_pred,
    VOL_BFNET_pred_default = pr_default$VOL_BFNET_pred,
    HT_top_default = pr_default$HT_top
  )
  if (i %% 5000 == 0) {
    elapsed <- (proc.time() - t_start)[3]
    rate <- i / elapsed
    cat(sprintf("  %d / %d (%.1f cond/sec)\n", i, nrow(acd_subset), rate))
  }
}
elapsed <- (proc.time() - t_start)[3]
cat(sprintf("Done. %.1f sec total (%.1f cond/sec)\n\n",
            elapsed, nrow(acd_subset) / elapsed))

# Combine projection results into validation data
pr_dt <- rbindlist(projection_results)
vd <- merge(acd_subset, pr_dt, by = c("PLT_CN_t1", "CONDID"), all.x = TRUE)
vd[, VARIANT := "ACD"]  # override regardless of source tag
vd[, calibration_source := calibration_source]
cat("Validation data assembled:", nrow(vd), "rows\n")
n_valid <- sum(!is.na(vd$BA_pred_calib) & !is.na(vd$BA_pred_default))
cat("  Non-NA prediction pairs:", n_valid, "/", nrow(vd), "\n\n")

# Save
fwrite(vd, file.path(inter_dir, "validation_data_acd_post.csv"))

if (n_valid == 0) {
  warning("Zero non-NA prediction pairs; downstream stats will be NA. ",
          "Check that calibration/output/variants/{acd,ne}/ has the posterior RDS files.",
          call. = FALSE)
}

# Compute statistics
cat("Computing ACD variant statistics...\n")
stats_list <- compute_variant_stats(vd, "ACD")
stats_dt <- as.data.table(stats_list)
fwrite(stats_dt,
        file.path(output_root, "manuscript_tables/fia_benchmark_results_acd_row.csv"))

# Append to fia_benchmark_results.csv if present and missing the ACD row
main_path <- file.path(output_root, "manuscript_tables/fia_benchmark_results.csv")
if (file.exists(main_path)) {
  main <- fread(main_path)
  if (!"ACD" %in% main$VARIANT) {
    common_cols <- intersect(names(main), names(stats_dt))
    new_row <- stats_dt[, ..common_cols]
    for (col in setdiff(names(main), names(new_row))) {
      new_row[[col]] <- NA
    }
    setcolorder(new_row, names(main))
    main <- rbind(main, new_row)
    fwrite(main, main_path)
    cat("Appended ACD row to", main_path, "\n")
  }
}

cat("\nDone. Headline:\n")
cat(sprintf("  ACD n=%d  BA R^2_calib=%.3f  RMSE=%.1f%%  bias=%.1f%%  [%s]\n",
            stats_dt$n_conditions,
            stats_dt$BA_r2_calib,
            stats_dt$BA_RMSE_pct_calib,
            stats_dt$BA_bias_pct_calib,
            calibration_source))
