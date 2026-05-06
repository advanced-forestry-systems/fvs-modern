## =============================================================================
## 60b_pivot_modifier_table.R
##
## Wide-format pivot of cross-model modifier alphas. Reads
## modifier_comparison.csv (output of 60_compare_modifiers_across_models.R)
## or rebuilds it from scratch by walking the modifier directories. Produces
## one row per modifier coefficient with columns model_lambda containing
## "mean [q5, q95]" for easy paste into deck/manuscript tables.
## =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(optparse)
  library(readr)
})

opts <- OptionParser(option_list = list(
  make_option("--base_dir", type = "character",
              default = "calibration/output/conus"),
  make_option("--out", type = "character",
              default = "calibration/output/conus/modifier_comparison_wide.csv")
)) |> parse_args()

key_alphas <- c("alpha_plant",
                "alpha_fire", "alpha_wind", "alpha_harvest",
                "alpha_insect", "alpha_disease",
                "alpha_cutting", "alpha_siteprep",
                "alpha_0", "sigma_L1", "sigma_resid")

candidate_dirs <- list.dirs(opts$base_dir, full.names = TRUE, recursive = TRUE)
modifier_dirs  <- candidate_dirs[grepl("/modifier($|/)", candidate_dirs)]

rows <- list()
for (d in modifier_dirs) {
  csvs <- list.files(d, pattern = "summary\\.csv$", full.names = TRUE)
  for (cs in csvs) {
    df <- tryCatch(read_csv(cs, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    setDT(df)
    if (!"variable" %in% names(df)) next

    model_name <- basename(dirname(d))
    if (model_name == "modifier") {
      model_name <- basename(dirname(dirname(d)))
    }
    fit_tag <- sub("_summary\\.csv$", "", basename(cs))
    lambda_match <- regmatches(fit_tag, regexpr("lambda[0-9]+", fit_tag))
    lambda <- ifelse(length(lambda_match) > 0, lambda_match, "lambda?")

    df[, model_label := paste0(model_name, "_", lambda)]
    keep <- df[variable %in% key_alphas, .(model_label, variable, mean, q5, q95, rhat, ess_bulk)]
    if (nrow(keep) > 0) rows[[length(rows) + 1L]] <- keep
  }
}

if (length(rows) == 0) stop("No modifier summaries found")
all <- rbindlist(rows, fill = TRUE)

all[, formatted := sprintf("%+0.3f [%+0.3f,%+0.3f]", mean, q5, q95)]
all[, rhat_flag := ifelse(is.na(rhat) | rhat < 1.05, "", " *")]
all[, formatted := paste0(formatted, rhat_flag)]

# Pivot to wide: coef on rows, model_lambda on cols
wide <- dcast(all, variable ~ model_label, value.var = "formatted")
setcolorder(wide, c("variable", sort(setdiff(names(wide), "variable"))))

# Order rows by intuitive importance
ord <- c("alpha_plant", "alpha_cutting", "alpha_siteprep",
         "alpha_insect", "alpha_disease", "alpha_fire",
         "alpha_wind", "alpha_harvest",
         "alpha_0", "sigma_L1", "sigma_resid")
wide <- wide[match(ord, variable, nomatch = 0)]

write_csv(wide, opts$out)
cat(sprintf("Wrote pivoted table to %s\n", opts$out))
cat("Columns:", paste(names(wide), collapse = "  "), "\n\n")

print(wide, nrows = 30)
cat("\n* indicates rhat >= 1.05 (potential convergence concern)\n")
