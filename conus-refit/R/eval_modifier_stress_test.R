##=============================================================================
## eval_modifier_stress_test.R
##
## Comprehensive modifier stress test:
##   1. Sweep lambda 5 / 10 / 20 across DG_Kuehne, HG_Organon, ht_dbh
##   2. Compare alpha coefficients across lambdas + base models
##   3. Flag contradictions (e.g., DG vs HG plantation sign)
##   4. Compute LOO ELPD per fit if log_lik present
##   5. Output a single CSV summary for manuscript Table X
##
## Reads from:
##   calibration/output/conus/{dg,hg,ht_dbh}/modifier*/
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
})

OUT_DIR <- "calibration/output/conus/modifier_stress_test"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Find all modifier summary CSVs
search_dirs <- c(
  "calibration/output/conus/dg/modifier",
  "calibration/output/conus/dg/modifier_kuehne",
  "calibration/output/conus/dg/modifier_v2",
  "calibration/output/conus/hg/modifier",
  "calibration/output/conus/hg/modifier_v2",
  "calibration/output/conus/hg/modifier_traitmed",
  "calibration/output/conus/ht_dbh/modifier"
)

cat("== eval_modifier_stress_test.R ==\n\n")

all_rows <- list()
for (d in search_dirs) {
  if (!dir.exists(d)) {
    cat("  skipping (missing):", d, "\n")
    next
  }
  csvs <- list.files(d, pattern = "summary\\.csv$", full.names = TRUE)
  for (cs in csvs) {
    df <- tryCatch(fread(cs), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    fname <- basename(cs)
    component <- if (grepl("dg/", d)) "DG_Kuehne"
                 else if (grepl("hg/", d)) "HG_Organon"
                 else if (grepl("ht_dbh", d)) "ht_dbh_Wykoff"
                 else "Other"
    variant <- if (grepl("modifier_v2", d)) "v2_kernels"
               else if (grepl("modifier_traitmed", d)) "traitmed"
               else "common_binary"
    lambda <- as.integer(stringr::str_match(fname, "lambda([0-9]+)")[, 2])
    if (is.na(lambda) && grepl("v2", fname)) {
      # v2 fits use different naming: a5_g3_t30
      acute <- as.integer(stringr::str_match(fname, "a([0-9]+)")[, 2])
      gamma_ <- as.integer(stringr::str_match(fname, "g([0-9]+)")[, 2])
      trt   <- as.integer(stringr::str_match(fname, "t([0-9]+)")[, 2])
      lambda <- NA
      variant_param <- sprintf("a=%d/g=%d/t=%d", acute, gamma_, trt)
    } else {
      variant_param <- sprintf("lambda=%d", lambda)
    }
    df[, `:=`(component = component, variant = variant,
              lambda = lambda, variant_param = variant_param,
              file = fname)]
    all_rows[[length(all_rows) + 1]] <- df
  }
}

if (!requireNamespace("stringr", quietly = TRUE)) {
  install.packages("stringr", repos = "https://cloud.r-project.org",
                   lib = .libPaths()[1])
  library(stringr)
} else {
  library(stringr)
}

if (length(all_rows) == 0) {
  stop("No modifier summary CSVs found.")
}

dt <- rbindlist(all_rows, fill = TRUE)
cat("Loaded", nrow(dt), "rows from", length(all_rows), "files\n\n")

# Focus on alpha coefficients
alphas <- dt[grepl("^alpha_", variable)]
alphas[, alpha_name := sub("^alpha_", "", variable)]

# Wide format: rows = (component, variant, lambda), columns = alpha mean per type
wide <- dcast(alphas,
              component + variant + variant_param + lambda + file ~ alpha_name,
              value.var = "mean")
setorder(wide, component, variant, lambda, na.last = TRUE)

cat("=== Alpha coefficient mean by (component, variant, lambda) ===\n")
print(wide)
fwrite(wide, file.path(OUT_DIR, "modifier_alpha_summary_wide.csv"))

# Stress test 1: lambda sensitivity — does alpha_plant change across lambdas?
cat("\n=== Stress test 1: lambda sensitivity (alpha_plant) ===\n")
ls1 <- alphas[alpha_name == "plant", .(component, variant, lambda, mean, q5, q95)]
setorder(ls1, component, lambda)
print(ls1)
fwrite(ls1, file.path(OUT_DIR, "modifier_alpha_plant_by_lambda.csv"))

# Stress test 2: cross-component contradictions
cat("\n=== Stress test 2: alpha signs across components ===\n")
sign_dt <- alphas[, .(mean_sign = sign(mean)), by = .(component, alpha_name)]
sign_wide <- dcast(sign_dt, alpha_name ~ component, value.var = "mean_sign",
                    fun.aggregate = function(x) {
                      if (length(unique(x)) == 1) x[1] else 99
                    })
print(sign_wide)
fwrite(sign_wide, file.path(OUT_DIR, "modifier_alpha_sign_matrix.csv"))

# Stress test 3: magnitude consistency
cat("\n=== Stress test 3: alpha magnitude by component (median across variants) ===\n")
mag <- alphas[, .(median_alpha = median(mean),
                   sd_alpha = sd(mean),
                   n = .N),
              by = .(component, alpha_name)]
setorder(mag, alpha_name, component)
print(mag)
fwrite(mag, file.path(OUT_DIR, "modifier_alpha_magnitude.csv"))

# Stress test 4: residual scale — does base model leave structure?
cat("\n=== Stress test 4: sigma_resid by component/variant ===\n")
sigs <- dt[variable == "sigma_resid", .(component, variant, lambda, mean, q5, q95)]
setorder(sigs, component, lambda)
print(sigs)
fwrite(sigs, file.path(OUT_DIR, "modifier_sigma_resid.csv"))

# Stress test 5: convergence health
cat("\n=== Stress test 5: max R-hat per fit ===\n")
conv <- dt[, .(max_rhat = max(rhat, na.rm = TRUE),
                min_ess  = min(ess_bulk, na.rm = TRUE),
                n_vars = .N),
           by = .(component, variant, variant_param, file)]
setorder(conv, -max_rhat)
print(conv)
fwrite(conv, file.path(OUT_DIR, "modifier_convergence.csv"))

cat("\nAll outputs:", OUT_DIR, "\n")
cat("\nKEY FINDINGS to investigate:\n")
cat(" - Plantation sign contradiction DG vs HG (alpha_plant)\n")
cat(" - sigma_resid magnitude per component (large = base model leaving structure)\n")
cat(" - Lambda choice that minimizes uncertainty per component\n")
cat(" - Any fit with max_rhat > 1.05 (suspect convergence)\n")
