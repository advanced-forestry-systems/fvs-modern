##=============================================================================
## eval_loo_v6_v7_v8.R
##
## LOO ELPD comparison across DG_Kuehne architecture variants v6, v7, v8.
## Reads each fit's log_lik from the saved fit object and produces
## loo_compare() with se_diff.
##
## CLI:
##   --v6_fit=PATH (path to v6 fit RDS)
##   --v7_fit=PATH (path to v7 fit RDS)
##   --v8_fit=PATH (path to v8 fit RDS)
##   --outdir=PATH (where to save comparison)
##
## Author: A. Weiskittel + Claude
## Date: 2026-05-15
##=============================================================================

library(data.table)
library(loo)

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

V6_FIT <- get_arg("v6_fit", "calibration/output/conus/dg_kue/v6/dg_kuehne_v6_100k_prod_fit.rds")
V7_FIT <- get_arg("v7_fit", "calibration/output/conus/dg_kue/v7/dg_kuehne_v7_100k_prod_fit.rds")
V8_FIT <- get_arg("v8_fit", "calibration/output/conus/dg_kue/v8/dg_kuehne_v8_100k_prod_fit.rds")
OUT_DIR <- get_arg("outdir", "calibration/output/conus/dg_kue/loo_compare")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("== eval_loo_v6_v7_v8.R ==\n")
cat("v6_fit:", V6_FIT, "(exists:", file.exists(V6_FIT), ")\n")
cat("v7_fit:", V7_FIT, "(exists:", file.exists(V7_FIT), ")\n")
cat("v8_fit:", V8_FIT, "(exists:", file.exists(V8_FIT), ")\n\n")

loo_one <- function(label, path) {
  if (!file.exists(path)) {
    cat("  ", label, ": file missing, skipping\n", sep = "")
    return(NULL)
  }
  cat("Loading ", label, " ...", sep = ""); flush.console()
  fit <- readRDS(path)
  ll <- fit$draws("log_lik", format = "draws_matrix")
  cat(" log_lik dims:", paste(dim(ll), collapse=" x "), "\n")
  l <- loo::loo(ll)
  cat("  ", label, " ELPD =", round(l$estimates["elpd_loo", "Estimate"], 1),
      " (SE", round(l$estimates["elpd_loo", "SE"], 1), ")\n", sep="")
  rm(fit, ll); gc()
  l
}

loos <- list(
  v6 = loo_one("v6", V6_FIT),
  v7 = loo_one("v7", V7_FIT),
  v8 = loo_one("v8", V8_FIT)
)
loos <- loos[!sapply(loos, is.null)]

if (length(loos) < 2) {
  cat("\nFewer than 2 fits available; cannot compare. Exiting.\n")
  quit(status = 0)
}

cat("\n=== loo_compare ===\n")
cmp <- loo::loo_compare(loos)
print(cmp)

saveRDS(list(loos = loos, comparison = cmp),
        file.path(OUT_DIR, "loo_compare_v6_v7_v8.rds"))
data.table::fwrite(as.data.frame(cmp),
                    file.path(OUT_DIR, "loo_compare_v6_v7_v8.csv"),
                    row.names = TRUE)
cat("\nSaved to:", OUT_DIR, "\n")
