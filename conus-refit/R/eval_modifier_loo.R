##=============================================================================
## eval_modifier_loo.R
## LOO ELPD comparison: modifier_common vs modifier_traitmed.
##=============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(loo)
})

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  m <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(m) == 0) return(default)
  sub(paste0("^--", name, "="), "", m[1])
}

COMMON_FIT   <- get_arg("common")
TRAITMED_FIT <- get_arg("traitmed")
OUT_FILE     <- get_arg("out", "modifier_loo_compare.rds")

stopifnot(!is.null(COMMON_FIT), file.exists(COMMON_FIT))
stopifnot(!is.null(TRAITMED_FIT), file.exists(TRAITMED_FIT))

cat("== eval_modifier_loo.R ==\n")
cat("  common  :", COMMON_FIT, "\n")
cat("  traitmed:", TRAITMED_FIT, "\n")
cat("  out     :", OUT_FILE, "\n\n")

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

loo_one <- function(label, path) {
  cat("Loading", label, "...\n"); flush.console()
  fit <- readRDS(path)
  ll <- fit$draws("log_lik", format = "draws_matrix")
  cat("  log_lik dims:", paste(dim(ll), collapse=" x "), "\n")
  l <- loo::loo(ll)
  cat(sprintf("  %s ELPD = %.1f (SE %.1f)\n",
              label,
              l$estimates["elpd_loo","Estimate"],
              l$estimates["elpd_loo","SE"]))
  rm(fit, ll); gc()
  l
}

l_common   <- loo_one("modifier_common",   COMMON_FIT)
l_traitmed <- loo_one("modifier_traitmed", TRAITMED_FIT)

cat("\n=== loo_compare ===\n")
cmp <- loo::loo_compare(list(common = l_common, traitmed = l_traitmed))
print(cmp)

elpd_diff <- abs(cmp[2, "elpd_diff"])
se_diff   <- cmp[2, "se_diff"]
n_se      <- elpd_diff / max(se_diff, 0.01)
winner    <- rownames(cmp)[1]
cat(sprintf("\nELPD difference: %.1f (SE %.1f), %.2f SE units\n",
            elpd_diff, se_diff, n_se))
if (n_se > 4) {
  cat(sprintf("DECISION: %s wins by > 4 SE (very strong evidence).\n", winner))
} else if (n_se > 2) {
  cat(sprintf("DECISION: %s wins by > 2 SE (strong evidence).\n", winner))
} else if (n_se > 1) {
  cat(sprintf("DECISION: %s wins by > 1 SE (moderate evidence).\n", winner))
} else {
  cat("DECISION: models tied within 1 SE.\n")
}

saveRDS(list(common = l_common, traitmed = l_traitmed, compare = cmp), OUT_FILE)
cat("\nSaved:", OUT_FILE, "\n")
