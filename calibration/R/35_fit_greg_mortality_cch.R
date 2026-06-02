#!/usr/bin/env Rscript
# =============================================================================
# fit_greg_mortality.R
# Greg Johnson's CONUS per-species mortality (Gompit on crown ratio + crown
# closure at tree tip), from "Mortality Equations for CONUS" (Johnson, Marshall,
# Weiskittel, 2026-05-26).
#
# Model (per species), annual survival hazard with an exposure offset so that
# variable FIA remeasurement interval lengths are integrated directly in the
# likelihood (the "integrative approach" in the paper):
#
#   eta  = b0 + b1 * (cr + 0.01)^b2 + b3 * cch^b4          (linear predictor)
#   H    = exp(eta)                                        (annual hazard)
#   P_surv(T) = exp(-H * T)                                (survive T years)
#   P_die(T)  = 1 - exp(-H * T)
#
# This reduces to the paper's Eq. (1) at T = 1:
#   P_surv = 1 - (1 - exp(-exp(eta)))  ... i.e. exp(-exp(eta)).
# cr  = crown ratio (fraction of total height)
# cch = crown closure at tree tip (fraction of acre)
#
# Fit per species by minimizing the negative log-likelihood with optim().
# b2 and b4 are exponents (nonlinear), so we fit on the natural scale with
# multiple starts and light bounds for stability.
#
# Input  : a CSV/RDS with columns  SPCD, cr, cch, alive (1/0), years
#          (alive = survived the interval; years = interval length)
# Output : <out>/greg_mortality_coefficients.csv  (SPCD, n, b0..b4, nll, conv)
#          <out>/greg_mortality_fit_summary.csv   (convergence + base-rate NLL)
#
# Usage:
#   Rscript fit_greg_mortality.R --data mort_cch.csv --out out/ --min-obs 5000
# =============================================================================

suppressWarnings(suppressMessages({
  ok <- requireNamespace("readr", quietly = TRUE)
}))

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) return(default); args[i + 1]
}
data_path <- getarg("--data")
out_dir   <- getarg("--out", "out")
min_obs   <- as.integer(getarg("--min-obs", "5000"))
stopifnot(!is.null(data_path))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_any <- function(p) {
  if (grepl("\\.rds$", p, ignore.case = TRUE)) return(as.data.frame(readRDS(p)))
  utils::read.csv(p, stringsAsFactors = FALSE)
}

# ---- negative log-likelihood for one species -------------------------------
# theta = c(b0, b1, b2, b3, b4). Guards keep the power bases valid and the
# hazard finite; returns a large penalty on non-finite evaluations so optim
# walks away from bad regions instead of erroring.
nll_species <- function(theta, cr, cch, alive, years) {
  b0 <- theta[1]; b1 <- theta[2]; b2 <- theta[3]; b3 <- theta[4]; b4 <- theta[5]
  cr_term  <- (cr + 0.01)^b2            # cr + 0.01 > 0 always
  cch_term <- ifelse(cch > 0, cch^b4, 0) # 0^b4 -> 0 for b4>0; guard cch==0
  eta <- b0 + b1 * cr_term + b3 * cch_term
  eta <- pmin(pmax(eta, -30), 30)        # clamp to keep exp() finite
  H   <- exp(eta)                        # annual hazard
  HT  <- H * years
  # log P(survive) = -H*T ; log P(die) = log(1 - exp(-H*T))
  log_surv <- -HT
  log_die  <- log1p(-exp(-HT))           # = log(1 - exp(-HT)), stable
  ll <- ifelse(alive == 1, log_surv, log_die)
  val <- -sum(ll)
  if (!is.finite(val)) return(1e12)
  val
}

# base-rate (intercept-only) NLL for the species, for the paper's comparison
nll_baserate <- function(cr, cch, alive, years) {
  f <- function(b0) {
    H <- exp(min(max(b0, -30), 30)); HT <- H * years
    -sum(ifelse(alive == 1, -HT, log1p(-exp(-HT))))
  }
  o <- optimize(f, c(-15, 5))
  o$objective
}

fit_one <- function(d) {
  cr <- d$cr; cch <- d$cch; alive <- d$alive; years <- d$years
  # multistart: vary b0 (hazard scale) and the exponents b2,b4
  starts <- list(
    c(-4, -0.5, -0.6,  0.0, 0.1),
    c(-3,  0.0, -0.8, -0.5, 0.3),
    c(-5, -1.0, -0.4,  0.2, 0.5),
    c(-4,  0.5, -0.6, -1.0, 1.0)
  )
  best <- NULL
  for (s in starts) {
    fit <- tryCatch(
      optim(s, nll_species, cr = cr, cch = cch, alive = alive, years = years,
            method = "Nelder-Mead",
            control = list(maxit = 2000, reltol = 1e-9)),
      error = function(e) NULL)
    if (!is.null(fit) && is.finite(fit$value)) {
      if (is.null(best) || fit$value < best$value) best <- fit
    }
  }
  best
}

# ---- main ------------------------------------------------------------------
dat <- read_any(data_path)
needed <- c("SPCD", "cr", "cch", "alive", "years")
miss <- setdiff(needed, names(dat))
if (length(miss)) stop("input missing columns: ", paste(miss, collapse = ", "))
dat <- dat[is.finite(dat$cr) & is.finite(dat$cch) &
           dat$cr > 0 & dat$cr <= 1 & dat$cch >= 0 &
           dat$years >= 1 & dat$years <= 20 &
           dat$alive %in% c(0, 1), , drop = FALSE]

sp_tab <- sort(table(dat$SPCD), decreasing = TRUE)
sp_keep <- as.integer(names(sp_tab[sp_tab >= min_obs]))
message(sprintf("species with >= %d obs: %d (of %d)",
                min_obs, length(sp_keep), length(sp_tab)))

rows <- list(); summ <- list()
for (sp in sp_keep) {
  d <- dat[dat$SPCD == sp, , drop = FALSE]
  fit <- fit_one(d)
  base_nll <- nll_baserate(d$cr, d$cch, d$alive, d$years)
  if (is.null(fit)) {
    message(sprintf("  SPCD %d (n=%d): FAILED", sp, nrow(d))); next
  }
  th <- fit$par
  rows[[as.character(sp)]] <- data.frame(
    SPCD = sp, n = nrow(d),
    b0 = th[1], b1 = th[2], b2 = th[3], b3 = th[4], b4 = th[5],
    nll = fit$value, nll_baserate = base_nll,
    improved = fit$value < base_nll, convergence = fit$convergence)
  message(sprintf("  SPCD %d (n=%d): nll=%.1f base=%.1f %s",
                  sp, nrow(d), fit$value, base_nll,
                  ifelse(fit$value < base_nll, "(improved)", "(NO improvement)")))
}

coef_df <- do.call(rbind, rows)
utils::write.csv(coef_df, file.path(out_dir, "greg_mortality_coefficients.csv"),
                 row.names = FALSE)
summary_df <- data.frame(
  n_species_fit = nrow(coef_df),
  n_improved = sum(coef_df$improved),
  median_b2 = median(coef_df$b2), median_b4 = median(coef_df$b4),
  total_nll = sum(coef_df$nll), total_nll_baserate = sum(coef_df$nll_baserate))
utils::write.csv(summary_df, file.path(out_dir, "greg_mortality_fit_summary.csv"),
                 row.names = FALSE)
message("wrote ", file.path(out_dir, "greg_mortality_coefficients.csv"))
print(summary_df)
