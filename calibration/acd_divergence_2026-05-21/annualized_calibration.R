## =============================================================================
## annualized_calibration.R
##
## Convert periodic calibration into ANNUAL species multipliers for the annual
## stepping Acadian model (AcadianGY.R / AcadianGYOneStand, cyclen = 1).
##
## Motivation
##   FVS calibration (and the FIA remeasurement it is fit against) is periodic:
##   a single factor scales a 5 or 10 year basal area increment. The Acadian R
##   model steps one year at a time, and it scales a LINEAR diameter increment
##   (dDBH.mult), not a basal area increment. Two mismatches follow:
##     1. scale     : basal area is proportional to diameter squared, so a basal
##                    area factor m is not a diameter factor m.
##     2. time base : a periodic factor applied once is not an annual factor
##                    applied L times, because growth compounds and competition
##                    evolves year to year.
##   This file provides both a closed form conversion (Tier 1) and a data driven
##   optimizer (Tier 2) that produce annual diameter and height multipliers the
##   model can use directly.
##
## Design
##   Pure base R so it runs on Cardinal without package loading. A self test at
##   the bottom recovers known multipliers from simulated data. Source the file
##   and the self test runs; set RUN_SELFTEST <- FALSE to suppress.
##
## Author: prepared for A. Weiskittel, 2026-05-21
## =============================================================================

## ---------------------------------------------------------------------------
## 0. constants and small helpers
## ---------------------------------------------------------------------------

## basal area of a stem, m2, for dbh in cm. Multiply by EXPF for per hectare.
K_BA <- 0.00007854

ba_stem <- function(dbh_cm) K_BA * dbh_cm^2

## per record basal area increment (m2 per ha if expf is trees per ha)
ba_increment <- function(dbh, ddbh, expf = 1) {
  expf * K_BA * ((dbh + ddbh)^2 - dbh^2)
}

## ---------------------------------------------------------------------------
## 1. TIER 1: closed form basal area factor to diameter factor
##
## For a single tree in a single year, find the diameter increment factor c so
## that the resulting basal area increment equals m times the unscaled basal
## area increment:
##     (dbh + c*ddbh)^2 - dbh^2 = m * [ (dbh+ddbh)^2 - dbh^2 ]
## This is a quadratic in c with the positive root taken.
## ---------------------------------------------------------------------------

dia_factor_from_ba_factor <- function(dbh, ddbh, m) {
  ## guard against zero growth
  ok <- ddbh > 0 & is.finite(ddbh)
  c_out <- rep(NA_real_, length(dbh))
  target <- m * ((dbh + ddbh)^2 - dbh^2)   # desired squared diameter increment
  ## c^2 ddbh^2 + 2 dbh ddbh c - target = 0
  a <- ddbh[ok]^2
  b <- 2 * dbh[ok] * ddbh[ok]
  cc <- -target[ok]
  disc <- b^2 - 4 * a * cc
  disc[disc < 0] <- 0
  c_out[ok] <- (-b + sqrt(disc)) / (2 * a)
  ## where there was no growth, the factor is undefined; return m as a neutral
  ## fallback so downstream code does not break
  c_out[!ok] <- m
  c_out
}

## ---------------------------------------------------------------------------
## 2. TIER 1 aggregate: one annual diameter factor that reproduces a target
##    PERIODIC basal area factor over L years for a whole tree list.
##
##    Because the annual increment is recomputed as the tree grows, the L year
##    behaviour is found by iteration, and the single annual factor c is solved
##    by one dimensional root finding so that
##        total_ba_increment(c) = M * total_ba_increment(c = 1)
##
##    ddbh_fn(dbh, trees) must return the UNSCALED annual diameter increment for
##    the current diameters. Supply your own (a thin wrapper over dDBH_fun from
##    AcadianGY.R is ideal). If omitted, a constant annual increment is assumed.
## ---------------------------------------------------------------------------

.total_ba_increment <- function(trees, c, L, ddbh_fn) {
  d <- trees$dbh
  ba0 <- sum(trees$expf * ba_stem(d))
  for (yr in seq_len(L)) {
    inc <- ddbh_fn(d, trees) * c
    inc[!is.finite(inc) | inc < 0] <- 0
    d <- d + inc
  }
  sum(trees$expf * ba_stem(d)) - ba0
}

annual_dia_factor_for_periodic_ba <- function(trees, M, L,
                                              ddbh_fn = NULL,
                                              interval = c(0.05, 12)) {
  if (is.null(ddbh_fn)) {
    ## degenerate constant increment assumption
    stopifnot(!is.null(trees$ddbh_annual))
    ddbh_fn <- function(d, tr) tr$ddbh_annual
  }
  base <- .total_ba_increment(trees, 1, L, ddbh_fn)
  if (base <= 0) return(M)            # nothing to scale, return neutral
  g <- function(c) .total_ba_increment(trees, c, L, ddbh_fn) - M * base
  ## bracket and solve
  lo <- interval[1]; hi <- interval[2]
  if (sign(g(lo)) == sign(g(hi))) {
    ## expand once if not bracketed
    hi <- hi * 4
    if (sign(g(lo)) == sign(g(hi))) return(M)  # fall back to naive factor
  }
  uniroot(g, lower = lo, upper = hi, tol = 1e-6)$root
}

## ---------------------------------------------------------------------------
## 3. TIER 2: fit annual multipliers against observed periodic remeasurement
##
##    This is the rigorous answer. Given remeasured plots (initial state, the
##    remeasurement interval L, and the observed end state), find the per species
##    annual diameter multiplier that, when the annual model is iterated L years,
##    best matches the observed periodic change.
##
##    project_fn is the pluggable annual stepper. Signature:
##        project_fn(state, c_by_spp, L) -> data.frame with at least $dbh
##    where state is a data.frame of one plot (or all plots) at year 0 with at
##    least columns dbh, sp, expf, plot, and c_by_spp is a named numeric vector
##    of annual diameter multipliers keyed by species code. The function must
##    apply c_by_spp[sp] to that species' annual diameter increment each year and
##    recompute competition between years. A thin wrapper over AcadianGYOneStand
##    does exactly this; a toy stepper is provided for the self test.
##
##    The objective is summed squared error in plot basal area per hectare at the
##    horizon, which is the quantity calibration usually targets. Switch to
##    diameter or volume by editing observed_metric / predicted_metric.
## ---------------------------------------------------------------------------

## per species basal area at the horizon, one value per species present
ba_by_species <- function(df) {
  tapply(df$expf * ba_stem(df$dbh), df$sp, sum)
}

## squared error between predicted and observed per species basal area for one
## plot, aligned over the union of species and treating absent species as zero.
## Calibrating to per species basal area (not total plot basal area) is both more
## correct, because a species multiplier should answer to that species' own
## growth, and far better identified, because it stops similar species from
## trading off against each other.
.sse_plot <- function(pred_df, obs_vec) {
  pred <- ba_by_species(pred_df)
  sps  <- union(names(pred), names(obs_vec))
  p <- setNames(rep(0, length(sps)), sps); p[names(pred)] <- pred
  o <- setNames(rep(0, length(sps)), sps); o[names(obs_vec)] <- obs_vec
  sum((p - o)^2)
}

## objective for a single species multiplier, holding the others fixed
.objective_one <- function(c_try, sp, c_fixed, plots, observed_sp, project_fn) {
  c_vec <- c_fixed
  c_vec[sp] <- c_try
  sse <- 0
  for (pid in names(plots)) {
    end <- project_fn(plots[[pid]], c_vec, attr(plots[[pid]], "L"))
    sse <- sse + .sse_plot(end, observed_sp[[pid]])
  }
  sse
}

fit_annual_multiplier_one_species <- function(sp, c_fixed, plots, observed_sp,
                                              project_fn, bounds = c(0.3, 3.0)) {
  opt <- optimize(.objective_one, interval = bounds,
                  sp = sp, c_fixed = c_fixed, plots = plots,
                  observed_sp = observed_sp, project_fn = project_fn)
  opt$minimum
}

## coordinate descent across species, because competition couples them
fit_annual_multipliers <- function(plots, observed_sp, project_fn,
                                   species = NULL, n_rounds = 6,
                                   bounds = c(0.3, 3.0), verbose = TRUE) {
  if (is.null(species)) {
    species <- sort(unique(unlist(lapply(plots, function(p) unique(p$sp)))))
  }
  c_fixed <- setNames(rep(1, length(species)), species)
  for (r in seq_len(n_rounds)) {
    for (sp in species) {
      c_fixed[sp] <- fit_annual_multiplier_one_species(
        sp, c_fixed, plots, observed_sp, project_fn, bounds)
    }
    if (verbose) {
      tot_sse <- sum(sapply(names(plots), function(pid) {
        end <- project_fn(plots[[pid]], c_fixed, attr(plots[[pid]], "L"))
        .sse_plot(end, observed_sp[[pid]])
      }))
      cat(sprintf("round %d  total SSE = %.5f\n", r, tot_sse))
      cat("   factors: ",
          paste(sprintf("%s=%.3f", names(c_fixed), c_fixed), collapse = "  "),
          "\n")
    }
  }
  c_fixed
}

## convenience: turn fitted ANNUAL diameter factors into the table the customRun
## bridge expects (one row per species). Keep height and mortality as separate
## columns; fit them with the same machinery on their own metrics.
build_annual_calib_table <- function(c_dia,
                                     c_ht = NULL,
                                     c_mort = NULL) {
  sp <- names(c_dia)
  data.frame(
    SP        = sp,
    dDBH.mult = as.numeric(c_dia[sp]),
    dHt.mult  = if (is.null(c_ht))   1 else as.numeric(c_ht[sp]),
    mort.mult = if (is.null(c_mort)) 1 else as.numeric(c_mort[sp]),
    stringsAsFactors = FALSE
  )
}

## =============================================================================
## SELF TEST
##   Build synthetic remeasured plots from a known toy annual diameter model with
##   known per species annual multipliers, then recover those multipliers. This
##   verifies the optimizer and the basal area conversion without needing FVS,
##   the DLL, or the real FIA pull.
## =============================================================================

if (!exists("RUN_SELFTEST")) RUN_SELFTEST <- TRUE

if (RUN_SELFTEST) {

  set.seed(42)

  ## --- toy annual diameter increment: depends on dbh and on stand competition
  ## (basal area in larger trees), so species are coupled, just like the real
  ## model. Units are cm per year. This stands in for dDBH_fun.
  toy_ddbh <- function(dbh, trees) {
    o <- order(-dbh)
    ba_cum <- cumsum(trees$expf[o] * ba_stem(dbh[o]))
    bal <- numeric(length(dbh))
    bal[o] <- ba_cum - trees$expf[o] * ba_stem(dbh[o])
    base_rate <- c(SW = 0.45, HW = 0.30)[trees$grp]
    inc <- base_rate * exp(-0.012 * bal) * (1 + 0.02 * pmin(dbh, 40))
    inc[!is.finite(inc)] <- 0
    inc
  }

  ## toy annual stepper that honours per species multipliers and recomputes
  ## competition each year, matching the project_fn contract
  toy_project <- function(state, c_by_spp, L) {
    d <- state$dbh
    for (yr in seq_len(L)) {
      tr <- state; tr$dbh <- d
      inc <- toy_ddbh(d, tr) * c_by_spp[state$sp]
      inc[!is.finite(inc) | inc < 0] <- 0
      d <- d + inc
    }
    out <- state; out$dbh <- d; out
  }

  ## build several plots
  make_plot <- function(pid, n = 40, L = 5) {
    sp  <- sample(c("RS", "BF", "RM", "YB"), n, replace = TRUE)
    grp <- ifelse(sp %in% c("RS", "BF"), "SW", "HW")
    df <- data.frame(plot = pid, sp = sp, grp = grp,
                     dbh = pmax(2, rgamma(n, 4, 0.25)),
                     expf = runif(n, 8, 25),
                     stringsAsFactors = FALSE)
    attr(df, "L") <- L
    df
  }

  plots0 <- lapply(sprintf("P%02d", 1:12), make_plot)
  names(plots0) <- sprintf("P%02d", 1:12)

  ## the TRUTH we will try to recover
  c_true <- c(RS = 1.20, BF = 0.85, RM = 1.05, YB = 0.70)

  ## generate observed per species end state basal area using the true factors
  observed_sp <- lapply(plots0, function(p)
    ba_by_species(toy_project(p, c_true, attr(p, "L"))))

  cat("\n=== Tier 2 optimizer: recover known annual multipliers ===\n")
  c_hat <- fit_annual_multipliers(plots0, observed_sp, toy_project,
                                  n_rounds = 6)

  cat("\nrecovery check (true vs fitted):\n")
  comp <- data.frame(species = names(c_true),
                     true = as.numeric(c_true),
                     fitted = as.numeric(c_hat[names(c_true)]))
  comp$abs_err <- abs(comp$true - comp$fitted)
  print(comp, row.names = FALSE)
  cat(sprintf("max absolute error = %.4f  (target < 0.05)\n", max(comp$abs_err)))
  stopifnot(max(comp$abs_err) < 0.05)
  cat("PASS: optimizer recovered the known multipliers.\n")

  cat("\n=== Tier 1 check: basal area factor to diameter factor ===\n")
  ## a tree growing 1 cm per year, dbh 20, target basal area factor 1.5
  d0 <- 20; dd <- 1; m <- 1.5
  cfac <- dia_factor_from_ba_factor(d0, dd, m)
  achieved <- ((d0 + cfac * dd)^2 - d0^2) / ((d0 + dd)^2 - d0^2)
  cat(sprintf("dbh=%.0f ddbh=%.1f target BA factor=%.2f -> dia factor=%.4f, achieved BA factor=%.4f\n",
              d0, dd, m, cfac, achieved))
  stopifnot(abs(achieved - m) < 1e-6)
  cat(sprintf("note: diameter factor %.4f differs from naive %.2f by %.1f%%\n",
              cfac, m, 100 * (m - cfac) / m))
  cat("PASS: closed form conversion is exact.\n")

  cat("\n=== Tier 1 aggregate: annual factor for a periodic BA target ===\n")
  trees <- plots0[["P01"]]
  ddbh_fn <- function(d, tr) toy_ddbh(d, tr)
  c_annual <- annual_dia_factor_for_periodic_ba(trees, M = 1.3, L = 5,
                                                ddbh_fn = ddbh_fn)
  ach <- .total_ba_increment(trees, c_annual, 5, ddbh_fn) /
         .total_ba_increment(trees, 1, 5, ddbh_fn)
  cat(sprintf("target periodic BA factor 1.30 over 5 yr -> annual dia factor %.4f (achieved %.4f)\n",
              c_annual, ach))
  stopifnot(abs(ach - 1.30) < 1e-3)
  cat("PASS: annual factor reproduces the periodic basal area target.\n")

  cat("\nALL SELF TESTS PASSED.\n")
}

## =============================================================================
## HOW TO USE WITH THE REAL MODEL (sketch, not run here)
##
##   source("AcadianGY.R")
##   # 1. wrap AcadianGYOneStand as a project_fn that applies c_by_spp to dDBH:
##   real_project <- function(state, c_by_spp, L) {
##     tree <- state
##     # attach a per record annual diameter multiplier from c_by_spp
##     for (yr in seq_len(L)) {
##       tree$YEAR <- yr
##       # inject the fitted factor by setting dDBH.mult before the call, or by
##       # multiplying dDBH inside a lightly patched AcadianGYOneStand
##       out <- AcadianGYOneStand(tree, stand = list(CSI = 12), ops = ops_annual)
##       tree <- out
##     }
##     tree
##   }
##   # 2. assemble remeasured FIA plots into the `plots` list (year 0 state) and
##   #    observed_ba (measured basal area per ha at remeasurement), then:
##   c_dia  <- fit_annual_multipliers(plots, observed_ba, real_project)
##   tab    <- build_annual_calib_table(c_dia)
##   write.csv(tab, "acd_annual_calibration.csv", row.names = FALSE)
##
##   The resulting dDBH.mult column is the ANNUAL diameter multiplier to feed the
##   bridge in place of the mis scaled FVS baimult. Repeat with a height metric
##   for dHt.mult and a survival metric for mort.mult.
## =============================================================================
