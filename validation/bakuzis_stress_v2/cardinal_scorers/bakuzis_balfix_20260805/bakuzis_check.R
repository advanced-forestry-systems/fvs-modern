#!/usr/bin/env Rscript
# bakuzis_check.R  --  automated Bakuzis-matrix consistency assessment
# Usage: Rscript bakuzis_check.R trajectories.csv [form_factor=0.40]
# Input CSV (long, one row per group x age) needs columns:
#   site (productivity class, numeric or ordered), age, HT, TPH, BAPH, QMD
#   optional: origin (a second grouping), VOL (else VOL = BAPH*HT*form_factor)
# Output: a PASS/FLAG report and a recommendations block to stdout.
# Base R only. Falsifies structurally flawed projections; does not confirm accuracy.

args <- commandArgs(trailingOnly = TRUE)
csv  <- if (length(args) >= 1) args[1] else stop("supply a trajectories CSV")
ff   <- if (length(args) >= 2) as.numeric(args[2]) else 0.40
d <- read.csv(csv, stringsAsFactors = FALSE)
if (!"origin" %in% names(d)) d$origin <- "all"
if (!"VOL"    %in% names(d)) d$VOL <- d$BAPH * d$HT * ff
d <- d[order(d$origin, d$site, d$age), ]
sites   <- sort(unique(d$site))
origins <- unique(d$origin)
flags <- character(0)
mono  <- function(x, inc = TRUE, tol = 1e-6) if (inc) all(diff(x) >= -tol) else all(diff(x) <= tol)
cat(strrep("=", 64), "\nBAKUZIS ASSESSMENT:", csv, "\n", strrep("=", 64), "\n", sep = "")

## 1 + 3 monotonicity (HT up, TPH down) -------------------------------------
cat("\n[1/3 Monotonicity] HT increasing, TPH decreasing with age:\n")
for (o in origins) for (s in sites) {
  g <- d[d$origin == o & d$site == s, ]
  okh <- mono(g$HT, TRUE); okt <- mono(g$TPH, FALSE)
  cat(sprintf("  %s site %s: HT_up=%s TPH_down=%s\n", o, s, okh, okt))
  if (!okh) flags <- c(flags, sprintf("HT not monotone up (%s site %s)", o, s))
  if (!okt) flags <- c(flags, sprintf("TPH not monotone down (%s site %s)", o, s))
}

## 1 + 7 site ordering (better site not overtaken in HT or VOL) -------------
cat("\n[1/7 Site ordering] HT and VOL ordered by site at each age:\n")
ageset <- sort(unique(d$age))
for (o in origins) {
  bh <- bv <- 0
  for (a in ageset) {
    h <- sapply(sites, function(s) { v <- d$HT [d$origin==o & d$site==s & d$age==a]; if(length(v)) v[1] else NA })
    v <- sapply(sites, function(s) { z <- d$VOL[d$origin==o & d$site==s & d$age==a]; if(length(z)) z[1] else NA })
    if (!any(is.na(h)) && is.unsorted(h)) bh <- bh + 1
    if (!any(is.na(v)) && is.unsorted(v)) bv <- bv + 1
  }
  cat(sprintf("  %s: HT inversions=%d/%d  VOL inversions=%d/%d\n", o, bh, length(ageset), bv, length(ageset)))
  if (bh) flags <- c(flags, sprintf("HT site-ordering inverted %dx (%s)", bh, o))
  if (bv) flags <- c(flags, sprintf("VOL site-ordering inverted %dx (%s) [check self-thinning/SDImax]", bv, o))
}

## 4 Reineke slope ----------------------------------------------------------
cat("\n[4 Reineke] log(TPH)~log(QMD) self-thinning slope (target ~ -1.605):\n")
for (o in origins) for (s in sites) {
  g <- d[d$origin==o & d$site==s, ]; g <- g[order(g$age), ]
  # self-thinning phase = points where QMD is genuinely increasing (exclude capped QMD)
  keep <- c(TRUE, diff(g$QMD) > 0.5)
  g <- g[keep, ]
  if (nrow(g) >= 3 && sd(log(g$QMD)) > 1e-3) {
    sl <- coef(lm(log(TPH) ~ log(QMD), g))[2]
    ok <- is.finite(sl) && sl > -2.2 && sl < -1.2
    cat(sprintf("  %s site %s: slope=%.2f %s\n", o, s, sl, ifelse(ok,"OK","FLAG")))
    if (!ok) flags <- c(flags, sprintf("Reineke slope %.2f off (%s site %s)", sl, o, s))
  } else {
    cat(sprintf("  %s site %s: not in a self-thinning regime (QMD capped/flat) - Reineke n/a\n", o, s))
  }
}

## 8 Eichhorn: volume-height site independence ------------------------------
cat("\n[8 Eichhorn] volume-height site independence (CV of VOL across sites at matched HT):\n")
for (o in origins) {
  hi <- max(sapply(sites, function(s) min(d$HT[d$origin==o & d$site==s])))
  hx <- min(sapply(sites, function(s) max(d$HT[d$origin==o & d$site==s])))
  if (hx > hi) {
    grid <- seq(hi, hx, length.out = 8)
    cv <- sapply(grid, function(hg) {
      vv <- sapply(sites, function(s) { g <- d[d$origin==o & d$site==s, ]; approx(g$HT, g$VOL, hg, rule=2)$y })
      sd(vv) / mean(vv) })
    m <- mean(cv); ok <- m < 0.15
    cat(sprintf("  %s: mean CV=%.1f%% %s\n", o, 100*m, ifelse(ok,"OK (site-independent)","FLAG (site-dependent V-H)")))
    if (!ok) flags <- c(flags, sprintf("Eichhorn: V-H site-dependent CV %.0f%% (%s)", 100*m, o))
  }
}

## summary ------------------------------------------------------------------
cat("\n", strrep("=", 64), "\nFLAGS: ", length(flags), "\n", sep = "")
if (length(flags)) for (f in flags) cat("  - ", f, "\n", sep = "") else cat("  none - projection passes all tested relations (falsification failed; not a confirmation of accuracy)\n")
cat(strrep("=", 64), "\nInterpret with references/bakuzis-method.md: map each flag to the responsible submodel and a refinement.\n", sep="")
