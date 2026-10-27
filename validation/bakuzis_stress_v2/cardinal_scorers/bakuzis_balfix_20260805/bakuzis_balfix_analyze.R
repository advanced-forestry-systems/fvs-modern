#!/usr/bin/env Rscript
# bakuzis_balfix_analyze.R
# Modified Bakuzis matrix (Leary 1997) on the FROZEN PTBALT-repaired FVS NE
# binary bin-balfix_20260804/FVSne (md5 b85e4c06995382425b01c756811d6d34),
# all FIVE model arms, plus the diagnostics addendum and a direct old-versus-new
# comparison against the pre-repair run in bakuzis_20260804/.
#
# Input : out/bakuzis_traj_raw.csv
# Output: traj_<arm>.csv, bakuzis_results.txt, bakuzis_matrix_table.csv,
#         bakuzis_addendum.txt, bakuzis_addendum_table.csv,
#         oldnew_reineke.csv, fig/*.png (300 dpi), error_log.txt

set.seed(20260805)
SEED <- 20260805

suppressPackageStartupMessages(library(ggplot2))
has_ggsci <- requireNamespace("ggsci", quietly = TRUE)

NEWDIR <- Sys.getenv("BAKDIR", "/fs/scratch/PUOM0008/crsfaaron/bakuzis_balfix_20260805")
OLDDIR <- Sys.getenv("BAKOLD", "/fs/scratch/PUOM0008/crsfaaron/bakuzis_20260804")
setwd(NEWDIR)
dir.create("fig", showWarnings = FALSE)

ERR <- file.path(NEWDIR, "error_log.txt")
cat(sprintf("bakuzis_balfix_analyze.R  seed=%d  started %s\n",
            SEED, format(Sys.time())), file = ERR)
logerr <- function(ctx, e) {
  cat(sprintf("[ERROR] %s :: %s\n", ctx, conditionMessage(e)),
      file = ERR, append = TRUE)
  invisible(NULL)
}
safe <- function(ctx, expr, fallback = NULL) {
  tryCatch(expr, error = function(e) { logerr(ctx, e); fallback })
}

# --- unit conversions (FVS imperial -> metric) -----------------------------
K_TPH  <- 2.4710538   # trees acre-1 -> trees ha-1
K_BAPH <- 0.2295684   # ft2 acre-1   -> m2 ha-1
K_QMD  <- 2.54        # in           -> cm
K_HT   <- 0.3048      # ft           -> m
K_VOL  <- 0.06997     # ft3 acre-1   -> m3 ha-1
PMSDIU <- 0.85; KSDI <- 0.5454154

ARMS <- c("fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate")
ARM_LAB <- c(fvs_base = "fvs_base", fvs_regional = "fvs_regional",
             organon = "organon (Greg)", conus_spdep = "conus_spdep (Greg)",
             conus_climate = "conus_climate (Greg)")
SITE_LAB <- c("S1 (SI 47)", "S2 (SI 61)", "S3 (SI 70)", "S4 (SI 83)")

prep <- function(path) {
  r <- read.csv(path, stringsAsFactors = FALSE)
  r <- r[is.finite(r$Tpa) & r$Tpa > 0 & is.finite(r$QMD) & r$QMD > 0, ]
  r$TPH  <- r$Tpa * K_TPH
  r$BAPH <- r$BA * K_BAPH
  r$QMDc <- r$QMD * K_QMD
  r$HTm  <- r$TopHt * K_HT
  r$VOLm <- r$TCuFt * K_VOL
  r$RD   <- ifelse(r$SDIMax > 0, r$SDI / r$SDIMax, NA)
  r$BAMAX_ft2 <- r$SDIMax * KSDI * PMSDIU     # vbase/sdical.f90 line 202
  r$BA_frac   <- r$BA / r$BAMAX_ft2
  r
}

raw <- prep("out/bakuzis_traj_raw.csv")
ARMS <- ARMS[ARMS %in% unique(raw$arm)]

# --- aggregate to a common stand-age grid ----------------------------------
AGE_GRID <- seq(70, 190, by = 10)
MIN_N    <- 20

interp_plot <- function(g, grid) {
  g <- g[order(g$Age), ]
  if (nrow(g) < 3) return(NULL)
  f <- function(y) approx(g$Age, y, grid, rule = 1)$y
  data.frame(age = grid, HT = f(g$HTm), TPH = f(g$TPH), BAPH = f(g$BAPH),
             QMD = f(g$QMDc), VOL = f(g$VOLm), SDI = f(g$SDI),
             SDIMAX = f(g$SDIMax), RD = f(g$RD))
}

aggregate_arm <- function(ra, a, write = TRUE) {
  parts <- list()
  for (p in unique(ra$PLOT)) {
    g <- ra[ra$PLOT == p, ]
    ip <- interp_plot(g, AGE_GRID)
    if (is.null(ip)) next
    ip$site <- g$site_cls[1]; ip$origin <- g$dens_cls[1]; ip$PLOT <- p
    parts[[length(parts) + 1]] <- ip
  }
  A <- do.call(rbind, parts)
  A <- A[!is.na(A$TPH), ]
  key <- paste(A$site, A$origin, A$age, sep = "|")
  num <- c("HT", "TPH", "BAPH", "QMD", "VOL", "SDI", "SDIMAX", "RD")
  M <- aggregate(A[num], by = list(key = key), FUN = function(x) mean(x, na.rm = TRUE))
  N <- aggregate(list(n = A$TPH), by = list(key = key), FUN = length)
  M <- merge(M, N, by = "key")
  sp <- do.call(rbind, strsplit(M$key, "|", fixed = TRUE))
  M$site <- as.integer(sp[, 1]); M$origin <- sp[, 2]; M$age <- as.numeric(sp[, 3])
  M$key <- NULL
  M <- M[M$n >= MIN_N, ]
  M <- M[order(M$origin, M$site, M$age), ]
  M$arm <- a
  if (write) write.csv(M[, c("site", "origin", "age", "HT", "TPH", "BAPH",
                             "QMD", "VOL", "SDI", "SDIMAX", "RD", "n")],
                       sprintf("traj_%s.csv", a), row.names = FALSE)
  M
}

agg_list <- list()
for (a in ARMS)
  agg_list[[a]] <- safe(paste("aggregate", a), aggregate_arm(raw[raw$arm == a, ], a))
AGG <- do.call(rbind, agg_list)
AGG$arm_f    <- factor(AGG$arm, levels = ARMS, labels = ARM_LAB[ARMS])
AGG$site_f   <- factor(AGG$site, levels = 1:4, labels = SITE_LAB)
AGG$origin_f <- factor(AGG$origin, levels = c("low", "high"),
                       labels = c("low initial density", "high initial density"))
write.csv(AGG, "out/bakuzis_agg_all.csv", row.names = FALSE)

# ===========================================================================
# Objective tests
# ===========================================================================
mono <- function(x, inc = TRUE, tol = 1e-6)
  if (inc) all(diff(x) >= -tol) else all(diff(x) <= tol)

# Reineke slope for one aggregated arm; returns per-group table
reineke_groups <- function(d) {
  rows <- list()
  for (o in unique(d$origin)) for (s in sort(unique(d$site))) {
    g <- d[d$origin == o & d$site == s, ]; g <- g[order(g$age), ]
    if (nrow(g) < 2) next
    keep <- c(TRUE, diff(g$QMD) > 0.2); g <- g[keep, ]
    if (nrow(g) < 3 || sd(log(g$QMD)) < 1e-3) next
    fit <- lm(log(TPH) ~ log(QMD), g)
    sl <- unname(coef(fit)[2]); ci <- confint(fit)[2, ]
    rows[[length(rows) + 1]] <- data.frame(
      site = s, origin = o, slope = sl, lo = ci[1], hi = ci[2],
      r2 = summary(fit)$r.squared, npt = nrow(g), stringsAsFactors = FALSE)
  }
  if (!length(rows)) return(NULL)
  do.call(rbind, rows)
}

sink("bakuzis_results.txt")
res <- list()
add <- function(...) res[[length(res) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

cat(strrep("=", 78), "\n")
cat("MODIFIED BAKUZIS MATRIX (Leary 1997) -- FVS NE, FROZEN PTBALT-REPAIRED BUILD\n")
cat("binary bin-balfix_20260804/FVSne  md5 b85e4c06995382425b01c756811d6d34\n")
cat("440 stratified plots ; 15 cycles x 10 yr = 150 yr ; five arms ; seed", SEED, "\n")
cat("Aggregated to age grid", min(AGE_GRID), "-", max(AGE_GRID),
    "by 10 yr, min n =", MIN_N, "\n")
cat(strrep("=", 78), "\n")

REIN <- list()
for (a in ARMS) {
  d <- AGG[AGG$arm == a, ]
  cat("\n\n############ ARM:", a, "############\n")
  sites <- sort(unique(d$site)); origins <- unique(d$origin)

  ## --- 1 Site curves ------------------------------------------------------
  cat("\n[1] SITE CURVES (HT-age)\n")
  bad_mono <- 0; inv_ht <- 0; nage <- 0
  for (o in origins) {
    for (s in sites) {
      g <- d[d$origin == o & d$site == s, ]
      if (!mono(g$HT, TRUE)) bad_mono <- bad_mono + 1
    }
    for (ag in sort(unique(d$age[d$origin == o]))) {
      h <- sapply(sites, function(s) {
        v <- d$HT[d$origin == o & d$site == s & d$age == ag]
        if (length(v)) v[1] else NA })
      if (!any(is.na(h))) { nage <- nage + 1; if (is.unsorted(h)) inv_ht <- inv_ht + 1 }
    }
  }
  gap <- mean(d$HT[d$site == 4], na.rm = TRUE) - mean(d$HT[d$site == 1], na.rm = TRUE)
  cat(sprintf("  HT non-monotone groups: %d/%d\n", bad_mono, length(sites) * length(origins)))
  cat(sprintf("  HT site-ordering inversions: %d/%d age points\n", inv_ht, nage))
  cat(sprintf("  mean HT separation S4-S1: %.2f m\n", gap))
  v1 <- (bad_mono == 0) && (inv_ht == 0)
  add(arm = a, relation = "1 Site curves (HT-age)",
      statistic = sprintf("inversions %d/%d; non-monotone %d; S4-S1 gap %.2f m",
                          inv_ht, nage, bad_mono, gap),
      decision = ifelse(v1, "PASS", "FLAG"))

  ## --- 2 Height-diameter --------------------------------------------------
  cat("\n[2] HEIGHT-DIAMETER\n")
  bad_hd <- 0; ngrp <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]; g <- g[order(g$QMD), ]
    ngrp <- ngrp + 1
    if (!mono(g$HT, TRUE)) bad_hd <- bad_hd + 1
  }
  inv_q <- 0; nq <- 0
  for (o in origins) for (ag in sort(unique(d$age[d$origin == o]))) {
    q <- sapply(sites, function(s) {
      v <- d$QMD[d$origin == o & d$site == s & d$age == ag]
      if (length(v)) v[1] else NA })
    if (!any(is.na(q))) { nq <- nq + 1; if (is.unsorted(q)) inv_q <- inv_q + 1 }
  }
  cat(sprintf("  HT non-monotone in QMD: %d/%d groups\n", bad_hd, ngrp))
  cat(sprintf("  QMD site-ordering inversions (REDPINE flaw test): %d/%d\n", inv_q, nq))
  v2 <- (bad_hd == 0) && (inv_q <= 0.2 * nq)
  add(arm = a, relation = "2 Height-diameter",
      statistic = sprintf("HT non-monotone %d/%d; QMD site inversions %d/%d",
                          bad_hd, ngrp, inv_q, nq),
      decision = ifelse(v2, "PASS", "FLAG"))

  ## --- 3 Sukachev ---------------------------------------------------------
  cat("\n[3] SUKACHEV (stems-age)\n")
  bad_t <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]
    if (!mono(g$TPH, FALSE)) bad_t <- bad_t + 1
  }
  loss <- sapply(sites, function(s) {
    g <- d[d$site == s, ]; a1 <- min(g$age); a2 <- max(g$age)
    1 - mean(g$TPH[g$age == a2]) / mean(g$TPH[g$age == a1]) })
  names(loss) <- paste0("S", sites)
  cat("  TPH non-monotone groups:", bad_t, "\n")
  cat("  fractional stem loss by site: ",
      paste(sprintf("%s=%.3f", names(loss), loss), collapse = "  "), "\n")
  ok_order <- loss[length(loss)] > loss[1]
  spread <- max(loss) - min(loss)
  cat(sprintf("  better-site-thins-faster: %s (spread %.3f)\n", ok_order, spread))
  v3 <- (bad_t == 0) && ok_order && spread > 0.03
  add(arm = a, relation = "3 Sukachev (stems-age)",
      statistic = sprintf("non-monotone %d; loss S1=%.3f S4=%.3f; spread %.3f",
                          bad_t, loss[1], loss[length(loss)], spread),
      decision = ifelse(v3, "PASS", "FLAG"))

  ## --- 4 Reineke ----------------------------------------------------------
  cat("\n[4] REINEKE (log TPH ~ log QMD), target -1.605, accept -2.2..-1.2\n")
  rg <- safe(paste("reineke", a), reineke_groups(d))
  if (!is.null(rg)) {
    rg$ok <- rg$slope > -2.2 & rg$slope < -1.2
    for (i in seq_len(nrow(rg)))
      cat(sprintf("  site %d %-4s: slope=%6.3f  [%6.3f,%6.3f]  R2=%.3f  %s\n",
                  rg$site[i], rg$origin[i], rg$slope[i], rg$lo[i], rg$hi[i],
                  rg$r2[i], ifelse(rg$ok[i], "OK", "FLAG")))
    msl <- mean(rg$slope); fl4 <- sum(!rg$ok); ng4 <- nrow(rg)
    rg$arm <- a; REIN[[a]] <- rg
  } else { msl <- NA; fl4 <- NA; ng4 <- 0 }
  cat(sprintf("  ARM MEAN REINEKE SLOPE = %.3f  (deviation from -1.605: %+.3f)\n",
              msl, msl - (-1.605)))
  cat(sprintf("  relative density SDI/SDImax: start %.3f, end %.3f, max %.3f\n",
              mean(d$RD[d$age == min(d$age)], na.rm = TRUE),
              mean(d$RD[d$age == max(d$age)], na.rm = TRUE),
              max(d$RD, na.rm = TRUE)))
  v4 <- isTRUE(fl4 == 0)
  add(arm = a, relation = "4 Reineke self-thinning",
      statistic = sprintf("mean slope %.3f (target -1.605); %d/%d groups outside -2.2..-1.2",
                          msl, fl4, ng4),
      decision = ifelse(v4, "PASS", "FLAG"))

  ## --- 5 Spacing percent --------------------------------------------------
  cat("\n[5] SPACING PERCENT (stems-height)\n")
  bad5 <- 0; ng5 <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]; g <- g[order(g$HT), ]
    ng5 <- ng5 + 1
    if (!mono(g$TPH, FALSE)) bad5 <- bad5 + 1
  }
  cat(sprintf("  reversals in TPH vs HT: %d/%d groups\n", bad5, ng5))
  v5 <- bad5 == 0
  add(arm = a, relation = "5 Spacing percent (stems-height)",
      statistic = sprintf("reversals %d/%d groups", bad5, ng5),
      decision = ifelse(v5, "PASS", "FLAG"))

  ## --- 6 Stocking guide ---------------------------------------------------
  cat("\n[6] STOCKING GUIDE (BA-stems)\n")
  hi <- max(d$BAPH, na.rm = TRUE)
  n_run <- sum(d$BAPH > 60, na.rm = TRUE)
  coll <- 0; ng6 <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]; ng6 <- ng6 + 1
    if (nrow(g) && max(g$BAPH) > 0 &&
        g$BAPH[nrow(g)] < 0.75 * max(g$BAPH)) coll <- coll + 1
  }
  cat(sprintf("  peak BAPH overall = %.1f m2/ha; points > 60 m2/ha = %d/%d\n",
              hi, n_run, nrow(d)))
  cat(sprintf("  groups collapsing >25%% below own BA peak: %d/%d\n", coll, ng6))
  v6 <- (n_run == 0) && (coll == 0)
  add(arm = a, relation = "6 Stocking guide (BA-stems)",
      statistic = sprintf("peak BA %.1f m2/ha; >60 m2/ha %d pts; BA collapse %d/%d groups",
                          hi, n_run, coll, ng6),
      decision = ifelse(v6, "PASS", "FLAG"))

  ## --- 7 Volume-age -------------------------------------------------------
  cat("\n[7] VOLUME-AGE (site ordering)\n")
  inv_v <- 0; nv <- 0; bad_vm <- 0
  for (o in origins) {
    for (s in sites) {
      g <- d[d$origin == o & d$site == s, ]
      if (!mono(g$VOL, TRUE)) bad_vm <- bad_vm + 1
    }
    for (ag in sort(unique(d$age[d$origin == o]))) {
      v <- sapply(sites, function(s) {
        z <- d$VOL[d$origin == o & d$site == s & d$age == ag]
        if (length(z)) z[1] else NA })
      if (!any(is.na(v))) { nv <- nv + 1; if (is.unsorted(v)) inv_v <- inv_v + 1 }
    }
  }
  cat(sprintf("  VOL site-ordering inversions: %d/%d age points\n", inv_v, nv))
  cat(sprintf("  VOL non-monotone groups: %d\n", bad_vm))
  v7 <- inv_v == 0
  add(arm = a, relation = "7 Volume-age",
      statistic = sprintf("site inversions %d/%d; declining-volume groups %d",
                          inv_v, nv, bad_vm),
      decision = ifelse(v7, "PASS", "FLAG"))

  ## --- 8 Eichhorn (aggregate view; raw view in the addendum) --------------
  cat("\n[8] EICHHORN (volume-height site independence), aggregate view\n")
  cvs <- c()
  for (o in origins) {
    hi0 <- max(sapply(sites, function(s) min(d$HT[d$origin == o & d$site == s])))
    hx0 <- min(sapply(sites, function(s) max(d$HT[d$origin == o & d$site == s])))
    if (hx0 > hi0) {
      grid <- seq(hi0, hx0, length.out = 8)
      cv <- sapply(grid, function(hg) {
        vv <- sapply(sites, function(s) {
          g <- d[d$origin == o & d$site == s, ]
          approx(g$HT, g$VOL, hg, rule = 2)$y })
        sd(vv) / mean(vv) })
      cat(sprintf("  %-4s: HT window %.1f-%.1f m, mean CV = %.1f%%\n",
                  o, hi0, hx0, 100 * mean(cv)))
      cvs <- c(cvs, mean(cv))
    }
  }
  mcv <- mean(cvs)
  cat(sprintf("  ARM MEAN Eichhorn CV = %.1f%% (threshold 15%%)\n", 100 * mcv))

  ## --- 9 Yield-density ----------------------------------------------------
  cat("\n[9] YIELD-DENSITY (volume-stems)\n")
  sl9 <- c()
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]
    if (nrow(g) < 3) next
    sl9 <- c(sl9, unname(coef(lm(log(VOL) ~ log(TPH), g))[2]))
  }
  cat(sprintf("  mean d log VOL / d log TPH = %.3f\n", mean(sl9)))
  add(arm = a, relation = "9 Yield-density (volume-stems)",
      statistic = sprintf("mean dlogV/dlogN = %.3f", mean(sl9)),
      decision = "SUPPORTING")
}

RES <- do.call(rbind, res)
cat("\n\n", strrep("=", 78), "\nSUMMARY MATRIX\n", strrep("=", 78), "\n", sep = "")
print(RES, right = FALSE)
sink()
write.csv(RES, "bakuzis_matrix_table.csv", row.names = FALSE)
REINEKE <- do.call(rbind, REIN)
write.csv(REINEKE, "reineke_groups_new.csv", row.names = FALSE)

# ===========================================================================
# ADDENDUM: raw Eichhorn, SDI gate, BAMAX ceiling -- all five arms
# ===========================================================================
sink("bakuzis_addendum.txt")
out <- list(); addx <- function(...) out[[length(out) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
cat(strrep("=", 78), "\nADDENDUM: Eichhorn (raw), SDI gate, BAMAX ceiling -- five arms\n",
    strrep("=", 78), "\n", sep = "")

cat("\n\n[8 EICHHORN, raw plot-level, 2 m height bins, min 20 plot-cycles per bin]\n")
eich_rows <- list()
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$HTm) & is.finite(raw$VOLm) & raw$VOLm > 0, ]
  d$hbin <- floor(d$HTm / 2) * 2 + 1
  tb <- aggregate(list(VOL = d$VOLm), by = list(hbin = d$hbin, site = d$site_cls), FUN = mean)
  nb <- aggregate(list(n = d$VOLm), by = list(hbin = d$hbin, site = d$site_cls), FUN = length)
  tb <- merge(tb, nb); tb <- tb[tb$n >= 20, ]
  ok <- names(which(table(tb$hbin) == length(unique(d$site_cls))))
  tb <- tb[as.character(tb$hbin) %in% ok, ]
  cv <- tapply(tb$VOL, tb$hbin, function(x) sd(x) / mean(x))
  mcv <- mean(cv, na.rm = TRUE)
  cat(sprintf("\n  %-14s bins with all 4 site classes: %d (HT %.0f-%.0f m)\n",
              a, length(cv), min(tb$hbin), max(tb$hbin)))
  cat(sprintf("  %-14s per-bin CV: %s\n", "",
              paste(sprintf("%.0fm=%.1f%%", as.numeric(names(cv)), 100 * cv), collapse = " ")))
  cat(sprintf("  %-14s MEAN Eichhorn CV = %.1f%%  -> %s\n", "", 100 * mcv,
              ifelse(mcv < 0.15, "PASS", "FLAG")))
  addx(arm = a, relation = "8 Eichhorn (volume-height)",
       statistic = sprintf("mean CV of VOL at matched HT = %.1f%% over %d common 2 m bins (threshold 15%%)",
                           100 * mcv, length(cv)),
       decision = ifelse(mcv < 0.15, "PASS", "FLAG"))
  tb$arm <- a; eich_rows[[a]] <- tb
}

cat("\n\n[SDI GATE] vls/morts.f90: density mortality engages only above PMSDIL = 55% of SDImax\n")
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$RD), ]
  f  <- mean(d$RD < 0.55, na.rm = TRUE)
  fl <- mean(d$RD[d$dens_cls == "low"]  < 0.55, na.rm = TRUE)
  fh <- mean(d$RD[d$dens_cls == "high"] < 0.55, na.rm = TRUE)
  bysite <- tapply(d$RD, d$site_cls, function(x) mean(x < 0.55, na.rm = TRUE))
  cat(sprintf("  %-14s below 55%% gate: all=%.1f%%  low-density=%.1f%%  high-density=%.1f%%\n",
              a, 100 * f, 100 * fl, 100 * fh))
  cat(sprintf("  %-14s by site class: %s\n", "",
              paste(sprintf("S%s=%.0f%%", names(bysite), 100 * bysite), collapse = "  ")))
  addx(arm = a, relation = "DIAG: SDI gate (55% SDImax)",
       statistic = sprintf("%.1f%% of plot-cycles below the gate (low-density %.1f%%, high-density %.1f%%)",
                           100 * f, 100 * fl, 100 * fh),
       decision = "DIAGNOSTIC")
}

cat("\n\n[BAMAX CEILING] vls/morts.f90 lines 750-786; BAMAX = SDImax * 0.5454154 * 0.85\n")
for (a in ARMS) {
  d <- raw[raw$arm == a & is.finite(raw$BA_frac), ]
  pinned <- mean(abs(d$BA - d$BAMAX_ft2) <= 1.0, na.rm = TRUE)
  over   <- mean(d$BA_frac > 0.99, na.rm = TRUE)
  cat(sprintf("  %-14s max BA/BAMAX = %.3f ; pinned within 1 ft2/ac = %.1f%% ; above 0.99 = %.1f%%\n",
              a, max(d$BA_frac, na.rm = TRUE), 100 * pinned, 100 * over))
  addx(arm = a, relation = "DIAG: BAMAX ceiling",
       statistic = sprintf("%.1f%% of plot-cycles pinned within 1 ft2 ac-1 of BAMAX (max BA/BAMAX %.3f)",
                           100 * pinned, max(d$BA_frac, na.rm = TRUE)),
       decision = "DIAGNOSTIC")
}
RESX <- do.call(rbind, out)
cat("\n\n", strrep("=", 78), "\n", sep = "")
print(RESX, right = FALSE)
sink()
write.csv(RESX, "bakuzis_addendum_table.csv", row.names = FALSE)

# ===========================================================================
# OLD versus NEW: the three arms that were run before the PTBALT repair
# ===========================================================================
OLDARMS <- c("fvs_base", "fvs_regional", "organon")
oldpath <- file.path(OLDDIR, "out", "bakuzis_traj_raw.csv")
OLDAGG <- NULL
if (file.exists(oldpath)) {
  oraw <- safe("read old raw", prep(oldpath))
  if (!is.null(oraw)) {
    ol <- list()
    for (a in OLDARMS) if (a %in% unique(oraw$arm))
      ol[[a]] <- safe(paste("old aggregate", a),
                      aggregate_arm(oraw[oraw$arm == a, ], a, write = FALSE))
    OLDAGG <- do.call(rbind, ol)
  }
} else logerr("old run", simpleError(paste("missing", oldpath)))

sink("bakuzis_oldnew.txt")
cat(strrep("=", 78), "\nOLD (bin-fixed, pre-PTBALT) versus NEW (bin-balfix, post-PTBALT)\n",
    strrep("=", 78), "\n", sep = "")
ON <- NULL
if (!is.null(OLDAGG)) {
  rows <- list()
  for (a in OLDARMS) {
    ro <- safe(paste("old reineke", a), reineke_groups(OLDAGG[OLDAGG$arm == a, ]))
    rn <- REIN[[a]]
    if (is.null(ro) || is.null(rn)) next
    m <- merge(ro[, c("site", "origin", "slope", "r2")],
               rn[, c("site", "origin", "slope", "r2")],
               by = c("site", "origin"), suffixes = c("_old", "_new"))
    m$arm <- a; rows[[a]] <- m
    cat(sprintf("\n%s\n", a))
    cat(sprintf("  mean Reineke slope  OLD %7.3f   NEW %7.3f   change %+.3f\n",
                mean(ro$slope), mean(rn$slope), mean(rn$slope) - mean(ro$slope)))
    cat(sprintf("  groups outside -2.2..-1.2  OLD %d/%d   NEW %d/%d\n",
                sum(ro$slope <= -2.2 | ro$slope >= -1.2), nrow(ro),
                sum(rn$slope <= -2.2 | rn$slope >= -1.2), nrow(rn)))
    for (i in seq_len(nrow(m)))
      cat(sprintf("    site %d %-4s: %7.3f -> %7.3f  (%+.3f)\n",
                  m$site[i], m$origin[i], m$slope_old[i], m$slope_new[i],
                  m$slope_new[i] - m$slope_old[i]))
  }
  ON <- do.call(rbind, rows)
  if (!is.null(ON)) write.csv(ON, "oldnew_reineke.csv", row.names = FALSE)

  cat("\n\n[Terminal-state comparison at age 190, high initial density]\n")
  for (a in OLDARMS) {
    o <- OLDAGG[OLDAGG$arm == a & OLDAGG$age == 190 & OLDAGG$origin == "high", ]
    n <- AGG[AGG$arm == a & AGG$age == 190 & AGG$origin == "high", ]
    if (!nrow(o) || !nrow(n)) next
    cat(sprintf("  %-14s QMD %5.1f -> %5.1f cm (%+.1f%%) ; TPH %6.1f -> %6.1f ; BA %5.1f -> %5.1f m2/ha ; HT %5.1f -> %5.1f m\n",
                a, mean(o$QMD), mean(n$QMD),
                100 * (mean(n$QMD) / mean(o$QMD) - 1),
                mean(o$TPH), mean(n$TPH), mean(o$BAPH), mean(n$BAPH),
                mean(o$HT), mean(n$HT)))
  }
} else cat("OLD run not available\n")
sink()

# ===========================================================================
# Figures
# ===========================================================================
pal_site <- if (has_ggsci) ggsci::pal_npg()(4) else
  c("#3B4992", "#008B45", "#EE7733", "#BB0021")
names(pal_site) <- SITE_LAB
pal_arm <- if (has_ggsci) ggsci::pal_npg()(5) else
  c("#3B4992", "#008B45", "#EE7733", "#BB0021", "#631879")

base_theme <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

gsave <- function(p, f, w, h) {
  safe(paste("figure", f), {
    png(file.path("fig", f), width = w, height = h, units = "in", res = 300)
    print(p); dev.off(); cat("wrote fig/", f, "\n", sep = "")
  })
}

D <- AGG[AGG$origin == "high", ]

# Fig 1: site curves, five arms
gsave(ggplot(AGG, aes(age, HT, colour = site_f, linetype = origin_f)) +
        geom_line(linewidth = 0.9) + facet_wrap(~arm_f, nrow = 1) +
        scale_colour_manual(values = pal_site, name = "Site class") +
        scale_linetype_discrete(name = NULL) +
        labs(x = "Stand age (yr)", y = expression(Dominant~height~(m))) +
        base_theme, "fig1_site_curves.png", 13, 4.2)

# Fig 2: Sukachev
gsave(ggplot(AGG, aes(age, TPH, colour = site_f, linetype = origin_f)) +
        geom_line(linewidth = 0.9) + facet_wrap(~arm_f, nrow = 1) +
        scale_colour_manual(values = pal_site, name = "Site class") +
        scale_linetype_discrete(name = NULL) +
        labs(x = "Stand age (yr)", y = expression(Stem~density~(trees~ha^-1))) +
        base_theme, "fig2_sukachev.png", 13, 4.2)

# Fig 3: five-arm Reineke panel
mkref <- function(dd) {
  a0 <- dd$TPH[which.max(dd$QMD)]; q0 <- max(dd$QMD)
  data.frame(QMD = c(10, 60), TPH = a0 * (c(10, 60) / q0)^(-1.605))
}
refl <- do.call(rbind, lapply(split(D, D$arm_f), function(dd) {
  if (!nrow(dd)) return(NULL); r <- mkref(dd); r$arm_f <- dd$arm_f[1]; r }))
gsave(ggplot(D, aes(QMD, TPH, colour = site_f)) +
        geom_line(linewidth = 0.9) + geom_point(size = 0.8) +
        geom_line(data = refl, aes(QMD, TPH), colour = "grey30",
                  linetype = "dashed", inherit.aes = FALSE) +
        facet_wrap(~arm_f, nrow = 1) +
        scale_x_log10() + scale_y_log10() +
        scale_colour_manual(values = pal_site, name = "Site class") +
        labs(x = expression(Quadratic~mean~diameter~(cm)),
             y = expression(Stem~density~(trees~ha^-1))) +
        base_theme, "fig3_reineke_5arm.png", 13, 4.2)

# Fig 4: stocking guide
gsave(ggplot(AGG, aes(TPH, BAPH, colour = site_f, linetype = origin_f)) +
        geom_line(linewidth = 0.9) + facet_wrap(~arm_f, nrow = 1) +
        geom_hline(yintercept = 60, colour = "grey40", linetype = "dotted") +
        scale_colour_manual(values = pal_site, name = "Site class") +
        scale_linetype_discrete(name = NULL) +
        labs(x = expression(Stem~density~(trees~ha^-1)),
             y = expression(Basal~area~(m^2~ha^-1))) +
        base_theme, "fig4_stocking.png", 13, 4.2)

# Fig 5: Eichhorn on raw binned data
E <- do.call(rbind, eich_rows)
E$site_f <- factor(E$site, levels = 1:4, labels = SITE_LAB)
E$arm_f  <- factor(E$arm, levels = ARMS, labels = ARM_LAB[ARMS])
gsave(ggplot(E, aes(hbin, VOL, colour = site_f)) +
        geom_line(linewidth = 0.9) + geom_point(size = 1.1) +
        facet_wrap(~arm_f, nrow = 1) +
        scale_colour_manual(values = pal_site, name = "Site class") +
        labs(x = expression(Dominant~height~(m)),
             y = expression(Gross~cubic~volume~(m^3~ha^-1))) +
        base_theme, "fig5_eichhorn.png", 13, 4.2)

# Fig 6: relative density
gsave(ggplot(AGG, aes(age, RD, colour = site_f, linetype = origin_f)) +
        geom_line(linewidth = 0.9) + facet_wrap(~arm_f, nrow = 1) +
        geom_hline(yintercept = c(0.35, 0.55), colour = "grey40", linetype = "dotted") +
        scale_colour_manual(values = pal_site, name = "Site class") +
        scale_linetype_discrete(name = NULL) +
        labs(x = "Stand age (yr)", y = "Relative density (SDI / SDImax)") +
        base_theme, "fig6_relative_density.png", 13, 4.2)

# Fig 7: BAMAX ceiling
B <- raw[is.finite(raw$BA_frac), ]
B$arm_f  <- factor(B$arm, levels = ARMS, labels = ARM_LAB[ARMS])
B$site_f <- factor(B$site_cls, levels = 1:4, labels = SITE_LAB)
aggB <- aggregate(list(BA_frac = B$BA_frac),
                  by = list(Age = B$Age, arm_f = B$arm_f, site_f = B$site_f), FUN = mean)
aggB <- aggB[aggB$Age <= 200, ]
gsave(ggplot(aggB, aes(Age, BA_frac, colour = site_f)) +
        geom_line(linewidth = 0.9) + facet_wrap(~arm_f, nrow = 1) +
        geom_hline(yintercept = 1, colour = "grey25", linetype = "dashed") +
        scale_colour_manual(values = pal_site, name = "Site class") +
        labs(x = "Stand age (yr)",
             y = "Basal area relative to the BAMAX ceiling") +
        base_theme, "fig7_bamax_ceiling.png", 13, 4.2)

# Fig 8: OLD versus NEW Reineke comparison
if (!is.null(ON)) {
  ON$arm_f <- factor(ON$arm, levels = OLDARMS, labels = ARM_LAB[OLDARMS])
  ON$grp <- paste0("S", ON$site, " ", ON$origin)
  ON$grp_f <- factor(ON$grp, levels = rev(sort(unique(ON$grp))))
  armmean <- aggregate(cbind(slope_old, slope_new) ~ arm_f, ON, mean)
  gsave(ggplot(ON) +
          annotate("rect", xmin = -2.2, xmax = -1.2, ymin = -Inf, ymax = Inf,
                   fill = "grey90", alpha = 0.7) +
          geom_vline(xintercept = -1.605, colour = "grey30", linetype = "dashed") +
          geom_segment(aes(x = slope_old, xend = slope_new,
                           y = grp_f, yend = grp_f),
                       colour = "grey55", linewidth = 0.7,
                       arrow = arrow(length = unit(0.06, "in"), type = "closed")) +
          geom_point(aes(slope_old, grp_f, shape = "pre-repair"),
                     colour = pal_arm[1], size = 2.4) +
          geom_point(aes(slope_new, grp_f, shape = "post-repair"),
                     colour = pal_arm[4], size = 2.4) +
          facet_wrap(~arm_f, nrow = 1) +
          scale_shape_manual(values = c(`pre-repair` = 1, `post-repair` = 16),
                             name = NULL,
                             breaks = c("pre-repair", "post-repair")) +
          labs(x = "Reineke slope, d log(TPH) / d log(QMD)",
               y = "Site class and initial density stratum") +
          base_theme, "fig8_oldnew_reineke.png", 11, 4.2)
  write.csv(armmean, "oldnew_armmeans.csv", row.names = FALSE)
}

cat("\nDONE\n")
cat(sprintf("finished %s\n", format(Sys.time())), file = ERR, append = TRUE)
