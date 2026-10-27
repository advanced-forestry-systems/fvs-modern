#!/usr/bin/env Rscript
# bakuzis_analyze.R
# Prepare group trajectories, run the modified Bakuzis matrix tests
# (Leary 1997) on each FVS NE model arm, and render figures.
#
# Input : out/bakuzis_traj_raw.csv  (from bakuzis_run.py)
# Output: traj_<arm>.csv, bakuzis_results.txt, bakuzis_matrix_table.csv,
#         fig/*.png  (300 dpi)

suppressPackageStartupMessages({
  library(ggplot2)
})
has_ggsci <- requireNamespace("ggsci", quietly = TRUE)

setwd(Sys.getenv("BAKDIR", "/fs/scratch/PUOM0008/crsfaaron/bakuzis_20260804"))
dir.create("fig", showWarnings = FALSE)

# --- unit conversions (FVS imperial -> metric) -----------------------------
K_TPH  <- 2.4710538   # trees/acre    -> trees/ha
K_BAPH <- 0.2295684   # ft2/acre      -> m2/ha
K_QMD  <- 2.54        # in            -> cm
K_HT   <- 0.3048      # ft            -> m
K_VOL  <- 0.06997     # ft3/acre      -> m3/ha

raw <- read.csv("out/bakuzis_traj_raw.csv", stringsAsFactors = FALSE)
raw <- raw[is.finite(raw$Tpa) & raw$Tpa > 0 & is.finite(raw$QMD) & raw$QMD > 0, ]
raw$TPH  <- raw$Tpa   * K_TPH
raw$BAPH <- raw$BA    * K_BAPH
raw$QMDc <- raw$QMD   * K_QMD
raw$HTm  <- raw$TopHt * K_HT
raw$VOLm <- raw$TCuFt * K_VOL
raw$RD   <- ifelse(raw$SDIMax > 0, raw$SDI / raw$SDIMax, NA)

ARMS <- c("fvs_base", "fvs_regional", "organon")
ARMS <- ARMS[ARMS %in% unique(raw$arm)]
ARM_LAB <- c(fvs_base = "fvs_base", fvs_regional = "fvs_regional",
             organon = "organon (Greg)")

# --- aggregate to a common stand-age grid ----------------------------------
# Each plot has its own AGE0; interpolate onto a shared grid, then average
# within (arm, site class, density class).
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

agg_list <- list()
for (a in ARMS) {
  ra <- raw[raw$arm == a, ]
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
  write.csv(M[, c("site", "origin", "age", "HT", "TPH", "BAPH", "QMD",
                  "VOL", "SDI", "SDIMAX", "RD", "n")],
            sprintf("traj_%s.csv", a), row.names = FALSE)
  agg_list[[a]] <- M
}
AGG <- do.call(rbind, agg_list)
AGG$arm_f  <- factor(AGG$arm, levels = ARMS, labels = ARM_LAB[ARMS])
AGG$site_f <- factor(AGG$site, levels = 1:4,
                     labels = c("S1 (SI 47)", "S2 (SI 61)",
                                "S3 (SI 70)", "S4 (SI 83)"))
AGG$origin_f <- factor(AGG$origin, levels = c("low", "high"),
                       labels = c("low initial density", "high initial density"))
write.csv(AGG, "out/bakuzis_agg_all.csv", row.names = FALSE)

# ===========================================================================
# Objective tests
# ===========================================================================
sink("bakuzis_results.txt")
res <- list()
add <- function(...) res[[length(res) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)

cat(strrep("=", 74), "\n")
cat("MODIFIED BAKUZIS MATRIX (Leary 1997) -- FVS NE, repaired binary\n")
cat("bin-fixed_20260804/FVSne ; 440 plots ; 15 cycles x 10 yr = 150 yr\n")
cat("Aggregated to age grid", min(AGE_GRID), "-", max(AGE_GRID),
    "by 10 yr, min n =", MIN_N, "\n")
cat(strrep("=", 74), "\n")

mono <- function(x, inc = TRUE, tol = 1e-6)
  if (inc) all(diff(x) >= -tol) else all(diff(x) <= tol)

for (a in ARMS) {
  d <- AGG[AGG$arm == a, ]
  cat("\n\n############ ARM:", a, "############\n")
  sites <- sort(unique(d$site)); origins <- unique(d$origin)

  ## --- 1 Site curves: HT-age monotone + site ordering ---------------------
  cat("\n[1] SITE CURVES (HT-age)\n")
  bad_mono <- 0; inv_ht <- 0; nage <- 0
  for (o in origins) {
    for (s in sites) {
      g <- d[d$origin == o & d$site == s, ]
      if (!mono(g$HT, TRUE)) bad_mono <- bad_mono + 1
    }
    ages <- sort(unique(d$age[d$origin == o]))
    for (ag in ages) {
      h <- sapply(sites, function(s) {
        v <- d$HT[d$origin == o & d$site == s & d$age == ag]
        if (length(v)) v[1] else NA })
      if (!any(is.na(h))) { nage <- nage + 1; if (is.unsorted(h)) inv_ht <- inv_ht + 1 }
    }
  }
  # site separation: mean HT gap between S4 and S1 at matched age
  gap <- mean(d$HT[d$site == 4] , na.rm = TRUE) - mean(d$HT[d$site == 1], na.rm = TRUE)
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
  # site inversion in QMD: is a poorer site producing a LARGER QMD at age?
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

  ## --- 3 Sukachev: stems-age ---------------------------------------------
  cat("\n[3] SUKACHEV (stems-age)\n")
  bad_t <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]
    if (!mono(g$TPH, FALSE)) bad_t <- bad_t + 1
  }
  # site sensitivity: fractional stem loss over the window by site class
  loss <- sapply(sites, function(s) {
    g <- d[d$site == s, ]
    a1 <- min(g$age); a2 <- max(g$age)
    1 - mean(g$TPH[g$age == a2]) / mean(g$TPH[g$age == a1]) })
  names(loss) <- paste0("S", sites)
  cat("  TPH non-monotone groups:", bad_t, "\n")
  cat("  fractional stem loss over window by site: ",
      paste(sprintf("%s=%.3f", names(loss), loss), collapse = "  "), "\n")
  # better sites should self-thin FASTER
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
  sl_all <- c(); fl4 <- 0; ng4 <- 0
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]; g <- g[order(g$age), ]
    keep <- c(TRUE, diff(g$QMD) > 0.2)
    g <- g[keep, ]
    if (nrow(g) < 3 || sd(log(g$QMD)) < 1e-3) { cat("  n/a\n"); next }
    fit <- lm(log(TPH) ~ log(QMD), g)
    sl <- unname(coef(fit)[2]); ci <- confint(fit)[2, ]
    ok <- sl > -2.2 && sl < -1.2
    ng4 <- ng4 + 1; if (!ok) fl4 <- fl4 + 1
    sl_all <- c(sl_all, sl)
    cat(sprintf("  site %d %-4s: slope=%6.3f  [%6.3f,%6.3f]  R2=%.3f  %s\n",
                s, o, sl, ci[1], ci[2], summary(fit)$r.squared,
                ifelse(ok, "OK", "FLAG")))
  }
  msl <- mean(sl_all)
  cat(sprintf("  ARM MEAN REINEKE SLOPE = %.3f  (deviation from -1.605: %+.3f)\n",
              msl, msl - (-1.605)))
  # relative density: is the stand ever near the self-thinning boundary?
  cat(sprintf("  relative density SDI/SDImax: start %.3f, end %.3f, max %.3f\n",
              mean(d$RD[d$age == min(d$age)], na.rm = TRUE),
              mean(d$RD[d$age == max(d$age)], na.rm = TRUE),
              max(d$RD, na.rm = TRUE)))
  v4 <- fl4 == 0
  add(arm = a, relation = "4 Reineke self-thinning",
      statistic = sprintf("mean slope %.3f (target -1.605); %d/%d groups outside -2.2..-1.2",
                          msl, fl4, ng4),
      decision = ifelse(v4, "PASS", "FLAG"))

  ## --- 5 Spacing percent (stems-height) ----------------------------------
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

  ## --- 6 Stocking guide (BA-stems) ---------------------------------------
  cat("\n[6] STOCKING GUIDE (BA-stems)\n")
  bamax <- tapply(d$BAPH, list(d$site, d$origin), max)
  baend <- tapply(d$BAPH, list(d$site, d$origin), function(x) x[length(x)])
  cat("  max BAPH by site x origin (m2/ha):\n"); print(round(bamax, 1))
  hi <- max(d$BAPH, na.rm = TRUE)
  # runaway test: any group exceeding a generous NE maximum (~60 m2/ha)
  n_run <- sum(d$BAPH > 60, na.rm = TRUE)
  # collapse test: BA falling more than 25% from its own peak
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
  cat(sprintf("  VOL non-monotone (declining standing volume): %d groups\n", bad_vm))
  v7 <- inv_v == 0
  add(arm = a, relation = "7 Volume-age",
      statistic = sprintf("site inversions %d/%d; declining-volume groups %d",
                          inv_v, nv, bad_vm),
      decision = ifelse(v7, "PASS", "FLAG"))

  ## --- 8 Eichhorn ---------------------------------------------------------
  cat("\n[8] EICHHORN (volume-height site independence)\n")
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
  v8 <- mcv < 0.15
  add(arm = a, relation = "8 Eichhorn (volume-height)",
      statistic = sprintf("mean CV of VOL at matched HT = %.1f%% (threshold 15%%)", 100 * mcv),
      decision = ifelse(v8, "PASS", "FLAG"))

  ## --- 9 Yield-density ----------------------------------------------------
  cat("\n[9] YIELD-DENSITY (volume-stems)\n")
  sl9 <- c()
  for (o in origins) for (s in sites) {
    g <- d[d$origin == o & d$site == s, ]
    if (nrow(g) < 3) next
    sl9 <- c(sl9, unname(coef(lm(log(VOL) ~ log(TPH), g))[2]))
  }
  cat(sprintf("  mean d log VOL / d log TPH = %.3f (self-thinning line ~ -0.5)\n",
              mean(sl9)))
  add(arm = a, relation = "9 Yield-density (volume-stems)",
      statistic = sprintf("mean dlogV/dlogN = %.3f", mean(sl9)),
      decision = "SUPPORTING")
}

RES <- do.call(rbind, res)
cat("\n\n", strrep("=", 74), "\nSUMMARY MATRIX\n", strrep("=", 74), "\n", sep = "")
print(RES, right = FALSE)
sink()
write.csv(RES, "bakuzis_matrix_table.csv", row.names = FALSE)

# ===========================================================================
# Figures
# ===========================================================================
pal_site <- if (has_ggsci) ggsci::pal_npg()(4) else
  c("#3B4992", "#008B45", "#EE7733", "#BB0021")
names(pal_site) <- levels(AGG$site_f)

base_theme <- theme_classic(base_size = 11) +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom",
        panel.grid.major = element_line(colour = "grey92", linewidth = 0.3))

D <- AGG[AGG$origin == "high", ]   # primary panel set: high initial density

gsave <- function(p, f, w, h) {
  png(file.path("fig", f), width = w, height = h, units = "in", res = 300)
  print(p); dev.off(); cat("wrote fig/", f, "\n", sep = "")
}

# Fig 1: site curves
p1 <- ggplot(AGG, aes(age, HT, colour = site_f, linetype = origin_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  scale_colour_manual(values = pal_site, name = "Site class") +
  scale_linetype_discrete(name = NULL) +
  labs(x = "Stand age (yr)", y = expression(Dominant~height~(m))) +
  base_theme
gsave(p1, "fig1_site_curves.png", 9, 4)

# Fig 2: Sukachev
p2 <- ggplot(AGG, aes(age, TPH, colour = site_f, linetype = origin_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  scale_colour_manual(values = pal_site, name = "Site class") +
  scale_linetype_discrete(name = NULL) +
  labs(x = "Stand age (yr)", y = expression(Stem~density~(trees~ha^-1))) +
  base_theme
gsave(p2, "fig2_sukachev.png", 9, 4)

# Fig 3: Reineke with the -1.605 reference line
ref <- data.frame(QMD = c(10, 60))
mkref <- function(dd) {
  a0 <- dd$TPH[which.max(dd$QMD)]; q0 <- max(dd$QMD)
  data.frame(QMD = ref$QMD, TPH = a0 * (ref$QMD / q0)^(-1.605))
}
refl <- do.call(rbind, lapply(split(D, list(D$arm_f)), function(dd) {
  r <- mkref(dd); r$arm_f <- dd$arm_f[1]; r }))
p3 <- ggplot(D, aes(QMD, TPH, colour = site_f)) +
  geom_line(linewidth = 0.9) + geom_point(size = 0.8) +
  geom_line(data = refl, aes(QMD, TPH), colour = "grey30",
            linetype = "dashed", inherit.aes = FALSE) +
  facet_wrap(~arm_f) +
  scale_x_log10() + scale_y_log10() +
  scale_colour_manual(values = pal_site, name = "Site class") +
  labs(x = expression(Quadratic~mean~diameter~(cm)),
       y = expression(Stem~density~(trees~ha^-1))) +
  base_theme
gsave(p3, "fig3_reineke.png", 9, 4)

# Fig 4: stocking guide
p4 <- ggplot(AGG, aes(TPH, BAPH, colour = site_f, linetype = origin_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  geom_hline(yintercept = 60, colour = "grey40", linetype = "dotted") +
  scale_colour_manual(values = pal_site, name = "Site class") +
  scale_linetype_discrete(name = NULL) +
  labs(x = expression(Stem~density~(trees~ha^-1)),
       y = expression(Basal~area~(m^2~ha^-1))) +
  base_theme
gsave(p4, "fig4_stocking.png", 9, 4)

# Fig 5: Eichhorn
p5 <- ggplot(AGG, aes(HT, VOL, colour = site_f, linetype = origin_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  scale_colour_manual(values = pal_site, name = "Site class") +
  scale_linetype_discrete(name = NULL) +
  labs(x = expression(Dominant~height~(m)),
       y = expression(Gross~cubic~volume~(m^3~ha^-1))) +
  base_theme
gsave(p5, "fig5_eichhorn.png", 9, 4)

# Fig 6: relative density trajectory (SDI / SDImax)
p6 <- ggplot(AGG, aes(age, RD, colour = site_f, linetype = origin_f)) +
  geom_line(linewidth = 0.9) + facet_wrap(~arm_f) +
  geom_hline(yintercept = c(0.35, 0.60), colour = "grey40", linetype = "dotted") +
  scale_colour_manual(values = pal_site, name = "Site class") +
  scale_linetype_discrete(name = NULL) +
  labs(x = "Stand age (yr)", y = "Relative density (SDI / SDImax)") +
  base_theme
gsave(p6, "fig6_relative_density.png", 9, 4)

cat("\nDONE\n")
