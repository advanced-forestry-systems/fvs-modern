#!/usr/bin/env Rscript
# Bakuzis-style biological-realism assessment of the four CONUS FVS arms (b1/b2 growth x
# native/gompit mortality) across a SITE_INDEX productivity gradient. Tests the law-like
# relations the emitted data supports: Reineke self-thinning (upper-boundary slope near -1.605,
# realized max density), site ordering (better sites not overtaken in BA/AGB), and yield-density
# (AGB rises with SDI). Eichhorn (volume vs mean height) is NOT testable: stand height is not
# emitted (CCH_MEAN empty) -> reported as FLAG-untestable, not pass.
# Memory bounded: incremental aggregation + 2% frontier sample. Seed set. Headless png.
suppressPackageStartupMessages({library(data.table)})
set.seed(20260629)
D    <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj"
SBV  <- "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant"
OUT  <- D
elog <- file.path(OUT, "error_log_bakuzis.txt")
logerr <- function(w, e) cat(sprintf("[%s] %s: %s\n", Sys.time(), w, conditionMessage(e)),
                             file = elog, append = TRUE)
arms <- list(
  b1_native = file.path(D, "out_conus_eq",          "*conus_b1_metrics.csv"),
  b2_native = file.path(D, "out_conus_eq",          "*conus_b2_metrics.csv"),
  b2_gompit = file.path(D, "out_conus_eq_gompit",   "*conus_b2_gompit_metrics.csv"),
  b1_gompit = file.path(D, "out_conus_eq_b1_gompit","*conus_b1_gompit_metrics.csv"))

## ---- per-stand SITE_INDEX map from all standinit files, + global tertiles ----
si_map <- rbindlist(lapply(Sys.glob(file.path(SBV, "standinit_*.csv")), function(f) tryCatch({
  d <- fread(f, select = c("STAND_CN","SITE_INDEX"),
             colClasses = list(character = "STAND_CN"))
  d[is.finite(SITE_INDEX) & SITE_INDEX > 0]
}, error = function(e) { logerr(paste("si", f), e); NULL })), fill = TRUE)
si_map <- unique(si_map, by = "STAND_CN")
qt <- quantile(si_map$SITE_INDEX, c(1/3, 2/3), na.rm = TRUE)
si_map[, site_class := fifelse(SITE_INDEX <= qt[1], "low",
                       fifelse(SITE_INDEX <= qt[2], "mid", "high"))]
setkey(si_map, STAND_CN)
cat(sprintf("SITE_INDEX tertile cutpoints: %.1f / %.1f ; n stands=%d\n", qt[1], qt[2], nrow(si_map)))

sel <- c("STAND_CN","PROJ_YEAR","BA_FT2AC","QMD_IN","TPH","AGB_TONS_AC")
add <- function(d) {
  d <- d[is.finite(BA_FT2AC) & is.finite(QMD_IN) & is.finite(TPH) & QMD_IN > 0 & TPH > 0]
  d[, TPA := TPH / 2.471][, SDI := TPA * (QMD_IN / 10)^1.605][]
}

agg_list <- list(); term_list <- list(); front_list <- list()
for (a in names(arms)) {
  fs <- Sys.glob(arms[[a]])
  AGG <- NULL; TERM <- NULL; FRONT <- NULL
  for (f in fs) tryCatch({
    d <- add(fread(f, select = sel, colClasses = list(character = "STAND_CN")))
    d <- d[si_map, on = "STAND_CN", nomatch = 0L]
    # (a) means by site_class x PROJ_YEAR
    g <- d[, .(BA = sum(BA_FT2AC), TPH = sum(TPH), QMD = sum(QMD_IN),
               SDI = sum(SDI), AGB = sum(AGB_TONS_AC[is.finite(AGB_TONS_AC)]),
               nAGB = sum(is.finite(AGB_TONS_AC)), n = .N),
           by = .(site_class, PROJ_YEAR)]
    AGG <- if (is.null(AGG)) g else rbindlist(list(AGG, g))
    # (b) terminal-year per-stand SDI (over-densification)
    tm <- d[d[, .I[PROJ_YEAR == max(PROJ_YEAR)], by = STAND_CN]$V1,
            .(STAND_CN, site_class, SDI)]
    TERM <- if (is.null(TERM)) tm else rbindlist(list(TERM, tm))
    # (c) 2% frontier sample of (QMD,TPH)
    s <- d[sample(.N, max(1L, as.integer(0.02 * .N))), .(QMD_IN, TPH)]
    FRONT <- if (is.null(FRONT)) s else rbindlist(list(FRONT, s))
  }, error = function(e) logerr(paste("read", a, f), e))
  # collapse AGG sums across files
  AGG <- AGG[, .(BA = sum(BA), TPH = sum(TPH), QMD = sum(QMD), SDI = sum(SDI),
                 AGB = sum(AGB), nAGB = sum(nAGB), n = sum(n)),
             by = .(site_class, PROJ_YEAR)]
  AGG[, `:=`(meanBA = BA/n, meanTPH = TPH/n, meanQMD = QMD/n, meanSDI = SDI/n,
             meanAGB = fifelse(nAGB>0, AGB/nAGB, NA_real_), arm = a)]
  agg_list[[a]] <- AGG; term_list[[a]] <- TERM[, arm := a]; front_list[[a]] <- FRONT[, arm := a]
  cat(sprintf("arm %s: %d site_class x year cells, %d terminal stands\n", a, nrow(AGG), nrow(TERM)))
  gc()
}
agg  <- rbindlist(agg_list, fill = TRUE); fwrite(agg,  file.path(OUT, "bakuzis_site_ordering.csv"))
term <- rbindlist(term_list, fill = TRUE)
front<- rbindlist(front_list, fill = TRUE)

## ---- TEST 1: Reineke self-thinning frontier slope + realized max density ----
reineke <- rbindlist(lapply(names(arms), function(a) {
  s <- front[arm == a & is.finite(QMD_IN) & QMD_IN > 1 & TPH > 0]
  s[, lq := log(QMD_IN)][, lt := log(TPH)]
  br <- seq(min(s$lq), max(s$lq), length.out = 21)
  s[, bin := cut(lq, br, include.lowest = TRUE)]
  fr <- s[, .(lq = mean(lq), lt95 = quantile(lt, 0.95)), by = bin][is.finite(lq)]
  sl <- tryCatch(coef(lm(lt95 ~ lq, fr))[2], error = function(e) NA_real_)
  tt <- term[arm == a]
  data.table(arm = a, frontier_slope = round(as.numeric(sl), 3),
             SDI_p99 = round(quantile(tt$SDI, 0.99, na.rm = TRUE)),
             SDI_p50 = round(median(tt$SDI, na.rm = TRUE)))
}))
fwrite(reineke, file.path(OUT, "bakuzis_reineke.csv")); print(reineke)

## ---- TEST 2: site ordering (high >= mid >= low in BA and AGB, per year) ----
ord <- dcast(agg, arm + PROJ_YEAR ~ site_class, value.var = "meanBA")
ord[, ba_ordered := (high >= mid) & (mid >= low)]
ordA <- dcast(agg, arm + PROJ_YEAR ~ site_class, value.var = "meanAGB")
ordA[, agb_ordered := (high >= mid) & (mid >= low)]
site_order <- merge(
  ord[, .(BA_inversions = sum(!ba_ordered, na.rm = TRUE), n_year = .N), by = arm],
  ordA[, .(AGB_inversions = sum(!agb_ordered, na.rm = TRUE)), by = arm], by = "arm")
fwrite(site_order, file.path(OUT, "bakuzis_site_order.csv")); print(site_order)

## ---- TEST 3: yield-density (AGB rises with SDI across site classes at mid-horizon) ----
mid_year <- agg[, .(yy = PROJ_YEAR[which.min(abs(PROJ_YEAR - 50))][1])]
yd <- agg[PROJ_YEAR == 50, .(arm, site_class, meanSDI, meanAGB)]
yd_test <- yd[, .(yd_monotone = cor(meanSDI, meanAGB, use = "complete.obs") > 0), by = arm]
fwrite(yd, file.path(OUT, "bakuzis_yield_density.csv"))

## ---- PASS/FLAG table ----
pf <- merge(reineke, site_order, by = "arm")
pf <- merge(pf, yd_test, by = "arm")
pf[, reineke_flag := fifelse(abs(frontier_slope + 1.605) <= 0.4, "PASS", "FLAG")]
pf[, siteorder_flag := fifelse(BA_inversions == 0, "PASS", "FLAG")]
pf[, yd_flag := fifelse(yd_monotone, "PASS", "FLAG")]
pf[, eichhorn_flag := "UNTESTABLE (no HT emitted)"]
fwrite(pf, file.path(OUT, "bakuzis_passflag.csv")); print(pf)

## ---- FIGURE: compact Bakuzis matrix (4 panels) colored by site class, lty by mortality ----
cols <- c(low = "#d95f02", mid = "#7570b3", high = "#1b9e77")
arm_lty <- c(b1_native = 2, b2_native = 2, b1_gompit = 1, b2_gompit = 1)
arm_col2<- c(b1_native = "#d62728", b2_native = "#e6ab02", b1_gompit = "#7b2d8e", b2_gompit = "#1f77b4")
tryCatch({
  png(file.path(OUT, "fvs_bakuzis_matrix.png"), width = 2000, height = 1500, res = 200)
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
  # P1 BA vs year by site class, gompit arms solid
  plot(NA, xlim = range(agg$PROJ_YEAR), ylim = c(0, max(agg$meanBA, na.rm=TRUE)),
       xlab="projection year", ylab="mean BA (ft2/ac)", main="BA vs age by site class")
  for (a in names(arms)) for (sc in names(cols)) { z <- agg[arm==a & site_class==sc][order(PROJ_YEAR)]
    lines(z$PROJ_YEAR, z$meanBA, col = cols[sc], lwd=1.6, lty = arm_lty[a]) }
  legend("topleft", c(names(cols),"native (dash)","gompit (solid)"),
         col=c(cols,"grey40","grey40"), lty=c(1,1,1,2,1), lwd=1.6, bty="n", cex=0.8)
  # P2 self-thinning frontier TPH vs QMD (log-log) per arm
  plot(NA, xlim=log(c(2,30)), ylim=log(c(20,2000)), xlab="log QMD (in)", ylab="log TPH",
       main="Reineke self-thinning frontier")
  for (a in names(arms)) { s <- front[arm==a & QMD_IN>1 & TPH>0]
    br <- seq(log(2), log(30), length.out=21); s[, b:=cut(log(QMD_IN),br)]
    fr <- s[, .(lq=mean(log(QMD_IN)), lt=quantile(log(TPH),0.95)), by=b][is.finite(lq)][order(lq)]
    lines(fr$lq, fr$lt, col=arm_col2[a], lwd=2, lty=arm_lty[a]) }
  abline(a = log(2000), b = -1.605, col="grey50", lty=3)  # reference -1.605
  legend("bottomleft", c("b1 nat","b2 nat","b1 gom","b2 gom","slope -1.605"),
         col=c(arm_col2,"grey50"), lty=c(2,2,1,1,3), lwd=2, bty="n", cex=0.8)
  # P3 mean SDI vs year (over-densification)
  plot(NA, xlim=range(agg$PROJ_YEAR), ylim=c(0, max(agg$meanSDI, na.rm=TRUE)),
       xlab="projection year", ylab="mean SDI", main="Stand density index vs age")
  for (a in names(arms)) { z <- agg[arm==a, .(s=mean(meanSDI,na.rm=TRUE)), by=PROJ_YEAR][order(PROJ_YEAR)]
    lines(z$PROJ_YEAR, z$s, col=arm_col2[a], lwd=2, lty=arm_lty[a]) }
  legend("topleft", names(arms), col=arm_col2, lty=arm_lty, lwd=2, bty="n", cex=0.8)
  # P4 yield-density AGB vs SDI at year 50
  yd50 <- agg[PROJ_YEAR==50]
  plot(yd50$meanSDI, yd50$meanAGB, col=arm_col2[yd50$arm], pch=19,
       xlab="mean SDI (year 50)", ylab="mean AGB (tons/ac)", main="Yield-density at year 50")
  legend("topleft", names(arms), col=arm_col2, pch=19, bty="n", cex=0.8)
  dev.off()
}, error = function(e) logerr("figure", e))
cat("bakuzis_assess done\n")
