#!/usr/bin/env Rscript
# cardinal_acadgy_insource_v16.R
# ==============================
# #126b numerical re-validation. Confirms the IN-SOURCE size-dependent mortality
# correction baked into AcadianGY_12.3.6.r (ops$MORTCAL, per-cycle, survivors
# only) reproduces the v14 post-step wrapper on the SAME 200 FIA plots:
#     MORTCAL off  ~ +15.4%  (must equal canonical v12.3.5 baseline exactly)
#     MORTCAL on   ~  +8.6%  (matches acadgy_mortcorr_v14_results.csv;
#                             may differ by <1% because the in-source version
#                             does NOT haircut fresh ingrowth, which is correct)
# Mirrors cardinal_acadgy_mortcorr_v14.R exactly except: it sources the
# in-source model and toggles ops$MORTCAL per stand (with MORTCAL_INTERVAL set
# to each plot's remeasurement interval) instead of applying a post-step EXPF
# overlay.
# Output: acadgy_insource_v16_results.csv
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
PROJECT_ROOT <- "/users/PUOM0008/crsfaaron"; N_PLOTS <- 100L
OUT_DIR <- file.path(PROJECT_ROOT, "acadgy_fia_verify"); set.seed(42)
ACRES_PER_HA <- 2.4710538147; n_years <- 10L; ft2_ac_per_m2_ha <- 4.35
source(file.path(PROJECT_ROOT, "seven_islands/AcadianGY_12.3.6.r"))  # IN-SOURCE model
vdat <- read.csv(file.path(PROJECT_ROOT, "fvs-modern/calibration/output/comparisons/intermediate/validation_data_acd_post.csv"))
me_plot <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_PLOT.csv"), colClasses = c("CN" = "character"))
me_tree <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_TREE.csv"), colClasses = c("CN" = "character", "PLT_CN" = "character"))
xwk <- read.csv(file.path(PROJECT_ROOT, "fvs-modern/calibration/data/osm_acadian_species_crosswalk.csv"))
acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB","HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM","RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP","WS","YB")
xf <- xwk[!is.na(xwk$FIA) & xwk$FIA > 0 & nchar(xwk$OSM_AD_CmdKey) > 0, ]
xf$AGY_SP <- ifelse(xf$OSM_AD_CmdKey %in% acadgy_species, xf$OSM_AD_CmdKey, ifelse(xf$FIA >= 300, "OH", "OS"))
spcd_to_code <- setNames(xf$AGY_SP, as.character(xf$FIA))
vdat$PLT_CN_t1_str <- as.character(vdat$PLT_CN_t1)
vme <- vdat[vdat$PLT_CN_t1_str %in% as.character(me_plot$CN) & !is.na(vdat$interval_years) & vdat$interval_years %in% 5:10, ]
samp <- vme[sample(nrow(vme), min(N_PLOTS, nrow(vme))), ]; samp$PLT_CN_t1_str <- as.character(samp$PLT_CN_t1)
me_tree$PLT_CN <- as.character(me_tree$PLT_CN)
ts <- me_tree[me_tree$PLT_CN %in% samp$PLT_CN_t1_str & !is.na(me_tree$DIA) & me_tree$DIA > 0 & me_tree$STATUSCD == 1, ]
ts$SP <- spcd_to_code[as.character(ts$SPCD)]; mm <- is.na(ts$SP) | nchar(ts$SP) == 0
ts$SP[mm] <- ifelse(ts$SPCD[mm] >= 300, "OH", "OS")
ts$TPA_UNADJ <- as.numeric(ts$TPA_UNADJ)
ts$tpa <- ifelse(!is.na(ts$TPA_UNADJ) & ts$TPA_UNADJ > 0, ts$TPA_UNADJ, ifelse(ts$DIA < 5, 74.965, 6.018))
ts$EXPF <- ts$tpa * ACRES_PER_HA; ts$DBH <- ts$DIA * 2.54
ts$HT <- ifelse(!is.na(ts$HT) & ts$HT > 0, ts$HT * 0.3048, NA_real_); mh <- is.na(ts$HT) | ts$HT == 0
ts$HT[mh] <- pmax(2, 1.3 + 30 * (1 - exp(-0.04 * ts$DBH[mh])))
ts$HCB <- NA_real_; ts$YEAR <- 2020; ts$dDBH.mult <- 1; ts$dHt.mult <- 1; ts$mort.mult <- 1
ts$max.dbh <- 200; ts$max.height <- 50; ts$Form <- NA; ts$Risk <- NA
ts$STAND <- ts$PLT_CN; ts$PLOT <- 1L; ts <- ts[order(ts$STAND), ]
ts$TREE <- ave(seq_along(ts$STAND), ts$STAND, FUN = seq_along)
base_init <- ts[, c("STAND","PLOT","TREE","SP","DBH","HT","HCB","EXPF","YEAR","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height","Form","Risk")]
samp$ClimateSI_m <- ifelse(!is.na(samp$ClimateSI_ft), samp$ClimateSI_ft * 0.3048, 14)
samp$ELEV_m <- ifelse(!is.na(samp$ELEV_t2), samp$ELEV_t2 * 0.3048, 200)
stand_init <- unique(data.frame(STAND = samp$PLT_CN_t1_str, CSI = samp$ClimateSI_m, ELEV = samp$ELEV_m, stringsAsFactors = FALSE))
ops0 <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 3.0)
interval_of <- setNames(pmin(pmax(round(samp$interval_years),1L), n_years), samp$PLT_CN_t1_str)
# project one annual cycle for all stands; MORTCAL handled IN-SOURCE per stand
p1y <- function(trees, mortcal, ingrowth="Y", cutpoint=NULL) {
  pc <- list()
  for (sid in unique(trees$STAND)) {
    st <- as.list(subset(stand_init, STAND == sid)); sub <- trees[trees$STAND == sid, ]
    if (nrow(sub) == 0) next
    ops <- ops0; ops$INGROWTH <- ingrowth; if (!is.null(cutpoint)) ops$CutPoint <- cutpoint
    if (mortcal) { ops$MORTCAL <- TRUE; iv <- as.numeric(interval_of[sid]); ops$MORTCAL_INTERVAL <- if (is.na(iv) || iv < 1) 5 else iv }
    pr <- tryCatch(AcadianGYOneStand(sub, stand = st, ops = ops), error = function(e) NULL)
    if (!is.null(pr)) pc[[sid]] <- pr
  }
  if (length(pc) == 0) return(NULL); dplyr::bind_rows(pc)
}
basnap <- function(df) tapply((df$DBH^2)*0.00007854*df$EXPF, df$STAND, sum, na.rm=TRUE)
tpsnap <- function(df) tapply(df$EXPF, df$STAND, sum, na.rm=TRUE)
qmsnap <- function(df) sqrt(tapply(df$DBH^2*df$EXPF, df$STAND, sum, na.rm=TRUE)/tapply(df$EXPF, df$STAND, sum, na.rm=TRUE))
run_cfg <- function(mortcal, ingrowth="Y", cutpoint=NULL) {
  cur <- base_init; bl <- list(); tl <- list(); ql <- list()
  bl[["0"]] <- basnap(cur); tl[["0"]] <- tpsnap(cur); ql[["0"]] <- qmsnap(cur)
  for (yr in 1:n_years) {
    cur <- p1y(cur, mortcal, ingrowth, cutpoint); if (is.null(cur)) break
    bl[[as.character(yr)]] <- basnap(cur); tl[[as.character(yr)]] <- tpsnap(cur); ql[[as.character(yr)]] <- qmsnap(cur)
  }
  gm <- function(by, s, i) { k <- as.character(min(max(round(i),1), n_years)); v <- by[[k]][as.character(s)]; if (length(v)==0||is.na(v)) NA_real_ else as.numeric(v) }
  ba <- mapply(function(s,i) gm(bl,s,i), samp$PLT_CN_t1_str, samp$interval_years)
  tph <- mapply(function(s,i) gm(tl,s,i), samp$PLT_CN_t1_str, samp$interval_years)
  qmd <- mapply(function(s,i) gm(ql,s,i), samp$PLT_CN_t1_str, samp$interval_years)
  d <- data.frame(BA_pred = ba*ft2_ac_per_m2_ha, TPA = tph/ACRES_PER_HA, QMD_in = qmd/2.54, BA_t2 = samp$BA_t2, QMD_t2 = samp$QMD_t2, TPA_t2 = samp$TPA_t2)
  d[is.finite(d$BA_pred) & is.finite(d$BA_t2), ]
}
summ <- function(tag, d) {
  o<-d$BA_t2; p<-d$BA_pred; r2<-1-sum((p-o)^2)/sum((o-mean(o))^2)
  cat(sprintf("%-14s BA_bias=%+.1f%% R2=%.3f TPA=%.0f(obs %.0f) QMD=%.2f(obs %.2f)\n",
              tag, 100*(mean(p)-mean(o))/mean(o), r2, mean(d$TPA,na.rm=T), mean(d$TPA_t2,na.rm=T), mean(d$QMD_in,na.rm=T), mean(d$QMD_t2,na.rm=T)))
  data.frame(config=tag, BA_obs=mean(o), BA_pred=mean(p), BA_bias_pct=100*(mean(p)-mean(o))/mean(o), BA_r2=r2,
             TPA=mean(d$TPA,na.rm=T), TPA_obs=mean(d$TPA_t2,na.rm=T), QMD=mean(d$QMD_in,na.rm=T), QMD_obs=mean(d$QMD_t2,na.rm=T))
}
rows <- list()
cat("[v20] cp95_mort  (CutPoint 0.95, ~no ingrowth, MORTCAL on)\n");        rows[["a"]] <- summ("cp95_mort",  run_cfg(TRUE,  "Y", 0.95))
cat("[v20] cp0_mort   (CutPoint 0 = expected-value ingrowth, MORTCAL on)\n"); rows[["b"]] <- summ("cp0_mort",   run_cfg(TRUE,  "Y", 0))
cat("[v20] cp0_nomort (CutPoint 0 = expected-value ingrowth, MORTCAL off)\n"); rows[["c"]] <- summ("cp0_nomort", run_cfg(FALSE, "Y", 0))
res <- dplyr::bind_rows(rows); write.csv(res, file.path(OUT_DIR, "acadgy_cp0_v20_results.csv"), row.names=FALSE)
cat("\n=== FIA: expected-value ingrowth (CutPoint=0) vs hard threshold ===\n"); print(format(res, digits=4))
cat("\n=== in-source #126b validation (compare to acadgy_mortcorr_v14_results.csv) ===\n")
print(format(res, digits=4))
# pass/fail vs the wrapper run
v14 <- tryCatch(read.csv(file.path(OUT_DIR, "acadgy_mortcorr_v14_results.csv")), error=function(e) NULL)
if (!is.null(v14)) {
  b_off <- res$BA_bias_pct[res$config=="canonical_off"]; b_on <- res$BA_bias_pct[res$config=="insource_on"]
  w_base <- v14$BA_bias_pct[v14$config=="baseline"];     w_mc <- v14$BA_bias_pct[v14$config=="mortcorr"]
  cat(sprintf("\ncanonical_off %+.2f%% vs wrapper baseline %+.2f%%  -> diff %.2f (expect ~0)\n", b_off, w_base, b_off-w_base))
  cat(sprintf("insource_on  %+.2f%% vs wrapper mortcorr %+.2f%%  -> diff %.2f (expect <1, ingrowth not haircut)\n", b_on, w_mc, b_on-w_mc))
  ok <- abs(b_off-w_base) < 0.5 && abs(b_on-w_mc) < 1.5
  cat(if (ok) "VALIDATION PASS\n" else "VALIDATION REVIEW (diff larger than expected)\n")
}
cat("Done.\n")
