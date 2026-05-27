## cardinal_magplot_insource_v17.R
## Cross-validate the IN-SOURCE 12.3.6 MORTCAL on Canadian CFI (MAGPlot NB).
## Completes the 2x2 matrix: FIA (v16) already done; this is MAGPlot in-source.
## Expectation: canonical_off ~ -0% (unbiased, == v15 baseline); insource_on
## over-corrects (negative BA bias) but LESS than the wrapper's -6.56% because the
## in-source correction haircuts survivors only and preserves ingrowth.
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
PROJECT_ROOT <- "/users/PUOM0008/crsfaaron"
IN <- "/users/PUOM0008/crsfaaron/magplot/verify"; OUT_DIR <- IN
n_years <- 10L

source(file.path(PROJECT_ROOT, "AcadianGY_12.3.8.r"))   # in-source, opt-in MORTCAL, bumped tag
cat("[v17-magplot] model:", AcadianVersionTag, "\n")

ti <- read.csv(file.path(IN, "magplot_tree_init.csv"), stringsAsFactors = FALSE)
pr <- read.csv(file.path(IN, "magplot_pairs.csv"), stringsAsFactors = FALSE)
ti$STAND <- as.character(ti$STAND); pr$STAND <- as.character(pr$STAND)
cat(sprintf("[v17-magplot] %d trees, %d MAGPlot NB pairs\n", nrow(ti), nrow(pr)))

tree_init <- data.frame(STAND=ti$STAND, PLOT=1L, TREE=ti$TREE, SP=ti$SP,
  DBH=ti$DBH, HT=ti$HT, HCB=NA_real_, EXPF=ti$EXPF, YEAR=ti$YEAR,
  dDBH.mult=1, dHt.mult=1, mort.mult=1, max.dbh=200, max.height=50,
  Form=NA, Risk=NA, stringsAsFactors=FALSE)
stand_init <- data.frame(STAND=pr$STAND, CSI=14.0,
  ELEV=ifelse(!is.na(pr$ELEV), pr$ELEV, 200), stringsAsFactors=FALSE)

p1y <- function(trees, ops) {
  pc <- list()
  for (sid in unique(trees$STAND)) {
    st <- as.list(subset(stand_init, STAND==sid)); sub <- trees[trees$STAND==sid, ]
    if (nrow(sub)==0) next
    prj <- tryCatch(AcadianGYOneStand(sub, stand=st, ops=ops), error=function(e) NULL)
    if (!is.null(prj)) pc[[sid]] <- prj
  }
  if (length(pc)==0) return(NULL); dplyr::bind_rows(pc)
}
ba_m2 <- function(df) tapply((df$DBH^2)*0.00007854*df$EXPF, df$STAND, sum, na.rm=TRUE)
tph   <- function(df) tapply(df$EXPF, df$STAND, sum, na.rm=TRUE)
qmd_cm<- function(df) sqrt(tapply(df$DBH^2*df$EXPF, df$STAND, sum, na.rm=TRUE)/tapply(df$EXPF, df$STAND, sum, na.rm=TRUE))

run_cfg <- function(mortcal, cutpoint=NULL) {
  ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=3.0)
  if (!is.null(cutpoint)) ops$CutPoint <- cutpoint
  if (mortcal) { ops$MORTCAL <- TRUE; ops$MORTCAL_INTERVAL <- 5 }
  cur <- tree_init; B<-list(); T<-list(); Q<-list()
  B[["0"]]<-ba_m2(cur); T[["0"]]<-tph(cur); Q[["0"]]<-qmd_cm(cur)
  for (yr in 1:n_years) {
    cur <- p1y(cur, ops); if (is.null(cur)) break
    B[[as.character(yr)]]<-ba_m2(cur); T[[as.character(yr)]]<-tph(cur); Q[[as.character(yr)]]<-qmd_cm(cur)
  }
  gm <- function(by,s,i){k<-as.character(min(max(round(i),1),n_years)); v<-by[[k]][as.character(s)]; if(length(v)==0||is.na(v)) NA_real_ else as.numeric(v)}
  data.frame(STAND=pr$STAND,
             BA=mapply(function(s,i) gm(B,s,i), pr$STAND, pr$interval_years),
             TPH=mapply(function(s,i) gm(T,s,i), pr$STAND, pr$interval_years),
             QMD=mapply(function(s,i) gm(Q,s,i), pr$STAND, pr$interval_years),
             BA_obs=pr$BA_t2_obs, TPH_obs=pr$TPH_t2_obs, QMD_obs=pr$QMD_t2_obs)
}
summ <- function(tag, d) {
  d <- d[is.finite(d$BA) & is.finite(d$BA_obs), ]
  bb <- 100*(mean(d$BA)-mean(d$BA_obs))/mean(d$BA_obs)
  r2 <- 1 - sum((d$BA-d$BA_obs)^2)/sum((d$BA_obs-mean(d$BA_obs))^2)
  cat(sprintf("%-14s n=%d BA_bias=%+.2f%% R2=%.3f QMD=%.2f(obs %.2f) TPH=%.0f(obs %.0f)\n",
              tag, nrow(d), bb, r2, mean(d$QMD,na.rm=T), mean(d$QMD_obs,na.rm=T), mean(d$TPH,na.rm=T), mean(d$TPH_obs,na.rm=T)))
  data.frame(config=tag, n=nrow(d), BA=mean(d$BA), BA_obs=mean(d$BA_obs), BA_bias_pct=bb, BA_r2=r2,
             QMD=mean(d$QMD,na.rm=T), QMD_obs=mean(d$QMD_obs,na.rm=T), TPH=mean(d$TPH,na.rm=T), TPH_obs=mean(d$TPH_obs,na.rm=T))
}
rows <- list()
cat("[v17-magplot] canonical_off  (12.3.8, MORTCAL off, default CutPoint)
")
rows[["a"]] <- summ("canonical_off",   run_cfg(FALSE, NULL))
cat("[v17-magplot] insource_on    (12.3.8, MORTCAL on, default CutPoint)
")
rows[["b"]] <- summ("insource_on",     run_cfg(TRUE,  NULL))
cat("[v17-magplot] insource_on_cp0(12.3.8, MORTCAL on, CutPoint=0 EV ingrowth)
")
rows[["c"]] <- summ("insource_on_cp0", run_cfg(TRUE,  0))
cat("[v17-magplot] ingrowth_only  (12.3.8, MORTCAL off, CutPoint=0 EV ingrowth)
")
rows[["d"]] <- summ("ingrowth_only",   run_cfg(FALSE, 0))
res <- dplyr::bind_rows(rows)
write.csv(res, file.path(OUT_DIR, "magplot_insource_v17_results.csv"), row.names=FALSE)
cat("
=== MAGPlot NB: does MORTCAL still over-correct once ingrowth flows? ===
")
print(format(res, digits=4))
cat("\nv15 (wrapper) for reference: baseline -0.04%, mortcorr -6.56%\n")
cat("Done.\n")
