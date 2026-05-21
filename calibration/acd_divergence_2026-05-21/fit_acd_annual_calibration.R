## fit_acd_annual_calibration.R
##
## Fit ANNUAL per species diameter increment multipliers for the Acadian model
## against real ACD FIA remeasurement, by iterating the Acadian dDBH_fun annually
## over each tree's remeasurement interval and matching observed periodic growth.
##
## Produces acd_annual_calibration.csv (SP, dDBH.mult, ...), the table the bridge
## patch reads in place of the mis scaled FVS baimult.
##
## Run on Cardinal: module load gcc/12.3.0 R/4.4.0; Rscript fit_acd_annual_calibration.R
suppressWarnings(suppressMessages({
  library(data.table); library(plyr); library(dplyr); library(purrr); library(tibble)
}))

ACD <- Sys.getenv("ACD_PATH",
  "~/upstream_fvs_check/ForestVegetationSimulator-Interface/fvsOL/inst/extdata/AcadianGY.R")
source(ACD)                          # defines dDBH_fun and ddbh.fun.spp
cat("sourced", AcadianVersionTag, "\n")

DG <- Sys.getenv("ACD_DG",
  "/users/PUOM0008/crsfaaron/fvs-conus/data/processed/acd/diameter_growth.csv")
CAP <- as.integer(Sys.getenv("FIT_CAP", "4000"))   # max trees per species (speed)

## SPCD (FIA numeric) to Acadian alpha code crosswalk for the species present
xwalk <- c(
  "12"="BF","91"="NS","94"="WS","95"="BS","97"="RS","71"="TA","105"="JP",
  "125"="RN","130"="SC","129"="WP","241"="WC","261"="EH",
  "315"="ST","316"="RM","317"="SV","318"="SM","319"="MM",
  "371"="YB","375"="PB","379"="GB",
  "531"="AB","541"="WA","543"="BA","544"="GA",
  "601"="BN","741"="BP","743"="BT","746"="QA","762"="BC",
  "802"="WO","833"="RO","837"="BO","951"="BW","972"="AE")

## climate site index from FVS site index (feet), same map the bridge uses
si_to_csi <- function(si_ft) {
  csi <- approxfun(c(0,8,14,20), c(0,8,12,14), rule = 2)(si_ft * 0.3048)
  ifelse(is.na(csi) | csi <= 0, 12, csi)
}

cat("reading", DG, "\n")
d <- fread(DG, select = c("DIA_t1","DIA_t2","CR_t1","BAL","SI",
                          "years_interval","SPCD"))
setnames(d, c("DIA_t1","DIA_t2","CR_t1","BAL","SI","years_interval","SPCD"),
            c("d1","d2","cr","bal","si","yi","spcd"))

## clean and derive
d <- d[is.finite(d1) & is.finite(d2) & d1 > 0 & yi >= 1 & yi <= 15]
d[, alpha := xwalk[as.character(spcd)]]
d <- d[!is.na(alpha)]
d[, `:=`(
  DBHcm   = d1 * 2.54,
  obs_cm  = (d2 - d1) * 2.54,                 # observed periodic diameter growth, cm
  CRp     = pmin(0.95, pmax(0.05, cr / 100)), # CR percent to proportion
  BALsw   = pmax(bal, 0),                     # use total BAL as larger softwood term
  BALhw   = 0,                                # split unavailable; b5 term small
  CSI     = si_to_csi(si),
  yi      = as.integer(round(yi)))]
d <- d[obs_cm > -1]                            # drop large negatives (measurement error)
cat("usable records:", nrow(d), " species:", length(unique(d$alpha)), "\n")

re <- as.data.frame(ddbh.fun.spp)              # species random effects table
get_re <- function(a) {
  r <- re[re$Spp == a, ]
  if (nrow(r) == 0) return(c(0,0,0,0))
  as.numeric(r[1, c("ddbh.b0.spp","ddbh.b2.spp","ddbh.b3.spp","ddbh.b4.spp")])
}

## predicted periodic growth (cm) for a species subset, given annual factor c
pred_periodic <- function(df, b, c) {
  DBH   <- df$DBHcm
  start <- DBH
  maxy  <- max(df$yi)
  for (yr in seq_len(maxy)) {
    inc <- c * dDBH_fun("X", pmax(DBH,1), df$CRp, df$BALsw, df$BALhw, df$CSI,
                        b[1], b[2], b[3], b[4])
    inc[!is.finite(inc) | inc < 0] <- 0
    act <- df$yi >= yr
    DBH[act] <- DBH[act] + inc[act]
  }
  DBH - start
}

set.seed(1)
species <- names(sort(table(d$alpha), decreasing = TRUE))
out <- vector("list", length(species))
for (i in seq_along(species)) {
  a   <- species[i]
  sub <- d[alpha == a]
  if (nrow(sub) > CAP) sub <- sub[sample(.N, CAP)]
  b   <- get_re(a)
  obs <- sub$obs_cm
  obj <- function(c) { p <- pred_periodic(sub, b, c); sum((p - obs)^2) }
  fit <- optimize(obj, c(0.2, 4.0))
  cfac <- fit$minimum
  p1 <- pred_periodic(sub, b, 1)               # default (no calibration)
  pc <- pred_periodic(sub, b, cfac)            # calibrated
  out[[i]] <- data.frame(
    SP          = a,
    n           = nrow(sub),
    dDBH.mult   = round(cfac, 4),
    mean_obs_cm = round(mean(obs), 3),
    mean_def_cm = round(mean(p1), 3),
    rmse_def_cm = round(sqrt(mean((p1 - obs)^2)), 3),
    rmse_cal_cm = round(sqrt(mean((pc - obs)^2)), 3),
    stringsAsFactors = FALSE)
  cat(sprintf("  %-3s n=%6d  mult=%.3f  obs=%.2f def=%.2f  rmse %.2f -> %.2f\n",
              a, nrow(sub), cfac, mean(obs), mean(p1),
              out[[i]]$rmse_def_cm, out[[i]]$rmse_cal_cm))
}
tab <- do.call(rbind, out)
tab <- tab[order(-tab$n), ]

## fill the remaining Acadian codes the bridge may ask for with neutral 1.0
all_codes <- sort(unique(as.character(ddbh.fun.spp$Spp)))
miss <- setdiff(all_codes, tab$SP)
if (length(miss)) {
  tab <- rbind(tab, data.frame(SP = miss, n = 0, dDBH.mult = 1,
    mean_obs_cm = NA, mean_def_cm = NA, rmse_def_cm = NA, rmse_cal_cm = NA))
}
tab$dHt.mult  <- 1     # height/mortality fit separately; neutral for now
tab$mort.mult <- 1

outfile <- Sys.getenv("OUT_CSV", "acd_annual_calibration.csv")
write.csv(tab[, c("SP","dDBH.mult","dHt.mult","mort.mult",
                  "n","mean_obs_cm","mean_def_cm","rmse_def_cm","rmse_cal_cm")],
          outfile, row.names = FALSE)

fit_rows <- tab[tab$n > 0, ]
cat(sprintf("\nwrote %s : %d fitted species, %d total rows\n",
            outfile, nrow(fit_rows), nrow(tab)))
cat(sprintf("mean RMSE default = %.3f cm, calibrated = %.3f cm (%.1f%% reduction)\n",
            mean(fit_rows$rmse_def_cm), mean(fit_rows$rmse_cal_cm),
            100*(mean(fit_rows$rmse_def_cm)-mean(fit_rows$rmse_cal_cm))/mean(fit_rows$rmse_def_cm)))
cat("multiplier range:", round(range(fit_rows$dDBH.mult),3), "\n")
cat("DONE\n")
