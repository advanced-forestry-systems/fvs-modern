## Layer 2 attribution: run the REAL AcadianGYOneStand under standalone vs
## bridge equivalent configurations and decompose the R vs FVS-ACD gap.
## Pure R model, no FVS DLL needed. Run on Cardinal with R/4.4.0.

suppressWarnings(suppressMessages({
  library(plyr); library(dplyr); library(purrr); library(tibble)
}))

ACD <- Sys.getenv("ACD_PATH",
  "~/upstream_fvs_check/ForestVegetationSimulator-Interface/fvsOL/inst/extdata/AcadianGY.R")
source(ACD)
cat("sourced", AcadianVersionTag, "\n")

K_BA <- 0.00007854   # m2 per stem for dbh in cm

## ---- representative mixed Acadian stand (spruce fir + northern hardwood) ----
set.seed(7)
n  <- 60
sp <- sample(c("RS","BF","RM","YB","WP","SM"), n, replace = TRUE,
             prob = c(0.28, 0.22, 0.20, 0.12, 0.08, 0.10))
tree0 <- data.frame(
  PLOT = 1L,
  TREE = seq_len(n),
  SP   = sp,
  DBH  = pmax(3, rgamma(n, shape = 4, rate = 0.22)),          # cm
  CR   = pmin(0.9, pmax(0.2, rnorm(n, 0.55, 0.12))),          # proportion
  EXPF = runif(n, 3, 12),                                     # trees per ha
  stringsAsFactors = FALSE)
tree0$HT  <- 1.37 + 0.85 * tree0$DBH^0.78                     # measured-style height (m)
tree0$HCB <- tree0$HT * (1 - tree0$CR)

## non binding size caps so driver F is off (isolate A and C)
tree0$max.dbh    <- 200
tree0$max.height <- 60

## injected FVS calibration multipliers (baimult / htgmult style)
fvs_mult <- c(RS = 1.18, BF = 1.05, RM = 0.92, YB = 0.88, WP = 1.10, SM = 0.95)

base_ops <- list(
  verbose = FALSE, INGROWTH = "N", MinDBH = 3.0, CutPoint = 0.95,
  SBW = NULL, usedHTCap = TRUE,
  rtnVars = c("YEAR","PLOT","TREE","SP","DBH","HT","HCB","CR","EXPF",
              "pHT","pHCB","Form","Risk",
              "dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height"))

## project the real model annually for `horizon` years under a config
project_real <- function(tree0, mult, CSI, horizon, ops = base_ops) {
  tr <- tree0
  tr$dDBH.mult <- as.numeric(mult[tr$SP]); tr$dDBH.mult[is.na(tr$dDBH.mult)] <- 1
  tr$dHt.mult  <- as.numeric(mult[tr$SP]); tr$dHt.mult[is.na(tr$dHt.mult)]  <- 1
  tr$mort.mult <- 1
  stand <- list(CSI = CSI, ELEV = 350)
  for (yr in seq_len(horizon)) {
    tr$YEAR <- yr
    ## reattach calib columns each year in case rtnVars dropped them
    if (is.null(tr$dDBH.mult)) {
      tr$dDBH.mult <- as.numeric(mult[tr$SP]); tr$dDBH.mult[is.na(tr$dDBH.mult)] <- 1
    }
    if (is.null(tr$dHt.mult)) {
      tr$dHt.mult <- as.numeric(mult[tr$SP]); tr$dHt.mult[is.na(tr$dHt.mult)] <- 1
    }
    if (is.null(tr$mort.mult)) tr$mort.mult <- 1
    if (is.null(tr$max.dbh))    tr$max.dbh    <- 200
    if (is.null(tr$max.height)) tr$max.height <- 60
    out <- AcadianGYOneStand(tr, stand = stand, ops = ops)
    tr  <- out
  }
  tr
}

ba_total <- function(tr) sum(tr$EXPF * K_BA * tr$DBH^2)
ba_by_sp <- function(tr) tapply(tr$EXPF * K_BA * tr$DBH^2, tr$SP, sum)

mult1 <- setNames(rep(1, length(unique(sp))), unique(sp))

cfgs <- list(
  "baseline standalone (mult=1, CSI=12)" = list(mult = mult1,    CSI = 12),
  "A: inject FVS multipliers"            = list(mult = fvs_mult, CSI = 12),
  "C: CSI 9 not 12"                      = list(mult = mult1,    CSI = 9),
  "C: CSI 15 not 12"                     = list(mult = mult1,    CSI = 15),
  "A+C: multipliers + CSI 9"             = list(mult = fvs_mult, CSI = 9)
)

cat(sprintf("\nReal AcadianGYOneStand attribution. Stand: %d trees, start BA = %.2f m2/ha\n",
            n, ba_total(tree0)))

for (horizon in c(5, 25)) {
  cat(sprintf("\n--- horizon = %d yr (real model) ---\n", horizon))
  base_ba <- NA
  rows <- lapply(names(cfgs), function(nm) {
    cf  <- cfgs[[nm]]
    end <- project_real(tree0, cf$mult, cf$CSI, horizon)
    ba_total(end)
  })
  names(rows) <- names(cfgs)
  base_ba <- rows[[1]]
  res <- data.frame(
    config = names(rows),
    ba = round(unlist(rows), 2),
    delta = round(unlist(rows) - base_ba, 2),
    pct = round(100 * (unlist(rows) - base_ba) / base_ba, 1),
    row.names = NULL)
  print(res, row.names = FALSE)
}

## per species at 25 yr for the multiplier driver
end_b <- project_real(tree0, mult1,    12, 25)
end_a <- project_real(tree0, fvs_mult, 12, 25)
bsp <- ba_by_sp(end_b); asp <- ba_by_sp(end_a)
sps <- sort(union(names(bsp), names(asp)))
comp <- data.frame(species = sps,
                   BA_standalone = round(as.numeric(bsp[sps]), 2),
                   BA_with_mult  = round(as.numeric(asp[sps]), 2),
                   pct = round(100*(as.numeric(asp[sps]) - as.numeric(bsp[sps])) /
                               as.numeric(bsp[sps]), 1))
cat("\n--- per species BA at 25 yr, multiplier driver (real model) ---\n")
print(comp, row.names = FALSE)
cat(sprintf("stand total moved %.1f%%, species moved up to %.1f%%\n",
            100*(sum(asp)-sum(bsp))/sum(bsp), max(abs(comp$pct), na.rm = TRUE)))
cat("\nLAYER 2 REAL ATTRIBUTION COMPLETE\n")
