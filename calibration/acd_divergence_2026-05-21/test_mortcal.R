## test_mortcal.R
## Test the updated FVS-ACD R code: AcadianGY_12.3.5_mortcal.r
## Checks (1) reversibility: MORTCAL!=TRUE is identical to canonical
##        (2) direction: MORTCAL=TRUE lowers EXPF/BA (removes over-projection)
## Pure R, no FVS DLL. Run on Cardinal with R/4.4.0.
suppressWarnings(suppressMessages({ library(plyr); library(dplyr); library(purrr); library(tibble) }))

MC <- Sys.getenv("MORTCAL_PATH", "/users/PUOM0008/crsfaaron/seven_islands/AcadianGY_12.3.5_mortcal.r")
source(MC)                      # sources canonical AcadianGY then defines AcadianGYOneStandMortCal
cat("sourced", AcadianVersionTag, "; MortCal defined:", exists("AcadianGYOneStandMortCal"), "\n")

K_BA <- 0.00007854
set.seed(7); n <- 60
sp <- sample(c("RS","BF","RM","YB","WP","SM"), n, replace=TRUE, prob=c(.28,.22,.20,.12,.08,.10))
tree0 <- data.frame(PLOT=1L, TREE=seq_len(n), SP=sp,
  DBH=pmax(3, rgamma(n,4,0.22)), CR=pmin(.9,pmax(.2,rnorm(n,.55,.12))),
  EXPF=runif(n,3,12), stringsAsFactors=FALSE)
tree0$HT <- 1.37 + 0.85*tree0$DBH^0.78; tree0$HCB <- tree0$HT*(1-tree0$CR)
tree0$max.dbh <- 200; tree0$max.height <- 60
tree0$dDBH.mult <- 1; tree0$dHt.mult <- 1; tree0$mort.mult <- 1

ops <- list(verbose=FALSE, INGROWTH="N", MinDBH=3.0, CutPoint=0.95, SBW=NULL,
  rtnVars=c("YEAR","PLOT","TREE","SP","DBH","HT","HCB","CR","EXPF","pHT","pHCB",
            "Form","Risk","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height"))
stand <- list(CSI=12, ELEV=350)

project <- function(fun, ops, horizon=25) {
  tr <- tree0
  for (yr in seq_len(horizon)) {
    tr$YEAR <- yr
    if (is.null(tr$dDBH.mult)) tr$dDBH.mult <- 1
    if (is.null(tr$dHt.mult))  tr$dHt.mult  <- 1
    if (is.null(tr$mort.mult)) tr$mort.mult <- 1
    if (is.null(tr$max.dbh))    tr$max.dbh    <- 200
    if (is.null(tr$max.height)) tr$max.height <- 60
    tr <- fun(tr, stand=stand, ops=ops)
  }
  tr
}
ba  <- function(tr) sum(tr$EXPF*K_BA*tr$DBH^2)
tph <- function(tr) sum(tr$EXPF)

for (H in c(5, 25)) {
  can <- project(AcadianGYOneStand,        ops,                              H)
  off <- project(AcadianGYOneStandMortCal, modifyList(ops,list(MORTCAL=FALSE)), H)
  on_ <- project(AcadianGYOneStandMortCal, modifyList(ops,list(MORTCAL=TRUE)),  H)
  cat(sprintf("\n--- horizon %d yr ---\n", H))
  cat(sprintf("canonical        : BA=%.2f  TPH=%.1f\n", ba(can), tph(can)))
  cat(sprintf("mortcal MORTCAL=F: BA=%.2f  TPH=%.1f\n", ba(off), tph(off)))
  cat(sprintf("mortcal MORTCAL=T: BA=%.2f  TPH=%.1f  (BA %+.1f%% vs canonical)\n",
              ba(on_), tph(on_), 100*(ba(on_)-ba(can))/ba(can)))
  ## reversibility: canonical == MORTCAL=FALSE
  rev_ok <- isTRUE(all.equal(ba(can), ba(off))) && isTRUE(all.equal(tph(can), tph(off)))
  cat("reversibility (canonical == MORTCAL=FALSE):", rev_ok, "\n")
  ## direction: MORTCAL=TRUE removes trees -> lower TPH and BA
  cat("direction (MORTCAL=TRUE lowers BA & TPH):",
      (ba(on_) < ba(can)) && (tph(on_) < tph(can)), "\n")
  if (H == 25) {
    stopifnot(rev_ok)
    stopifnot(ba(on_) < ba(can), tph(on_) < tph(can))
  }
}
cat("\nALL MORTCAL TESTS PASSED\n")
