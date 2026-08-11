## test_12_3_6.R
## Validate AcadianGY 12.3.6 (in-source opt-in MORTCAL) on Cardinal:
##   A. equivalence: 12.3.6 with MORTCAL off == canonical 12.3.5
##   B. direction:   12.3.6 with MORTCAL=TRUE lowers survivor EXPF (BA/TPH)
## Pure R, no FVS DLL. Run with module load gcc/12.3.0 R/4.4.0.
suppressWarnings(suppressMessages({ library(plyr); library(dplyr); library(purrr); library(tibble) }))

V125 <- Sys.getenv("ACD125", "/users/PUOM0008/crsfaaron/seven_islands/AcadianGY_12.3.5.r")
V126 <- Sys.getenv("ACD126", "/users/PUOM0008/crsfaaron/AcadianGY_12.3.6_uploaded.r")
K_BA <- 0.00007854

make_stand <- function() {
  set.seed(7); n <- 60
  sp <- sample(c("RS","BF","RM","YB","WP","SM"), n, replace=TRUE, prob=c(.28,.22,.20,.12,.08,.10))
  d <- data.frame(PLOT=1L, TREE=seq_len(n), SP=sp,
    DBH=pmax(3, rgamma(n,4,0.22)), CR=pmin(.9,pmax(.2,rnorm(n,.55,.12))),
    EXPF=runif(n,3,12), stringsAsFactors=FALSE)
  d$HT <- 1.37 + 0.85*d$DBH^0.78; d$HCB <- d$HT*(1-d$CR)
  d$max.dbh <- 200; d$max.height <- 60
  d$dDBH.mult <- 1; d$dHt.mult <- 1; d$mort.mult <- 1
  d
}
OPS <- function(extra=list()) modifyList(list(verbose=FALSE, INGROWTH="N", MinDBH=3.0,
  CutPoint=0.95, SBW=NULL,
  rtnVars=c("YEAR","PLOT","TREE","SP","DBH","HT","HCB","CR","EXPF","pHT","pHCB",
            "Form","Risk","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height")), extra)
STAND <- list(CSI=12, ELEV=350)
ba  <- function(tr) sum(tr$EXPF*K_BA*tr$DBH^2)
tph <- function(tr) sum(tr$EXPF)

project <- function(ops, horizon=25) {
  tr <- make_stand()
  for (yr in seq_len(horizon)) {
    tr$YEAR <- yr
    for (cc in c("dDBH.mult","dHt.mult","mort.mult")) if (is.null(tr[[cc]])) tr[[cc]] <- 1
    if (is.null(tr$max.dbh)) tr$max.dbh <- 200
    if (is.null(tr$max.height)) tr$max.height <- 60
    tr <- AcadianGYOneStand(tr, stand=STAND, ops=ops)
  }
  tr
}

## --- canonical 12.3.5 ---
source(V125); cat("sourced canonical:", AcadianVersionTag, "from", basename(V125), "\n")
can <- project(OPS(), 25); ba125 <- ba(can); tph125 <- tph(can)

## --- 12.3.6 (overwrites functions) ---
source(V126); cat("sourced 12.3.6:", AcadianVersionTag, "from", basename(V126), "\n")
off <- project(OPS(), 25)                        # MORTCAL not set
on_ <- project(OPS(list(MORTCAL=TRUE)), 25)      # MORTCAL on

cat(sprintf("\n12.3.5 canonical      : BA=%.3f TPH=%.2f\n", ba125, tph125))
cat(sprintf("12.3.6 MORTCAL off    : BA=%.3f TPH=%.2f\n", ba(off), tph(off)))
cat(sprintf("12.3.6 MORTCAL=TRUE   : BA=%.3f TPH=%.2f  (BA %+.1f%% vs off)\n",
            ba(on_), tph(on_), 100*(ba(on_)-ba(off))/ba(off)))

equiv <- isTRUE(all.equal(ba125, ba(off))) && isTRUE(all.equal(tph125, tph(off)))
cat("\nA. equivalence (12.3.6 MORTCAL off == canonical 12.3.5):", equiv, "\n")
cat("B. direction  (12.3.6 MORTCAL=TRUE lowers BA & TPH):",
    (ba(on_) < ba(off)) && (tph(on_) < tph(off)), "\n")
stopifnot(equiv, ba(on_) < ba(off), tph(on_) < tph(off))
cat("\nALL 12.3.6 TESTS PASSED\n")
