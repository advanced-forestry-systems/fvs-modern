## test_acd_comprehensive.R
## One consolidated PASS/FAIL pass over the FVS-ACD 12.3.6 work (R side).
## T1 version tag; T2 equivalence (off == 12.3.5); T3 direction (on lowers BA);
## T6 calibration table + bridge helper. (T7 patch + T8 Fortran are run separately.)
## Run: module load gcc/12.3.0 R/4.4.0; Rscript test_acd_comprehensive.R
suppressWarnings(suppressMessages({ library(plyr); library(dplyr); library(purrr); library(tibble) }))

V125 <- Sys.getenv("ACD125", "/users/PUOM0008/crsfaaron/seven_islands/AcadianGY_12.3.5.r")
V126 <- Sys.getenv("ACD126", "/users/PUOM0008/crsfaaron/AcadianGY_12.3.6.r")
WT   <- Sys.getenv("WTDIR",  "/users/PUOM0008/crsfaaron/fvs-modern-acdbridge/calibration/acd_divergence_2026-05-21")
CSV  <- file.path(WT, "acd_annual_calibration.csv")
HELP <- file.path(WT, "make_acd_calib_from_table.R")
K_BA <- 0.00007854

results <- list()
chk <- function(name, ok) { results[[name]] <<- isTRUE(ok); cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name)) }

make_stand <- function() {
  set.seed(7); n <- 60
  sp <- sample(c("RS","BF","RM","YB","WP","SM"), n, replace=TRUE, prob=c(.28,.22,.20,.12,.08,.10))
  d <- data.frame(PLOT=1L, TREE=seq_len(n), SP=sp, DBH=pmax(3, rgamma(n,4,0.22)),
    CR=pmin(.9,pmax(.2,rnorm(n,.55,.12))), EXPF=runif(n,3,12), stringsAsFactors=FALSE)
  d$HT <- 1.37 + 0.85*d$DBH^0.78; d$HCB <- d$HT*(1-d$CR)
  d$max.dbh <- 200; d$max.height <- 60; d$dDBH.mult <- 1; d$dHt.mult <- 1; d$mort.mult <- 1; d
}
OPS <- function(extra=list()) modifyList(list(verbose=FALSE, INGROWTH="N", MinDBH=3.0, CutPoint=0.95,
  SBW=NULL, rtnVars=c("YEAR","PLOT","TREE","SP","DBH","HT","HCB","CR","EXPF","pHT","pHCB",
  "Form","Risk","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height")), extra)
STAND <- list(CSI=12, ELEV=350)
ba <- function(tr) sum(tr$EXPF*K_BA*tr$DBH^2); tph <- function(tr) sum(tr$EXPF)
project <- function(ops, horizon=25) {
  tr <- make_stand()
  for (yr in seq_len(horizon)) { tr$YEAR <- yr
    for (cc in c("dDBH.mult","dHt.mult","mort.mult")) if (is.null(tr[[cc]])) tr[[cc]] <- 1
    if (is.null(tr$max.dbh)) tr$max.dbh <- 200; if (is.null(tr$max.height)) tr$max.height <- 60
    tr <- AcadianGYOneStand(tr, stand=STAND, ops=ops) }
  tr
}

cat("=== FVS-ACD 12.3.6 comprehensive test ===\n")
source(V125); ba125 <- ba(project(OPS())); tph125 <- tph(project(OPS()))
source(V126)
chk("T1 version tag is AcadianV12.3.6", identical(AcadianVersionTag, "AcadianV12.3.6"))
off <- project(OPS()); on_ <- project(OPS(list(MORTCAL=TRUE)))
chk("T2 equivalence: 12.3.6 MORTCAL off == canonical 12.3.5",
    isTRUE(all.equal(ba125, ba(off))) && isTRUE(all.equal(tph125, tph(off))))
chk("T3 direction: 12.3.6 MORTCAL=TRUE lowers BA & TPH", ba(on_) < ba(off) && tph(on_) < tph(off))

RUN_SELFTEST <- FALSE; source(HELP)
sc <- data.frame(Species=c("BF","RM","RS","RO","ZZ"), max.dbh=c(115,161,113,221,100),
                 max.height=c(30,30,35,32,30), stringsAsFactors=FALSE)
tab <- make_acd_calib_from_table(CSV, sc)
chk("T6 calibration table + helper builds valid calib.spp",
    all(c("SP","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height") %in% names(tab)) &&
    all(is.finite(tab$dDBH.mult)) && tab$dDBH.mult[tab$SP=="ZZ"] == 1 && tab$max.dbh[tab$SP=="RO"] == 221)

cat(sprintf("\nfigures: 12.3.5 BA=%.3f | 12.3.6 off BA=%.3f | on BA=%.3f (%+.1f%%)\n",
            ba125, ba(off), ba(on_), 100*(ba(on_)-ba(off))/ba(off)))
np <- sum(unlist(results)); nt <- length(results)
cat(sprintf("\nR-SIDE: %d/%d PASS\n", np, nt))
if (np == nt) cat("R-SIDE ALL PASS\n") else { cat("R-SIDE FAILURES:\n"); print(names(results)[!unlist(results)]) }
