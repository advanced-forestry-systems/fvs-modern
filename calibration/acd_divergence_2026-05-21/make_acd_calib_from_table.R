## make_acd_calib_from_table.R
##
## Drop-in helper for the Acadian bridge. Builds the calib.spp data frame that
## AcadianGYOneStand expects (SP, dDBH.mult, dHt.mult, mort.mult, max.dbh,
## max.height) by reading the fitted Acadian annual calibration table instead of
## the FVS-NE baimult/htgmult/mortmult. Size caps still come from FVS.
##
## How to wire into make_fvs_calib() (AcadianGY.R): replace the block that does
##   calib.fvs = fvsGetSpeciesAttrs(c("baimult","htgmult","mortmult", ...))
##   ... rename(dDBH.mult=baimult, dHt.mult=htgmult, mort.mult=mortmult) ...
## with
##   caps = fvsGetSpeciesAttrs(c("maxdbh","maxht"))          # keep FVS size caps
##   size_caps = data.frame(Species = spcodes_alpha,
##                          max.dbh = caps$maxdbh * INtoCM,
##                          max.height = caps$maxht * FTtoM)
##   calib.fvs = make_acd_calib_from_table(csv, size_caps)
##
## Pure base R so it is portable and testable without the FVS DLL.

make_acd_calib_from_table <- function(csv_path, size_caps) {
  ## size_caps: data.frame with columns Species (Acadian/FVS alpha), max.dbh (cm),
  ##            max.height (m). Provides the size caps and the full species list.
  stopifnot(all(c("Species","max.dbh","max.height") %in% names(size_caps)))
  cal <- utils::read.csv(csv_path, stringsAsFactors = FALSE)
  cal <- cal[, c("SP","dDBH.mult","dHt.mult","mort.mult")]

  out <- merge(size_caps, cal, by.x = "Species", by.y = "SP", all.x = TRUE)
  out$dDBH.mult[is.na(out$dDBH.mult)] <- 1   # species absent from the table -> neutral
  out$dHt.mult[is.na(out$dHt.mult)]   <- 1
  out$mort.mult[is.na(out$mort.mult)] <- 1

  data.frame(SP         = out$Species,
             dDBH.mult  = out$dDBH.mult,
             dHt.mult   = out$dHt.mult,
             mort.mult  = out$mort.mult,
             max.dbh    = out$max.dbh,
             max.height = out$max.height,
             stringsAsFactors = FALSE)
}

## ---------------------------------------------------------------------------
## self test (base R; uses the shipped calibration table next to this file)
## ---------------------------------------------------------------------------
if (!exists("RUN_SELFTEST")) RUN_SELFTEST <- TRUE
if (RUN_SELFTEST) {
  csv <- Sys.getenv("ACD_CAL_CSV", "acd_annual_calibration.csv")
  if (!file.exists(csv)) {
    cat("self test skipped: calibration csv not found at", csv, "\n")
  } else {
    ## mock FVS size caps for a few species plus one not in the table (ZZ)
    size_caps <- data.frame(
      Species    = c("BF","RM","RS","RO","ZZ"),
      max.dbh    = c(115, 161, 113, 221, 100),   # cm
      max.height = c(30,  30,  35,  32,  30),    # m
      stringsAsFactors = FALSE)
    tab <- make_acd_calib_from_table(csv, size_caps)
    print(tab, row.names = FALSE)

    ## checks
    cal <- read.csv(csv, stringsAsFactors = FALSE)
    bf  <- cal$dDBH.mult[cal$SP == "BF"]
    ro  <- cal$dDBH.mult[cal$SP == "RO"]
    stopifnot(nrow(tab) == 5)
    stopifnot(abs(tab$dDBH.mult[tab$SP == "BF"] - bf) < 1e-9)
    stopifnot(abs(tab$dDBH.mult[tab$SP == "RO"] - ro) < 1e-9)
    stopifnot(tab$dDBH.mult[tab$SP == "ZZ"] == 1)        # unknown species neutral
    stopifnot(tab$dHt.mult[tab$SP == "ZZ"] == 1)         # unknown species neutral
    stopifnot(all(is.finite(tab$dHt.mult) & tab$dHt.mult > 0))
    stopifnot(tab$max.dbh[tab$SP == "RO"] == 221)        # FVS cap preserved
    cat("\nPASS: helper maps fitted dDBH.mult and dHt.mult, defaults unknown",
        "species to 1, and preserves FVS size caps.\n")
  }
}
