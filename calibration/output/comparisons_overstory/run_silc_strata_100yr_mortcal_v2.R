#!/usr/bin/env Rscript
# run_silc_strata_100yr_mortcal_v2.R
# =====================================================================
# 100-year AcadianGY MORTCAL=TRUE projection on the SILC byStrata
# stands. Uses the GrownDB year-2023 snapshot as the starting tree
# list (already correctly scaled by FVS/AGM upstream) rather than the
# raw TreeInit (whose Tree_Count semantics are variable radius and
# don't reduce to a simple per-tree EXPF). This guarantees the
# starting state matches the existing default trajectory exactly.
# =====================================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
source("~/AcadianGY_12.3.9.r")

ACRES_PER_HA <- 2.4710538147

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
                     "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM",
                     "RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP",
                     "WS","YB")

prep_stand_from_growndb <- function(td_stand, stand_id, inv_year) {
  sp <- as.character(td_stand$Species)
  sp[!(sp %in% acadgy_species)] <- "OH"
  DBH_cm <- td_stand$DBH * 2.54
  HT_m   <- td_stand$Ht * 0.3048
  HT_m[is.na(HT_m) | HT_m == 0] <-
    pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[is.na(HT_m) | HT_m == 0])))
  data.frame(
    STAND = stand_id, PLOT = 1L,
    TREE = seq_len(nrow(td_stand)),
    SP = sp, DBH = DBH_cm, HT = HT_m, HCB = NA_real_,
    EXPF = td_stand$TPA * ACRES_PER_HA,  # GrownDB TPA is per-acre, AGM wants per-ha
    YEAR = inv_year,
    dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
    max.dbh = 200, max.height = 50,
    Form = NA, Risk = NA, stringsAsFactors = FALSE
  )
}

summarize_stand <- function(tree, stand_id, year) {
  live <- tree[!is.na(tree$DBH) & tree$DBH >= 11.43, ]  # >= 4.5"
  DBH_in <- live$DBH / 2.54
  HT_ft  <- live$HT / 0.3048
  EXPF_ac <- live$EXPF / ACRES_PER_HA
  BA  <- sum(0.005454 * DBH_in^2 * EXPF_ac)
  TPA <- sum(EXPF_ac)
  QMD <- if (TPA > 0) sqrt((BA / TPA) / 0.005454) else NA_real_
  tcuft <- 0.0025 * DBH_in^2 * HT_ft
  Cords <- sum(tcuft * 0.90 * EXPF_ac, na.rm=TRUE) / 79
  ok9 <- DBH_in >= 9.0
  BdFt <- sum(0.01 * DBH_in[ok9]^2 * HT_ft[ok9] * EXPF_ac[ok9], na.rm=TRUE)
  data.frame(StandID = stand_id, Year = year,
             TPA = TPA, BA = BA, QMD = QMD,
             Cords = Cords, NetCords = Cords, BdFt = BdFt)
}

WD <- getwd()
si <- read.csv(file.path(WD, "Acadian_Matrix_StandInit_2023.csv"), stringsAsFactors = FALSE)
gr <- read.csv(file.path(WD, "GrownDB_byStrata_ALL.csv"), stringsAsFactors = FALSE)
gr2023 <- gr[gr$Year == 2023 & is.finite(gr$DBH) & is.finite(gr$TPA) & gr$TPA > 0, ]
stand_ids <- sort(unique(gr2023$StandID))
cat(sprintf("Running 100-yr MORTCAL=TRUE on %d stands from GrownDB 2023 snapshot\n", length(stand_ids)))

out_rows <- list()
START_YR <- 2023
N_YR <- 100
for (sid in stand_ids) {
  cat(sprintf("\n--- Stand %s ---\n", sid))
  td <- gr2023[gr2023$StandID == sid, ]
  if (nrow(td) == 0) next
  cur <- prep_stand_from_growndb(td, sid, START_YR)
  ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 5.0,
              MORTCAL = TRUE, MORTCAL_INTERVAL = 5)
  stand <- list(
    CSI = as.numeric(si$ClimateSiteIndexMeters[si$STAND_ID == sid][1]),
    ELEV = as.numeric(si$ElevationMeters[si$STAND_ID == sid][1])
  )
  out_rows[[length(out_rows)+1]] <- summarize_stand(cur, sid, START_YR)
  ok <- TRUE
  for (yr_off in seq_len(N_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops),
                    error = function(e) { message("ERR yr ", yr_off, ": ", e$message); NULL })
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt
    cur$YEAR <- START_YR + yr_off
    if (yr_off %% 5 == 0) {
      out_rows[[length(out_rows)+1]] <- summarize_stand(cur, sid, START_YR + yr_off)
      last <- tail(out_rows, 1)[[1]]
      cat(sprintf("  yr+%d: BA=%.1f TPA=%.0f QMD=%.1f Cords=%.1f\n",
                  yr_off, last$BA, last$TPA, last$QMD, last$Cords))
    }
  }
  if (!ok) cat(sprintf("Stand %s failed at yr %d (partial data kept)\n", sid, yr_off))
}
out <- do.call(rbind, out_rows)
write.csv(out, file.path(WD, "silc_strata_100yr_mortcal_trajectories.csv"), row.names=FALSE)
cat(sprintf("\nDone. Wrote %d rows for %d stands.\n", nrow(out), length(unique(out$StandID))))
