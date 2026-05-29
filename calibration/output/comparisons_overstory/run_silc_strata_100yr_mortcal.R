#!/usr/bin/env Rscript
# run_silc_strata_100yr_mortcal.R
# =====================================================================
# 100-year AcadianGY (FVS-ACD) projection on the 11 SILC byStrata
# stands with the #126b in-source MORTCAL=TRUE size-dependent mortality
# correction enabled. Outputs annual stand-level trajectory rolled up
# to 5-yr cycles for comparison with the existing default GrownDB.
# =====================================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
source("~/AcadianGY_12.3.9.r")

ACRES_PER_HA  <- 2.4710538147
FT2_AC_PER_M2_HA <- 4.35

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
                     "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM",
                     "RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP",
                     "WS","YB")
silc_to_agy <- c("AB"="AB","AS"="OH","BA"="BA","BC"="BC","BF"="BF","BP"="BP",
                 "BS"="BS","BT"="BT","EC"="EC","EH"="EH","GA"="GA","GB"="GB",
                 "HH"="HH","HW"="HW","JP"="JP","NS"="NS","NC"="WC","PB"="PB",
                 "PC"="PC","PR"="PR","QA"="QA","RB"="RB","RM"="RM","RN"="RN",
                 "RO"="RO","RP"="RP","RS"="RS","SB"="SB","SM"="SM","ST"="ST",
                 "SW"="SW","TA"="TA","WA"="WA","WC"="WC","WP"="WP","WS"="WS",
                 "YB"="YB")

prep_stand_trees <- function(td_stand, stand_id, inv_year) {
  sp <- silc_to_agy[td_stand$Species]
  sp[is.na(sp) | !(sp %in% acadgy_species)] <- "OH"
  # Diameter in TreeInit is cm, HT in m, Tree_Count is EXPF per hectare
  DBH_cm <- td_stand$Diameter
  HT_m   <- td_stand$HT
  HT_m[is.na(HT_m) | HT_m == 0] <-
    pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[is.na(HT_m) | HT_m == 0])))
  data.frame(
    STAND = stand_id, PLOT = 1L,
    TREE = seq_len(nrow(td_stand)),
    SP = sp, DBH = DBH_cm, HT = HT_m, HCB = NA_real_,
    EXPF = td_stand$Tree_Count,  # already per-ha
    YEAR = inv_year,
    dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
    max.dbh = 200, max.height = 50,
    Form = NA, Risk = NA, stringsAsFactors = FALSE
  )
}

# Stand-level summary at a given year
summarize_stand <- function(tree, stand_id, year) {
  live <- tree[!is.na(tree$DBH) & tree$DBH >= 11.43, ]  # >= 4.5"
  DBH_in <- live$DBH / 2.54
  HT_ft  <- live$HT / 0.3048
  EXPF_ac <- live$EXPF / ACRES_PER_HA
  BA  <- sum(0.005454 * DBH_in^2 * EXPF_ac)
  TPA <- sum(EXPF_ac)
  QMD <- if (TPA > 0) sqrt((BA / TPA) / 0.005454) else NA_real_
  # Merch volume (4.5" DBH, 1ft stump to 4" top, simple Honer)
  tcuft <- 0.0025 * DBH_in^2 * HT_ft
  Cords <- sum(tcuft * 0.90 * EXPF_ac, na.rm=TRUE) / 79
  # NetCords excluding cull (none modeled here = same)
  NetCords <- Cords
  # BdFt Intl 1/4 (>= 9")
  ok9 <- DBH_in >= 9.0
  BdFt <- sum(0.01 * DBH_in[ok9]^2 * HT_ft[ok9] * EXPF_ac[ok9], na.rm=TRUE)
  data.frame(StandID = stand_id, Year = year,
             TPA = TPA, BA = BA, QMD = QMD,
             Cords = Cords, NetCords = NetCords, BdFt = BdFt)
}

WD <- getwd()
si <- read.csv(file.path(WD, "Acadian_Matrix_StandInit_2023.csv"), stringsAsFactors = FALSE)
ti <- read.csv(file.path(WD, "Acadian_Matrix_TreeInit_2023.csv"), stringsAsFactors = FALSE)
stand_ids <- sort(unique(ti$STAND_ID))
cat(sprintf("Running 100-yr MORTCAL=TRUE on %d stands\n", length(stand_ids)))

out_rows <- list()
START_YR <- 2023
N_YR <- 100  # project 2023 -> 2123
for (sid in stand_ids) {
  cat(sprintf("\n--- Stand %s ---\n", sid))
  td <- ti[ti$STAND_ID == sid, ]
  if (nrow(td) == 0) next
  cur <- prep_stand_trees(td, sid, START_YR)
  ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 5.0,
              MORTCAL = TRUE, MORTCAL_INTERVAL = 5)
  stand <- list(CSI = as.numeric(si$ClimateSiteIndexMeters[si$STAND_ID == sid][1]),
                ELEV = as.numeric(si$ElevationMeters[si$STAND_ID == sid][1]))
  # year 0 snapshot
  out_rows[[length(out_rows)+1]] <- summarize_stand(cur, sid, START_YR)
  ok <- TRUE
  for (yr_off in seq_len(N_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops),
                    error = function(e) { message("ERR yr ", yr_off, ": ", e$message); NULL })
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt
    cur$YEAR <- START_YR + yr_off
    # snapshot at every 5-yr boundary
    if (yr_off %% 5 == 0) {
      out_rows[[length(out_rows)+1]] <- summarize_stand(cur, sid, START_YR + yr_off)
      cat(sprintf("  yr+%d (%d): BA=%.1f TPA=%.0f QMD=%.1f Cords=%.1f\n",
                  yr_off, START_YR + yr_off,
                  tail(out_rows, 1)[[1]]$BA, tail(out_rows, 1)[[1]]$TPA,
                  tail(out_rows, 1)[[1]]$QMD, tail(out_rows, 1)[[1]]$Cords))
    }
  }
  if (!ok) cat(sprintf("Stand %s failed at yr %d, partial data captured\n", sid, yr_off))
}
out <- do.call(rbind, out_rows)
write.csv(out, file.path(WD, "silc_strata_100yr_mortcal_trajectories.csv"), row.names=FALSE)
cat(sprintf("\nDone. Wrote %d rows for %d stands.\n", nrow(out), length(unique(out$StandID))))
