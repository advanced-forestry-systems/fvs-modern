#!/usr/bin/env Rscript
# run_cfi_backfill_100yr_mortcal.R
# 100-yr AGM MORTCAL on CFI plots assigned to empty strata cells.
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
source("~/AcadianGY_12.3.9.r")
ACRES_PER_HA <- 2.4710538147
CSI_M <- 12; ELEV_M <- 305

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
                     "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM",
                     "RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP",
                     "WS","YB")
cfi_to_agy <- c("BF"="BF","RS"="RS","PB"="PB","YB"="YB","RM"="RM","SM"="SM",
                 "WA"="WA","RO"="RO","BS"="BS","WS"="WS","JP"="JP","RP"="RP",
                 "WP"="WP","EH"="EH","HM"="EH","CE"="WC","WC"="WC","BE"="BC",
                 "QA"="QA","GB"="GB","BC"="BC","BA"="BA","AB"="AB","ST"="ST",
                 "TA"="TA")
SPCD_to_AGY <- c("12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
                  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
                  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
                  "541"="WA","934"="GB")

prep_tree <- function(td, year_prev) {
  td$SP <- ifelse(td$SP %in% names(cfi_to_agy), cfi_to_agy[td$SP], NA)
  miss <- is.na(td$SP) | nchar(td$SP) == 0
  td$SP[miss] <- SPCD_to_AGY[as.character(td$SPCD[miss])]
  td$SP[is.na(td$SP) | !(td$SP %in% acadgy_species)] <-
    ifelse(td$SPCD[is.na(td$SP) | !(td$SP %in% acadgy_species)] >= 300, "OH", "OS")
  DBH_cm <- td$DBH * 2.54
  HT_m <- ifelse(!is.na(td$HT) & td$HT > 0, td$HT * 0.3048, NA_real_)
  miss_h <- is.na(HT_m) | HT_m == 0
  HT_m[miss_h] <- pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[miss_h])))
  data.frame(STAND = sprintf("CFI_%04d", td$PLOT[1]), PLOT = 1L,
             TREE = seq_len(nrow(td)), SP = td$SP, DBH = DBH_cm, HT = HT_m,
             HCB = NA_real_, EXPF = td$EXPF * ACRES_PER_HA, YEAR = year_prev,
             dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
             max.dbh = 200, max.height = 50, Form = NA, Risk = NA,
             stringsAsFactors = FALSE)
}

summarize_stand <- function(tree, stand_id, year) {
  live <- tree[!is.na(tree$DBH) & tree$DBH >= 11.43, ]
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
  data.frame(StandID = stand_id, Year = year, TPA = TPA, BA = BA, QMD = QMD,
             Cords = Cords, NetCords = Cords, BdFt = BdFt)
}

WD <- getwd()
manifest <- read.csv(file.path(WD, "cfi_backfill_manifest.csv"), stringsAsFactors=FALSE)
START_YR <- 2000; N_YR <- 100
out_rows <- list()
for (i in seq_len(nrow(manifest))) {
  pr <- manifest[i, ]
  tf <- file.path(WD, pr$tree_list_file)
  if (!file.exists(tf)) { cat("missing", tf, "\n"); next }
  td <- read.csv(tf, stringsAsFactors=FALSE)
  if (nrow(td) == 0) next
  sid <- sprintf("CFI_%04d", pr$PLOT)
  cat(sprintf("\n--- %s (%s / %s) ---\n", sid, pr$forest_type, pr$density_class))
  cur <- prep_tree(td, START_YR)
  ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=5.0, MORTCAL=TRUE, MORTCAL_INTERVAL=5)
  stand <- list(CSI=CSI_M, ELEV=ELEV_M)
  out_rows[[length(out_rows)+1]] <- cbind(summarize_stand(cur, sid, START_YR),
                                          forest_type=pr$forest_type, density_class=pr$density_class)
  ok <- TRUE
  for (yr_off in seq_len(N_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand=stand, ops=ops), error=function(e) NULL)
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt; cur$YEAR <- START_YR + yr_off
    if (yr_off %% 5 == 0) {
      r <- summarize_stand(cur, sid, START_YR + yr_off)
      out_rows[[length(out_rows)+1]] <- cbind(r, forest_type=pr$forest_type, density_class=pr$density_class)
      cat(sprintf("  yr+%d: BA=%.1f TPA=%.0f Cords=%.1f\n", yr_off, r$BA, r$TPA, r$Cords))
    }
  }
  if (!ok) cat(sprintf("%s failed at yr %d (partial kept)\n", sid, yr_off))
}
out <- do.call(rbind, out_rows)
write.csv(out, file.path(WD, "silc_cfi_backfill_100yr_mortcal_trajectories.csv"), row.names=FALSE)
cat(sprintf("\nDone. Wrote %d rows for %d stands.\n", nrow(out), length(unique(out$StandID))))
