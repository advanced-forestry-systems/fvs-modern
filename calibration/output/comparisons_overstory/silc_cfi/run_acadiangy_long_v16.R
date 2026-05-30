#!/usr/bin/env Rscript
# run_acadiangy_long_v16.R - AcadianGY 12.3.9 MORTCAL=TRUE on the
# v16 long horizon CFI pair set (n=73 plots).
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

# Species lookup by common name (since v16 TREE has COMMON_NAME not SP code)
COMMON_TO_AGY <- c("Balsam fir"="BF","Red spruce"="RS","Paper birch"="PB","Yellow birch"="YB",
                   "Red maple"="RM","Sugar maple"="SM","White ash"="WA","Northern red oak"="RO",
                   "Black spruce"="BS","White spruce"="WS","Jack pine"="JP","Red pine"="RP",
                   "White pine"="WP","Eastern hemlock"="EH","Northern white cedar"="WC",
                   "Black cherry"="BC","Quaking aspen"="QA","Eastern white pine"="WP",
                   "Gray birch"="GB","Bigtooth aspen"="BA","American beech"="AB",
                   "Striped maple"="ST","Tamarack"="TA")

prep_tree <- function(td, year_prev) {
  # Try SP code first, then COMMON_NAME, then SPCD numeric
  sp <- rep(NA_character_, nrow(td))
  if ("SP" %in% names(td)) {
    sp <- ifelse(td$SP %in% names(cfi_to_agy), cfi_to_agy[td$SP], NA)
  }
  if ("COMMON_NAME" %in% names(td)) {
    miss <- is.na(sp)
    sp[miss] <- COMMON_TO_AGY[td$COMMON_NAME[miss]]
  }
  miss <- is.na(sp) & !is.na(td$SPCD)
  sp[miss] <- SPCD_to_AGY[as.character(td$SPCD[miss])]
  sp[is.na(sp) | !(sp %in% acadgy_species)] <-
    ifelse(td$SPCD[is.na(sp) | !(sp %in% acadgy_species)] >= 300, "OH", "OS")
  DBH_cm <- td$DBH * 2.54
  HT_m <- ifelse(!is.na(td$HT) & td$HT > 0, td$HT * 0.3048, NA_real_)
  miss_h <- is.na(HT_m) | HT_m == 0
  HT_m[miss_h] <- pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[miss_h])))
  data.frame(STAND = sprintf("CFI_%04d", td$PLOT[1]), PLOT = 1L,
             TREE = seq_len(nrow(td)), SP = sp, DBH = DBH_cm, HT = HT_m,
             HCB = NA_real_, EXPF = td$EXPF * ACRES_PER_HA,
             YEAR = year_prev,
             dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
             max.dbh = 200, max.height = 50, Form = NA, Risk = NA,
             stringsAsFactors = FALSE)
}

vol_intl_bdft_imp <- function(d, h, e) {
  ok <- d >= 9.0 & is.finite(d) & is.finite(h)
  sum(0.01 * d[ok]^2 * h[ok] * e[ok], na.rm=TRUE)
}

WD <- getwd()
manifest <- read.csv(file.path(WD, "silc_cfi_longhorizon_pairs_v16.csv"), stringsAsFactors=FALSE)
cat(sprintf("Running AcadianGY MORTCAL=TRUE on %d v16 pairs\n", nrow(manifest)))
out_rows <- list()
fail <- 0
for (i in seq_len(nrow(manifest))) {
  pr <- manifest[i, ]
  tf <- file.path(WD, pr$tree_list_file)
  if (!file.exists(tf)) { fail <- fail + 1; next }
  td <- read.csv(tf, stringsAsFactors=FALSE)
  if (nrow(td) == 0) { fail <- fail + 1; next }
  cur <- prep_tree(td, pr$YEAR_PREV)
  ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=5.0,
              MORTCAL=TRUE, MORTCAL_INTERVAL=5)
  stand <- list(CSI=CSI_M, ELEV=ELEV_M)
  ok <- TRUE
  for (yr in seq_len(pr$PERIOD_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand=stand, ops=ops),
                    error=function(e) NULL)
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt; cur$YEAR <- pr$YEAR_PREV + yr
  }
  if (!ok) { fail <- fail + 1; next }
  live <- cur[!is.na(cur$DBH) & cur$DBH >= 11.43, ]
  if (nrow(live) == 0) { fail <- fail + 1; next }
  DBH_in <- live$DBH / 2.54
  HT_ft  <- live$HT / 0.3048
  EXPF_ac<- live$EXPF / ACRES_PER_HA
  BA  <- sum(0.005454 * DBH_in^2 * EXPF_ac)
  TPA <- sum(EXPF_ac)
  QMD <- if (TPA > 0) sqrt((BA / TPA) / 0.005454) else NA
  tcuft <- 0.0025 * DBH_in^2 * HT_ft
  Cords <- sum(tcuft * 0.90 * EXPF_ac) / 79
  BdFt_intl <- vol_intl_bdft_imp(DBH_in, HT_ft, EXPF_ac)
  out_rows[[i]] <- data.frame(
    PLOT=pr$PLOT, YEAR_PREV=pr$YEAR_PREV, YEAR_CURR=pr$YEAR_CURR,
    PERIOD_YR=pr$PERIOD_YR,
    BA_PRED_ft2ac=BA, TPA_PRED=TPA, QMD_PRED_in=QMD,
    Cords_PRED_ac=Cords, BdFt_intl_PRED_ac=BdFt_intl,
    BA_OBS_PREV=pr$BA_PREV_FT2AC, BA_OBS_CURR=pr$BA_CURR_FT2AC
  )
  if (i %% 10 == 0) cat(sprintf("  pair %d/%d done\n", i, nrow(manifest)))
}
out <- do.call(rbind, out_rows)
write.csv(out, file.path(WD, "silc_cfi_long_agy_mortcal_v16_results.csv"), row.names=FALSE)
cat(sprintf("\nDone. Wrote %d pairs (%d failed)\n", nrow(out), fail))
