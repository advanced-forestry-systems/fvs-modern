#!/usr/bin/env Rscript
# run_acadiangy_long.R - Cardinal driver for long-horizon CFI pairs
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
source("~/AcadianGY_12.3.9.r")
ACRES_PER_HA <- 2.4710538147
FT2_AC_PER_M2_HA <- 4.35
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
    ifelse(td$SPCD[is.na(td$SP) | !(td$SP %in% acadgy_species)] >= 300,
           "OH", "OS")
  DBH_cm <- td$DBH * 2.54
  HT_m   <- ifelse(!is.na(td$HT) & td$HT > 0, td$HT * 0.3048, NA_real_)
  miss_h <- is.na(HT_m) | HT_m == 0
  HT_m[miss_h] <- pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[miss_h])))
  data.frame(STAND = sprintf("CFI_%04d", td$PLOT[1]), PLOT = 1L,
             TREE = seq_len(nrow(td)), SP = td$SP, DBH = DBH_cm, HT = HT_m,
             HCB = NA_real_, EXPF = td$EXPF * ACRES_PER_HA,
             YEAR = year_prev, dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
             max.dbh = 200, max.height = 50, Form = NA, Risk = NA,
             stringsAsFactors = FALSE)
}

vol_intl_bdft_imp <- function(d, h, e) {
  ok <- d >= 9.0 & is.finite(d) & is.finite(h)
  sum(0.01 * d[ok]^2 * h[ok] * e[ok], na.rm=TRUE)
}

WD <- getwd()
manifest <- read.csv(file.path(WD, "silc_cfi_longhorizon_pairs.csv"),
                     stringsAsFactors=FALSE)
out_rows <- list()
tl_rows <- list()
for (i in seq_len(nrow(manifest))) {
  pr <- manifest[i, ]
  tf <- file.path(WD, pr$tree_list_file)
  if (!file.exists(tf)) next
  td <- read.csv(tf, stringsAsFactors=FALSE)
  if (nrow(td) == 0) next
  cur <- prep_tree(td, pr$YEAR_PREV)
  ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=5.0)
  stand <- list(CSI=CSI_M, ELEV=ELEV_M)
  ok <- TRUE
  for (yr in seq_len(pr$PERIOD_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand=stand, ops=ops),
                    error=function(e) NULL)
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt; cur$YEAR <- pr$YEAR_PREV + yr
  }
  if (!ok) next
  live <- cur[!is.na(cur$DBH) & cur$DBH >= 11.43, ]
  DBH_in <- live$DBH / 2.54
  HT_ft  <- live$HT  / 0.3048
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
  # Tree list for species comp
  tl_rows[[i]] <- data.frame(
    PLOT=pr$PLOT, YEAR_PREV=pr$YEAR_PREV, YEAR_CURR=pr$YEAR_CURR,
    tree_idx=seq_len(nrow(live)), SP=live$SP,
    DBH_cm=live$DBH, HT_m=live$HT, EXPF_ha=live$EXPF
  )
  cat(sprintf("  pair %d/%d  plot %d  %d-yr  pred BA=%.1f  cords=%.1f  BdFt=%.0f\n",
              i, nrow(manifest), pr$PLOT, pr$PERIOD_YR, BA, Cords, BdFt_intl))
}
write.csv(do.call(rbind, out_rows),
          file.path(WD, "silc_cfi_long_agy_results.csv"), row.names=FALSE)
write.csv(do.call(rbind, tl_rows),
          file.path(WD, "silc_cfi_long_agy_treelist.csv"), row.names=FALSE)
cat("\nDone\n")
