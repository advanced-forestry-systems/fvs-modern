#!/usr/bin/env Rscript
# run_acadiangy_treelist.R
# =====================================================================
# Lightweight driver: write predicted year_curr tree list per CFI pair
# from AcadianGY 12.3.9 (in metric units, with species). Output is
# silc_cfi_acadiangy_treelist.csv (long form: PLOT,YEAR_PREV,YEAR_CURR,
# tree_id,SP,DBH_cm,HT_m,EXPF_ha).
# =====================================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })

source("~/AcadianGY_12.3.9.r")
ACRES_PER_HA <- 2.4710538147
CSI_M  <- 12
ELEV_M <- 305

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH",
                     "GA","GB","HH","HW","JP","NS","OH","OS","PB","PC",
                     "PR","QA","RB","RM","RN","RO","RP","RS","SB","SM",
                     "ST","SW","TA","WA","WC","WP","WS","YB")
cfi_to_agy <- c("BF"="BF","RS"="RS","PB"="PB","YB"="YB","RM"="RM",
                 "SM"="SM","WA"="WA","RO"="RO","BS"="BS","WS"="WS",
                 "JP"="JP","RP"="RP","WP"="WP","EH"="EH","HM"="EH",
                 "CE"="WC","WC"="WC","BE"="BC","QA"="QA","GB"="GB",
                 "BC"="BC","BA"="BA","AB"="AB","ST"="ST","TA"="TA")
SPCD_to_AGY <- c("12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
                  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
                  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
                  "541"="WA","934"="GB")

WD <- getwd()
manifest <- read.csv(file.path(WD, "silc_cfi_pair_summary.csv"),
                     stringsAsFactors = FALSE)

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

out_rows <- list()
for (i in seq_len(nrow(manifest))) {
  pr <- manifest[i, ]
  tf <- file.path(WD, pr$tree_list_file)
  if (!file.exists(tf)) next
  td <- read.csv(tf, stringsAsFactors = FALSE)
  if (nrow(td) == 0) next
  cur <- prep_tree(td, pr$YEAR_PREV)
  ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 5.0)
  stand <- list(CSI = CSI_M, ELEV = ELEV_M)
  ok <- TRUE
  for (yr in seq_len(pr$PERIOD_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops),
                    error = function(e) NULL)
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt; cur$YEAR <- pr$YEAR_PREV + yr
  }
  if (!ok) next
  # Filter to live trees with DBH >= 4.5 in == 11.43 cm
  live <- cur[!is.na(cur$DBH) & cur$DBH >= 11.43, ]
  if (nrow(live) == 0) next
  out_rows[[i]] <- data.frame(
    PLOT = pr$PLOT, YEAR_PREV = pr$YEAR_PREV, YEAR_CURR = pr$YEAR_CURR,
    tree_idx = seq_len(nrow(live)),
    SP = live$SP, DBH_cm = live$DBH, HT_m = live$HT, EXPF_ha = live$EXPF
  )
  cat(sprintf("  pair %d/%d  plot %d  %d trees written\n",
              i, nrow(manifest), pr$PLOT, nrow(live)))
}
out <- do.call(rbind, out_rows)
write.csv(out, file.path(WD, "silc_cfi_acadiangy_treelist.csv"),
          row.names = FALSE)
cat(sprintf("\nWrote silc_cfi_acadiangy_treelist.csv: %d tree rows across %d pairs\n",
            nrow(out), length(unique(paste(out$PLOT, out$YEAR_PREV)))))
