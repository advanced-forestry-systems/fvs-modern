#!/usr/bin/env Rscript
# run_acadiangy_on_cfi_v2.R
# =====================================================================
# Cardinal driver for the SILC CFI four-model benchmark.
# Calls AcadianGY_12.3.9.r (current in-source mortality-corrected,
# ingrowth-fixed) on each of the 24 CFI prediction pairs, then writes
# silc_cfi_acadiangy_pred.csv with per-pair predicted vs observed.
#
# Tree-list contract follows acadgy_fia_verify wrapper:
#   STAND, PLOT, TREE, SP, DBH (cm), HT (m), HCB, EXPF (trees/ha),
#   YEAR, dDBH.mult=1, dHt.mult=1, mort.mult=1, max.dbh=200,
#   max.height=50, Form=NA, Risk=NA
# Projects one annual cycle at a time for PERIOD_YR years per pair.
# =====================================================================
.libPaths(c("~/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))
suppressMessages({ library(dplyr); library(plyr); library(purrr) })

ACRES_PER_HA <- 2.4710538147
FT2_AC_PER_M2_HA <- 4.35
ACADGY_VERSION <- "AcadianGY_12.3.9.r"

source(file.path("~", ACADGY_VERSION))
cat(sprintf("Sourced %s -- AcadianVersionTag=%s\n",
            ACADGY_VERSION,
            if (exists("AcadianVersionTag")) AcadianVersionTag else "unknown"))

WD <- getwd()
manifest <- read.csv(file.path(WD, "silc_cfi_pair_summary.csv"),
                     stringsAsFactors = FALSE)

# CFI site: Davistown, Maine. CSI default 12 m (NE softwood typical),
# ELEV ~ 305 m (1000 ft).
CSI_M  <- 12
ELEV_M <- 305

# Species crosswalk: SILC CFI uses 2-3 char COMMON_NAME or SPCD.
# AcadianGY accepts these 38 species codes:
acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH",
                     "GA","GB","HH","HW","JP","NS","OH","OS","PB","PC",
                     "PR","QA","RB","RM","RN","RO","RP","RS","SB","SM",
                     "ST","SW","TA","WA","WC","WP","WS","YB")

# Map SILC CFI COMMON_NAME or SPCD to AcadianGY species code
# CFI uses 2-char codes (BF, RS, PB, YB, etc.) which already align.
# WC = white cedar (SILC code CE -> AcadianGY WC)
cfi_to_agy <- c(
  "BF"="BF", "RS"="RS", "PB"="PB", "YB"="YB", "RM"="RM",
  "SM"="SM", "WA"="WA", "RO"="RO", "BS"="BS", "WS"="WS",
  "JP"="JP", "RP"="RP", "WP"="WP", "EH"="EH", "HM"="EH",  # HM->EH hemlock
  "CE"="WC", "WC"="WC",                                     # cedar
  "BE"="BC", "QA"="QA", "GB"="GB", "BC"="BC",
  "BA"="BA", "AB"="AB", "ST"="ST", "TA"="TA"
)
SPCD_to_AGY <- c(
  "12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
  "541"="WA","934"="GB"
)

prep_tree_metric <- function(td, year_prev) {
  # td columns from pair_<P>_<Yprev>_tree.csv (imperial input)
  # Build SP token from COMMON_NAME first, else SPCD
  td$SP <- ifelse(td$SP %in% names(cfi_to_agy), cfi_to_agy[td$SP], NA)
  miss <- is.na(td$SP) | nchar(td$SP) == 0
  td$SP[miss] <- SPCD_to_AGY[as.character(td$SPCD[miss])]
  td$SP[is.na(td$SP) | !(td$SP %in% acadgy_species)] <-
    ifelse(td$SPCD[is.na(td$SP) | !(td$SP %in% acadgy_species)] >= 300,
           "OH", "OS")  # hardwoods >= 300, softwoods < 300

  # Convert to metric
  DBH_cm   <- td$DBH * 2.54      # in -> cm
  HT_m     <- ifelse(!is.na(td$HT) & td$HT > 0, td$HT * 0.3048, NA_real_)
  # Fill missing heights with a simple Acadian height model from DBH
  miss_h <- is.na(HT_m) | HT_m == 0
  HT_m[miss_h] <- pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[miss_h])))
  # EXPF: input EXPF=5 trees/ac -> trees/ha
  EXPF_ha  <- td$EXPF * ACRES_PER_HA

  out <- data.frame(
    STAND    = sprintf("CFI_%04d", td$PLOT[1]),
    PLOT     = 1L,
    TREE     = seq_len(nrow(td)),
    SP       = td$SP,
    DBH      = DBH_cm,
    HT       = HT_m,
    HCB      = NA_real_,
    EXPF     = EXPF_ha,
    YEAR     = year_prev,
    dDBH.mult= 1,
    dHt.mult = 1,
    mort.mult= 1,
    max.dbh  = 200,
    max.height = 50,
    Form     = NA,
    Risk     = NA,
    stringsAsFactors = FALSE
  )
  out
}

stand_metrics_metric <- function(df) {
  # given a metric tree list, return BA(m2/ha), TPA(/ha), QMD(cm)
  ba <- sum(0.00007854 * df$DBH^2 * df$EXPF, na.rm = TRUE)
  tp <- sum(df$EXPF, na.rm = TRUE)
  qm <- if (tp > 0) sqrt(sum(df$DBH^2 * df$EXPF, na.rm=TRUE) / tp) else NA
  list(BA_m2ha = ba, TPA_ha = tp, QMD_cm = qm)
}

run_one_pair <- function(p_row) {
  tf <- file.path(WD, p_row$tree_list_file)
  if (!file.exists(tf)) {
    cat(sprintf("  skip: %s missing\n", p_row$tree_list_file))
    return(NULL)
  }
  td <- read.csv(tf, stringsAsFactors = FALSE)
  if (nrow(td) == 0) return(NULL)
  cur <- prep_tree_metric(td, p_row$YEAR_PREV)
  if (nrow(cur) == 0) return(NULL)

  ops <- list(verbose = FALSE,
              INGROWTH = "Y",
              MinDBH = 5.0)        # cm  (4.5 in ~= 11.4 cm; choose 5 for conservative ingrowth)
  stand <- list(CSI = CSI_M, ELEV = ELEV_M)

  ok <- TRUE
  for (yr in seq_len(p_row$PERIOD_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops),
                    error = function(e) {
                      cat(sprintf("  AcadianGY err pair %d %d->yr%d: %s\n",
                                  p_row$PLOT, p_row$YEAR_PREV, yr, e$message))
                      NULL
                    })
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt
    cur$YEAR <- p_row$YEAR_PREV + yr
  }
  if (!ok) return(NULL)

  # Aggregate live trees with DBH >= 4.5 in == 11.43 cm
  live <- cur[!is.na(cur$DBH) & cur$DBH >= 11.43, ]
  sm_m <- stand_metrics_metric(live)
  data.frame(
    PLOT       = p_row$PLOT,
    YEAR_PREV  = p_row$YEAR_PREV,
    YEAR_CURR  = p_row$YEAR_CURR,
    PERIOD_YR  = p_row$PERIOD_YR,
    BA_PRED_ft2ac = sm_m$BA_m2ha * FT2_AC_PER_M2_HA,
    TPA_PRED   = sm_m$TPA_ha / ACRES_PER_HA,
    QMD_PRED_in= sm_m$QMD_cm / 2.54,
    BA_OBS_PREV= p_row$BA_PREV_FT2AC,
    BA_OBS_CURR= p_row$BA_CURR_FT2AC,
    PAI_NET_OBS  = p_row$PAI_BA_NET_FT2ACY,
    PAI_NET_PRED = (sm_m$BA_m2ha * FT2_AC_PER_M2_HA -
                    p_row$BA_PREV_FT2AC) / p_row$PERIOD_YR
  )
}

cat(sprintf("\nRunning AcadianGY on %d CFI pairs ...\n", nrow(manifest)))
res <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  cat(sprintf("  pair %d/%d  plot %d  %d->%d\n", i, nrow(manifest),
              manifest$PLOT[i], manifest$YEAR_PREV[i], manifest$YEAR_CURR[i]))
  res[[i]] <- run_one_pair(manifest[i, ])
}
out <- do.call(rbind, res)
write.csv(out, file.path(WD, "silc_cfi_acadiangy_pred.csv"),
          row.names = FALSE)

cat("\n=== AcadianGY CFI scorecard ===\n")
ok <- !is.na(out$BA_PRED_ft2ac) & !is.na(out$BA_OBS_CURR)
cat(sprintf("  n pairs OK   : %d / %d\n", sum(ok), nrow(out)))
if (sum(ok) > 0) {
  bias_pct <- 100*(mean(out$BA_PRED_ft2ac[ok])/mean(out$BA_OBS_CURR[ok]) - 1)
  rmse     <- sqrt(mean((out$BA_PRED_ft2ac[ok] - out$BA_OBS_CURR[ok])^2))
  cat(sprintf("  BA bias-of-means : %+.2f%%\n", bias_pct))
  cat(sprintf("  BA RMSE          : %.2f ft^2/ac\n", rmse))
  cat(sprintf("  Net PAI bias     : %+.3f vs obs mean %.3f ft^2/ac/yr\n",
              mean(out$PAI_NET_PRED[ok] - out$PAI_NET_OBS[ok]),
              mean(out$PAI_NET_OBS[ok])))
}
