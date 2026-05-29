#!/usr/bin/env Rscript
# run_acadiangy_on_cfi_v3.R
# =====================================================================
# Adds merchantable volume (cords + Scribner BdFt) to v2 driver.
#
# Per-tree volume formulas (uniform across observed and predicted so
# the CFI scorecard is apples-to-apples):
#   tcuft_imperial = 0.0025 * DBH_in^2 * HT_ft           total cuft
#   merch_cuft = tcuft * 0.90 if DBH_in >= 4.5           pulp + saw
#              = 0           otherwise
#   sawtimber_cuft = tcuft * 0.55 if DBH_in >= 9.0       Scribner core
#                  = 0          otherwise
#   cords    = sum(merch_cuft * EXPF) / 79               1 cord = 79 cuft
#   bdft_msf = sum(sawtimber_cuft * EXPF) * 6.0 / 1000   approx Scribner
#
# Outputs:
#   silc_cfi_acadiangy_pred_v3.csv   per-pair predicted BA + cords + BdFt
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

CSI_M  <- 12
ELEV_M <- 305

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH",
                     "GA","GB","HH","HW","JP","NS","OH","OS","PB","PC",
                     "PR","QA","RB","RM","RN","RO","RP","RS","SB","SM",
                     "ST","SW","TA","WA","WC","WP","WS","YB")
cfi_to_agy <- c(
  "BF"="BF", "RS"="RS", "PB"="PB", "YB"="YB", "RM"="RM",
  "SM"="SM", "WA"="WA", "RO"="RO", "BS"="BS", "WS"="WS",
  "JP"="JP", "RP"="RP", "WP"="WP", "EH"="EH", "HM"="EH",
  "CE"="WC", "WC"="WC", "BE"="BC", "QA"="QA", "GB"="GB",
  "BC"="BC", "BA"="BA", "AB"="AB", "ST"="ST", "TA"="TA"
)
SPCD_to_AGY <- c(
  "12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
  "541"="WA","934"="GB"
)

prep_tree_metric <- function(td, year_prev) {
  td$SP <- ifelse(td$SP %in% names(cfi_to_agy), cfi_to_agy[td$SP], NA)
  miss <- is.na(td$SP) | nchar(td$SP) == 0
  td$SP[miss] <- SPCD_to_AGY[as.character(td$SPCD[miss])]
  td$SP[is.na(td$SP) | !(td$SP %in% acadgy_species)] <-
    ifelse(td$SPCD[is.na(td$SP) | !(td$SP %in% acadgy_species)] >= 300,
           "OH", "OS")
  DBH_cm   <- td$DBH * 2.54
  HT_m     <- ifelse(!is.na(td$HT) & td$HT > 0, td$HT * 0.3048, NA_real_)
  miss_h <- is.na(HT_m) | HT_m == 0
  HT_m[miss_h] <- pmax(2, 1.3 + 25 * (1 - exp(-0.04 * DBH_cm[miss_h])))
  EXPF_ha  <- td$EXPF * ACRES_PER_HA
  data.frame(
    STAND = sprintf("CFI_%04d", td$PLOT[1]), PLOT = 1L,
    TREE = seq_len(nrow(td)), SP = td$SP,
    DBH = DBH_cm, HT = HT_m, HCB = NA_real_, EXPF = EXPF_ha,
    YEAR = year_prev, dDBH.mult = 1, dHt.mult = 1, mort.mult = 1,
    max.dbh = 200, max.height = 50, Form = NA, Risk = NA,
    stringsAsFactors = FALSE
  )
}

# Stand-level metrics in imperial units
stand_metrics_imperial <- function(df_metric) {
  # df_metric has DBH (cm), HT (m), EXPF (trees/ha)
  DBH_in <- df_metric$DBH / 2.54
  HT_ft  <- df_metric$HT / 0.3048
  EXPF_ac<- df_metric$EXPF / ACRES_PER_HA
  live <- DBH_in >= 4.5 & !is.na(DBH_in)
  d <- DBH_in[live]; h <- HT_ft[live]; e <- EXPF_ac[live]
  if (length(d) == 0)
    return(list(BA_ft2ac = 0, TPA_ac = 0, QMD_in = NA,
                Cords_ac = 0, BdFt_ac = 0))
  BA  <- sum(0.005454 * d^2 * e)
  TPA <- sum(e)
  QMD <- if (TPA > 0) sqrt((BA / TPA) / 0.005454) else NA
  # Volume per tree
  tcuft <- 0.0025 * d^2 * h
  merch_cuft <- tcuft * 0.90
  sawd <- d >= 9.0
  saw_cuft <- ifelse(sawd, tcuft * 0.55, 0)
  Cords <- sum(merch_cuft * e) / 79
  BdFt  <- sum(saw_cuft  * e) * 6.0   # ~6 BdFt per cuft (Scribner)
  list(BA_ft2ac = BA, TPA_ac = TPA, QMD_in = QMD,
       Cords_ac = Cords, BdFt_ac = BdFt)
}

# Observed volume from the input CSV (year_prev tree list)
observed_metrics <- function(p_row) {
  tf <- file.path(WD, p_row$tree_list_file)
  if (!file.exists(tf)) return(NULL)
  td <- read.csv(tf, stringsAsFactors = FALSE)
  if (nrow(td) == 0) return(NULL)
  # td DBH in inches, HT in ft, EXPF trees/ac per the pair_list emitter
  d <- td$DBH; h <- td$HT; e <- td$EXPF
  # Where HT is missing, impute (use same model as prep_tree_metric)
  miss_h <- is.na(h) | h == 0
  HT_imp_ft <- pmax(6, 4.27 + 82 * (1 - exp(-0.04 * (d * 2.54))))
  h[miss_h] <- HT_imp_ft[miss_h] / 1  # already ft
  live <- d >= 4.5 & !is.na(d)
  d <- d[live]; h <- h[live]; e <- e[live]
  tcuft <- 0.0025 * d^2 * h
  merch_cuft <- tcuft * 0.90
  saw <- d >= 9.0
  saw_cuft <- ifelse(saw, tcuft * 0.55, 0)
  list(
    BA_obs_prev    = sum(0.005454 * d^2 * e),
    Cords_obs_prev = sum(merch_cuft * e) / 79,
    BdFt_obs_prev  = sum(saw_cuft  * e) * 6.0
  )
}

run_one_pair <- function(p_row) {
  tf <- file.path(WD, p_row$tree_list_file)
  if (!file.exists(tf)) return(NULL)
  td <- read.csv(tf, stringsAsFactors = FALSE)
  if (nrow(td) == 0) return(NULL)
  cur <- prep_tree_metric(td, p_row$YEAR_PREV)
  obs_prev <- observed_metrics(p_row)

  ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 5.0)
  stand <- list(CSI = CSI_M, ELEV = ELEV_M)
  ok <- TRUE
  for (yr in seq_len(p_row$PERIOD_YR)) {
    nxt <- tryCatch(AcadianGYOneStand(cur, stand = stand, ops = ops),
                    error = function(e) NULL)
    if (is.null(nxt)) { ok <- FALSE; break }
    cur <- nxt; cur$YEAR <- p_row$YEAR_PREV + yr
  }
  if (!ok) return(NULL)
  pm <- stand_metrics_imperial(cur)
  data.frame(
    PLOT       = p_row$PLOT,
    YEAR_PREV  = p_row$YEAR_PREV,
    YEAR_CURR  = p_row$YEAR_CURR,
    PERIOD_YR  = p_row$PERIOD_YR,
    BA_PRED_ft2ac    = pm$BA_ft2ac,
    TPA_PRED         = pm$TPA_ac,
    QMD_PRED_in      = pm$QMD_in,
    Cords_PRED_ac    = pm$Cords_ac,
    BdFt_PRED_ac     = pm$BdFt_ac,
    BA_OBS_PREV      = p_row$BA_PREV_FT2AC,
    BA_OBS_CURR      = p_row$BA_CURR_FT2AC,
    Cords_OBS_PREV   = if(is.null(obs_prev)) NA else obs_prev$Cords_obs_prev,
    BdFt_OBS_PREV    = if(is.null(obs_prev)) NA else obs_prev$BdFt_obs_prev,
    PAI_NET_OBS      = p_row$PAI_BA_NET_FT2ACY
  )
}

cat(sprintf("\nRunning AcadianGY (v3, merch vol) on %d CFI pairs ...\n",
            nrow(manifest)))
res <- vector("list", nrow(manifest))
for (i in seq_len(nrow(manifest))) {
  cat(sprintf("  pair %d/%d  plot %d  %d->%d\n", i, nrow(manifest),
              manifest$PLOT[i], manifest$YEAR_PREV[i], manifest$YEAR_CURR[i]))
  res[[i]] <- run_one_pair(manifest[i, ])
}
out <- do.call(rbind, res)
write.csv(out, file.path(WD, "silc_cfi_acadiangy_pred_v3.csv"),
          row.names = FALSE)

cat("\n=== AcadianGY CFI scorecard (BA + merch vol) ===\n")
ok <- !is.na(out$BA_PRED_ft2ac) & !is.na(out$BA_OBS_CURR)
cat(sprintf("  n pairs OK : %d / %d\n", sum(ok), nrow(out)))
cat(sprintf("  Pred BA mean   : %.1f ft^2/ac\n",   mean(out$BA_PRED_ft2ac[ok])))
cat(sprintf("  Obs BA mean    : %.1f ft^2/ac\n",   mean(out$BA_OBS_CURR[ok])))
cat(sprintf("  Pred Cords mean: %.2f cords/ac\n",  mean(out$Cords_PRED_ac[ok])))
cat(sprintf("  Pred BdFt mean : %.0f bd ft/ac\n",  mean(out$BdFt_PRED_ac[ok])))
