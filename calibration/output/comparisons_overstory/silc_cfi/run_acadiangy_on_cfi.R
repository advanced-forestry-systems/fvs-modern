#!/usr/bin/env Rscript
# run_acadiangy_on_cfi.R
# =====================================================================
# Cardinal-side driver. Loops over the 24 SILC CFI prediction pairs,
# runs AcadianGY_12.3.5.r per pair, captures predicted year_curr
# BA / TPA / QMD / NetPAI components, writes one row per pair to
# silc_cfi_acadiangy_pred.csv.
#
# Pair-list manifest: silc_cfi_pair_summary.csv
# Per-pair tree lists: pair_input/pair_<PLOT>_<YEAR_PREV>_tree.csv
# AcadianGY source: AcadianGY_12.3.5.r (sourced here)
# =====================================================================
suppressPackageStartupMessages({
  library(plyr)
  library(dplyr)
  library(purrr)
})

WD <- getwd()
source(file.path(WD, "AcadianGY_12.3.5.r"))

# CFI site context — Davistown, Maine
CSI_DEFAULT  <- 12     # AcadianGY CSI default (m)
ELEV_DEFAULT <- 350    # elevation m (Davistown ~ 1000 ft = 305 m)

manifest <- read.csv(file.path(WD, "silc_cfi_pair_summary.csv"),
                     stringsAsFactors = FALSE)

# AcadianGY species code -> abbreviation map (script uses 2-char SP)
spcd_to_sp <- c(
  "12"  = "BF",  "97" = "RS",  "375"= "PB",  "371" = "YB",
  "316" = "RM",  "241"= "CE",  "261"= "HM",  "95"  = "BS",
  "91"  = "WS",  "105"= "JP",  "129"= "WP",  "318" = "SM",
  "531" = "BE",  "746"= "QA",  "833"= "RO",  "541" = "WA",
  "934" = "GB"   # any missing -> map to nearest analog in run_one()
)

run_one <- function(p_row) {
  treef <- file.path(WD, p_row$tree_list_file)
  if (!file.exists(treef)) return(data.frame())
  td <- read.csv(treef, stringsAsFactors = FALSE)
  if (nrow(td) == 0) return(data.frame())

  # Map SPCD -> AcadianGY SP token if needed
  td$SP <- ifelse(is.na(td$SP) | td$SP == "",
                  spcd_to_sp[as.character(td$SPCD)], td$SP)
  td$SP[is.na(td$SP)] <- "BF"  # fallback

  # AcadianGY needs YEAR, PLOT, TREE, SP, DBH, HT, HCB, EXPF, Form, Risk
  td <- td[, c("STAND","YEAR","PLOT","TREE","SP","DBH","HT","HCB","EXPF","Form","Risk")]

  # AcadianGY runs at 1-yr cycle; project PERIOD_YR cycles
  ops <- list(verbose = FALSE,
              INGROWTH = "Y",
              MinDBH = 4.5,
              CutPoint = 0.95)
  stand <- list(CSI = CSI_DEFAULT, ELEV = ELEV_DEFAULT)
  cur <- td
  for (yr in seq_len(p_row$PERIOD_YR)) {
    cur <- tryCatch(AcadianGYOneStand(tree = cur, stand = stand, ops = ops),
                    error = function(e) {
                      cat(sprintf("  AcadianGY error pair %d %d: %s\n",
                                  p_row$PLOT, p_row$YEAR_PREV, e$message))
                      NULL
                    })
    if (is.null(cur)) return(data.frame())
    cur$YEAR <- p_row$YEAR_PREV + yr
  }

  # Aggregate to stand-level
  live <- cur[!is.na(cur$DBH) & cur$DBH >= 4.5, ]
  TPA_pred <- sum(live$EXPF, na.rm = TRUE)
  BA_pred  <- sum(0.005454 * live$DBH^2 * live$EXPF, na.rm = TRUE)
  QMD_pred <- if (TPA_pred > 0) sqrt((BA_pred / TPA_pred) / 0.005454) else NA

  data.frame(
    PLOT       = p_row$PLOT,
    YEAR_PREV  = p_row$YEAR_PREV,
    YEAR_CURR  = p_row$YEAR_CURR,
    PERIOD_YR  = p_row$PERIOD_YR,
    BA_PRED    = BA_pred,
    TPA_PRED   = TPA_pred,
    QMD_PRED   = QMD_pred,
    BA_OBS_PREV= p_row$BA_PREV_FT2AC,
    BA_OBS_CURR= p_row$BA_CURR_FT2AC,
    PAI_NET_PRED = (BA_pred - p_row$BA_PREV_FT2AC) / p_row$PERIOD_YR,
    PAI_NET_OBS  = p_row$PAI_BA_NET_FT2ACY
  )
}

cat(sprintf("Running AcadianGY on %d pairs...\n", nrow(manifest)))
out <- do.call(rbind, lapply(seq_len(nrow(manifest)), function(i) {
  cat(sprintf("  pair %d/%d: plot %d  %d->%d\n",
              i, nrow(manifest), manifest$PLOT[i],
              manifest$YEAR_PREV[i], manifest$YEAR_CURR[i]))
  run_one(manifest[i, ])
}))

write.csv(out, file.path(WD, "silc_cfi_acadiangy_pred.csv"),
          row.names = FALSE)

cat("\n=== AcadianGY scorecard (BA at year_curr) ===\n")
ok <- !is.na(out$BA_PRED) & !is.na(out$BA_OBS_CURR)
if (sum(ok) > 0) {
  bias_pct <- 100*(mean(out$BA_PRED[ok])/mean(out$BA_OBS_CURR[ok]) - 1)
  rmse     <- sqrt(mean((out$BA_PRED[ok] - out$BA_OBS_CURR[ok])^2))
  cat(sprintf("  n=%d  BA bias %+.2f%%  RMSE %.2f ft^2/ac\n",
              sum(ok), bias_pct, rmse))
}
