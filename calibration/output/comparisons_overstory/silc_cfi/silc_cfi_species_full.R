#!/usr/bin/env Rscript
# silc_cfi_species_full.R
# =====================================================================
# Multi-model species composition scorecard:
#   * observed (CFI TREE.csv)
#   * AGM / AcadianGY (silc_cfi_acadiangy_treelist.csv)
#   * OSM-ACD (TreeListProjections.csv at matched year_curr)
#   * FVS-NE / FVS-ACD calibrated (silc_cfi_fvs_treelist.csv)
#
# Output: % share of stand basal area in each of four groups
#   (Cedar, Comm SW, Other SW, Hardwood). Per-pair contributions
#   averaged across the 17 routine-growth pairs.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

# --- Species group mapping (matches the earlier scorecard) ---
sg_softwood_comm  <- c("RS","BF","BS","WS","NS","EH","HM")
sg_softwood_other <- c("JP","RP","WP","TL","TA")
sg_cedar          <- c("WC","CE")
sg_hardwood       <- c("PB","WB","YB","RM","SM","BE","BC","QA","RO","WA","GB","AB")
species_group <- function(sp) {
  ifelse(sp %in% sg_cedar, "Cedar",
    ifelse(sp %in% sg_softwood_comm, "Comm SW",
      ifelse(sp %in% sg_softwood_other, "Other SW",
        ifelse(sp %in% sg_hardwood, "Hardwood", "Other"))))
}
spcd_to_code <- c("12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
                  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
                  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
                  "541"="WA","934"="GB")
EXPF_CFI <- 5.0; ACRES_PER_HA <- 2.4710538147

# --- Load ---
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))
tr <- read.csv(file.path(od, "TREE.csv"))
ag_tl <- read.csv(file.path(od, "silc_cfi_acadiangy_treelist.csv"),
                  stringsAsFactors = FALSE)
osm_t <- read.csv(file.path(od, "silc_cfi_TreeListProjections.csv"))
fvs_tl <- tryCatch(
  read.csv(file.path(od, "silc_cfi_fvs_treelist.csv"), stringsAsFactors = FALSE),
  error = function(e) NULL)

# Identify routine-growth pairs
m_naive <- read.csv(file.path(od, "silc_cfi_acadiangy_pred_v3.csv"))
m_naive$establishment <- m_naive$BA_OBS_PREV < 10 |
                         abs(m_naive$PAI_NET_OBS) > 5
core_keys <- m_naive[!m_naive$establishment, c("PLOT","YEAR_PREV","YEAR_CURR")]

# --- compute shares per (pair, model) ---
share_from_tree <- function(td, sp_col, expf_col, dbh_col, ht_col = NULL) {
  td <- td[is.finite(td[[dbh_col]]) & td[[dbh_col]] >= 4.5, ]
  if (nrow(td) == 0) return(NULL)
  td$BA  <- 0.005454 * td[[dbh_col]]^2 * td[[expf_col]]
  td$grp <- species_group(td[[sp_col]])
  agg <- tapply(td$BA, td$grp, sum, na.rm = TRUE)
  # Renormalize within the 4 identified groups (exclude "Other" from denom)
  identified <- sum(c(agg["Cedar"], agg["Comm SW"],
                       agg["Other SW"], agg["Hardwood"]), na.rm = TRUE)
  if (identified <= 0) return(NULL)
  data.frame(
    Cedar    = as.numeric(agg["Cedar"])    %||% 0,
    `Comm SW`  = as.numeric(agg["Comm SW"])  %||% 0,
    `Other SW` = as.numeric(agg["Other SW"]) %||% 0,
    Hardwood = as.numeric(agg["Hardwood"]) %||% 0,
    Total    = identified,
    check.names = FALSE
  )
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

# --- Observed shares ---
obs_share <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y & tr$STATUSCD == 1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  if (nrow(t) == 0) return(NULL)
  t$sp_code <- spcd_to_code[as.character(t$SPCD)]
  t$sp_code[is.na(t$sp_code)] <- "OH"
  t$BA  <- 0.005454 * t$DIA_IN^2 * EXPF_CFI
  t$grp <- species_group(t$sp_code)
  agg <- tapply(t$BA, t$grp, sum, na.rm = TRUE)
  identified <- sum(c(agg["Cedar"], agg["Comm SW"],
                       agg["Other SW"], agg["Hardwood"]), na.rm = TRUE)
  if (identified <= 0) return(NULL)
  data.frame(Cedar = as.numeric(agg["Cedar"]) %||% 0,
             `Comm SW` = as.numeric(agg["Comm SW"]) %||% 0,
             `Other SW`= as.numeric(agg["Other SW"]) %||% 0,
             Hardwood = as.numeric(agg["Hardwood"]) %||% 0,
             Total = identified, check.names = FALSE)
}

# --- AGM shares ---
ag_tl$DBH_in <- ag_tl$DBH_cm / 2.54
ag_tl$EXPF_ac<- ag_tl$EXPF_ha / ACRES_PER_HA
ag_share <- function(p, yp) {
  sub <- ag_tl[ag_tl$PLOT == p & ag_tl$YEAR_PREV == yp, ]
  share_from_tree(sub, "SP", "EXPF_ac", "DBH_in")
}

# --- OSM shares ---
osm_t$DBH_in <- osm_t$DBH / 2.54
osm_t$EXPF_ac<- osm_t$Stems / ACRES_PER_HA
sid_to_pair <- ps[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR")]
sid_to_pair$SurveyID <- seq_len(nrow(sid_to_pair))
osm_share <- function(p, yp, yc, period_yr) {
  pp <- sid_to_pair[sid_to_pair$PLOT == p & sid_to_pair$YEAR_PREV == yp, ]
  if (nrow(pp) == 0) return(NULL)
  sub <- osm_t[osm_t$SurveyID == pp$SurveyID, ]
  if (nrow(sub) == 0) return(NULL)
  y0 <- min(sub$Year); sub$yr_off <- sub$Year - y0
  died_ok <- is.na(sub$Died) | as.numeric(sub$Died) == 0
  cut_ok  <- is.na(sub$Cut)  | sub$Cut %in% c("False","false","FALSE",FALSE)
  sub_t <- sub[sub$yr_off == period_yr & died_ok & cut_ok &
               is.finite(sub$DBH_in) & sub$DBH_in >= 4.5, ]
  share_from_tree(sub_t, "Species", "EXPF_ac", "DBH_in")
}

# --- FVS shares (per variant/config) ---
fvs_share <- function(p, yp, variant, config) {
  if (is.null(fvs_tl)) return(NULL)
  sub <- fvs_tl[fvs_tl$PLOT == p & fvs_tl$YEAR_PREV == yp &
                fvs_tl$variant == variant & fvs_tl$config == config, ]
  share_from_tree(sub, "SP", "EXPF_ac", "DBH_in")
}

# --- Assemble: one row per pair, columns per model ---
out <- core_keys
out$obs_Cedar <- out$obs_CommSW <- out$obs_OtherSW <- out$obs_Hardwood <- NA
out$agm_Cedar <- out$agm_CommSW <- out$agm_OtherSW <- out$agm_Hardwood <- NA
out$osm_Cedar <- out$osm_CommSW <- out$osm_OtherSW <- out$osm_Hardwood <- NA
out$fvsne_Cedar <- out$fvsne_CommSW <- out$fvsne_OtherSW <- out$fvsne_Hardwood <- NA
out$fvsacd_Cedar <- out$fvsacd_CommSW <- out$fvsacd_OtherSW <- out$fvsacd_Hardwood <- NA

fill_row <- function(r, prefix, sh) {
  if (is.null(sh) || sh$Total <= 0) return(r)
  r[[paste0(prefix, "_Cedar")]]    <- 100 * sh$Cedar    / sh$Total
  r[[paste0(prefix, "_CommSW")]]   <- 100 * sh$`Comm SW`  / sh$Total
  r[[paste0(prefix, "_OtherSW")]]  <- 100 * sh$`Other SW` / sh$Total
  r[[paste0(prefix, "_Hardwood")]] <- 100 * sh$Hardwood / sh$Total
  r
}

for (i in seq_len(nrow(out))) {
  r <- out[i, ]
  p <- r$PLOT; yp <- r$YEAR_PREV; yc <- r$YEAR_CURR
  period_yr <- yc - yp
  out[i, ] <- fill_row(r, "obs",    obs_share(p, yc))
  out[i, ] <- fill_row(out[i, ], "agm",    ag_share(p, yp))
  out[i, ] <- fill_row(out[i, ], "osm",    osm_share(p, yp, yc, period_yr))
  out[i, ] <- fill_row(out[i, ], "fvsne",  fvs_share(p, yp, "NE", "calibrated"))
  out[i, ] <- fill_row(out[i, ], "fvsacd", fvs_share(p, yp, "ACD","calibrated"))
}
write.csv(out, file.path(od, "silc_cfi_species_full.csv"), row.names = FALSE)

# Mean across pairs, four groups, five models
groups <- c("Cedar","CommSW","OtherSW","Hardwood")
models_all <- list(observed = "obs", `AGM (AcadianGY)` = "agm",
                    `OSM-ACD` = "osm", `FVS-NE cal` = "fvsne",
                    `FVS-ACD cal` = "fvsacd")
# Drop models with no data (FVS treelist deferred)
models <- models_all
for (mn in names(models_all)) {
  pref <- models_all[[mn]]
  if (all(is.na(out[[paste0(pref, "_Cedar")]]))) {
    models[[mn]] <- NULL
  }
}
M <- matrix(0, nrow = length(models), ncol = length(groups),
            dimnames = list(names(models), groups))
for (mn in names(models)) {
  pref <- models[[mn]]
  for (g in groups) {
    col <- paste0(pref, "_", g)
    M[mn, g] <- mean(out[[col]], na.rm = TRUE)
  }
}
cat("=== Mean species-group BA share (%, n=17 routine-growth) ===\n")
print(round(M, 1))

# === Figure: grouped bar (4 groups x 5 models) ===
type_order <- c("Cedar","CommSW","OtherSW","Hardwood")
cols_grp <- c("#7b3294","#1A3D28","#5e4fa2","#a6dba0")
png(file.path(od, "silc_cfi_species_full.png"),
    width = 2400, height = 1100, res = 165)
par(mar = c(4.6, 4.8, 3.6, 8), mgp = c(2.8, 0.6, 0), xpd = NA)
bp <- barplot(t(M), beside = TRUE, col = cols_grp, border = NA, las = 1,
              ylab = "% of stand basal area",
              main = "Species composition by model -- mean across 17 routine-growth CFI pairs",
              cex.main = 1.1, font.main = 2)
legend(par("usr")[2] + 0.3, par("usr")[4],
       legend = c("Cedar","Comm SW","Other SW","Hardwood"),
       fill = cols_grp, border = NA, bty = "n",
       cex = 1.0, title = "Species group")
# Bias annotations
for (j in seq_along(groups)) {
  bias <- M[-1, j] - M[1, j]
  for (i in seq_along(bias)) {
    text(bp[j, i + 1], M[i + 1, j] + 1.5,
         sprintf("%+.1f pp", bias[i]), cex = 0.66, font = 2,
         col = ifelse(abs(bias[i]) > 5, "#aa3333", "#444"))
  }
}
dev.off()
cat("\nwrote silc_cfi_species_full.png\n")
