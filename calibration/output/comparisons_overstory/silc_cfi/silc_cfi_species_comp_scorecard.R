#!/usr/bin/env Rscript
# silc_cfi_species_comp_scorecard.R
# =====================================================================
# Species composition scorecard. Compute proportion of total BA by
# species group at year_curr for:
#   * observed (from TREE.csv)
#   * OSM-ACD predicted (TreeListProjections.csv)
#
# AGM and FVS species composition need per-tree output enabled in
# their drivers (defer). Documented here as a gap.
# =====================================================================
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

# Species group definitions (SILC operational)
sg_softwood_comm  <- c("RS","BF","BS","WS","NS","EH","HM")   # SF group + hemlock
sg_softwood_other <- c("JP","RP","WP","TL","TA")             # pine + larch
sg_cedar          <- c("WC","CE")
sg_hardwood       <- c("PB","WB","YB","RM","SM","BE","BC","QA","RO","WA","GB","AB")

species_group <- function(sp) {
  ifelse(sp %in% sg_cedar,           "Cedar",
    ifelse(sp %in% sg_softwood_comm, "Comm SW",
      ifelse(sp %in% sg_softwood_other, "Other SW",
        ifelse(sp %in% sg_hardwood,  "Hardwood",
                                     "Other"))))
}

# SPCD -> 2-char code (CFI uses both)
spcd_to_code <- c("12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
                  "241"="WC","261"="EH","95"="BS","91"="WS","105"="JP",
                  "129"="WP","318"="SM","531"="BC","746"="QA","833"="RO",
                  "541"="WA","934"="GB")

# === Observed composition at year_curr ===
tr <- read.csv(file.path(od, "TREE.csv"))
ps <- read.csv(file.path(od, "silc_cfi_pair_summary.csv"))

obs_comp <- function(p, y) {
  t <- tr[tr$PLOT == p & tr$MEASYEAR == y & tr$STATUSCD == 1 &
          is.finite(tr$DIA_IN) & tr$DIA_IN >= 4.5, ]
  if (nrow(t) == 0) return(NULL)
  # COMMON_NAME -> 2-char via FIA SPCD lookup
  t$sp_code <- spcd_to_code[as.character(t$SPCD)]
  t$sp_code[is.na(t$sp_code)] <- "OH"
  t$BA  <- 0.005454 * t$DIA_IN^2 * 5   # EXPF = 5
  t$grp <- species_group(t$sp_code)
  agg <- tapply(t$BA, t$grp, sum, na.rm=TRUE)
  total <- sum(agg, na.rm=TRUE)
  if (total <= 0) return(NULL)
  data.frame(PLOT=p, MEASYEAR=y,
             Cedar    = as.numeric(agg["Cedar"])    %||% 0,
             CommSW   = as.numeric(agg["Comm SW"])  %||% 0,
             OtherSW  = as.numeric(agg["Other SW"]) %||% 0,
             Hardwood = as.numeric(agg["Hardwood"]) %||% 0,
             Total    = total)
}

`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

obs_rows <- list()
for (i in seq_len(nrow(ps))) {
  pr <- ps[i, ]
  oc <- obs_comp(pr$PLOT, pr$YEAR_CURR)
  if (!is.null(oc)) {
    oc$YEAR_PREV <- pr$YEAR_PREV
    obs_rows[[i]] <- oc
  }
}
obs <- do.call(rbind, obs_rows)
obs$pct_Cedar    <- 100 * obs$Cedar    / obs$Total
obs$pct_CommSW   <- 100 * obs$CommSW   / obs$Total
obs$pct_OtherSW  <- 100 * obs$OtherSW  / obs$Total
obs$pct_Hardwood <- 100 * obs$Hardwood / obs$Total

# === OSM predicted composition at year_curr ===
osm_tree <- read.csv(file.path(od, "silc_cfi_TreeListProjections.csv"))
ACRES_PER_HA <- 2.4710538147
osm_tree$DBH_in <- osm_tree$DBH / 2.54
osm_tree$EXPF_ac<- osm_tree$Stems / ACRES_PER_HA
sid_to_pair <- ps[, c("PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR")]
sid_to_pair$SurveyID <- seq_len(nrow(sid_to_pair))

osm_comp <- function(pp) {
  sub <- osm_tree[osm_tree$SurveyID == pp$SurveyID, ]
  if (nrow(sub) == 0) return(NULL)
  y0 <- min(sub$Year); sub$yr_off <- sub$Year - y0
  died_ok <- is.na(sub$Died) | as.numeric(sub$Died) == 0
  cut_ok  <- is.na(sub$Cut)  | sub$Cut %in% c("False","false","FALSE",FALSE)
  sub_t <- sub[sub$yr_off == pp$PERIOD_YR & died_ok & cut_ok &
               is.finite(sub$DBH_in) & sub$DBH_in >= 4.5, ]
  if (nrow(sub_t) == 0) return(NULL)
  sub_t$BA  <- 0.005454 * sub_t$DBH_in^2 * sub_t$EXPF_ac
  sub_t$grp <- species_group(sub_t$Species)
  agg <- tapply(sub_t$BA, sub_t$grp, sum, na.rm=TRUE)
  total <- sum(agg, na.rm=TRUE)
  if (total <= 0) return(NULL)
  data.frame(PLOT=pp$PLOT, YEAR_PREV=pp$YEAR_PREV, YEAR_CURR=pp$YEAR_CURR,
             osm_Cedar    = as.numeric(agg["Cedar"])    %||% 0,
             osm_CommSW   = as.numeric(agg["Comm SW"])  %||% 0,
             osm_OtherSW  = as.numeric(agg["Other SW"]) %||% 0,
             osm_Hardwood = as.numeric(agg["Hardwood"]) %||% 0,
             osm_Total    = total)
}
osm_rows <- list()
for (i in seq_len(nrow(sid_to_pair))) {
  oc <- osm_comp(sid_to_pair[i, ])
  if (!is.null(oc)) osm_rows[[i]] <- oc
}
osm_comp_df <- do.call(rbind, osm_rows)
osm_comp_df$osm_pct_Cedar    <- 100 * osm_comp_df$osm_Cedar    / osm_comp_df$osm_Total
osm_comp_df$osm_pct_CommSW   <- 100 * osm_comp_df$osm_CommSW   / osm_comp_df$osm_Total
osm_comp_df$osm_pct_OtherSW  <- 100 * osm_comp_df$osm_OtherSW  / osm_comp_df$osm_Total
osm_comp_df$osm_pct_Hardwood <- 100 * osm_comp_df$osm_Hardwood / osm_comp_df$osm_Total

# === Compare ===
m <- merge(obs, osm_comp_df,
           by.x = c("PLOT","YEAR_PREV","MEASYEAR"),
           by.y = c("PLOT","YEAR_PREV","YEAR_CURR"),
           all.x = TRUE)

bias_pct <- function(p, o) 100*(mean(p, na.rm=TRUE)/mean(o, na.rm=TRUE) - 1)

cat("=== Observed species composition (mean across 24 plot-year_curr records) ===\n")
cat(sprintf("  Cedar       %.1f%%\n", mean(obs$pct_Cedar)))
cat(sprintf("  Comm SW     %.1f%%\n", mean(obs$pct_CommSW)))
cat(sprintf("  Other SW    %.1f%%\n", mean(obs$pct_OtherSW)))
cat(sprintf("  Hardwood    %.1f%%\n", mean(obs$pct_Hardwood)))

cat("\n=== OSM-ACD predicted vs observed (% BA share, bias points) ===\n")
ok <- !is.na(m$osm_Total)
cat(sprintf("  Cedar       obs %.1f%%   pred %.1f%%   bias %+.1f pp\n",
            mean(m$pct_Cedar[ok]), mean(m$osm_pct_Cedar[ok]),
            mean(m$osm_pct_Cedar[ok]) - mean(m$pct_Cedar[ok])))
cat(sprintf("  Comm SW     obs %.1f%%   pred %.1f%%   bias %+.1f pp\n",
            mean(m$pct_CommSW[ok]), mean(m$osm_pct_CommSW[ok]),
            mean(m$osm_pct_CommSW[ok]) - mean(m$pct_CommSW[ok])))
cat(sprintf("  Other SW    obs %.1f%%   pred %.1f%%   bias %+.1f pp\n",
            mean(m$pct_OtherSW[ok]), mean(m$osm_pct_OtherSW[ok]),
            mean(m$osm_pct_OtherSW[ok]) - mean(m$pct_OtherSW[ok])))
cat(sprintf("  Hardwood    obs %.1f%%   pred %.1f%%   bias %+.1f pp\n",
            mean(m$pct_Hardwood[ok]), mean(m$osm_pct_Hardwood[ok]),
            mean(m$osm_pct_Hardwood[ok]) - mean(m$pct_Hardwood[ok])))

write.csv(m, file.path(od, "silc_cfi_species_comp.csv"), row.names = FALSE)

# === Stacked bar figure ===
type_order <- c("Cedar","Comm SW","Other SW","Hardwood")
cols <- c("#7b3294", "#1A3D28", "#5e4fa2", "#a6dba0")

# Mean observed composition through time
yrs <- sort(unique(obs$MEASYEAR))
M <- matrix(0, nrow=4, ncol=length(yrs),
            dimnames=list(type_order, as.character(yrs)))
for (y in yrs) {
  sub <- obs[obs$MEASYEAR == y, ]
  M["Cedar",    as.character(y)] <- mean(sub$pct_Cedar)
  M["Comm SW",  as.character(y)] <- mean(sub$pct_CommSW)
  M["Other SW", as.character(y)] <- mean(sub$pct_OtherSW)
  M["Hardwood", as.character(y)] <- mean(sub$pct_Hardwood)
}

png(file.path(od, "silc_cfi_species_composition_deck.png"),
    width = 2400, height = 1000, res = 165)
par(mfrow = c(1, 2), mar = c(4.8, 4.8, 3.6, 7), mgp = c(2.8, 0.6, 0), xpd=NA)

bp <- barplot(M, col=cols, border=NA, las=1,
              xlab="Measurement year",
              ylab="% of stand basal area",
              main="Observed CFI species composition over time",
              cex.main=1.1, font.main=2)
legend(par("usr")[2] + 0.3, par("usr")[4],
       legend = rev(type_order), fill = rev(cols),
       border = NA, bty = "n", cex = 1.0, title = "Species group")

# Observed vs OSM-ACD predicted at year_curr (mean across pairs)
M2 <- rbind(
  observed = c(Cedar = mean(m$pct_Cedar[ok]),
               `Comm SW` = mean(m$pct_CommSW[ok]),
               `Other SW` = mean(m$pct_OtherSW[ok]),
               Hardwood = mean(m$pct_Hardwood[ok])),
  `OSM-ACD pred` = c(Cedar = mean(m$osm_pct_Cedar[ok]),
                     `Comm SW` = mean(m$osm_pct_CommSW[ok]),
                     `Other SW` = mean(m$osm_pct_OtherSW[ok]),
                     Hardwood = mean(m$osm_pct_Hardwood[ok]))
)
bp <- barplot(t(M2), beside = TRUE, col = cols, border = NA, las = 1,
              ylab = "% of stand basal area",
              main = "Observed vs OSM-ACD species share (mean)",
              cex.main = 1.1, font.main = 2)
legend(par("usr")[2] + 0.3, par("usr")[4],
       legend = type_order, fill = cols, border = NA, bty = "n",
       cex = 1.0, title = "Species group")
dev.off()
cat("\nwrote silc_cfi_species_composition_deck.png\n")
cat("wrote silc_cfi_species_comp.csv\n")
