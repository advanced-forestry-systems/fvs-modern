#!/usr/bin/env Rscript
# silc_cfi_overview.R   (base R port of Aaron's tidyverse silc_cfi_stand_metrics.r)
# =====================================================================
# Stand-level overview of the SILC CFI (Davistown / Town 13 Tract 2,
# northern Maine, STATECD 23), 10 fixed-area 1/5-acre plots, EXPF = 5,
# measured 1981, 1986, 1990, 1995, 2000.
#
# Produces:
#   * silc_cfi_plot_summary.csv   per-plot per-year stand metrics
#   * silc_cfi_pai_summary.csv    PAI components summary by interval
#   * silc_cfi_plot_strata_map.csv  CFI plot -> (forest_type, density_class)
#                                    using top-2 species share + REL_DENSITY
#   * silc_cfi_ba_trajectory.png  BA trajectory per plot 1981-2000
#   * silc_cfi_sdmd.png           Reineke SDMD with all plot-year points
#   * silc_cfi_pai_components.png stacked PAI bar (surv / ingrowth / mort)
#   * silc_cfi_species_composition.png stacked BA by top 6 species over time
#   * silc_cfi_dia_growth_by_species.png annual DBH growth boxplots
# =====================================================================
EXPF <- 5.0; REF_DBH <- 10.0; SDI_MAX <- 450
od <- "/sessions/friendly-compassionate-rubin/mnt/repos--fvs-modern/calibration/output/comparisons_overstory/silc_cfi"

sm  <- read.csv(file.path(od, "STAND_METRICS.csv"))
pai <- read.csv(file.path(od, "STAND_PAI.csv"))
sp  <- read.csv(file.path(od, "STAND_SPECIES.csv"))
tr  <- read.csv(file.path(od, "TREE.csv"))
grm <- read.csv(file.path(od, "GRM.csv"))

# ----- 1. per-plot summary CSV ---------------------------------------
summ <- sm[sm$METRICS_RELIABLE == "Y", ]
summ <- summ[, c("PLOT", "MEASYEAR", "N_LIVE", "TPA", "TPA_ALL",
                 "BA_FT2_AC", "QMD_IN", "SDI", "REL_DENSITY",
                 "CURTIS_RD", "DBH_MEAN_IN", "LOREYS_HT_FT")]
write.csv(summ, file.path(od, "silc_cfi_plot_summary.csv"),
          row.names = FALSE)

# ----- 2. PAI summary CSV --------------------------------------------
pai_ok <- pai[!is.na(pai$PAI_BA_NET_FT2ACY), ]
write.csv(pai_ok, file.path(od, "silc_cfi_pai_summary.csv"),
          row.names = FALSE)
cat("=== PAI summary (mean +/- sd, ft^2/ac/yr) ===\n")
for (col in c("PAI_BA_SURV_FT2ACY","PAI_BA_MORT_FT2ACY",
              "PAI_BA_INGR_FT2ACY","PAI_BA_NET_FT2ACY")) {
  v <- pai_ok[[col]]
  cat(sprintf("  %-22s  %6.3f +/- %5.3f   (n=%d)\n",
              col, mean(v, na.rm=TRUE), sd(v, na.rm=TRUE), sum(!is.na(v))))
}

# ----- 3. CFI plot -> 5x2 strata mapping -----------------------------
# Use year-2000 composition and REL_DENSITY to assign each plot
sp_lookup <- read.csv(file.path(od, "SPP_LOOKUP.csv"))
# Group species into the SILC matrix categories
# Cedar      = CE (white cedar)
# Hardwood   = PB, YB, RM, SM, BE, WA, RO, BC, BTA, QA, BA, AB (broadleaf)
# Mixedwood  = mixed >= 30% SW and >= 30% HW per plot (composition rule)
# Comm SW    = RS, BS, WS, NS, BF, EH (red/black/white/Norway spruce,
#              balsam fir, eastern hemlock = HM is treated as comm SW)
# Other SW   = JP, RP, WP, NS (jack/red/white pine; pines outside spruce-fir)
hw_set <- c("PB","YB","RM","SM","BE","WA","RO","BTA","QA","BA","AB",
            "BC","BTA","BTA","WB","GB")
sw_comm_set <- c("RS","BS","WS","BF","HM")
sw_other_set <- c("JP","RP","WP")
cedar_set <- c("CE")

assign_type_density <- function(plot_id, year, sp_df, sm_df) {
  s <- sp_df[sp_df$PLOT == plot_id & sp_df$MEASYEAR == year, ]
  if (nrow(s) == 0) return(c(forest_type=NA, density_class=NA,
                              ba_pct_cedar=NA, ba_pct_hw=NA,
                              ba_pct_sw_comm=NA, ba_pct_sw_other=NA))
  bp <- function(set) sum(s$BA_PCT[s$SPECIES %in% set], na.rm=TRUE)
  pct_c  <- bp(cedar_set)
  pct_h  <- bp(hw_set)
  pct_sc <- bp(sw_comm_set)
  pct_so <- bp(sw_other_set)
  # classify forest type using a Mixedwood-only-when-both-present rule:
  # Mixedwood when both HW and SW components exceed 30% (Acadian
  # HS/SH definition; both leading groups must each carry >=30% BA).
  # Below 30% on the minor group, the dominant group wins.
  # Unclassifiable when neither group reaches 50% AND total < 50%
  # (covers the unknown-species cases like CFI plot 1103).
  ft <- if (pct_c >= 30)                                      "Cedar"
        else if (pct_so >= 50)                                "Other Softwood"
        else if (pct_h + pct_sc + pct_c + pct_so < 50)        "Unclassifiable"
        else if (pct_h >= 30 && pct_sc >= 30)                 "Mixedwood"
        else if (pct_h >= pct_sc)                             "Hardwood"
        else                                                  "Commercial Softwood"
  rd <- sm_df$REL_DENSITY[sm_df$PLOT==plot_id & sm_df$MEASYEAR==year]
  rd <- if (length(rd) == 0 || !is.finite(rd)) NA_real_ else rd
  # Use the CFI sample median (~0.28) as the A+B vs C+D cut, matching
  # SILC's "A+B = upper half of stocking" intent in their sample.
  # SILC Matrix BA-stocking lines are not directly comparable to the
  # Long 1985 SDI_max = 450 reference, so we adopt a sample-relative
  # threshold for the CFI evaluation.
  dc <- if (is.na(rd))      NA_character_
        else if (rd >= 0.28) "A+B (high)"
        else                 "C+D (low)"
  c(forest_type = ft, density_class = dc,
    ba_pct_cedar = round(pct_c,1), ba_pct_hw = round(pct_h,1),
    ba_pct_sw_comm = round(pct_sc,1), ba_pct_sw_other = round(pct_so,1))
}

plot_map <- do.call(rbind, lapply(sort(unique(sm$PLOT)), function(p) {
  out <- assign_type_density(p, 2000, sp, sm)
  c(PLOT = p, MEASYEAR = 2000, out)
}))
plot_map <- as.data.frame(plot_map, stringsAsFactors = FALSE)
plot_map$BA_FT2_AC_2000 <- sapply(plot_map$PLOT, function(p) {
  v <- sm$BA_FT2_AC[sm$PLOT == p & sm$MEASYEAR == 2000]
  if (length(v) == 0) NA_real_ else v
})
plot_map$REL_DENSITY_2000 <- sapply(plot_map$PLOT, function(p) {
  v <- sm$REL_DENSITY[sm$PLOT == p & sm$MEASYEAR == 2000]
  if (length(v) == 0) NA_real_ else v
})
write.csv(plot_map, file.path(od, "silc_cfi_plot_strata_map.csv"),
          row.names = FALSE)
cat("\n=== CFI plot -> 5x2 strata assignment (year 2000) ===\n")
print(plot_map[, c("PLOT","forest_type","density_class",
                    "BA_FT2_AC_2000","REL_DENSITY_2000",
                    "ba_pct_cedar","ba_pct_hw","ba_pct_sw_comm","ba_pct_sw_other")],
      row.names = FALSE)

# ----- 4. BA trajectory figure ---------------------------------------
plots <- sort(unique(summ$PLOT))
cols  <- rainbow(length(plots), s=0.7, v=0.85)
names(cols) <- as.character(plots)
png(file.path(od, "silc_cfi_ba_trajectory.png"),
    width = 2200, height = 1100, res = 170)
par(mar=c(4.5,4.5,3.4,1.0), mgp=c(2.7,0.6,0))
plot(NA, xlim=c(1980, 2001), ylim=c(0, max(summ$BA_FT2_AC)*1.1),
     xlab="Measurement year", ylab="Basal area (ft^2/ac)",
     main="SILC CFI: per-plot BA trajectory 1981-2000  (reliable plots only)",
     las=1, font.main=2)
for (p in plots) {
  s <- summ[summ$PLOT == p, ]
  s <- s[order(s$MEASYEAR), ]
  lines(s$MEASYEAR, s$BA_FT2_AC, col=cols[as.character(p)], lwd=2)
  points(s$MEASYEAR, s$BA_FT2_AC, pch=19, col=cols[as.character(p)], cex=1.0)
}
legend("topleft", legend=paste("Plot", plots),
       col=cols, lwd=2, pch=19, bty="n", ncol=2, cex=0.8)
mtext("Davistown, Maine | Plot design: 1/5-acre fixed area, EXPF=5",
      side=1, line=3.3, cex=0.78, col="#444")
dev.off()

# ----- 5. Reineke SDMD -----------------------------------------------
reineke_tpa <- function(x, pct=1.0) SDI_MAX * pct * (REF_DBH / x)^1.605
png(file.path(od, "silc_cfi_sdmd.png"),
    width = 1700, height = 1500, res = 175)
par(mar=c(4.5,4.5,3.5,1.0), mgp=c(2.7,0.6,0))
xs <- exp(seq(log(4), log(20), length.out=200))
plot(NA, log="xy", xlim=c(4,20), ylim=c(5,700),
     xlab="QMD (in)", ylab="TPA",
     main="SILC CFI: stand density management diagram (Reineke)",
     las=1, font.main=2)
lines(xs, reineke_tpa(xs, 1.00), lwd=2,   col="black")
lines(xs, reineke_tpa(xs, 0.60), lwd=1.5, col="grey40", lty=2)
lines(xs, reineke_tpa(xs, 0.35), lwd=1.5, col="grey65", lty=3)
text(15, reineke_tpa(15, 1.00)*1.1, "Max density (SDI=450)", cex=0.78, col="black")
text(15, reineke_tpa(15, 0.60)*1.1, "60% RD", cex=0.78, col="grey40")
text(15, reineke_tpa(15, 0.35)*1.1, "35% RD", cex=0.78, col="grey65")
sm2 <- summ[summ$TPA > 0 & summ$QMD_IN > 0, ]
for (p in plots) {
  s <- sm2[sm2$PLOT == p, ]
  s <- s[order(s$MEASYEAR), ]
  if (nrow(s) >= 2) lines(s$QMD_IN, s$TPA,
                          col=cols[as.character(p)], lwd=1.5, lty=1)
  points(s$QMD_IN, s$TPA, col=cols[as.character(p)], pch=19, cex=1.2)
}
legend("bottomleft", legend=paste("Plot", plots),
       col=cols, pch=19, lwd=1.5, bty="n", ncol=2, cex=0.8)
mtext("SDI_max = 450 (Long 1985 NE softwood reference)",
      side=1, line=3.3, cex=0.78, col="#444")
dev.off()

# ----- 6. PAI components figure --------------------------------------
png(file.path(od, "silc_cfi_pai_components.png"),
    width = 2200, height = 1000, res = 170)
par(mar=c(7.0,4.5,3.4,1.0), mgp=c(2.8,0.6,0))
pai_ok$LABEL <- sprintf("P%d\n%d-%d", pai_ok$PLOT,
                        pai_ok$YEAR_PREV, pai_ok$YEAR_CURR)
M <- rbind(
  surv  = pai_ok$PAI_BA_SURV_FT2ACY,
  ingr  = pai_ok$PAI_BA_INGR_FT2ACY,
  mort  = -pai_ok$PAI_BA_MORT_FT2ACY
)
colnames(M) <- pai_ok$LABEL
cols_pai <- c("#2E86AB", "#A8D5BA", "#E84855")
bp <- barplot(M, col=cols_pai, border=NA,
              ylab="PAI BA (ft^2/ac/yr)", las=2, cex.names=0.62)
abline(h=0, lwd=1.2)
# net (sum of surv + ingr - mort) overlay
net <- pai_ok$PAI_BA_NET_FT2ACY
points(bp, net, pch=18, col="black", cex=1.4)
lines(bp, net, col="black", lwd=1)
legend("topleft", legend=c("Survivor growth","Ingrowth","Mortality loss",
                            "Net PAI"),
       fill=c(cols_pai, NA), border=NA,
       pch=c(NA,NA,NA,18), col=c(rep(NA,3), "black"),
       lwd=c(NA,NA,NA,1), bty="n", cex=0.85)
title(main="SILC CFI: PAI in basal area by component, 1981-2000",
      cex.main=1.05, font.main=2)
dev.off()

# ----- 7. Species composition stacked bar ----------------------------
top_spp <- aggregate(BA_FT2_AC ~ SPECIES, data=sp, FUN=sum)
top_spp <- top_spp[order(-top_spp$BA_FT2_AC), ]
top6 <- head(top_spp$SPECIES, 6)
sp$SPP_G <- ifelse(sp$SPECIES %in% top6, sp$SPECIES, "Other")
mn <- aggregate(BA_FT2_AC ~ MEASYEAR + SPP_G, data=sp,
                FUN=function(x) mean(x, na.rm=TRUE))
yrs <- sort(unique(mn$MEASYEAR))
gs  <- c(top6, "Other")
Mc <- matrix(0, nrow=length(gs), ncol=length(yrs),
             dimnames=list(gs, as.character(yrs)))
for (i in seq_len(nrow(mn))) {
  Mc[mn$SPP_G[i], as.character(mn$MEASYEAR[i])] <- mn$BA_FT2_AC[i]
}
png(file.path(od, "silc_cfi_species_composition.png"),
    width = 1900, height = 1100, res = 170)
par(mar=c(4.5, 4.5, 3.4, 7), mgp=c(2.7,0.6,0), xpd=NA)
cols_sp <- c("#66c2a5","#fc8d62","#8da0cb","#e78ac3","#a6d854",
              "#ffd92f","#cccccc")
bp2 <- barplot(Mc, col=cols_sp, border=NA, las=1,
               xlab="Measurement year", ylab="Mean BA (ft^2/ac)",
               main="SILC CFI: mean BA by species over time (reliable plots)",
               cex.main=1.05, font.main=2)
legend(par("usr")[2] + 0.1, par("usr")[4],
       legend=gs, fill=cols_sp, border=NA, bty="n",
       cex=0.9, title="Species")
dev.off()

# ----- 8. Annual DBH growth by species (boxplot) ---------------------
spcd_abbr <- c("12"="BF","97"="RS","375"="PB","371"="YB","316"="RM",
                "241"="CE","261"="HM","95"="BS","91"="WS","105"="JP",
                "129"="WP","318"="SM","531"="BE","746"="QA","833"="RO")
g <- grm[grm$GRM_CLASS=="GROWTH" & !is.na(grm$ANN_DIA_GROWTH_IN) &
        grm$ANN_DIA_GROWTH_IN > 0, ]
trkey <- tr[, c("CN", "SPCD")]
g <- merge(g, trkey, by.x="CURR_TRE_CN", by.y="CN", all.x=TRUE)
g$SP <- spcd_abbr[as.character(g$SPCD)]
g$SP[is.na(g$SP)] <- "Other"
g <- g[g$SP != "Other", ]
ord <- names(sort(tapply(g$ANN_DIA_GROWTH_IN, g$SP, median)))
g$SP <- factor(g$SP, levels=ord)
png(file.path(od, "silc_cfi_dia_growth_by_species.png"),
    width = 1900, height = 1100, res = 170)
par(mar=c(4.5,4.5,3.5,1.0), mgp=c(2.7,0.6,0))
boxplot(ANN_DIA_GROWTH_IN ~ SP, data=g, col="#2E86AB",
        outline=TRUE, outpch=20, outcex=0.4,
        ylab="Annual diameter increment (in/yr)", xlab="Species",
        main="SILC CFI: annual DBH growth by species (paired remeasurement)",
        cex.main=1.05, font.main=2, las=1)
abline(h=mean(g$ANN_DIA_GROWTH_IN, na.rm=TRUE), col="red", lwd=1.5, lty=2)
text(0.6, mean(g$ANN_DIA_GROWTH_IN, na.rm=TRUE),
     "  overall mean", col="red", pos=4, cex=0.8)
dev.off()

cat("\n=== Stand summary (reliable plot-years, mean across plots) ===\n")
s1 <- aggregate(cbind(TPA, BA_FT2_AC, QMD_IN, REL_DENSITY) ~ MEASYEAR,
                data = summ, FUN = mean)
print(round(s1, 2), row.names=FALSE)

cat("\nDone. Artifacts in", od, "\n")
