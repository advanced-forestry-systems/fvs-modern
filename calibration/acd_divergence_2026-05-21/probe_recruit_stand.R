## probe_recruit_stand.R - confirm what STAND value the recruits inherit after
## one cycle on a single FIA plot. Hypothesis: recruits get STAND=1 (not the
## parent PLT_CN), which is why the v21 multi-cycle harness silently drops them
## on cycle 2 (no matching stand_init row -> AcadianGYOneStand error -> NULL).
suppressMessages({ library(dplyr); library(plyr); library(purrr) })
source("/users/PUOM0008/crsfaaron/AcadianGY_12.3.7.r")
cat("model:", AcadianVersionTag, "\n")

# Re-use the v20 harness builders just enough to get one real FIA stand + trees
PROJECT_ROOT <- "/users/PUOM0008/crsfaaron"; set.seed(42); ACRES_PER_HA <- 2.4710538147
vdat <- read.csv(file.path(PROJECT_ROOT, "fvs-modern/calibration/output/comparisons/intermediate/validation_data_acd_post.csv"))
me_plot <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_PLOT.csv"), colClasses = c("CN" = "character"))
me_tree <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_TREE.csv"), colClasses = c("CN" = "character", "PLT_CN" = "character"))
xwk <- read.csv(file.path(PROJECT_ROOT, "fvs-modern/calibration/data/osm_acadian_species_crosswalk.csv"))
acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB","HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM","RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP","WS","YB")
xf <- xwk[!is.na(xwk$FIA) & xwk$FIA > 0 & nchar(xwk$OSM_AD_CmdKey) > 0, ]
xf$AGY_SP <- ifelse(xf$OSM_AD_CmdKey %in% acadgy_species, xf$OSM_AD_CmdKey, ifelse(xf$FIA >= 300, "OH", "OS"))
spcd_to_code <- setNames(xf$AGY_SP, as.character(xf$FIA))
vdat$PLT_CN_t1_str <- as.character(vdat$PLT_CN_t1)
vme <- vdat[vdat$PLT_CN_t1_str %in% as.character(me_plot$CN) & !is.na(vdat$interval_years) & vdat$interval_years %in% 5:10, ]
samp <- vme[sample(nrow(vme), 1), ]; sid <- as.character(samp$PLT_CN_t1)
cat("using plot:", sid, "\n")
ts <- me_tree[me_tree$PLT_CN == sid & !is.na(me_tree$DIA) & me_tree$DIA > 0 & me_tree$STATUSCD == 1, ]
ts$SP <- spcd_to_code[as.character(ts$SPCD)]; mm <- is.na(ts$SP) | nchar(ts$SP) == 0; ts$SP[mm] <- ifelse(ts$SPCD[mm] >= 300, "OH", "OS")
ts$TPA_UNADJ <- as.numeric(ts$TPA_UNADJ)
ts$tpa <- ifelse(!is.na(ts$TPA_UNADJ) & ts$TPA_UNADJ > 0, ts$TPA_UNADJ, ifelse(ts$DIA < 5, 74.965, 6.018))
ts$EXPF <- ts$tpa * ACRES_PER_HA; ts$DBH <- ts$DIA * 2.54
ts$HT <- ifelse(!is.na(ts$HT) & ts$HT > 0, ts$HT * 0.3048, pmax(2, 1.3 + 30 * (1 - exp(-0.04 * ts$DBH))))
ts$HCB <- NA_real_; ts$YEAR <- 2020; ts$dDBH.mult <- 1; ts$dHt.mult <- 1; ts$mort.mult <- 1
ts$max.dbh <- 200; ts$max.height <- 50; ts$Form <- NA; ts$Risk <- NA
ts$STAND <- sid; ts$PLOT <- 1L; ts <- ts[order(ts$STAND), ]; ts$TREE <- seq_len(nrow(ts))
sub <- ts[, c("STAND","PLOT","TREE","SP","DBH","HT","HCB","EXPF","YEAR","dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height","Form","Risk")]

st <- list(STAND = sid, CSI = 14, ELEV = 200)
ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 3.0, CutPoint = 0, MORTCAL = TRUE, MORTCAL_INTERVAL = 5)
out <- AcadianGYOneStand(sub, stand = st, ops = ops)
cat("nrow_in:", nrow(sub), " nrow_out:", nrow(out), " delta:", nrow(out)-nrow(sub), "\n")
cat("unique STAND values in OUTPUT frame:\n"); print(table(out$STAND))
cat("\nrecruits (DBH near 3 cm):\n")
recruits <- out[out$DBH <= 3.5 & !out$TREE %in% sub$TREE, ]
print(head(recruits[, c("STAND","PLOT","TREE","SP","DBH","HT","EXPF")], 12))
