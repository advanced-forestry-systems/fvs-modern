## run_ingrowth_probe.R - 10 FIA stands, ONE annual cycle, CutPoint=0 +
## MORTCAL on, against AcadianGY_12.3.7_dbg.r. The instrumented model prints
## NA counts for every Ingrowth.FUN input covariate per stand. This pins down
## which covariate is silently dropping IPH to NA on FIA stands.
#!/usr/bin/env Rscript
suppressMessages({ library(dplyr); library(plyr); library(purrr) })

PROJECT_ROOT <- "/users/PUOM0008/crsfaaron"; N_PLOTS <- 10L; set.seed(42)
OUT_DIR <- file.path(PROJECT_ROOT, "acadgy_fia_verify", "ingrowth_probe_out")
dir.create(OUT_DIR, showWarnings = FALSE)
ACRES_PER_HA <- 2.4710538147

source("/users/PUOM0008/crsfaaron/AcadianGY_12.3.7_dbg.r")
cat("[probe] model:", AcadianVersionTag, "\n")

## Stand selection mirrors v20 (validation_data_acd_post + ME_PLOT/TREE)
vdat <- read.csv(file.path(PROJECT_ROOT,
  "fvs-modern/calibration/output/comparisons/intermediate/validation_data_acd_post.csv"))
me_plot <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_PLOT.csv"),
  colClasses = c("CN" = "character"))
me_tree <- read.csv(file.path(PROJECT_ROOT, "fia_data/ME_TREE.csv"),
  colClasses = c("CN" = "character", "PLT_CN" = "character"))
xwk <- read.csv(file.path(PROJECT_ROOT,
  "fvs-modern/calibration/data/osm_acadian_species_crosswalk.csv"))

acadgy_species <- c("AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
  "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM","RN","RO","RP",
  "RS","SB","SM","ST","SW","TA","WA","WC","WP","WS","YB")
xf <- xwk[!is.na(xwk$FIA) & xwk$FIA > 0 & nchar(xwk$OSM_AD_CmdKey) > 0, ]
xf$AGY_SP <- ifelse(xf$OSM_AD_CmdKey %in% acadgy_species, xf$OSM_AD_CmdKey,
  ifelse(xf$FIA >= 300, "OH", "OS"))
spcd_to_code <- setNames(xf$AGY_SP, as.character(xf$FIA))

vdat$PLT_CN_t1_str <- as.character(vdat$PLT_CN_t1)
vme <- vdat[vdat$PLT_CN_t1_str %in% as.character(me_plot$CN) &
            !is.na(vdat$interval_years) & vdat$interval_years %in% 5:10, ]
samp <- vme[sample(nrow(vme), min(N_PLOTS, nrow(vme))), ]
samp$PLT_CN_t1_str <- as.character(samp$PLT_CN_t1)

me_tree$PLT_CN <- as.character(me_tree$PLT_CN)
ts <- me_tree[me_tree$PLT_CN %in% samp$PLT_CN_t1_str &
              !is.na(me_tree$DIA) & me_tree$DIA > 0 & me_tree$STATUSCD == 1, ]
ts$SP <- spcd_to_code[as.character(ts$SPCD)]
mm <- is.na(ts$SP) | nchar(ts$SP) == 0
ts$SP[mm] <- ifelse(ts$SPCD[mm] >= 300, "OH", "OS")
ts$TPA_UNADJ <- as.numeric(ts$TPA_UNADJ)
ts$tpa <- ifelse(!is.na(ts$TPA_UNADJ) & ts$TPA_UNADJ > 0, ts$TPA_UNADJ,
  ifelse(ts$DIA < 5, 74.965, 6.018))
ts$EXPF <- ts$tpa * ACRES_PER_HA
ts$DBH <- ts$DIA * 2.54
ts$HT <- ifelse(!is.na(ts$HT) & ts$HT > 0, ts$HT * 0.3048, NA_real_)
mh <- is.na(ts$HT) | ts$HT == 0
ts$HT[mh] <- pmax(2, 1.3 + 30 * (1 - exp(-0.04 * ts$DBH[mh])))
ts$HCB <- NA_real_; ts$YEAR <- 2020
ts$dDBH.mult <- 1; ts$dHt.mult <- 1; ts$mort.mult <- 1
ts$max.dbh <- 200; ts$max.height <- 50; ts$Form <- NA; ts$Risk <- NA
ts$STAND <- ts$PLT_CN; ts$PLOT <- 1L
ts <- ts[order(ts$STAND), ]
ts$TREE <- ave(seq_along(ts$STAND), ts$STAND, FUN = seq_along)

cols <- c("STAND","PLOT","TREE","SP","DBH","HT","HCB","EXPF","YEAR",
  "dDBH.mult","dHt.mult","mort.mult","max.dbh","max.height","Form","Risk")
base_init <- ts[, cols]

samp$ClimateSI_m <- ifelse(!is.na(samp$ClimateSI_ft), samp$ClimateSI_ft * 0.3048, 14)
samp$ELEV_m <- ifelse(!is.na(samp$ELEV_t2), samp$ELEV_t2 * 0.3048, 200)
stand_init <- unique(data.frame(STAND = samp$PLT_CN_t1_str,
  CSI = samp$ClimateSI_m, ELEV = samp$ELEV_m, stringsAsFactors = FALSE))

ops <- list(verbose = FALSE, INGROWTH = "Y", MinDBH = 3.0, CutPoint = 0,
  MORTCAL = TRUE, MORTCAL_INTERVAL = 5)

## Run one annual cycle per plot; the [INGRWDBG] cat() inside AcadianGYOneStand
## will fire each call. Send everything to a log file so we can grep.
log_path <- file.path(OUT_DIR, "ingrowth_probe.log")
sink(log_path, split = TRUE)
cat(sprintf("[probe] plots=%d  trees=%d\n", length(unique(base_init$STAND)),
            nrow(base_init)))
for (sid in unique(base_init$STAND)) {
  st <- as.list(subset(stand_init, STAND == sid))
  sub <- base_init[base_init$STAND == sid, ]
  if (nrow(sub) == 0) next
  cat(sprintf("[probe] --- plot %s ntree=%d ---\n", sid, nrow(sub)))
  out <- tryCatch(AcadianGYOneStand(sub, stand = st, ops = ops),
                  error = function(e) e)
  if (inherits(out, "error")) {
    cat("[probe] ERR:", conditionMessage(out), "\n")
  } else {
    cat(sprintf("[probe] in=%d  out=%d  added=%d\n",
                nrow(sub), nrow(out), nrow(out) - nrow(sub)))
  }
}
sink()
cat("[probe] log written:", log_path, "\n")
