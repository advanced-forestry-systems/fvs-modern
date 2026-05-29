## v39_magplot_guardrail.R
## Empirical guardrail test: apply v37 density correction to Canadian MAGPlot
## even though we said not to. Expected: bias goes from +0.4% (clean) to
## something worse. Confirms the "do not apply" recommendation is correct.
## Reuses the v17 MAGPlot harness with ingrowth_only config (best Canadian).
suppressMessages({ library(dplyr); library(plyr); library(purrr) })

PROJECT_ROOT <- "/users/PUOM0008/crsfaaron"
IN <- "/users/PUOM0008/crsfaaron/magplot/verify"
OUT_DIR <- IN
n_years <- 10L

source(file.path(PROJECT_ROOT, "AcadianGY_12.3.9.r"))
source(file.path(PROJECT_ROOT, "acadgy_fia_verify/apply_density_correction.R"))
cat("[v39-magplot] model:", AcadianVersionTag, "\n")
cat(sprintf("[v39-magplot] v37 coefficients: a=%.4f b=%.6f cap=(%.0f,%.0f) fit on n=%d\n",
            ACD_DENSITY_CORRECTION$a, ACD_DENSITY_CORRECTION$b,
            ACD_DENSITY_CORRECTION$lower_cap, ACD_DENSITY_CORRECTION$upper_cap,
            ACD_DENSITY_CORRECTION$n))

ti <- read.csv(file.path(IN, "magplot_tree_init.csv"), stringsAsFactors=FALSE)
pr <- read.csv(file.path(IN, "magplot_pairs.csv"), stringsAsFactors=FALSE)
ti$STAND <- as.character(ti$STAND); pr$STAND <- as.character(pr$STAND)

# Use 50-stand subset for speed
set.seed(2030)
sample_stands <- sample(unique(pr$STAND), min(50, length(unique(pr$STAND))))
ti <- ti[ti$STAND %in% sample_stands, ]
pr <- pr[pr$STAND %in% sample_stands, ]
cat(sprintf("[v39-magplot] %d trees, %d MAGPlot NB pairs (subset of 50 stands)\n", nrow(ti), nrow(pr)))

tree_init <- data.frame(STAND=ti$STAND, PLOT=1L, TREE=ti$TREE, SP=ti$SP,
  DBH=ti$DBH, HT=ti$HT, HCB=NA_real_, EXPF=ti$EXPF, YEAR=ti$YEAR,
  dDBH.mult=1, dHt.mult=1, mort.mult=1, max.dbh=200, max.height=50,
  Form=NA, Risk=NA, stringsAsFactors=FALSE)
stand_init <- data.frame(STAND=pr$STAND, CSI=14.0,
  ELEV=ifelse(!is.na(pr$ELEV), pr$ELEV, 200), stringsAsFactors=FALSE)

p1y <- function(trees, ops) {
  pc <- list()
  for (sid in unique(trees$STAND)) {
    st <- as.list(subset(stand_init, STAND==sid)); sub <- trees[trees$STAND==sid, ]
    if (nrow(sub)==0) next
    prj <- tryCatch(AcadianGYOneStand(sub, stand=st, ops=ops), error=function(e) NULL)
    if (!is.null(prj)) pc[[sid]] <- prj
  }
  if (length(pc)==0) return(NULL); dplyr::bind_rows(pc)
}
ba_m2 <- function(df) tapply((df$DBH^2)*0.00007854*df$EXPF, df$STAND, sum, na.rm=TRUE)

# Run optimal Canadian config: MORTCAL off + CutPoint = 0
ops <- list(verbose=FALSE, INGROWTH="Y", MinDBH=3.0, CutPoint=0)
cur <- tree_init
for (yr in 1:n_years) {
  cur <- p1y(cur, ops); if (is.null(cur)) break
}
ba_pred_m2ha <- ba_m2(cur)
stands <- names(ba_pred_m2ha)

# MAGPlot BA_t1 and BA_t2 are in m^2/ha. Pull obs and BA_t1 per stand.
obs <- setNames(pr$BA_t2_obs, pr$STAND)[stands]
ba_t1_m2ha <- setNames(pr$BA_t1_obs, pr$STAND)[stands]

# v37 was fit on ft^2/ac. Convert MAGPlot to ft^2/ac, apply correction, convert back.
M2HA_TO_FT2AC <- 4.356
ba_pred_ft2ac <- as.numeric(ba_pred_m2ha) * M2HA_TO_FT2AC
ba_t1_ft2ac   <- as.numeric(ba_t1_m2ha)   * M2HA_TO_FT2AC

ba_pred_corr_ft2ac <- apply_density_correction(ba_pred_ft2ac, ba_t1_ft2ac)
ba_pred_corr_m2ha  <- ba_pred_corr_ft2ac / M2HA_TO_FT2AC

bias_pct <- function(o, p) 100*(mean(p, na.rm=TRUE) - mean(o, na.rm=TRUE))/mean(o, na.rm=TRUE)
r2 <- function(o, p) {
  ok <- is.finite(o) & is.finite(p)
  1 - sum((p[ok]-o[ok])^2)/sum((o[ok]-mean(o[ok]))^2)
}

cat("\n=== v39 MAGPlot guardrail test ===\n")
cat(sprintf("Sample: %d Canadian (NB) stands, MORTCAL off + CutPoint=0 (production Canadian)\n",
            length(stands)))
cat(sprintf("BA range: obs %.1f-%.1f m^2/ha (mean %.1f)\n",
            min(obs, na.rm=TRUE), max(obs, na.rm=TRUE), mean(obs, na.rm=TRUE)))
cat(sprintf("BA_t1 range: %.1f-%.1f m^2/ha (mean %.1f)\n",
            min(ba_t1_m2ha, na.rm=TRUE), max(ba_t1_m2ha, na.rm=TRUE), mean(ba_t1_m2ha, na.rm=TRUE)))

cat("\n=== Effect of applying v37 correction to Canadian (should be WORSE) ===\n")
cat(sprintf("Uncorrected (production Canadian):  BA bias = %+.2f%%  R^2 = %.4f\n",
            bias_pct(obs, ba_pred_m2ha), r2(obs, ba_pred_m2ha)))
cat(sprintf("v37 correction applied (WRONG):     BA bias = %+.2f%%  R^2 = %.4f\n",
            bias_pct(obs, ba_pred_corr_m2ha), r2(obs, ba_pred_corr_m2ha)))

delta_bias <- bias_pct(obs, ba_pred_corr_m2ha) - bias_pct(obs, ba_pred_m2ha)
delta_r2 <- r2(obs, ba_pred_corr_m2ha) - r2(obs, ba_pred_m2ha)
cat(sprintf("\nDelta: bias %+.2f pp, R^2 %+.4f\n", delta_bias, delta_r2))
if (abs(delta_bias) > 1 || delta_r2 < -0.01) {
  cat("RESULT: v37 correction degrades Canadian fit. Guardrail confirmed.\n")
} else {
  cat("RESULT: v37 correction does not meaningfully change Canadian fit. Unexpected.\n")
}

write.csv(data.frame(STAND=stands, BA_t1_m2ha=ba_t1_m2ha, BA_obs=as.numeric(obs),
                     BA_pred_raw=as.numeric(ba_pred_m2ha),
                     BA_pred_v37corr=as.numeric(ba_pred_corr_m2ha)),
          file.path(OUT_DIR, "magplot_guardrail_v39_perstand.csv"), row.names=FALSE)
cat("\nDone. Per-stand CSV at magplot_guardrail_v39_perstand.csv\n")
