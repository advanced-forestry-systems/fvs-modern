#!/usr/bin/env Rscript
# =============================================================================
# 75_fit_stand_topheight.R -- STAND-LEVEL García/GADA TOP-HEIGHT TRANSITION.
#
# gada_refit.r (cspi_v3) fits the STATIC Cieszewski-Bailey base-age-invariant
# Chapman-Richards site form  H = b1*(1-exp(-b2*t))^b3  (site via b1), giving a
# per-species (b2,b3) and per-plot SI50. That is a SITE index, NOT a state-space
# H2|H1 transition. This script fits the minimal Bayesian STATE-SPACE transition
# of stand TOP height (mean HT of the tallest ~top_n_per_ha stems) over the
# remeasurement interval, so the tree height-growth constraint in
# stand_constraint.py::stand_constrain_topheight has a fitted H(t) TARGET:
#
#   r = ln(H2 / H1)              (log height-growth ratio, >= 0)
#   r ~ Normal(eta_r, sigma)
#   eta_r = b0 + b_lnH1*ln(H1) + b_lnyr*ln(YEARS) + b_rd*rd + b_lnqmd*ln_qmd
#           + b_bgi*bgi + (1 | EPA_L1)
#   H2_target = H1 * exp(eta_r)
#
# This is the García state-space form: next state H2 is a function of the current
# state H1 and the interval. ln(H1) and ln(YEARS) capture the Chapman-Richards
# deceleration; rd/qmd/bgi carry density and site. The base-age invariance is
# preserved implicitly (no absolute age needed, which is unreliable in FIA).
#
# STAND TOP HEIGHT: the pairs RDS carries HT40_1 (start top height, tallest 40
# stems basis) but NO t2 stand top-height column. We RECOMPUTE both endpoints
# from the tree HT1 / HT2 records: per stand-remeasurement, top height = mean HT
# of the TPA-accumulated tallest `top_n_per_ha` stems/ha (matches
# stand_constrain_topheight). t2 uses SURVIVORS' HT2.
#
# Usage:
#   Rscript 75_fit_stand_topheight.R --n_sub=120000 --top_n=100 \
#     --pairs=<rds> --out_dir=<dir>
#
# Author: A. Weiskittel + Claude (OODA autopilot)  Date: 2026-07-02
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(jsonlite) })
set.seed(20260702)
ELOG <- "error_log.txt"
elog <- function(msg) cat(sprintf("[%s] %s\n", format(Sys.time()), msg), file = ELOG, append = TRUE)
trycatch_run <- function(expr, what) tryCatch(expr, error = function(e) { elog(paste(what, ":", conditionMessage(e))); NULL })

args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d = NULL) { m <- grep(paste0("^--", n, "="), args, value = TRUE); if (length(m)) sub(paste0("^--", n, "="), "", m[1]) else d }
PAIRS <- ga("pairs", "/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds")
OUT_DIR <- ga("out_dir", "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus/stand_topheight")
N_SUB <- as.integer(ga("n_sub", "120000"))
TOP_N_PER_HA <- as.numeric(ga("top_n", "100"))   # tallest stems/ha defining top height
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("== 75_fit_stand_topheight.R ==\n top_n_per_ha:", TOP_N_PER_HA, " n_sub:", N_SUB, "\n")

d <- trycatch_run(as.data.table(readRDS(PAIRS)), "read pairs")
if (is.null(d)) quit(status = 1)
cat(" pairs rows:", nrow(d), "\n")

# ---- row filter: SURVIVORS with valid start/end tree heights -----------------
# Top height is defined by the DOMINANT cohort; both endpoints need HT1 & HT2.
# Restrict to survivors (TREESTATUS2==1) with finite positive heights so the same
# stems define H1 and H2 (a state-space transition of the surviving canopy).
if (!"TREESTATUS1" %in% names(d)) d[, TREESTATUS1 := 1L]
need <- c("HT1", "HT2", "TPH_UNADJ1", "YEARS", "PLT_CN_cond1", "CONDID_cond1",
          "INVYR1", "INVYR2", "EPA_L1_CODE")
miss <- setdiff(need, names(d))
if (length(miss)) { elog(paste("75 missing cols:", paste(miss, collapse=","))); cat(" MISSING COLS:", paste(miss, collapse=", "), "\n") }
# TPH weight: prefer per-tree TPH_UNADJ1 (stems/ha); fall back to TPA1 * factor
if (!"TPH_UNADJ1" %in% names(d) && "TPA1" %in% names(d)) d[, TPH_UNADJ1 := TPA1 * 2.4710538147]
d <- d[TREESTATUS1 == 1 & !is.na(TREESTATUS2) & TREESTATUS2 == 1 &
       is.finite(HT1) & HT1 > 1.3 & is.finite(HT2) & HT2 > 1.3 &
       is.finite(TPH_UNADJ1) & TPH_UNADJ1 > 0 &
       is.finite(YEARS) & YEARS >= 1 & YEARS <= 20]
cat(" survivor tree rows with heights:", nrow(d), "\n")

d[, stand_key := paste(PLT_CN_cond1, CONDID_cond1, INVYR1, INVYR2, sep = "_")]

# ---- stand TOP height at t1 and t2 (tallest top_n_per_ha stems/ha) -----------
# For each stand key, rank surviving trees by HT (t1 for H1, t2 for H2),
# accumulate TPH until top_n_per_ha stems included (fractional boundary stem),
# and take the TPH-weighted mean height. Mirrors stand_constrain_topheight.
top_ht <- function(ht, w, need) {
  o <- order(-ht); w <- w[o]; ht <- ht[o]
  cum <- cumsum(w)
  if (cum[length(cum)] <= need) return(sum(w * ht) / sum(w))
  full <- cum <= need
  wt <- ifelse(full, w, 0)
  k <- which(!full)[1]
  prev <- if (k > 1) cum[k - 1] else 0
  wt[k] <- need - prev
  sum(wt * ht) / sum(wt)
}

agg <- d[, {
  H1 <- top_ht(HT1, TPH_UNADJ1, TOP_N_PER_HA)
  H2 <- top_ht(HT2, TPH_UNADJ1, TOP_N_PER_HA)   # survivors' end heights, same weights
  .(YEARS = YEARS[1], H1 = H1, H2 = H2,
    qmd  = if ("QMD1" %in% names(.SD)) weighted.mean(QMD1, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    ba   = if ("BA1"  %in% names(.SD)) weighted.mean(BA1,  TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    sdi_add = if ("sdi_additive1" %in% names(.SD)) weighted.mean(sdi_additive1, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    sdimax  = if ("SDImax_brms" %in% names(.SD)) weighted.mean(SDImax_brms, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    rd_row  = if ("rd_add" %in% names(.SD)) weighted.mean(rd_add, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    bgi  = if ("bgi" %in% names(.SD)) weighted.mean(bgi, TPH_UNADJ1, na.rm = TRUE) else NA_real_,
    L1   = EPA_L1_CODE[1])
}, by = stand_key,
   .SDcols = intersect(c("QMD1","BA1","sdi_additive1","SDImax_brms","rd_add","bgi"), names(d))]
rm(d); gc()
cat(" stand-remeasurement top-height records:", nrow(agg), "\n")

# ---- covariates + response ---------------------------------------------------
agg[, rd := ifelse(is.finite(sdi_add) & is.finite(sdimax) & sdimax > 0, sdi_add / sdimax, rd_row)]
agg <- agg[is.finite(H1) & H1 > 1.3 & is.finite(H2) & H2 > 1.3 & is.finite(YEARS)]
# top height should not DROP for survivors; clip tiny negative measurement noise.
agg[H2 < H1, H2 := H1]
agg[, r := log(H2 / H1)]                          # log height-growth ratio (>= 0)
agg[, ln_h1 := log(H1)]
agg[, ln_years := log(YEARS)]
if (!"qmd" %in% names(agg)) agg[, qmd := NA_real_]
agg[, ln_qmd := ifelse(is.finite(qmd) & qmd > 0, log(qmd), NA_real_)]
bgi_med <- suppressWarnings(median(agg$bgi[is.finite(agg$bgi)], na.rm = TRUE))
if (!is.finite(bgi_med)) bgi_med <- 0
agg[!is.finite(bgi), bgi := bgi_med]
agg[!is.finite(rd), rd := suppressWarnings(median(agg$rd[is.finite(agg$rd)], na.rm = TRUE))]
agg[!is.finite(ln_qmd), ln_qmd := suppressWarnings(median(agg$ln_qmd[is.finite(agg$ln_qmd)], na.rm = TRUE))]
agg[, L1 := as.character(L1)]
agg <- agg[is.finite(r) & r >= 0 & r < 1.5 & !is.na(L1) & L1 != "" &
           is.finite(ln_h1) & is.finite(ln_years) & is.finite(rd) & is.finite(ln_qmd)]
cat(" modeling-ready records:", nrow(agg), "\n")
cat(" mean log height-growth ratio r:", round(mean(agg$r), 5),
    " mean H1:", round(mean(agg$H1), 2), " mean H2:", round(mean(agg$H2), 2), "\n")

# ---- subsample ---------------------------------------------------------------
n_use <- min(nrow(agg), N_SUB)
ds <- agg[sample(.N, n_use)]
cat(sprintf(" fit sample: %d of %d\n", nrow(ds), nrow(agg)))
rm(agg); gc()

# ---- Bayesian fit (brms): Gaussian on the log height-growth ratio ------------
ok <- trycatch_run({ suppressPackageStartupMessages(library(brms)); TRUE }, "load brms")
if (is.null(ok)) quit(status = 1)

form <- bf(r ~ ln_h1 + ln_years + rd + ln_qmd + bgi + (1 | L1))
priors <- c(set_prior("normal(0,1)", class = "b"),
            set_prior("normal(0,0.5)", class = "sd"),
            set_prior("normal(0,0.5)", class = "sigma"))
fit <- trycatch_run(
  brm(form, data = ds, family = gaussian(), prior = priors,
      chains = 4, iter = 800, warmup = 400, cores = 4, seed = 20260702,
      refresh = 100, control = list(adapt_delta = 0.9)),
  "brms stand-topheight fit")
if (is.null(fit)) quit(status = 1)

fx <- fixef(fit)
re <- trycatch_run(ranef(fit)$L1[, , "Intercept"], "ranef L1")
sm <- as.data.frame(fx); sm$param <- rownames(sm)
fwrite(as.data.table(sm), file.path(OUT_DIR, "stand_topheight_fixed.csv"))

fe_names <- rownames(fx)
name_map <- c("Intercept"="Intercept","ln_h1"="ln_h1","ln_years"="ln_years",
              "rd"="rd","ln_qmd"="ln_qmd","bgi"="bgi")
fixed_list <- setNames(
  lapply(fe_names, function(nm) list(mean = unname(fx[nm,"Estimate"]),
                                     sd = unname(fx[nm,"Est.Error"]),
                                     q2.5 = unname(fx[nm,"Q2.5"]),
                                     q97.5 = unname(fx[nm,"Q97.5"]))),
  fe_names)
re_table <- if (!is.null(re)) list(level = rownames(re), mean = unname(re[,"Estimate"]),
                                   sd = unname(re[,"Est.Error"])) else NULL
sd_L1 <- trycatch_run({ vc <- summary(fit)$random$L1; unname(vc["sd(Intercept)","Estimate"]) }, "sd L1")

bundle <- list(
  model = "stand-level García/GADA state-space TOP-HEIGHT transition",
  scale = "Gaussian on r = ln(H2/H1) (log height-growth ratio); H2 = H1*exp(eta_r)",
  response = "r = ln(H2/H1); H1,H2 = TPH-weighted mean height of tallest top_n_per_ha stems/ha (survivors)",
  stand_key = "PLT_CN_cond1 x CONDID_cond1 x INVYR1 x INVYR2 (FIA plot-condition-remeasurement)",
  top_n_per_ha = TOP_N_PER_HA,
  transition = "H2_target = H1 * exp(b0 + b_lnH1*ln(H1) + b_lnyr*ln(YEARS) + b_rd*rd + b_lnqmd*ln_qmd + b_bgi*bgi + z_L1)",
  gada_static_form = "Cieszewski-Bailey base-age-invariant Chapman-Richards H=b1*(1-exp(-b2*t))^b3 (gada_refit.r; site via b1); this bundle is the STATE-SPACE H2|H1 transition, not the static site form",
  covariates = list(
    ln_h1 = "log start stand top height (m)",
    ln_years = "log remeasurement interval (years) -- Chapman-Richards deceleration",
    rd = "stand SDI / SDIMAX (sdi_additive1 / SDImax_brms)",
    ln_qmd = "log stand QMD (cm)",
    bgi = "stand mean bgi (climate/site driver)"),
  n_stand_records = nrow(ds),
  fixed_effects = fixed_list,
  sd_L1 = sd_L1,
  re_L1 = re_table,
  seed = 20260702,
  notes = "H2 target feeds stand_constraint.py::stand_constrain_topheight to scale tree height growth so stand top height tracks this García/GADA trajectory."
)
write_json(bundle, file.path(OUT_DIR, "stand_topheight_bundle.json"),
           auto_unbox = TRUE, pretty = TRUE, digits = 8, null = "null")

cat("\nStand top-height transition fixed effects (log height-growth ratio):\n")
print(round(fx[, c("Estimate","Est.Error","Q2.5","Q97.5"), drop = FALSE], 5))
cat("\nsd(1|L1):", round(sd_L1, 5), "\n")
gc(); cat("\nDONE_STAND_TOPHEIGHT\n")
