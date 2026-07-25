
## ============================================================================
## stress_test_components.R
##
## Per-tree stress test for all 5 FVS NE component equations.
## Verifies biological plausibility across realistic input ranges and
## compares FVS native (periodic) vs Greg/Aaron (annualized) equations.
##
## Key species (NE): balsam fir (SPCD=11), red spruce (SPCD=97),
##   eastern white pine (SPCD=129), sugar maple (SPCD=318),
##   yellow birch (SPCD=371)
##
## ANNUALIZATION NOTE (flagged throughout):
##   Greg's DG and HG equations are ANNUAL (in/yr, ft/yr).
##   FVS default DG/HG are 5-YEAR periodic.
##   If FVS calls Greg's equations once per 5-year cycle and uses the
##   single annual increment, the 5-year projection equals 1 annual
##   increment (1/5th expected growth). This script runs Greg equations
##   ANNUALLY (×1/yr) in the multi-year simulation, which is the
##   correct application for annualized forms.
##
## Model forms (verified against validate_gregdghg.R):
##   HT-DBH : HT = 4.5 + B1*(1-exp(-B2*DBH))^B3  (Chapman-Richards, ft/in)
##   DG      : exp(B0 + B1*ln((DBH+1)^2/(CR*HT+1)^B3) +
##              B2*BAL^B4/ln(DBH+2.7) + B5*ELEV_ft + B6*EMT_C)  [in/yr]
##   HG      : B0*B1*B2*CR^B3 * exp(-B1*HT - B4*CCFL - B8*CCH^0.5
##              - B5*ELEV + B6*sqrt(TD) + B7*EMT) * (1-exp(-B1*HT))^(B2-1)
##              where B0 = asymptotic max height  [ft/yr]
##   CRW     : DHTLC = CL*(1 - exp(b0 + b1*DHTANN + b2*CCH))  [Fortran form]
##              CL=HT*CR (crown length ft), DHTANN=annual ht increment (ft/yr),
##              CCH=crown-closure fraction (0-1, same scale as DG/HG/mort),
##              output=annual
##              change in ht-to-live-crown (ft/yr; positive = recession upward)
##   MORT    : eta = b0 + b1*DBH_cm + b2*DBH_cm^2 + b3*CR
##              + b4*BGI_m + b5*ln(BAL_m2ha+5) + b6*sqrt(BA_m2ha*RD)
##              P_survive_1yr = exp(-exp(eta))
##
## Author: A. Weiskittel / G. Johnson (stress test scaffold)
## Date  : 2026-07-11
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

## ── Paths ──────────────────────────────────────────────────────────────────
COEF_DIR_HTDBH <- "/fs/scratch/PUOM0008/crsfaaron/wt-htdbh-ser/config"
COEF_DIR_NE    <- "/fs/scratch/PUOM0008/crsfaaron/wt-ne-dg/config"
OUT_DIR        <- "/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/stress_test"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## ── Species table ──────────────────────────────────────────────────────────
SPECIES <- data.frame(
  SPCD      = c(11,   97,   129,  318,  371),
  JSP       = c("BF", "RS", "WP", "SM", "YB"),
  NAME      = c("Balsam fir", "Red spruce", "Eastern white pine",
                "Sugar maple", "Yellow birch"),
  MAX_HT_FT = c(80,   85,   140,  100,  95),
  stringsAsFactors = FALSE
)

## ── NE environmental defaults ───────────────────────────────────────────────
NE_ELEV_FT <- 1000   # elevation in feet (Maine mid-elevation)
NE_EMT_C   <- -15    # extreme minimum temperature (°C), Maine winters
NE_TD_F    <- 22     # temperature differential (°F seasonal range)
NE_CCFL    <- 80     # crown competition factor from larger trees (%)
NE_CCH     <- 0.35   # canopy cover fraction (0-1), used by HG
NE_CCH_CRW <- 0.35   # crown-closure fraction (0-1), used by CRW (same CCH scale as Greg DG/HG/mort)
NE_BA      <- 80     # stand BA ft2/ac
NE_RD      <- 0.60   # relative density (Curtis RD)
NE_BGI     <- 0.45   # basal growth index (m, site productivity proxy)

## Unit converters
FT2AC_TO_M2HA <- 0.22957
IN_TO_CM      <- 2.54

cat("Loading coefficient files...\n")

co_htdbh <- read.csv(file.path(COEF_DIR_HTDBH, "conus_htdbh_coefficients.csv"))
co_dg    <- read.csv(file.path(COEF_DIR_NE, "greg_dg_coefficients.csv"),
                     comment.char = "#")
co_hg    <- read.csv(file.path(COEF_DIR_NE, "greg_hg_coefficients.csv"))
co_hgcn  <- read.csv(file.path(COEF_DIR_NE, "greg_hg_coefficients_compound_none.csv"))
co_crw   <- read.csv(file.path(COEF_DIR_NE, "greg_crown_change_coefficients.csv"))
co_mort  <- read.csv(file.path(COEF_DIR_NE, "greg_mortality_coefficients_size_bgi.csv"))

## Filter to key species
sp5 <- SPECIES$SPCD
co_htdbh <- co_htdbh[co_htdbh$SPCD %in% sp5, ]
co_dg    <- co_dg[co_dg$SPCD    %in% sp5, ]
co_hg    <- co_hg[co_hg$SPCD    %in% sp5, ]
co_hgcn  <- co_hgcn[co_hgcn$SPCD %in% sp5, ]
co_crw   <- co_crw[co_crw$SPCD  %in% sp5, ]
co_mort  <- co_mort[co_mort$SPCD %in% sp5, ]

cat(sprintf("Loaded: htdbh=%d  dg=%d  hg=%d  hgcn=%d  crw=%d  mort=%d species\n",
            nrow(co_htdbh), nrow(co_dg), nrow(co_hg),
            nrow(co_hgcn), nrow(co_crw), nrow(co_mort)))

## ── FVS NE approximate default HT-DBH (reference line) ────────────────────
## Chapman-Richards approximation for NE species using typical asymptote/
## rate values. These are representative reference values, not official
## FVS NE table coefficients.
fvs_ne_htdbh <- data.frame(
  SPCD   = c(11,    97,    129,   318,   371),
  a      = c(75,    80,    138,   98,    92),   # asymptote ft
  b      = c(0.100, 0.110, 0.080, 0.090, 0.090),# rate
  c      = c(1.20,  1.15,  1.25,  1.20,  1.18)  # shape
)

## ── Section 1: HT-DBH stress test ──────────────────────────────────────────
cat("\n=== Section 1: HT-DBH ===\n")

pred_htdbh_aaron <- function(coef_row, dbh) {
  B1 <- coef_row$B1; B2 <- coef_row$B2; B3 <- coef_row$B3
  bh <- 4.5
  ## CONUS best-AIC HT-DBH forms are dispatched per species by model_id
  ## (Marshall fits; a single form does not fit all species, and the sign is
  ## already carried in the fitted coefficients). Forms (bh = 4.5 ft):
  ##   1: HT = bh + exp(B1 + B2*DBH^-1)
  ##   2: HT = bh + exp(B1 + B2*DBH^B3)
  ##   3: HT = bh + B1*(1 - exp(B2*DBH))^B3
  ##   4: HT = bh + exp(B1 + B2/(DBH+1))
  ##   5: HT = bh + exp(B1 + B2/(DBH+B3))
  ##   6: HT = bh + B1*(1 - exp(B2*DBH^B3))
  ht <- switch(as.character(coef_row$model_id),
    "1" = bh + exp(B1 + B2 * dbh^(-1.0)),
    "2" = bh + exp(B1 + B2 * dbh^B3),
    "3" = bh + B1 * (1 - exp(B2 * dbh))^B3,
    "4" = bh + exp(B1 + B2 / (dbh + 1.0)),
    "5" = bh + exp(B1 + B2 / (dbh + B3)),
    "6" = bh + B1 * (1 - exp(B2 * dbh^B3)),
    rep(NA_real_, length(dbh)))
  pmax(ht, bh)
}

pred_htdbh_fvs <- function(ref_row, dbh) {
  pmax(4.5 + ref_row$a * (1 - exp(-ref_row$b * dbh))^ref_row$c, 4.5)
}

dbh_seq <- seq(1, 40, by = 0.5)

ht_preds <- lapply(sp5, function(sp) {
  cr <- co_htdbh[co_htdbh$SPCD == sp, ]
  rv <- fvs_ne_htdbh[fvs_ne_htdbh$SPCD == sp, ]
  nm <- SPECIES$NAME[SPECIES$SPCD == sp]
  ht_a <- pred_htdbh_aaron(cr, dbh_seq)
  ht_f <- pred_htdbh_fvs(rv, dbh_seq)
  ## Monotonicity check
  mono_a <- all(diff(ht_a) >= 0)
  mono_f <- all(diff(ht_f) >= 0)
  asym_a <- max(ht_a)
  cat(sprintf("  %s (SPCD %d): Aaron asymptote=%.0f ft | mono=%s | FVS ref=%.0f ft\n",
              nm, sp, asym_a, mono_a, max(ht_f)))
  data.frame(DBH=dbh_seq, HT_aaron=ht_a, HT_fvs=ht_f, SPCD=sp, NAME=nm)
}) %>% bind_rows()

p1 <- ht_preds %>%
  pivot_longer(c(HT_aaron, HT_fvs), names_to="Model", values_to="HT") %>%
  mutate(Model = recode(Model,
    HT_aaron = "Aaron CONUS (best-AIC HT-DBH)",
    HT_fvs   = "FVS NE reference")) %>%
  ggplot(aes(x=DBH, y=HT, color=Model, linetype=Model)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~NAME, scales="free_y", ncol=3) +
  scale_color_manual(values=c("steelblue","firebrick")) +
  scale_linetype_manual(values=c("solid","dashed")) +
  geom_hline(yintercept=4.5, linetype="dotted", color="grey50", linewidth=0.4) +
  labs(title="HT-DBH stress test: Aaron CONUS vs FVS NE reference",
       subtitle="CONUS best-AIC HT-DBH forms dispatched per species by model_id (Curtis-Arney and Wykoff family)",
       x="DBH (in)", y="Total height (ft)", color=NULL, linetype=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", strip.text=element_text(size=9))

ggsave(file.path(OUT_DIR, "stress_ht_dbh.png"), p1,
       width=10, height=6.5, dpi=150)
cat("  -> stress_ht_dbh.png saved\n")

## ── Section 2: DG stress test ──────────────────────────────────────────────
cat("\n=== Section 2: DG ===\n")
## ANNUALIZATION FLAG: Greg's DG coefficients are fit to ANNUAL diameter
## growth (in/yr). FVS defaults produce 5-year periodic DG. If FVS calls
## Greg's equation once per cycle, output is 1 annual increment applied
## as if it were a 5-year total — a 5x underestimate.

pred_dg_greg <- function(cr_row, dbh, cr, ht, bal,
                         elev=NE_ELEV_FT, emt=NE_EMT_C) {
  B0 <- cr_row$B0; B1 <- cr_row$B1; B2 <- cr_row$B2
  B3 <- cr_row$B3; B4 <- cr_row$B4; B5 <- cr_row$B5; B6 <- cr_row$B6
  cr  <- pmax(cr, 1e-4)
  ht  <- pmax(ht, 1)
  bal <- pmax(bal, 0)
  z <- B0 + B1 * log((dbh+1)^2 / (cr*ht+1)^B3) +
       B2 * bal^B4 / log(dbh+2.7) +
       B5 * elev + B6 * emt
  z <- pmin(pmax(z, -30), 5)
  pmax(exp(z), 0)
}

dbh_grid <- c(2, 4, 8, 12, 20)
bal_grid <- c(0, 25, 50, 100, 150)
si_grid  <- c(50, 70, 90)
cr_grid  <- c(0.3, 0.5, 0.7)

# Reference plot: annual DG vs DBH at BAL=25 and BAL=100, CR=0.5, typical HT
dg_preds <- lapply(sp5, function(sp) {
  cr_row <- co_dg[co_dg$SPCD == sp, ]
  nm <- SPECIES$NAME[SPECIES$SPCD == sp]
  lapply(c(25, 100), function(bal) {
    ht_vals <- pmax(4.5 + 50*(1-exp(-0.09*dbh_seq)), 5)  # rough HT for context
    dg <- pred_dg_greg(cr_row, dbh_seq, cr=0.5, ht=ht_vals, bal=bal)
    data.frame(DBH=dbh_seq, DG=dg, BAL=as.character(bal), SPCD=sp, NAME=nm)
  }) %>% bind_rows()
}) %>% bind_rows()

# Check range
dg_range <- dg_preds %>% group_by(NAME) %>%
  summarise(dg_min=min(DG), dg_max=max(DG), dg_med=median(DG), .groups="drop")
cat("  DG range (in/yr, CR=0.5, NE elev/EMT defaults):\n")
print(dg_range)
cat("  Expected: 0.05-0.40 in/yr for typical NE conditions\n")

## BAL is a numeric level stored as character; set explicit factor levels so the
## colour/label mapping is not reordered alphabetically ("100" before "25"),
## which previously swapped the two curves' legend labels.
p2 <- dg_preds %>%
  mutate(BAL = factor(BAL, levels = c("25", "100"))) %>%
  ggplot(aes(x=DBH, y=DG, color=BAL, group=BAL)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~NAME, ncol=3, scales="free_y") +
  scale_color_manual(values=c("25"="steelblue", "100"="firebrick"),
                     labels=c("25"="BAL = 25 ft² ac⁻¹", "100"="BAL = 100 ft² ac⁻¹")) +
  geom_hline(yintercept=c(0.05, 0.40), linetype="dashed",
             color="grey40", linewidth=0.35) +
  labs(title="DG stress test: Greg ORGANON-PC form (ANNUAL, in/yr)",
       subtitle=paste0("ANNUALIZATION FLAG: If FVS calls this once per 5-yr cycle,\n",
         "result = 1 annual increment (1/5 expected 5-yr growth). CR=0.5, ELEV=",
         NE_ELEV_FT," ft, EMT=",NE_EMT_C,"°C"),
       x="DBH (in)", y="Annual DG (in yr⁻¹)",
       color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", strip.text=element_text(size=9))

ggsave(file.path(OUT_DIR, "stress_dg.png"), p2,
       width=10, height=6.5, dpi=150)
cat("  -> stress_dg.png saved\n")

## ── Section 3: HG stress test ──────────────────────────────────────────────
cat("\n=== Section 3: HG ===\n")
## ANNUALIZATION FLAG: same as DG — Greg HG coefficients are annual (ft/yr).

pred_hg <- function(cr_row, ht, cr,
                    ccfl=NE_CCFL, cch=NE_CCH,
                    elev=NE_ELEV_FT, td=NE_TD_F, emt=NE_EMT_C) {
  B0 <- cr_row$B0; B1 <- cr_row$B1; B2 <- cr_row$B2; B3 <- cr_row$B3
  B4 <- cr_row$B4; B5 <- cr_row$B5; B6 <- cr_row$B6
  B7 <- cr_row$B7; B8 <- cr_row$B8
  crp  <- pmax(cr, 1e-4)
  cchp <- pmax(cch, 0)
  tdp  <- pmax(td, 0)
  hg <- B0 * B1 * B2 * crp^B3 *
        exp(-B1*ht - B4*ccfl - B8*sqrt(cchp) - B5*elev +
            B6*sqrt(tdp) + B7*emt) *
        (1 - exp(-B1*ht))^(B2 - 1)
  pmax(hg, 0)
}

ht_seq <- seq(5, 120, by=2)

hg_preds <- lapply(sp5, function(sp) {
  cr_hg <- co_hg[co_hg$SPCD == sp, ]
  cr_cn <- co_hgcn[co_hgcn$SPCD == sp, ]
  nm <- SPECIES$NAME[SPECIES$SPCD == sp]
  lapply(c(0.3, 0.6), function(cr_val) {
    hg_org <- pred_hg(cr_hg, ht_seq, cr=cr_val)
    hg_cn  <- if (nrow(cr_cn) > 0) pred_hg(cr_cn, ht_seq, cr=cr_val) else rep(NA_real_, length(ht_seq))
    data.frame(HT=ht_seq, HG_organon=hg_org, HG_cn=hg_cn,
               CR=as.character(cr_val), SPCD=sp, NAME=nm)
  }) %>% bind_rows()
}) %>% bind_rows()

cat("  HG range (ft/yr):\n")
hg_range <- hg_preds %>% filter(!is.na(HG_organon)) %>%
  group_by(NAME) %>%
  summarise(min_o=round(min(HG_organon),3), max_o=round(max(HG_organon),3),
            min_cn=round(min(HG_cn, na.rm=TRUE),3),
            max_cn=round(max(HG_cn, na.rm=TRUE),3), .groups="drop")
print(hg_range)
cat("  Expected: 0.5-4 ft/yr for young trees, declining with size\n")

p3 <- hg_preds %>%
  pivot_longer(c(HG_organon, HG_cn), names_to="Model", values_to="HG") %>%
  mutate(Model = recode(Model,
    HG_organon = "Greg ORGANON",
    HG_cn      = "Aaron compound_none")) %>%
  filter(!is.na(HG)) %>%
  ggplot(aes(x=HT, y=HG, color=Model, linetype=as.character(CR))) +
  geom_line(linewidth=0.9) +
  facet_wrap(~NAME, ncol=3, scales="free_y") +
  scale_color_manual(values=c("darkgreen","steelblue")) +
  labs(title="HG stress test: Greg ORGANON vs Aaron compound_none (ANNUAL, ft/yr)",
       subtitle=paste0("ANNUALIZATION FLAG: annual equation. CCFL=",NE_CCFL,
         ", CCH=",NE_CCH,", ELEV=",NE_ELEV_FT," ft, EMT=",NE_EMT_C,"°C"),
       x="Height (ft)", y="Annual HG (ft yr⁻¹)",
       color="Model", linetype="Crown ratio") +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", strip.text=element_text(size=9))

ggsave(file.path(OUT_DIR, "stress_hg.png"), p3,
       width=10, height=6.5, dpi=150)
cat("  -> stress_hg.png saved\n")

## ── Section 4: Crown recession stress test ─────────────────────────────────
cat("\n=== Section 4: Crown recession (CRW) ===\n")
## Correct Fortran form (greghghg.f90 CROWNDRIVER):
##   DHTLC = CL * (1 - exp(b0 + b1*DHTANN + b2*CCH))
##   CL     = HT * CR                  (crown length, ft)
##   DHTANN = annual height increment   (ft/yr)
##   CCH    = crown-closure fraction (0-1), the same CCH the deployed Greg
##            DG/HG/mortality forms use (documented crown closure at tip, 0-1),
##            NOT a 50-150 CCF index; feeding CCH on a 0-100 scale drives the
##            exponential argument strongly negative and makes DHTLC approach the
##            whole crown length (tens of ft/yr), which is unphysical.
##   DHTLC  = annual change in ht-to-live-crown (ft/yr; +ve = recession upward)
## Previous form (b0 + b1*CR + b2*BAL_m2ha) was inferred from stan model and
## is WRONG; replaced here with the Fortran functional form.

pred_crw <- function(crw_row, ht, cr, dhtann, cch) {
  cl  <- ht * cr
  arg <- crw_row$b0 + crw_row$b1 * dhtann + crw_row$b2 * cch
  cl * (1 - exp(arg))
}

crw_preds <- lapply(sp5, function(sp) {
  crw_row <- co_crw[co_crw$SPCD == sp, ]
  nm      <- SPECIES$NAME[SPECIES$SPCD == sp]
  expand.grid(
    DBH    = seq(2, 24, by = 2),
    CR     = seq(0.2, 0.8, by = 0.1),
    DHTANN = c(0.5, 1.0, 2.0, 3.0),
    CCH    = c(0.3, 0.6, 0.9)
  ) %>%
    mutate(
      HT     = 4.5 + exp(log(80) + log(1 - exp(-0.09 * DBH))),  # rough ht-dbh
      HT     = pmax(HT, 5),
      CL     = HT * CR,
      ARG    = crw_row$b0 + crw_row$b1 * DHTANN + crw_row$b2 * CCH,
      DHTLC  = CL * (1 - exp(ARG)),       # annual delta in ht-to-live-crown (ft/yr)
      HLC    = HT * (1 - CR),             # current ht to live crown
      new_HLC = pmax(0, HLC + DHTLC),
      new_CR  = pmax(0.05, (HT - new_HLC) / HT),
      SPCD = sp, NAME = nm
    )
}) %>% bind_rows()

cat("  Annual DHTLC range (ft/yr; +ve = recession upward):\n")
crw_range <- crw_preds %>% group_by(NAME) %>%
  summarise(min_dhtlc = round(min(DHTLC), 3),
            max_dhtlc = round(max(DHTLC), 3),
            mean_dhtlc = round(mean(DHTLC), 3), .groups = "drop")
print(crw_range)
cat("  Expected: positive DHTLC (upward recession) under high CCH / fast growth\n")

p4 <- crw_preds %>%
  filter(DBH %in% c(4, 10, 20), CCH == 0.6) %>%
  mutate(DBH_label = paste0("DBH = ", DBH, " in")) %>%
  ggplot(aes(x = CR, y = DHTLC, color = as.factor(DHTANN), group = as.factor(DHTANN))) +
  geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  facet_grid(NAME ~ DBH_label, scales = "free_y") +
  scale_color_manual(values = c("steelblue","darkorange","firebrick","purple4"),
                     labels = c("0.5","1.0","2.0","3.0")) +
  labs(title = "Crown recession stress test: Fortran CRW form (annual DHTLC, ft yr⁻¹)",
       subtitle = paste0("DHTLC = CL×(1−exp(b0+b1×DHTANN+b2×CCH)); CCH=0.6, CL=HT×CR\n",
         "Positive DHTLC = crown base recedes upward; shown at three DBH classes"),
       x = "Initial crown ratio", y = "Annual DHTLC (ft yr⁻¹)",
       color = "DHTANN\n(ft yr⁻¹)") +
  theme_bw(base_size = 10) +
  theme(legend.position = "bottom", strip.text = element_text(size = 8))

ggsave(file.path(OUT_DIR, "stress_crw.png"), p4,
       width = 12, height = 8, dpi = 150)
cat("  -> stress_crw.png saved\n")

## ── Section 5: Mortality stress test ───────────────────────────────────────
cat("\n=== Section 5: Mortality ===\n")
## Gompertz-based annual survival: P_surv = exp(-exp(eta))
## eta = b0 + b1*DBH_cm + b2*DBH_cm^2 + b3*CR + b4*BGI_m
##       + b5*ln(BAL_m2ha+5) + b6*sqrt(BA_m2ha*RD)

pred_surv <- function(mort_row, dbh_in, cr=0.5,
                      bgi=NE_BGI,
                      ba_m2ha=NE_BA*FT2AC_TO_M2HA,
                      bal_m2ha, rd=NE_RD) {
  dbh_cm   <- dbh_in * IN_TO_CM
  ln_bal   <- log(bal_m2ha + 5)
  sqrt_bard <- sqrt(ba_m2ha * rd)
  eta <- mort_row$b0 +
         mort_row$b1 * dbh_cm +
         mort_row$b2 * dbh_cm^2 +
         mort_row$b3 * cr +
         mort_row$b4 * bgi +
         mort_row$b5 * ln_bal +
         mort_row$b6 * sqrt_bard
  exp(-exp(eta))
}

dbh_seq_mort <- seq(2, 30, by=0.5)

mort_preds <- lapply(sp5, function(sp) {
  m_row <- co_mort[co_mort$SPCD == sp, ]
  nm <- SPECIES$NAME[SPECIES$SPCD == sp]
  lapply(c(0, 50, 100, 150), function(bal) {
    bal_m <- bal * FT2AC_TO_M2HA
    s <- pred_surv(m_row, dbh_seq_mort, bal_m2ha=bal_m)
    data.frame(DBH=dbh_seq_mort, P_surv=s,
               BAL=as.character(bal), SPCD=sp, NAME=nm)
  }) %>% bind_rows()
}) %>% bind_rows()

cat("  Annual survival range:\n")
ms <- mort_preds %>% group_by(NAME) %>%
  summarise(min_s=round(min(P_surv),4), max_s=round(max(P_surv),4), .groups="drop")
print(ms)
cat("  Expected: >0.97 at large DBH low density, lower for small/suppressed trees\n")

p5 <- mort_preds %>%
  ggplot(aes(x=DBH, y=1-P_surv, color=BAL, group=BAL)) +
  geom_line(linewidth=0.9) +
  facet_wrap(~NAME, ncol=3, scales="free_y") +
  scale_color_manual(values=c("steelblue","darkorange","firebrick","purple"),
                     labels=c("BAL=0","BAL=50","BAL=100","BAL=150")) +
  ## Annual mortality probabilities here are very small (order 1e-4 and below),
  ## so percent formatting rounded every tick to 0.0%. Scientific tick labels
  ## keep the (free_y) per-species curves legible.
  scale_y_continuous(labels=function(x) formatC(x, format="e", digits=1)) +
  labs(title="Mortality stress test: Greg GOMPIT (size+BGI, annual)",
       subtitle=paste0("P_survive = exp(-exp(η)); CR=0.5, BGI=",NE_BGI,
         " m, BA=",NE_BA," ft²/ac"),
       x="DBH (in)", y="Annual mortality probability",
       color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom", strip.text=element_text(size=9))

ggsave(file.path(OUT_DIR, "stress_mortality.png"), p5,
       width=10, height=6.5, dpi=150)
cat("  -> stress_mortality.png saved\n")

## ── Section 6: 50-year single-tree simulation ─────────────────────────────
cat("\n=== Section 6: 50-year simulation ===\n")
## Representative tree: red spruce, DBH=6 in, HT=40 ft, CR=0.5, SI=60, BAL=80 ft2/ac
## Models:
##   fvs_base       - FVS NE default periodic (5-yr) DG/HG scaled annually
##   fvs_regional   - Same with 1.0x calibration multiplier (placeholder)
##   organon        - Greg annual DG + Greg ORGANON HG + GOMPIT survival
##   conus_spdep    - Greg annual DG + Aaron compound_none HG + GOMPIT survival

sim_sp <- 97   # red spruce
nm_sp  <- "Red spruce"
cat(sprintf("  Running 50-year simulation for %s (SPCD=%d)\n", nm_sp, sim_sp))

## Starting tree state
init <- list(DBH=6.0, HT=40.0, CR=0.5, BAL=80.0,
             BA=NE_BA, RD=NE_RD)

## FVS NE default annual DG/HG approximations for RS
## Typical 5-yr periodic DG at DBH=6, BAL=80 ~ 0.35 in/5yr = 0.070 in/yr
## Typical 5-yr periodic HG at HT=40, CR=0.5 ~ 1.5 ft/5yr = 0.30 ft/yr
fvs_dg_base_rate <- 0.070   # in/yr at starting conditions
fvs_hg_base_rate <- 0.30    # ft/yr at starting conditions

## Simple FVS-type taper functions (decline with DBH/HT)
fvs_dg_ann <- function(dbh, bal) {
  fvs_dg_base_rate * exp(-0.02*(dbh-6)) * exp(-0.003*(bal-80))
}
fvs_hg_ann <- function(ht) {
  fvs_hg_base_rate * exp(-0.015*(ht-40))
}

## Greg coefficient rows for RS
dg_rs  <- co_dg[co_dg$SPCD == sim_sp, ]
hg_rs  <- co_hg[co_hg$SPCD == sim_sp, ]
hgcn_rs <- co_hgcn[co_hgcn$SPCD == sim_sp, ]
mort_rs <- co_mort[co_mort$SPCD == sim_sp, ]

## Simulation runner
run_sim <- function(label, dg_fn, hg_fn, mort_fn=NULL,
                    years=50, init_state=init) {
  s <- init_state
  traj <- data.frame(Year=0, DBH=s$DBH, HT=s$HT, CR=s$CR,
                     BAL=s$BAL, Model=label)
  for (yr in 1:years) {
    ## Growth increments
    dg <- dg_fn(s)
    hg <- hg_fn(s)
    ## Crown recession: Fortran form DHTLC = CL*(1-exp(b0+b1*DHTANN+b2*CCH))
    ## Use hg as DHTANN (annual height increment, ft/yr); CCH from NE defaults
    crw_row <- co_crw[co_crw$SPCD == sim_sp, ]
    dhtlc  <- pred_crw(crw_row, s$HT, s$CR, dhtann = hg, cch = NE_CCH_CRW)
    ## Survival check
    alive <- TRUE
    if (!is.null(mort_fn)) {
      p_s <- mort_fn(s)
      alive <- runif(1) < p_s
    }
    if (!alive) {
      cat(sprintf("    %s: tree died at year %d\n", label, yr))
      break
    }
    s$DBH <- max(s$DBH + dg, s$DBH)
    s$HT  <- max(s$HT  + hg, s$HT)
    ## Update crown ratio via ht-to-live-crown change
    hlc_old <- s$HT - hg   # pre-growth HT, approximately
    hlc_old <- hlc_old * (1 - s$CR)
    new_hlc  <- max(0, hlc_old + dhtlc)
    s$CR  <- max(0.05, min((s$HT - new_hlc) / s$HT, 1.0))
    ## BAL declines slightly as tree grows relative to stand
    s$BAL <- max(s$BAL - 0.3, 20)
    traj <- rbind(traj, data.frame(Year=yr, DBH=s$DBH, HT=s$HT,
                                   CR=s$CR, BAL=s$BAL, Model=label))
  }
  traj
}

set.seed(42)
traj_list <- list(
  run_sim("FVS base",
    dg_fn   = function(s) fvs_dg_ann(s$DBH, s$BAL),
    hg_fn   = function(s) fvs_hg_ann(s$HT)),
  run_sim("FVS regional (1.0×)",
    dg_fn   = function(s) 1.0 * fvs_dg_ann(s$DBH, s$BAL),
    hg_fn   = function(s) 1.0 * fvs_hg_ann(s$HT)),
  run_sim("ORGANON (Greg annual)",
    dg_fn   = function(s) pred_dg_greg(dg_rs, s$DBH, s$CR, s$HT, s$BAL),
    hg_fn   = function(s) pred_hg(hg_rs, s$HT, s$CR),
    mort_fn = function(s) pred_surv(mort_rs, s$DBH,
                                    bal_m2ha=s$BAL*FT2AC_TO_M2HA)),
  run_sim("CONUS sp-dep (Aaron HG + Greg DG)",
    dg_fn   = function(s) pred_dg_greg(dg_rs, s$DBH, s$CR, s$HT, s$BAL),
    hg_fn   = function(s) pred_hg(hgcn_rs, s$HT, s$CR),
    mort_fn = function(s) pred_surv(mort_rs, s$DBH,
                                    bal_m2ha=s$BAL*FT2AC_TO_M2HA))
)

all_traj <- bind_rows(traj_list)
model_colors <- c("FVS base"="grey40",
                  "FVS regional (1.0×)"="grey60",
                  "ORGANON (Greg annual)"="steelblue",
                  "CONUS sp-dep (Aaron HG + Greg DG)"="darkorange")

p_dbh <- ggplot(all_traj, aes(x=Year, y=DBH, color=Model)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=model_colors) +
  labs(title=sprintf("%s (SPCD=%d): 50-yr DBH trajectory", nm_sp, sim_sp),
       subtitle="Start: DBH=6 in, HT=40 ft, CR=0.5, BAL=80 ft²/ac\nANNUALIZATION FLAG: FVS base uses periodic rate converted annually; Greg uses native annual rate",
       x="Year", y="DBH (in)", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

p_ht <- ggplot(all_traj, aes(x=Year, y=HT, color=Model)) +
  geom_line(linewidth=1.0) +
  scale_color_manual(values=model_colors) +
  labs(title=sprintf("%s (SPCD=%d): 50-yr height trajectory", nm_sp, sim_sp),
       x="Year", y="Total height (ft)", color=NULL) +
  theme_bw(base_size=11) +
  theme(legend.position="bottom")

p6 <- p_dbh / p_ht
ggsave(file.path(OUT_DIR, "stress_multiyear_sim.png"), p6,
       width=10, height=9, dpi=150)
cat("  -> stress_multiyear_sim.png saved\n")

## ── Summary ────────────────────────────────────────────────────────────────
cat("\n==================================================================\n")
cat("Stress test complete. Output files in:\n")
cat(paste0("  ", OUT_DIR, "\n"))
cat("  stress_ht_dbh.png\n  stress_dg.png\n  stress_hg.png\n")
cat("  stress_crw.png\n  stress_mortality.png\n  stress_multiyear_sim.png\n")
cat("\nKey flags to resolve before production:\n")
cat("  1. ANNUALIZATION: Confirm FVS calls Greg DG/HG equations how many\n")
cat("     times per 5-year cycle (1x vs 5x) — determines whether\n")
cat("     Greg annual rates match FVS periodic output or 1/5 of it.\n")
cat("  2. CRW FORM: Now uses Fortran CROWNDRIVER form DHTLC=CL*(1-exp(b0+b1*DHTANN+b2*CCH)).\n")
cat("     CCH is a crown-closure fraction (0-1), same scale as DG/HG/mort; NE_CCH_CRW=0.35.\n")
cat("  3. DG CLIMATE VARS: CSV header says PC1/PC2 but validate script\n")
cat("     uses ELEV_ft + EMT_C directly. Confirm coefficient interpretation.\n")
cat("==================================================================\n")
