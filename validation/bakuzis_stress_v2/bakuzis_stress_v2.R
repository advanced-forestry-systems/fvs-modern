#!/usr/bin/env Rscript
## bakuzis_stress_v2.R
## Step-aware Bakuzis stress test for FVS stand projections (fvs-modern engine
## harness, fvs-conus projector v4, fvs-calibration synthetic gate, or any long table).
##
## Written 2026-10-27 against the weiskittel-toolkit bakuzis-assessment skill and its
## 2026-09-07 defect record (rate-as-difference screens, silent n/a, hardcoded verdicts).
## Every relation reports EVALUATED, NOT EVALUATED and FLAGGED counts. No verdict is
## hardcoded. An axis that cannot be tested is reported as NOT EVALUATED, never as PASS.
##
## Usage:
##   Rscript bakuzis_stress_v2.R --input=<csv> --source=engine|projector|gate|long \
##       [--out=<dir>] [--qmd_rate=0.10] [--sdimax=<csv: group,sdimax>] \
##       [--obs_ba_p99=<m2/ha>] [--obs_ba_p999=<m2/ha>] [--ff=0.40] [--label=<tag>]
##
## Sources and the columns each adapter expects
##   engine    fvs-modern validation/bakuzis_engine_verdict engine_bakuzis_100yr.csv
##             (year, age, tpa, tcuft, atba, attopht, species_group, site_class,
##              density_class, scenario; imperial)
##   projector fvs-modern validation/conus_projector *_metrics.csv
##             (STAND_CN, PROJ_YEAR, VARIANT, CONFIG, BA_FT2AC, QMD_IN, TPH, HT_DOM_M)
##             site = CONFIG is not a site axis; pass --site_col to name one if present.
##   gate      fvs-calibration conus_integration_pipeline bakuzis_trajectories.csv
##             (arm, spname, SI, age, TPA, BA, QMD, TOPHT, VOL; imperial)
##   long      already metric: arm, group, site, density, traj, age, TPH, BAPH, QMD, HT, [VOL]
##
## Output (in --out): <label>_scorecard_by_trajectory.csv, <label>_scorecard_summary.csv,
##   <label>_reineke.csv, <label>_long_metric.csv, <label>_bakuzis_matrix.png, and the
##   report on stdout. All units metric: TPH trees ha-1, BAPH m2 ha-1, QMD cm, HT m, VOL m3 ha-1.

suppressPackageStartupMessages({ library(data.table) })
args <- commandArgs(trailingOnly = TRUE)
ga <- function(n, d = NULL) { m <- grep(paste0("^--", n, "="), args, value = TRUE)
  if (!length(m)) return(d); sub(paste0("^--", n, "="), "", m[1]) }
IN     <- ga("input");  if (is.null(IN)) stop("--input= is required")
SRC    <- ga("source", "long")
OUT    <- ga("out", dirname(IN))
QRATE  <- as.numeric(ga("qmd_rate", "0.10"))     # cm yr-1, self-thinning screen
FF     <- as.numeric(ga("ff", "0.40"))           # form factor when VOL is absent
LABEL  <- ga("label", paste0(SRC, "_", format(Sys.Date(), "%Y%m%d")))
SDIMAX <- ga("sdimax"); OBS99 <- as.numeric(ga("obs_ba_p99", "NA")); OBS999 <- as.numeric(ga("obs_ba_p999", "NA"))
SITECOL <- ga("site_col")
OBSFILE <- ga("obs_ba_file")   # optional csv: group, ba_p99, ba_p999 (m2 ha-1) overriding the scalar thresholds per group
SHOCK  <- as.numeric(ga("shock", "0.25"))        # first-step stem or BA loss fraction that flags
MAXTRAJ <- as.integer(ga("max_traj_per_group", "0"))   # >0: random sample of trajectories per arm x group (seed 20261027)
WINDOW <- c(-2.2, -1.2)                          # Reineke acceptance window (skill default)
RD_GATE <- as.numeric(ga("rd_gate", "0.55"))     # SDI/SDImax above which the self-thinning limb is fitted (when SDImax is supplied)
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- adapters
FT2AC_M2HA <- 0.2295684; TPA_TPH <- 2.4710538; FT_M <- 0.3048; CUFTAC_M3HA <- 0.0699729
qmd_cm <- function(BAPH, TPH) ifelse(TPH > 0, sqrt(BAPH / TPH * 40000 / pi), NA_real_)
read_engine <- function(f) {
  d <- fread(f)
  d[, .(arm = "engine", group = species_group, site = site_class, density = density_class,
        traj = sprintf("%s|%s|%s", species_group, site_class, density_class),
        age = as.numeric(age), TPH = tpa * TPA_TPH, BAPH = atba * FT2AC_M2HA,
        HT = attopht * FT_M, VOL = tcuft * CUFTAC_M3HA)][, QMD := qmd_cm(BAPH, TPH)][]
}
read_projector <- function(f) {
  ## f may be one file or a comma-separated list (one per variant); STAND_CN is read as character so a
  ## 64-bit CN written as a double is never coerced. HT_DOM_M is optional (older metrics files lack it);
  ## without it HT and VOL are NA and every height or volume relation reports NOT EVALUATED.
  files <- strsplit(f, ",")[[1]]
  d <- rbindlist(lapply(files, function(x) fread(x, colClasses = list(character = "STAND_CN"))), fill = TRUE)
  site <- if (!is.null(SITECOL) && SITECOL %in% names(d)) d[[SITECOL]] else "unspecified"
  if (!"HT_DOM_M" %in% names(d)) d[, HT_DOM_M := NA_real_]
  d[, .(arm = CONFIG, group = VARIANT, site = site, density = "obs",
        traj = paste(VARIANT, STAND_CN, sep = ":"), age = as.numeric(PROJ_YEAR),
        TPH = TPH, BAPH = BA_FT2AC * FT2AC_M2HA, HT = HT_DOM_M,
        VOL = ifelse(is.finite(HT_DOM_M), BA_FT2AC * FT2AC_M2HA * HT_DOM_M * FF, NA_real_))][, QMD := qmd_cm(BAPH, TPH)][]
}
read_gate <- function(f) {
  d <- fread(f); nm <- names(d)
  pick <- function(cands) { h <- cands[cands %in% nm]; if (!length(h)) NA_character_ else h[1] }
  cTPA <- pick(c("TPA","tpa")); cBA <- pick(c("BA","ba","atba")); cHT <- pick(c("TOPHT","topht","HT","attopht"))
  cV <- pick(c("VOL","vol","MCUFT","tcuft")); cAge <- pick(c("age","AGE","year"))
  cArm <- pick(c("arm","ARM")); cSp <- pick(c("spname","species","spp")); cSI <- pick(c("SI","si","site"))
  if (any(is.na(c(cTPA, cBA, cHT, cAge)))) stop("gate adapter: could not find TPA/BA/HT/age columns in ", paste(nm, collapse = ","))
  d[, .(arm = if (is.na(cArm)) "gate" else get(cArm), group = if (is.na(cSp)) "all" else get(cSp),
        site = if (is.na(cSI)) "unspecified" else get(cSI), density = "fixed",
        traj = paste(if (is.na(cArm)) "gate" else get(cArm), if (is.na(cSp)) "all" else get(cSp), if (is.na(cSI)) "" else get(cSI), sep = "|"),
        age = as.numeric(get(cAge)), TPH = get(cTPA) * TPA_TPH, BAPH = get(cBA) * FT2AC_M2HA,
        HT = get(cHT) * FT_M, VOL = if (is.na(cV)) NA_real_ else get(cV) * CUFTAC_M3HA)][, QMD := qmd_cm(BAPH, TPH)][]
}
## uncert: fvs-modern calibration/output/bakuzis/bakuzis_uncertainty_summary_long.csv (posterior medians, real FIA stands)
read_uncert <- function(f) {
  d <- fread(f)[config %in% c("default", "calibrated_map") & n_draws > 0 & is.finite(median) & variable %in% c("tpa","atba","tcuft","attopht")]
  w <- dcast(d, variant + config + scenario + species_group + site_class + density_class + year ~ variable, value.var = "median")
  w <- w[is.finite(tpa) & is.finite(atba) & is.finite(attopht) & tpa > 0]
  w[, .(arm = paste(variant, config, sep = "_"), group = species_group, site = site_class, density = density_class,
        traj = sprintf("%s|%s|%s", species_group, site_class, density_class), age = as.numeric(year - min(year)),
        TPH = tpa * TPA_TPH, BAPH = atba * FT2AC_M2HA, HT = attopht * FT_M, VOL = tcuft * CUFTAC_M3HA)][, QMD := qmd_cm(BAPH, TPH)][]
}
read_long <- function(f) { d <- fread(f)
  need <- c("arm","group","site","density","traj","age","TPH","BAPH","QMD","HT")
  miss <- setdiff(need, names(d)); if (length(miss)) stop("long adapter: missing ", paste(miss, collapse = ","))
  if (!"VOL" %in% names(d)) d[, VOL := NA_real_]; d[, age := as.numeric(age)]; d }
d <- switch(SRC, engine = read_engine(IN), projector = read_projector(IN), gate = read_gate(IN), long = read_long(IN), uncert = read_uncert(IN),
            stop("unknown --source"))
d[(!is.finite(VOL) | VOL <= 0) & is.finite(HT), VOL := BAPH * HT * FF]
d[, site := as.character(site)][, density := as.character(density)]
if (MAXTRAJ > 0) { set.seed(20261027)
  keep_tr <- d[, .(traj = unique(traj)), by = .(arm, group)][, .(traj = if (.N > MAXTRAJ) sample(traj, MAXTRAJ) else traj), by = .(arm, group)]
  d <- d[traj %in% keep_tr$traj]; cat(sprintf("subsampled to at most %d trajectories per arm x group: %d trajectories kept\n", MAXTRAJ, uniqueN(d$traj))) }
setorder(d, arm, group, site, density, traj, age)
d[, SDI := ifelse(TPH > 0 & QMD > 0, TPH * (QMD / 25.4)^1.605, NA_real_)]   # Reineke SDI, metric form
fwrite(d, file.path(OUT, paste0(LABEL, "_long_metric.csv")))

## site ordering helper: site levels sorted by a numeric if parseable, else Low<Medium<High, else as given
site_rank <- function(s) { u <- unique(s); n <- suppressWarnings(as.numeric(u))
  if (all(is.finite(n))) return(match(s, u[order(n)]))
  lv <- c("Low","Medium","High","Poor","Med","Good"); if (all(u %in% lv)) return(match(s, c("Low","Poor","Medium","Med","High","Good")[c("Low","Poor","Medium","Med","High","Good") %in% u]))
  match(s, sort(u)) }
d[, site_r := site_rank(site)]
dens_rank <- function(s) { u <- unique(s); n <- suppressWarnings(as.numeric(u))
  if (all(is.finite(n))) return(match(s, u[order(n)])); lv <- c("Low","Medium","High"); if (all(u %in% lv)) return(match(s, lv[lv %in% u])); match(s, sort(u)) }
d[, dens_r := dens_rank(density)]

sdimax_map <- NULL
obs_map <- NULL
if (!is.null(OBSFILE) && file.exists(OBSFILE)) { ob <- fread(OBSFILE); obs_map <- ob[, .(group = as.character(group), ba_p99, ba_p999)] }
if (!is.null(SDIMAX) && file.exists(SDIMAX)) { sm <- fread(SDIMAX); sdimax_map <- setNames(sm$sdimax, sm$group) }
## Threshold key check (added 2026-10-28 after a silent miss: thresholds keyed by variant against data grouped by
## species group left T4 NOT EVALUATED without a word). A threshold file whose groups match none of the data groups
## is almost certainly keyed wrong, so say so loudly; a partial match is reported group by group.
if (!is.null(obs_map) || !is.null(sdimax_map)) {
  dg <- unique(as.character(d$group)); tg <- unique(c(if (!is.null(obs_map)) obs_map$group, names(sdimax_map)))
  if (!any(dg %in% tg)) cat("WARNING: no threshold group matches any data group. data groups:", paste(dg, collapse = ","),
                            " threshold groups:", paste(tg, collapse = ","), ". T4 will report NOT EVALUATED everywhere.\n")
  else if (!all(dg %in% tg)) cat("NOTE: data groups without thresholds (T4 not evaluated there):", paste(setdiff(dg, tg), collapse = ","), "\n")
}

rows <- list(); add <- function(key, ...) rows[[length(rows) + 1]] <<- as.data.table(c(key, list(...)))
step_of <- function(a) { s <- stats::median(diff(sort(unique(a)))); if (!is.finite(s) || s <= 0) NA_real_ else s }
cat(strrep("=", 72), "\nBAKUZIS STRESS TEST v2  ", LABEL, "\n  input: ", IN, "\n  source: ", SRC, "   qmd_rate: ", QRATE,
    " cm yr-1   form factor: ", FF, "\n", strrep("=", 72), "\n", sep = "")

## ---------------------------------------------------------------- per-trajectory tests
for (tr in unique(d$traj)) for (am in unique(d[traj == tr, arm])) {
  g <- d[traj == tr & arm == am]; setorder(g, age); if (nrow(g) < 2) next
  stp <- step_of(g$age); key <- list(arm = am, group = g$group[1], site = g$site[1], density = g$density[1], traj = tr, step_yr = stp)
  ## T0 initial-cycle shock: a stand should not lose a quarter of its stems or basal area in the first step
  st_loss <- 1 - g$TPH[2] / g$TPH[1]; ba_loss <- 1 - g$BAPH[2] / g$BAPH[1]
  add(key, relation = "T0 initial-step shock", status = ifelse(st_loss > SHOCK | ba_loss > SHOCK, "FLAG", "PASS"),
      value = max(st_loss, ba_loss), detail = sprintf("first %g yr: stems %+.0f%%, BA %+.0f%%", stp, -100 * st_loss, -100 * ba_loss))
  ## T1 monotone development (HT up; QMD up; TPH down unless ingrowth present)
  tol <- 1e-6
  ht_na <- !all(is.finite(g$HT)); ht_ok <- if (ht_na) NA else all(diff(g$HT) >= -0.05)   # 5 cm tolerance for top-height rounding
  qmd_ok <- all(diff(g$QMD) >= -0.05)
  tph_up <- sum(diff(g$TPH) > 0.5)           # steps with stem gain = ingrowth or a defect; reported, not flagged alone
  add(key, relation = "T1 monotone HT_dom", status = ifelse(ht_na, "NOT EVALUATED", ifelse(ht_ok, "PASS", "FLAG")), value = if (ht_na) NA_real_ else min(diff(g$HT)), detail = if (ht_na) "no dominant height in this source" else sprintf("min step dHT %.2f m", min(diff(g$HT))))
  add(key, relation = "T1 monotone QMD", status = ifelse(qmd_ok, "PASS", ifelse(tph_up > 0, "NOT EVALUATED", "FLAG")), value = min(diff(g$QMD)),
      detail = ifelse(qmd_ok, "monotone", ifelse(tph_up > 0, sprintf("QMD dips with %d stem-gain steps (ingrowth present); QMD monotonicity not diagnostic", tph_up), sprintf("QMD declines %.2f cm with no ingrowth", min(diff(g$QMD))))))
  ## T9 stand-level inversions: BA or VOL falls > 10% in one step (Paper 2 counted these as inversions)
  ba_inv <- sum(diff(g$BAPH) / head(g$BAPH, -1) < -0.10, na.rm = TRUE); v_inv <- sum(diff(g$VOL) / head(g$VOL, -1) < -0.10, na.rm = TRUE)
  add(key, relation = "T9 BA inversions (>10% step drop)", status = ifelse(ba_inv == 0, "PASS", "FLAG"), value = ba_inv, detail = sprintf("%d of %d steps", ba_inv, nrow(g) - 1))
  add(key, relation = "T9 VOL inversions (>10% step drop)", status = ifelse(!all(is.finite(g$VOL)), "NOT EVALUATED", ifelse(v_inv == 0, "PASS", "FLAG")), value = v_inv, detail = ifelse(!all(is.finite(g$VOL)), "no volume in this source", sprintf("%d of %d steps", v_inv, nrow(g) - 1)))
  ## T3 Reineke, step aware (defect record 2026-09-07): screen is a RATE times this trajectory's own step.
  ## Added 2026-10-27: a relative-density gate. Reineke describes the self-thinning limb, so when SDImax is
  ## known only points at SDI/SDImax >= RD_GATE are fitted; a trajectory that never reaches the zone is
  ## NOT APPLICABLE (with its peak RD printed), not FLAG and not PASS. Without SDImax the stem-loss screen stands.
  smx <- if (!is.null(sdimax_map) && g$group[1] %in% names(sdimax_map)) sdimax_map[[g$group[1]]] else NA_real_
  rd <- if (is.finite(smx)) g$SDI / smx else rep(NA_real_, nrow(g))
  keep <- c(TRUE, diff(g$QMD) > QRATE * stp)
  if (is.finite(smx)) keep <- keep & rd >= RD_GATE
  gk <- g[keep]
  loses <- (max(g$TPH) - min(g$TPH)) / max(g$TPH) > 0.01
  in_zone <- if (is.finite(smx)) sum(rd >= RD_GATE, na.rm = TRUE) >= 3 else loses
  if (nrow(gk) >= 3 && sd(log(gk$QMD)) > 1e-3 && in_zone) {
    sl <- coef(lm(log(TPH) ~ log(QMD), gk))[2]
    ok <- is.finite(sl) && sl > WINDOW[1] && sl < WINDOW[2]
    add(key, relation = "T3 Reineke slope (trajectory OLS)", status = ifelse(ok, "PASS", "FLAG"), value = sl,
        detail = sprintf("slope %.2f on %d of %d points (window %g to %g%s)", sl, nrow(gk), nrow(g), WINDOW[1], WINDOW[2],
                         ifelse(is.finite(smx), sprintf("; RD gate %.2f, peak RD %.2f", RD_GATE, max(rd, na.rm = TRUE)), "")))
  } else if (is.finite(smx) && !in_zone) {
    add(key, relation = "T3 Reineke slope (trajectory OLS)", status = "NOT APPLICABLE", value = NA_real_,
        detail = sprintf("never in the self-thinning zone: peak SDI/SDImax %.2f < gate %.2f (stems %+.0f%% over the run)", max(rd, na.rm = TRUE), RD_GATE, -100 * (max(g$TPH) - min(g$TPH)) / max(g$TPH)))
  } else if (loses) {
    add(key, relation = "T3 Reineke slope (trajectory OLS)", status = "NOT EVALUATED", value = NA_real_,
        detail = sprintf("stand loses %.0f%% of stems but screen left %d of %d points; check qmd_rate against the %g yr step", 100 * (max(g$TPH) - min(g$TPH)) / max(g$TPH), nrow(gk), nrow(g), stp))
  } else add(key, relation = "T3 Reineke slope (trajectory OLS)", status = "NOT APPLICABLE", value = NA_real_, detail = "stem count flat; no self-thinning phase")
  ## T4 magnitude: peak SDI against SDImax and BA against observed percentiles
  pk <- max(g$SDI, na.rm = TRUE); smx <- if (!is.null(sdimax_map) && g$group[1] %in% names(sdimax_map)) sdimax_map[[g$group[1]]] else NA_real_
  add(key, relation = "T4 peak SDI / SDImax", status = ifelse(is.na(smx), "NOT EVALUATED", ifelse(pk / smx > 1.05, "FLAG", "PASS")), value = ifelse(is.na(smx), NA_real_, pk / smx),
      detail = ifelse(is.na(smx), "no SDImax supplied for this group (--sdimax=)", sprintf("peak SDI %.0f vs SDImax %.0f (ratio %.2f)", pk, smx, pk / smx)))
  bmx <- max(g$BAPH, na.rm = TRUE)
  OBS99 <- if (!is.null(obs_map) && g$group[1] %in% obs_map$group) obs_map[group == g$group[1], ba_p99][1] else as.numeric(ga("obs_ba_p99", "NA"))
  OBS999 <- if (!is.null(obs_map) && g$group[1] %in% obs_map$group) obs_map[group == g$group[1], ba_p999][1] else as.numeric(ga("obs_ba_p999", "NA"))
  add(key, relation = "T4 max BA vs observed p99", status = ifelse(is.na(OBS99), "NOT EVALUATED", ifelse(bmx > OBS99, "FLAG", "PASS")), value = bmx,
      detail = ifelse(is.na(OBS99), sprintf("max BA %.1f m2 ha-1; no observed p99 supplied (--obs_ba_p99=)", bmx), sprintf("max BA %.1f vs observed p99 %.1f m2 ha-1", bmx, OBS99)))
  add(key, relation = "T4 max BA vs observed p99.9", status = ifelse(is.na(OBS999), "NOT EVALUATED", ifelse(bmx > OBS999, "FLAG", "PASS")), value = bmx,
      detail = ifelse(is.na(OBS999), sprintf("max BA %.1f m2 ha-1; no observed p99.9 supplied (--obs_ba_p999=)", bmx), sprintf("max BA %.1f vs observed p99.9 %.1f m2 ha-1", bmx, OBS999)))
}

## ---------------------------------------------------------------- cross-trajectory tests (per arm x group)
for (am in unique(d$arm)) for (gr in unique(d[arm == am, group])) {
  h <- d[arm == am & group == gr]; key <- list(arm = am, group = gr, site = "all", density = "all", traj = "cross", step_yr = step_of(h$age))
  ages <- sort(unique(h$age)); fin <- max(ages); inert <- FALSE
  ## site axis inert? identical trajectories across site levels mean the lever did nothing
  if (uniqueN(h$site) >= 2) {
    sp <- h[age == fin, .(ht = mean(HT, na.rm = TRUE), v = mean(VOL, na.rm = TRUE)), by = .(site, site_r)]
    inert <- isTRUE((max(sp$ht) - min(sp$ht)) / mean(sp$ht) < 0.01 & (max(sp$v) - min(sp$v)) / mean(sp$v) < 0.01)
    if (inert) {
      add(key, relation = "T2 site ordering HT", status = "NOT EVALUATED", value = NA_real_, detail = "site axis inert: final HT and VOL identical (<1%) across site levels; the site lever did not reach the engine")
      add(key, relation = "T2 site ordering VOL", status = "NOT EVALUATED", value = NA_real_, detail = "site axis inert")
    } else {
      ## within each density level, at each age, ordering of HT and VOL by site rank
      inv <- h[, { o <- order(site_r); list(ht_inv = is.unsorted(HT[o]), v_inv = is.unsorted(VOL[o]), n = .N) }, by = .(density, age)][n >= 2]
      add(key, relation = "T2 site ordering HT", status = ifelse(sum(inv$ht_inv) == 0, "PASS", "FLAG"), value = sum(inv$ht_inv), detail = sprintf("%d of %d density x age cells inverted", sum(inv$ht_inv), nrow(inv)))
      add(key, relation = "T2 site ordering VOL", status = ifelse(sum(inv$v_inv) == 0, "PASS", "FLAG"), value = sum(inv$v_inv), detail = sprintf("%d of %d density x age cells inverted", sum(inv$v_inv), nrow(inv)))
    }
  } else { add(key, relation = "T2 site ordering HT", status = "NOT EVALUATED", value = NA_real_, detail = "single site level"); add(key, relation = "T2 site ordering VOL", status = "NOT EVALUATED", value = NA_real_, detail = "single site level") }
  ## T5 Eichhorn: VOL at matched dominant HT should not depend on site (CV across sites < 15%)
  if (uniqueN(h$site) >= 2 && !inert) {
    cvs <- h[, { ss <- split(.SD, site); lo <- max(sapply(ss, function(z) min(z$HT))); hi <- min(sapply(ss, function(z) max(z$HT)))
      if (hi > lo + 0.5) { grid <- seq(lo, hi, length.out = 8)
        cv <- sapply(grid, function(hg) { vv <- sapply(ss, function(z) approx(z$HT, z$VOL, hg, rule = 2, ties = mean)$y); sd(vv) / mean(vv) }); list(cv = mean(cv), n = length(grid)) } else list(cv = NA_real_, n = 0L) }, by = density]
    m <- mean(cvs$cv, na.rm = TRUE)
    add(key, relation = "T5 Eichhorn VOL~HT site independence", status = ifelse(!is.finite(m), "NOT EVALUATED", ifelse(m < 0.15, "PASS", "FLAG")), value = m,
        detail = ifelse(is.finite(m), sprintf("mean CV of VOL across sites at matched HT %.1f%%", 100 * m), "no overlapping HT range across sites"))
  } else add(key, relation = "T5 Eichhorn VOL~HT site independence", status = "NOT EVALUATED", value = NA_real_, detail = "site axis inert or single level")
  ## T6 Sukachev: per-tree volume at the final age falls with initial density within site
  if (uniqueN(h$density) >= 2) {
    sk <- h[age == fin, .(ptv = VOL / TPH, dens_r), by = .(site, density)][, .(ok = cor(dens_r, ptv, method = "spearman") < 0, n = .N), by = site][n >= 2]
    add(key, relation = "T6 Sukachev per-tree VOL falls with density", status = ifelse(nrow(sk) == 0, "NOT EVALUATED", ifelse(all(sk$ok), "PASS", "FLAG")), value = mean(sk$ok),
        detail = sprintf("%d of %d site levels ordered", sum(sk$ok), nrow(sk)))
    ## T7 yield-density: final BA should not fall with initial density by more than 10%
    yd <- h[age == fin, .(ba = mean(BAPH)), by = .(site, density, dens_r)][order(site, dens_r)][, .(drop = min(diff(ba) / head(ba, -1)), n = .N), by = site][n >= 2]
    add(key, relation = "T7 yield-density final BA vs initial density", status = ifelse(nrow(yd) == 0, "NOT EVALUATED", ifelse(all(yd$drop > -0.10), "PASS", "FLAG")), value = min(yd$drop),
        detail = sprintf("worst relative BA drop from lower to higher initial density %.0f%%", 100 * min(yd$drop)))
  } else { add(key, relation = "T6 Sukachev per-tree VOL falls with density", status = "NOT EVALUATED", value = NA_real_, detail = "single density level"); add(key, relation = "T7 yield-density final BA vs initial density", status = "NOT EVALUATED", value = NA_real_, detail = "single density level") }
  ## T8 pinned plateau: final BA identical (within 1%) across scenarios that started differently is a cap, not a biology
  fb <- h[age == fin, .(ba = BAPH), by = traj]; b0 <- h[age == min(ages), .(b0 = BAPH), by = traj]; fb <- merge(fb, b0, by = "traj")
  if (nrow(fb) >= 3) { cvf <- sd(fb$ba) / mean(fb$ba); cv0 <- sd(fb$b0) / mean(fb$b0)
    add(key, relation = "T8 pinned BA plateau", status = ifelse(cvf < 0.03 & cv0 > 0.10, "FLAG", "PASS"), value = cvf,
        detail = sprintf("final BA CV %.1f%% across %d trajectories that started with CV %.0f%%", 100 * cvf, nrow(fb), 100 * cv0)) }
  ## T3b pooled self-thinning frontier slope by quantile regression (tau 0.95) on the self-thinning limb
  smx_g <- if (!is.null(sdimax_map) && gr %in% names(sdimax_map)) sdimax_map[[gr]] else NA_real_
  lim <- h[, { stp <- step_of(age); k <- c(TRUE, diff(QMD) > QRATE * stp) & c(FALSE, diff(TPH) < 0)
               if (is.finite(smx_g)) k <- k & (SDI / smx_g >= RD_GATE); .SD[k] }, by = traj]
  if (nrow(lim) >= 10 && requireNamespace("quantreg", quietly = TRUE)) {
    fr <- tryCatch(coef(quantreg::rq(log(TPH) ~ log(QMD), tau = 0.95, data = lim))[2], error = function(e) NA_real_)
    add(key, relation = "T3b frontier slope (rq tau 0.95, pooled limb)", status = ifelse(!is.finite(fr), "NOT EVALUATED", ifelse(fr > WINDOW[1] & fr < WINDOW[2], "PASS", "FLAG")), value = fr,
        detail = sprintf("slope %.2f on %d limb points from %d trajectories", fr, nrow(lim), uniqueN(lim$traj)))
  } else add(key, relation = "T3b frontier slope (rq tau 0.95, pooled limb)", status = ifelse(is.finite(smx_g) && nrow(lim) < 10, "NOT APPLICABLE", "NOT EVALUATED"), value = NA_real_,
           detail = sprintf("%d limb points%s; need >= 10 and quantreg", nrow(lim), ifelse(is.finite(smx_g), sprintf(" at SDI/SDImax >= %.2f (peak RD %.2f)", RD_GATE, max(h$SDI, na.rm = TRUE) / smx_g), "")))
}

sc <- rbindlist(rows, fill = TRUE)
fwrite(sc, file.path(OUT, paste0(LABEL, "_scorecard_by_trajectory.csv")))
sm <- sc[, .(evaluated = sum(status %in% c("PASS","FLAG")), not_evaluated = sum(status == "NOT EVALUATED"),
             not_applicable = sum(status == "NOT APPLICABLE"), flagged = sum(status == "FLAG"),
             verdict = { e <- sum(status %in% c("PASS","FLAG")); f <- sum(status == "FLAG"); na <- sum(status == "NOT EVALUATED")
               if (e == 0) "NOT EVALUATED" else if (f == 0 && na == 0) "PASS" else if (f == 0) "PASS (partial)" else if (f == e) "FLAG (all)" else "FLAG" }),
         by = .(arm, relation)][order(arm, relation)]
fwrite(sm, file.path(OUT, paste0(LABEL, "_scorecard_summary.csv")))
fwrite(sc[grepl("^T3", relation)], file.path(OUT, paste0(LABEL, "_reineke.csv")))

cat("\nSUMMARY (evaluated / not evaluated / flagged), median age step ", stats::median(sc$step_yr, na.rm = TRUE), " yr\n", sep = "")
print(sm, nrows = 200)
cat("\nFLAGGED DETAIL\n"); fl <- sc[status == "FLAG"][order(arm, relation, group, site, density)]
if (nrow(fl)) print(fl[, .(arm, relation, group, site, density, detail)], nrows = 500) else cat("  none\n")
cat("\nNOT EVALUATED DETAIL\n"); ne <- sc[status == "NOT EVALUATED"]
if (nrow(ne)) print(ne[, .(arm, relation, group, site, density, detail)], nrows = 500) else cat("  none\n")
cat("\nReineke axis: ", sm[grepl("^T3 ", relation), sum(evaluated)], " trajectories evaluated, ", sm[grepl("^T3 ", relation), sum(not_evaluated)],
    " not evaluated, ", sm[grepl("^T3 ", relation), sum(flagged)], " flagged. A zero flag count rests only on the evaluated count.\n", sep = "")
cat("This test falsifies; a pass is not a confirmation of accuracy.\n")

## ---------------------------------------------------------------- matrix figure (lower triangle, Leary form)
if (requireNamespace("ggplot2", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggplot2))
  vars <- c("age","HT","TPH","BAPH","QMD","VOL")
  lab  <- c(age = "Age (yr)", HT = "Dominant height (m)", TPH = "Stems (ha^-1)", BAPH = "Basal area (m^2 ha^-1)", QMD = "QMD (cm)", VOL = "Volume (m^3 ha^-1)")
  pairs <- list(); for (i in 2:length(vars)) for (j in 1:(i - 1)) pairs[[length(pairs) + 1]] <- c(vars[j], vars[i])
  pd <- rbindlist(lapply(pairs, function(p) d[, .(arm, group, site, density, traj, x = get(p[1]), y = get(p[2]), xv = p[1], yv = p[2])]))
  pd[, xv := factor(xv, vars)][, yv := factor(yv, vars)]
  safe <- function(z) gsub("[^A-Za-z0-9_]", "_", z)
  for (am in unique(pd$arm)) {
    for (gr in unique(pd[arm == am, group])) {
      pp <- pd[arm == am & group == gr]
      g <- ggplot(pp, aes(x, y, group = traj, colour = site, linetype = density)) + geom_path(linewidth = 0.4, alpha = 0.9) +
        facet_grid(yv ~ xv, scales = "free", switch = "both", labeller = labeller(xv = lab, yv = lab)) +
        labs(title = sprintf("Modified Bakuzis matrix: %s, %s (%s)", am, gr, LABEL), x = NULL, y = NULL, colour = "Site", linetype = "Initial density") +
        theme_bw(base_size = 9) + theme(panel.grid = element_blank(), strip.background = element_blank(), strip.placement = "outside",
                                        legend.position = "bottom", axis.text = element_text(size = 6))
      ggsave(file.path(OUT, sprintf("%s_%s_%s_bakuzis_matrix.png", LABEL, safe(am), safe(gr))), g, width = 10, height = 9, dpi = 300)
    }
    ## the key relation (TPH vs QMD, log-log) by group so the frontier is visible
    dd <- d[arm == am]
    g2 <- ggplot(dd, aes(QMD, TPH, group = traj, colour = site, linetype = density)) + geom_path(linewidth = 0.4) +
      scale_x_log10() + scale_y_log10() + facet_wrap(~ group) +
      labs(x = "QMD (cm)", y = expression(Stems~(ha^-1)), colour = "Site", linetype = "Initial density",
           title = sprintf("Self-thinning limb by group: %s", am), subtitle = "dotted: Reineke slope -1.605 through the maximum SDI reached in this arm") +
      geom_abline(slope = -1.605, intercept = log10(max(dd$SDI, na.rm = TRUE)) + 1.605 * log10(25.4), linetype = "dotted", colour = "grey40") +
      theme_bw(base_size = 9) + theme(panel.grid = element_blank(), strip.background = element_blank(), legend.position = "bottom")
    ggsave(file.path(OUT, sprintf("%s_%s_selfthinning.png", LABEL, safe(am))), g2, width = 9, height = 7, dpi = 300)
  }
}
cat("Wrote outputs to ", OUT, "\n", sep = "")
