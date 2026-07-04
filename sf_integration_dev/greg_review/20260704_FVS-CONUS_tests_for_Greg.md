# FVS-CONUS three-arm + stand-constraint tests — for Greg's recommendations

Date: 2026-07-04. Prepared for the Monday call. All work on OSC Cardinal; three parallel
test tracks. Nothing merged to production main; prototypes are on branch/isolated dirs for review.

## TL;DR
- **Mortality form: keep GOMPIT** (your native hook). The A/B is decisive — the exponential
  alternative annihilates every stand with your coefficient signs. Confirmed, not just assumed.
- **Two real bugs found in the stand-constraint layer, both prototyped fixed:** the top-height
  transition mean-reverts and inverts site order (fixed with a GADA envelope), and the SDIMAX
  density cap overshoots ~4% per step (fixed by re-capping after growth).
- **One engine blocker for the team:** a pre-existing unit-17 SQLite DB-input segfault blocks
  all in-engine multi-cycle projection (independent of the mortality form).

## 1. gompmort mortality form — GOMPIT confirmed (decisive)
The reconciliation had to pick between two `gompmort.f90` survival forms. The engine builds
clean either way (NE lib, 553 objects, `gregdghg.f90` in the source list; Python loads it,
`sf_engine_hook_test.py` prints MECHANISM_OK). Driving BOTH forms with the identical ETA from
your 92 fitted species coefficients (CR x CCH grid, FINTL = 5 yr):

| metric | GOMPIT `S = 1 - exp(-exp(ETA))` | EXP-HAZARD `HZ = exp(ETA); S = exp(-HZ*t)` |
|---|---|---|
| median 5-yr period survival | 0.928 | 7.0e-10 |
| median annual mortality | 1.48% | 98.52% |
| 100-yr cohort survivorship | 22.6% | 8.3e-182 |
| tree-conditions with S < 0.01 in one cycle | 0.0% | 100.0% |

The two forms read ETA with OPPOSITE sign. Your b0 are positive (~2.5), so ETA is positive;
GOMPIT reads that as high survival (realistic ~1.5%/yr background mortality), the exponential
hazard reads it as high log-hazard and kills every stand in a single cycle. **Recommendation:
keep the GOMPIT form** — it is the correct link for your cloglog/gompit fit; the exponential
form would need a sign-flipped, refit coefficient set to be usable at all.
Figure: `gompmort_ab_offline.png`; table: `gompmort_ab_table.csv`.
QUESTION FOR GREG: confirm GOMPIT is the intended engine mortality link, and confirm the b0
sign convention so we lock it in.

## 2. Engine build + the unit-17 blocker (team item)
The NE variant library builds cleanly from the reconciled source and the injection mechanism
fires (MECHANISM_OK). But the full 36-scenario engine Bakuzis pass could not run: both the
gompit and exponential libs **segfault identically at the first scenario**, at `fvs2py/_base.py:756`
inside `fvs.run()` right after "OPEN FAILED FOR 17" — the SQLite DB-input (DSNIN) projection
path fails to open Fortran unit 17 and then walks uninitialized memory. This is the same
DB-input plumbing that crashed the perseus keyfile path; it is independent of the mortality
form and blocks any engine-level multi-cycle projection.
QUESTION FOR GREG: is the unit-17 / DSNIN DB-input path known-broken, and is a keyfile-inventory
(inline TREEDATA) path the sanctioned route, or should we fix the unit-17 open? This gates the
definitive in-engine realism test and the fvs2py in-process story.

## 3. Top-height constraint — real bug, prototype fix (GADA envelope)
Off-engine, the García/GADA top-height H2|H1 state-space transition
`H2 = H1*exp(0.423 - 0.155*lnH1 + 0.059*lnyr - 0.011*rd - 0.031*lnqmd - 0.011*bgi)` behaves
correctly for a single FIA remeasurement but is a **contraction map** under repeated projection:
the negative ln(H1) slope drives H toward a step-dependent fixed point (~9.7 m for a 5-yr step),
so top height rises to ~year 25 then DECLINES; the negative bgi slope also inverts site order.
Dominant height should not shrink in an even-aged stand — see `topht_before_after.png` (left).
Prototype fix: wrap the fitted transition in a monotone base-age-invariant GADA Chapman-Richards
site envelope `H = b1*(1-exp(-b2*A))^b3` (refit b2 = 0.0358, b3 = 1.590 for the NE spruce-fir
group from NA_SITREE, n = 197,122), with a per-site asymptote and a monotone ratchet. After the
fix top height is monotone, asymptotes in the realistic 20-26 m band, and is correctly site-ordered
(`topht_before_after.png`, right).
QUESTION FOR GREG: preferred path — (a) re-constrain the fit so the composed multi-step map
converges to a site asymptote, or (b) keep the fitted transition for local signal but wrap it in
the GADA envelope tied to the real cspi_v7 site table (the prototype uses placeholder per-species
asymptotes)? Also: the stems N2|N1 transition uses the same negative-ln-slope form and should be
checked for the same pathology.

## 4. SDIMAX density cap — real ~4% overshoot, prototype fix
With the calibrated NE Reineke SDIMAX (~480 imperial) the cap is applied BEFORE the growth step,
so QMD growth (SDI ~ QMD^1.605) drifts end-of-step SDI back to ~498 (3.8% over). Not a proxy
artifact — measured on the projector's own internal summation SDI. Fix: re-apply the Reineke cap
once more AFTER growth+survival each step (one pass suffices; rescaling TPA does not change QMD),
which pins end-of-step SDI to 480.0 exactly with an essentially unchanged trajectory. Low stakes,
clean fix.

## 5. Site ordering is actually correct (on the right axis)
The earlier "site ordering violated" flag came from varying bgi (a secondary climate covariate with
a small negative slope) — the wrong axis. Varying the TRUE site anchor (initial dominant/top height)
with bgi held fixed, the constrained projection orders correctly: top height at year 100
17.1 < 18.8 < 23.8 m and yield proxy 1320 < 1453 < 1841 (low < mid < high). BA is near-flat across
sites because all are density-limited at the same SDIMAX ceiling — site correctly expresses through
height and volume, not standing BA at max density (Eichhorn-consistent). See
`track_c_true_site_ordering.png`.

## Consolidated questions for Greg
1. Confirm GOMPIT as the engine mortality link and the b0 sign convention.
2. Unit-17 / DSNIN DB-input segfault: known issue? keyfile-inventory the sanctioned path, or fix the open?
3. Top-height: re-constrain the fit to asymptote, or adopt the GADA-envelope wrap (with the real cspi_v7 table)?
4. Should the stems N2|N1 transition get the same GADA/monotonicity treatment?
5. Your climate-only (ELEV-free) DG form — when it lands we re-fit and re-land the ORGANON arm via 62c.

## Artifacts
Figures (this folder): gompmort_ab_offline.png, topht_before_after.png, track_c_true_site_ordering.png.
Prototype code + verdicts on Cardinal: /fs/scratch/PUOM0008/crsfaaron/{wt-engine/ab_engine, track_topht, track_site};
also committed to branch analysis/greg-review (sf_integration_dev/greg_review). Production main untouched;
PR #92 (reconciliation) remains the open gate.
