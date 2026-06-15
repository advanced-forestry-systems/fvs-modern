# Stress test of the whole program (2026-06-15)

A deliberately adversarial review of every major claim and of what is actually implemented, before
anything is finalized or sent to the FVS team. Each item: the claim, the evidence, the weakness, and a
confidence rating. Bugs found this session are flagged, because they show the harness needs vigilance.

## A. Maximum-SDI thread

**A1. Species-weighted maximum is +28% biased, near-zero plot skill.** Evidence: 95,206 plots against
the brms maximum. Weakness: the reference is the brms estimate, which is itself a model of an
unobservable quantity, so this is partly circular (acknowledged). The bias sign and rough magnitude are
robust, but "0.02 bias-corrected R2" is reference-dependent. Confidence: medium-high on direction, lower
on the exact number.

**A2. Localized maximum predicts observed self-thinning ~85% better (deviance 0.107 vs 0.058).** This is
the non-circular test and the strongest single result. Weakness: the absolute deviance explained is low
(0.05 to 0.12) because annual plot-level mortality is very noisy, so "85% better" is a ratio of two
small numbers; it is a real relative improvement but not a high-skill model. Confidence: high on
direction and sign, medium on magnitude framing. Action: always report it as relative improvement on a
weak signal, never imply high absolute skill.

**A3. All 20 variants: level-calibrated localized beats native in 19.** Weakness, and a bug found this
session: the sweep hard-coded the sample seed, so the four OR/WA variants drew one pathological sample
(RMSE inflated to ~98%); I reran those four with a good seed (~47-50%), but the other 16 variants still
rest on a single fixed-seed sample of ~80-140 plots each. The within-variant native-vs-calibrated
comparison is valid (same plots), but cross-variant magnitudes are noisy and not seed-consistent.
Confidence: medium on "calibrated >= native in most variants," low on any specific variant's magnitude.
Action: rerun all variants with multiple seeds before any quantitative claim is published; for the
briefing, state it qualitatively. **Quantified (this session):** a three-seed check on Pacific Northwest
gave native RMSE of 49, 68, and 62 percent across seeds 7, 11, 23, an 18-point spread, which directly
confirms the per-variant magnitudes are seed-sensitive and must be read qualitatively until a
multi-seed average is computed. The within-variant native-versus-calibrated direction was stable across
those seeds (calibration neutral-to-helpful in each), so the directional conclusion holds; only the
magnitudes are noisy.

**A4. Joint fit: per-variant level doubles self-thinning R2 (0.021 -> 0.046).** Weakness: both numbers
are tiny, and the data-estimated levels correlate only 0.12 with the engine levels, i.e. the levels are
poorly identified. The qualitative conclusion (level must be co-estimated with mortality) holds; the
level values themselves are uncertain. Confidence: medium on the conclusion, low on the level estimates.

**A5. The engine sweep measures density (TPH) error, which conflates the SDIMAX effect with FVS
mortality.** This is a real interpretive limit: "calibrated SDIMAX beats native" partly reflects
interaction with the already-tuned native mortality, not the maximum alone. It is consistent with A4's
conclusion but means the sweep is not a clean isolation of the maximum. Confidence: note as a caveat.

## B. Implementation status (audited in the repo today)

**B1. Component calibration is only partially populated.** Mortality, crown, and height-diameter
multipliers are real for all 25 variants; diameter growth is real for only ~7 variants (NE, PN, SN, CR
and most others are 1.0); height growth is 1.0 for every variant. So "Step 1 fully implemented across
all variants" is false: it is complete for mortality and two others, partial for diameter growth, and
absent for height growth. Filling the gap needs the per-variant component fits, which are compute-gated.
Confidence: high (audited).

**B2. The revised max SDI is not in production.** The config loader emits the native species-weighted
SDIMAX; the localized maximum lives only in a standalone module and the benchmark path. Confidence:
high.

**B3. Species-free equations are not merged.** PR #70 is open; the bundles are banked; the injection is
a prototype. Confidence: high.

**B4. Mortality multipliers at the clip bounds: checked, not a systemic problem.** Across variants only
0 to 4 percent of species hit the 0.10 or 10.0 bound (acd 0, ne 0, cr 0, cs 0, sn 0, ls 3 percent, ak 4
percent, on 1 percent), and those are rare species with sparse data. So the mortality calibration is not
broadly saturating; the concern is resolved. Confidence: resolved (audited).

## C. Species-free and Greg comparisons

**C1. Our DG competitive with Greg on Douglas-fir.** Weakness: the two predictors ran on non-identical
test sets from different pipelines, so the RMSE comparison is only approximate; the R2 gap is partly a
difference in observed variance. Confidence: medium; needs one identical held-out set before publishing.

**C2. Species-free within ~22% with zero DF data.** One species, one component. Encouraging but not
general. Confidence: medium as a proof of concept, low as a general claim.

**C3. HCB injection benchmark: pure species-free beats per-species on held-out NE.** Real and
well-calibrated, but one component, one region, and the in-sample LOO favors the species-dependent form.
Confidence: medium; extend to the other components (running) before generalizing.

## D. Benchmark harness

**D1. The brk_dbh fix is validated (year-0 reproduces observed t1).** Solid. **D2. But two harness bugs
surfaced this session (the brk_dbh tree-drop earlier, the fixed-seed sample now), so latent bugs are
plausible.** Every quantitative result should be reproduced with an independent script before
publication. Confidence: the harness is usable but not yet publication-hardened.

**D3. "Unified beats regional NE/ACD on basal area (-0.6% vs +12-13%)."** One region, modest n, and it
depends on what "unified/calibrated" config was used (mortality calibrated, growth largely 1.0). It
should be re-verified with the seed-consistency lesson and a clear statement of which multipliers were
active. Confidence: medium; do not headline without re-verification.

## E. External corroboration

**E1. Batista et al. 2026.** Already corrected: tree-ring growth is biased high (survivorship), and we
calibrated DG on FIA remeasurement, so it is directional support, not validation. No residual
overstatement after today's edit. Confidence: handled.

## F. The injection prototype

**F1. sf_injector.py is wiring-complete but untested end to end, and the fvsTreeAttr call signature is
inferred from the API, not verified against a live run.** It should be treated as a design artifact
until a shadow-mode run on Cardinal confirms it loads and the attribute names are correct. Confidence:
design only; do not describe as working.

## Prioritized fixes before finalizing or sending

1. Reword the FVS-team materials so every claim is stated at its true strength: max-SDI self-thinning
   improvement as a relative gain on a weak signal; the all-variant result as qualitative; calibration
   coverage as partial (mortality complete, growth partial/absent); species-free and Greg as promising
   not proven; injection as designed not running.
2. Rerun the all-variant level sweep with multiple seeds for seed-consistent magnitudes (compute-gated).
3. Run the Greg DG comparison on one identical held-out set.
4. Reproduce the headline NE basal-area result with an independent script and a stated config.
5. Inspect the clipped mortality multipliers.
6. Shadow-mode run of the injector before calling it working.

## Net

The program's direction is well supported and the max-SDI self-thinning result is the most robust
single finding, but several headline numbers rest on weak absolute signals, small or non-identical
samples, or partially populated configs, and the harness has had two bugs. Nothing here is fatal; the
fix is to state every claim at its true strength and to reproduce the quantitative headlines before they
go to the FVS team. The materials should be updated to this honest footing now, and the quantitative
re-verifications queued for when the cluster clears.
