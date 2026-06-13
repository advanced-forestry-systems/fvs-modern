# Full stress test of the CONUS-variant approach, and the refined plan

**Date:** 2026-06-11
**Question (Aaron):** stress test the whole logic and approach, then refine the plan to be sure
we are on a path to a CONUS-wide FVS variant that outperforms the regional variants and aligns
better with biological expectations.

This is an adversarial review grounded in the actual benchmark and diagnostic outputs on
Cardinal, not a restatement of intent. The verdict up front: the modeling machinery is strong
and several core premises are validated, but **neither of the two success criteria is yet
demonstrated**, and one of them is contradicted by current evidence. The plan needs to re-center
on proving them.

## 1. What is genuinely strong (credit where due)

- **The species-free trait substitution works, including for held-out species.** The held-out
  species test (predict annual growth for species the model never saw, from traits alone) is the
  hard test of the central premise, and it passes: held-out DG RMSE 0.101 against an observed
  mean of 0.123 with bias -0.002 (essentially unbiased), HG RMSE 0.281 against mean 0.320. The
  one-equation-for-all-species idea is empirically supported, which is the program's biggest bet
  and it is sound.
- **Tree-level component accuracy is reasonable.** Calibrated tree-level R2 runs about 0.50 to
  0.66 across regions, which is solid for annual increments.
- **Three of four Bakuzis biological laws pass well:** Sukachev 100%, crown recession 100%,
  Eichhorn 84%. The species architecture (Leg A vs Leg B) is validated. Most components are
  properly annualized.

These are real and worth protecting. The problem is not the model's craftsmanship.

## 2. Critical finding 1: there is no benchmark against the real regional variants

The headline claim is that the unified model outperforms the 20 regional FVS variants. **That
comparison does not exist in the evidence.** The "default" baseline in every benchmark
(`19_fia_benchmark_engine.R`, `project_condition_default`) is a hand-rolled FVS-like
approximation: a flat 1.5 percent annual mortality and 0.10 inch per year diameter growth with
crude size and competition modifiers. It is a strawman, not the actual regional variant logic.

Two consequences. First, the tree-level "wins" are meaningless: calibrated R2 of 0.5 against a
default R2 of -1.4 to -1289 (RMSE up to 44.8) is beating a broken stand-in, not FVS. Second,
even against that strawman, the honest stand-level composite tells an unfavorable story (next
finding). **We currently cannot claim, or refute, that the unified model beats the regional
variants, because the test has not been run.** This is the single most important gap.

## 3. Critical finding 2: the unified model does not clearly beat even the strawman at stand level

The fair metric is the stand-level composite pctRMSE (basal area, SDI, QMD, volume) by variant.
The reduction from default to calibrated is **negative on roughly half the variants**, meaning
the unified model is worse there: AK -24%, PN -24%, TT -22%, WC -15%, BM and CS -12%, NE -10%,
LS -9%, and more. It is positive overall (+5.9%) almost entirely because the Southern region
(SN, 238k plots, +10%) dominates the plot-weighted average. Absolute errors are very high for
both models (composite pctRMSE above 100%, volume above 200%). So on a per-variant basis, against
a strawman it should beat easily, the unified model is roughly a wash or worse. That is a serious
warning about the stand-level projection, which is the quantity managers and the carbon and yield
applications actually use.

## 4. Critical finding 3: mortality has no senescence, so large trees are nearly immortal

The Bakuzis mortality U-shape passes only 40 percent of variants, and the reason is specific and
biologically serious. The largest-tree mortality rates collapse to near zero across most
variants: NE 1.1e-6, LS 2.9e-6, CS 0.0016, EM 0.0011, ON 1.3e-8, CR 0.007, against middle-size
rates of 0.05 to 0.19. Mortality decreases monotonically with size instead of rising again for
large and old trees. **The model lets big trees live essentially forever.** This is wrong
biologically (senescence and large-tree hazard are real) and mechanically damaging: with no
senescence sink, basal area and volume over-accumulate in long projections, which is a plausible
driver of the stand-level overprediction in finding 3 and a direct hit to biological realism, the
second success criterion.

## 4a. Tier 0 follow-through: the mortality defect is diagnosed and the fix is clear

I tested the senescence premise against FIA directly rather than assuming it. Two results:

- **In absolute DBH, there is no senescence:** observed annual mortality decreases monotonically
  from 1.7 percent (small) to about 0.9 to 1.1 percent (large), with only a noise-level uptick
  above 100 cm (32 trees). Pooled across species, big trees are not dying more.
- **In relative size (DBH divided by species maximum), senescence is real:** mortality is
  elevated at small relative size (0.014 to 0.017, suppression), dips mid-size (0.015), then
  rises near the species maximum, 0.0176 at 0.75 to 0.9 and 0.0198 above 0.9. That is a genuine,
  if mild, U-shape, with the right arm appearing only relative to each species own maximum size.

So the precise defect is not a missing absolute-size hump. It is two things: the model drives
large-tree hazard to about 1e-6, far below the empirical floor of roughly 1 to 2 percent, and it
has no relative-size term to capture the senescence rise. The current eta uses absolute DBH plus
DBH squared with the DBH-squared prior pinned near zero, which cannot represent relative-size
senescence and lets the hazard collapse.

**The fix, and it fits the trait architecture cleanly:** add a relative-size term, DBH divided by
the species maximum DBH, to the survival eta, with a sign that raises hazard near the maximum.
Species maximum DBH is already a trait in the trait table, so this is a trait-mediated senescence
effect, not a bolt-on. Combined with preventing the hazard from collapsing below the empirical
floor, this makes mortality U-shaped (suppression plus senescence) and biologically correct. This
is a concrete, low-risk model change, ready to implement in the gompit survival Stan model and
re-check against the Bakuzis U-shape.

## 4b. Tier 0 results: senescence fix validated, regional benchmark scoped

**The senescence fix is validated.** A cloglog mortality model with exposure, base (absolute
size) versus plus a relative-size term, on 400k trees: the relsize-squared coefficient is
overwhelmingly significant (p = 5e-58, delta AIC 1604), and the predicted mortality by
relative-size class is decisive. The base model predicts mortality declining toward the species
maximum (0.0163 down to 0.0111), the exact defect, while the plus-relsize model predicts it
rising (to 0.0226), matching the observed 0.0200. Adding DBH over species-maximum DBH (a
quadratic, trait-mediated since the maximum is a trait) recovers the observed U-shape that the
current absolute-size form inverts. Production step: add relsz and relsz-squared to the gompit
survival eta, refit, re-run the Bakuzis U-shape, and confirm long-horizon basal area and volume
stop over-accumulating.

**The real regional-variant benchmark is proven executable (de-risked).** A smoke test ran the
actual Alaska variant end to end: `from fvs2py import FVS; FVS("lib/FVSak.so").load_keyfile(...).run()`
under Python 3.12 returned a stand-metrics DataFrame (year, age, tpa, cubic-foot volume, mortality,
forest type, size and stocking class) with realistic dynamics, trees per acre declining 100 to 61
from self-thinning while volume rose 109 to 155 from growth. So the real variants execute on
Cardinal today; the only reason the benchmark uses a strawman is that the R engine was wired to one.
The remaining work is purely assembly: `fia_stand_generator.py` already builds FIA stand keyfiles,
so the harness generates keyfiles for spatially held-out plots, runs each real variant forward to
the remeasurement, extracts BA, QMD, density, and volume, and compares to observed and to the
unified model on the identical plots. This is now a wiring task on proven pieces.

**The real regional-variant benchmark is feasible.** The pieces to run the actual variants exist:
`deployment/fvs2py` is a Python package that drives the real FVS Fortran engine (built from
`src-converted`), with per-variant keyfiles (ACD.key, AK.key, BC.key and the rest) and
`run_variant.sh` submit scripts. The current benchmark does not use them; it uses the hand-rolled
`project_condition_default`. The harness to build: for held-out FIA plots, construct the FVS
keyfile and tree list, project forward to the remeasurement under each real variant via fvs2py,
extract stand metrics (BA, QMD, density, volume), and compare to observed and to the unified model
on the identical plots, with spatial blocking. This is well-defined engineering on existing
pieces, not a research unknown, and it is the definitive test of success criterion 1.

## 4c. NE/ACD head-to-head pilot is built and running (first real numbers)

The pilot harness now runs the real NE variant, the real ACD variant, and the unified model on
the same real Northeast FIA stands, via fvs2py and the existing FIA stand generator. First result,
8 NE medium-site, medium-density stands projected 10 years (mean of the projected stand):

| model | mean BA (m2/ha) | mean TPH | mean volume (m3/ha) |
|---|---:|---:|---:|
| real NE variant | 27.1 | 3222 | 150.7 |
| real ACD variant | 27.4 | 3509 | 145.6 |
| unified | 21.2 | 2268 | 126.1 |

Two things stand out. First, the real NE and ACD variants are nearly identical (BA 27.1 vs 27.4,
volume 150.7 vs 145.6), which is the expected sanity check for two neighboring eastern variants.
Second, the unified model diverges clearly: it carries far fewer stems (TPH 2268 vs 3200 plus) and
accumulates less basal area and volume, meaning it removes more trees over the interval than the
real variants do. Whether that is better or worse is exactly the open question, and it is decided
only by the observed remeasurement, which is the immediate next step. Notably, the direction is the
opposite of the senescence concern: here the real variants are the ones retaining a lot of small
stems, and the unified model thins harder.

This is the milestone that matters: the definitive success test now produces real numbers on real
plots. The remaining work to make it a verdict is to attach the observed t2 stand state for these
plots (from the FIA remeasurement) and compute, by variant and overall with spatial blocking,
which model is closest on basal area, density, QMD, and volume. One of the eight stands returned an
empty projection for all three models (a single-plot tree-load edge case) and needs a look. The
pilot script is `calibration/python/ne_acd_pilot.py`.

## 4d. First verdict against observed: no clear win yet (one bug to fix first)

The benchmark now closes the loop: real NE, real ACD, and the unified model projected on 30 real
NE plots forward by each plot's true remeasurement interval, compared to the observed t2 stand
state. Result (percent RMSE / percent bias, n = 30):

| model | BA %RMSE | BA bias | TPH %RMSE | TPH bias | QMD %RMSE | QMD bias |
|---|---:|---:|---:|---:|---:|---:|
| real NE | 63.0 | -10.0 | 384 | +290 | 77 | -59 |
| real ACD | 62.8 | -0.4 | 348 | +255 | 74 | -54 |
| unified | 63.9 | -12.9 | 374 | +271 | 77 | -58 |

Two honest readings, one of them a caution about the harness and one a real signal.

First, the caution. Density (TPH) is off by roughly a factor of four for **all three** models
(+255 to +290 percent bias), and QMD is correspondingly too small (-54 to -59 percent). Because
the error is nearly identical across the real variants and the unified model, it is not a model
difference, it is a harness problem: the initial tree expansion or the FVS plot-design fields
(inv_plot_size, basal_area_factor) are scaling the loaded density wrong, so every projection starts
over-dense. That must be fixed before the TPH and QMD numbers mean anything.

Second, the signal that does survive the caution. On basal area, which is less sensitive to the
count scaling, all three models are essentially tied at about 63 percent RMSE, and **the unified
model is not better than the real variants, it is marginally worse** (BA bias -12.9 percent versus
the real ACD variant at -0.4 percent, which is nearly unbiased). This is the first direct
measurement against the actual success criterion, and the early answer is the one the stress test
warned about: the unified model does not yet outperform the real NE and ACD variants. It is
consistent with the composite pctRMSE finding that the unified model trails on roughly half the
variants.

This is exactly why the benchmark had to be built. The program now has, for the first time, a
direct and honest measurement against real regional variants, and it says the central claim is not
yet met. That makes the Tier 0 fixes (senescence, and whatever else the diagnosis surfaces) and the
density-scaling bug the immediate priorities, with the benchmark as the scoreboard. The pieces:
`stand_level/export_ne_bench.R` (exports the plots) and `calibration/python/ne_acd_observed.py`
(runs the three models against observed).

## 4e. The benchmark scaling issue is diagnosed (FIA plot expansion), verdict pending the fix

Adding an initial-state check (FVS year 0 versus observed t1) localized the problem to input
loading, not projection: at the original setting FVS started every stand about 4x too dense
(TPH year-0 bias +329%). Setting num_plots to 4 (the FIA subplot count) fixed density almost
exactly (TPH year-0 bias +7%, first plot FVS 6899 vs observed 6938), and on the density metric the
unified model then came out slightly best (TPH RMSE 48.8% vs 53.3%, smallest bias). But the same
divisor drove basal area 4x too low (BA year-0 bias -79%).

That contradiction is the real diagnosis: the FIA per-tree expansion is size-dependent (microplot
small stems carry a large expansion, subplot large trees a small one), so the four times too many
trees at num_plots=1 are mostly small and inflate the count without much affecting basal area,
which is why BA looked roughly right while TPH was 4x high. A single global plot count cannot
reconcile both. The correct fix is to honor the FIA microplot and subplot expansion per tree when
building the FVS tree list (the standard FIA to FVS plot-design handling), not a flat divisor.

So the verdict is not final: the benchmark pipeline is complete and correct in structure (real NE,
real ACD, unified, projected to observed), but the FIA-to-FVS expansion must be wired per tree
before the basal area and density numbers can be trusted together. This is a bounded, well
understood calibration task, not a flaw in the comparison itself. The honest interim reading is
unchanged: across the configurations tried, the unified model is at best tied and not yet clearly
beating the real NE and ACD variants, with a mildly encouraging density result once scaling is
roughly right. The next step is the per-tree FIA expansion fix, then a clean re-run, then iterate
the unified model (starting with the senescence fix) against this scoreboard.

## 4f. Root cause isolated: a unit inconsistency in the pairs, and the clean path

A direct input check (does the exported tree list reproduce the observed t1 stand metrics, before
FVS runs) localized the basal-area problem exactly. On plots where the tree count reproduces the
observed density exactly (for example plot 34-1-29: input 631 equals observed 631 trees per
hectare), the basal area is still off by about 6.5x, which is 2.54 squared. That is unambiguous: the
per-tree DBH in the remeasurement pairs and the plot-level basal area in the same file are not on the
same unit (a centimeter versus inch inconsistency). Density reproduces, basal area is off by exactly
the inch-to-centimeter factor squared, so the pairs-derived FVS input is unreliable for basal area.

This explains the whole saga. The original pilot, which sourced its stands from `fia_stand_generator`
(the full FIA tree list, in the documented FIA units of inches and feet), produced sensible basal
areas of 11 to 70 m2/ha. The observed comparison broke only when I reconstructed the input from the
matched-pair tree records, which carry the unit inconsistency. So the input source is the issue, not
the harness or the models.

The clean path is now definite: build the benchmark input from `fia_stand_generator` (full FIA tree
lists with correct units, already proven in the pilot to give sensible stands) and obtain the
observed t2 by matching those plots to their FIA remeasurement, rather than reconstructing stands from
the pairs. With that, the year-0 stand will reproduce observed t1 in both density and basal area, and
the projected-versus-observed verdict becomes trustworthy. The senescence fix and the rest then
iterate against a clean scoreboard.

This is the right place to pause the benchmark thread: the pipeline is built and the two blockers are
now precisely understood (the FIA per-tree expansion handling and this unit source), both bounded
data-engineering fixes rather than open questions. The interim scientific reading stands unchanged and
honest: the unified model is not yet shown to beat the real NE and ACD variants, and proving or
disproving that is exactly what the cleaned benchmark will do.

## 4g. Clean all-FIA benchmark built; a robust signal emerges (unified over-thins)

Rebuilt the benchmark entirely from FIA, using FIA's own remeasurement linkage (PREV_PLT_CN) so
both the input tree list and the observed t2 come from FIA in consistent units, sidestepping the
pairs unit inconsistency and the control-number precision problem. This is the correct architecture
and it runs: 20 remeasured NE plots, t1 trees in, projected to observed t2.

The absolute scaling still has one bounded issue: FVS's single plot-expansion setting cannot
reproduce FIA's variable-radius design (small trees on the 6.8 ft microplot, large trees on the
24 ft subplot). With inv_plot_size at its default of 6.0 the density inflates ~5x; setting it to
1.0 brings density close (TPH year-0 bias -15%) but then basal area loads about a quarter low.
Neither setting gets both right, which is the signature of the unrepresented microplot/subplot
expansion. This is a known FVS-to-FIA calibration task, and it hits all three models identically,
so it does not bias the comparison between them, only the absolute level.

And the comparison, run with density roughly scaled, gives a consistent and now robust signal:

| model | TPH %RMSE | TPH bias | BA bias |
|---|---:|---:|---:|
| real NE | 30.7 | -12.3 | -57.6 |
| real ACD | 32.8 | -14.5 | -56.5 |
| unified | 75.5 | -44.7 | -68.7 |

The unified model **over-thins**. On density it is far worse than the real variants (TPH bias
-45% versus -13%, RMSE 76% versus 31%), removing far too many trees, and it carries the largest
basal-area shortfall. Across every benchmark configuration tried this session, the conclusion has
been the same: the unified model does not beat the real NE and ACD variants, and the specific
failure mode is now clear, it kills too many trees over the interval. That is a mortality
calibration problem, and it points the model-improvement work precisely: the unified mortality is
too aggressive overall, on top of the separate large-tree senescence gap diagnosed earlier (which
is the opposite-direction error, too few large-tree deaths). Net, the mortality model removes the
wrong trees, too many small and medium stems and too few large ones.

So the extended session lands here: the benchmark is built on the right foundation, the remaining
absolute-scaling fix is bounded and understood, and the scientific verdict is consistent and
actionable. The unified model's path to beating the regional variants runs through fixing its
mortality, both the over-thinning and the senescence shape, then re-scoring against this benchmark.

## 4h. Capstone finding: the engine only applies multipliers, not the species-free equations

Following the over-thinning to its source produced the most important strategic finding of the
session. The "unified" config the engine runs (`config/calibrated/ne.json`) is not the trait-driven
species-free model at all. It is a set of **per-species multipliers** (MORTMULT, BAIMULT, SDIMAX)
applied on top of the native NE variant's own equations. The keyword generator emits only those
multiplier keywords. The separate species-free preview config (`ne.sf_preview.json`, the one with
the `categories_conus_sf` blocks) generates the **same kind** of keywords, and its mortality
multipliers are byte-for-byte identical to the calibrated config (mean 1.16, median 1.00, 19 percent
above 1.5x, max 6.7x). Through the engine, the species-free model and the multiplier calibration are
the same thing.

Two consequences, both central to the program:

1. **The CONUS species-free architecture is not yet realized in the runnable engine.** The goal,
   one trait-driven equation set replacing the 20 regional variants, requires injecting the actual
   species-free growth, mortality, crown, and height equations as the engine's functions. The
   current integration cannot do that; it can only scale the native regional equations by
   per-species factors. So what runs today is per-variant multiplier tuning of the old variants,
   not the unified model. Testing and realizing the unified variant requires completing the
   equation-injection integration, which is exactly the incomplete PR #70 work.
2. **The multiplier calibration over-thins because it raises mortality.** The mean mortality
   multiplier is 1.16 (16 percent higher than the native variant) with a long tail to 6.7x. That
   is the direct cause of the benchmark over-thinning, and it explains why the "unified" run loses
   to the native variants: it is the native variants plus too much mortality.

So the extended session converges on a clear, honest, and strategically important picture. The
benchmark scoreboard is built and trustworthy for comparison. The thing it has been scoring is the
multiplier calibration, which over-thins and does not beat the native variants. The actual
species-free CONUS model cannot yet be scored because it is not injected into the engine as
equations. The critical path to a CONUS variant that beats the regionals is therefore: (a) complete
the species-free equation injection into FVS so the trait-driven forms drive the engine, not
multipliers; (b) fix the mortality, both the over-aggressive overall level and the senescence shape;
then (c) score the real thing against this benchmark. The benchmark, the senescence fix, and the
integration scaffold built this session are exactly the pieces that critical path needs.

## 4i. Over-thinning root cause found by controlled experiment: SDIMAX set too low

A clean controlled experiment settled the over-thinning. Dampening the per-species mortality
multipliers had zero effect (bit-identical results), ruling them out. Removing the SDIMAX keyword
(the calibrated maximum stand density index) made the unified config snap to the native NE variant
almost exactly:

| run | TPH bias | BA bias |
|---|---:|---:|
| unified, as calibrated | -44.7% | -68.7% |
| unified, SDIMAX removed | -12.7% | -57.8% |
| native NE variant | -12.3% | -57.6% |

So the entire over-thinning was caused by the calibrated SDIMAX being too low: FVS self-thins the
stand down to that maximum density, and a too-low maximum forces excess density-dependent mortality.
This is a calibration bug in the SDImax values, not a flaw in the demographic equations, and it has a
direct fix (recalibrate or rescale SDImax).

Two important corollaries close the loop on the whole benchmark investigation:

1. **Once SDIMAX is corrected, the multiplier-based "unified" config reproduces the native variant**
   (the mortality and growth multipliers are near-neutral in net effect here). So the current
   engine integration, fixed, does not beat the native variants, it matches them. That is consistent
   with the capstone finding: the multiplier path cannot improve on the native equations, only scale
   them. Beating the regionals requires injecting the species-free equations themselves.
2. **The residual basal-area bias (about -57% for every model) is the FIA microplot/subplot
   plot-expansion issue**, affecting all models identically, separate from SDIMAX.

The benchmark investigation is therefore complete and fully explained: the harness is correct, the
absolute BA level needs the FIA variable-radius expansion fix, the over-thinning was a too-low
SDIMAX (now isolated and fixable), and the multiplier integration at best matches the native
variants. The path to a CONUS variant that beats them is unchanged and now very precisely motivated:
fix SDImax, complete the species-free equation injection so the trait-driven forms drive the engine,
fix the mortality senescence shape, fix the FIA expansion in the harness, then re-score. Every one of
those is a bounded, identified task, and the scoreboard to measure them on is built.

## 5. Secondary risks

- **Height measurement error is systemic.** We showed top-height growth is unfittable from FIA
  because the estimator noise exceeds the increment. HG, HT-DBH, HCB, and CR all depend on tree
  height, which carries the same error, so crown and height dynamics may be weakly identified.
- **A few regions fail Eichhorn:** growth does not increase with site in CR, SN, TT, UT (the
  growth-site response is flat or inverted there).
- **Validation honesty:** CSPI spatial CV is weak (0.42), and the eventual regional-variant
  comparison must use spatial blocking, not random holdout, or it will flatter the unified model.
- **Closure is untested:** no end-to-end forward stand simulation has been shown to reproduce
  realistic basal area, density, composition, and volume trajectories. Full-projection path
  invariance (all components composed) is also untested even though components are annualized.
- The CR2 annualization comparison job failed at launch (4 seconds) and needs a quick debug.

## 6. The honest bottom line

We have spent recent effort on real but secondary refinements: CSPI version selection, the
species blend, García-style annualization, stand-level prototypes. These are sound. But they
polish a model whose core value proposition is undemonstrated: there is no test against real
regional variants, the stand-level accuracy is not clearly better than a strawman, and mortality
realism is broken in a way that biases long-horizon yield. The risk is an elegant, well-engineered
model that does not actually win on the two things that define success. The plan must pivot to
prove those two things first.

## 7. Refined, re-prioritized plan

**Tier 0, the make-or-break (do these before more refinement):**

1. **Build the real regional-variant benchmark.** Replace the strawman default with the actual
   FVS regional variant logic (the fvs-modern engine runs it) and run the unified model against
   each real variant on held-out FIA plots, stand level (BA, QMD, density, volume), by region,
   with spatial blocking. This is the definition of success criterion 1. Report per-variant and
   overall. If the unified model loses per-variant, that is the central problem to solve, not a
   detail.
2. **Fix mortality senescence.** Add a large-tree or old-age hazard so mortality is U-shaped
   (raise the right arm of the U). Re-run the Bakuzis U-shape and confirm the long-horizon basal
   area and volume stop over-accumulating. This likely also improves finding 3.

**Tier 1, diagnose and close:**

3. **Diagnose the per-variant stand-level underperformance** once senescence is fixed: is the
   residual gap from the components, initialization, horizon, or the benchmark itself (the
   above-100% pctRMSE suggests the benchmark or horizon also needs scrutiny).
4. **Fix the Eichhorn-failing regions** (CR, SN, TT, UT growth-site response).
5. **Demonstrate closure:** a forward stand simulation reproducing FIA basal area, density,
   composition, and volume trajectories, with the stand-level density and basal-area constraints
   turned on, and a full-projection path-invariance check.

**Tier 2, the refinements already in progress (keep, but downstream of Tier 0):**

6. Height measurement-error refit of HG and HT-DBH; derive top height from tree-level HG.
7. The species blend, CSPI v7 site term (pending the running ΔLOO), annualized CR2, the coupled
   stand state-space module, and the engine integration and uncertainty work.

**The reframed success test, stated once:** the unified CONUS variant succeeds when, on held-out
FIA plots with spatial blocking, it matches or beats each real regional FVS variant on stand-level
basal area, density, QMD, and volume, and passes all four Bakuzis laws including a U-shaped
mortality. Everything else is in service of that, and that test has not yet been run.

## 8. What is strong enough to keep building on

To be clear this is not a teardown: the trait-substitution premise is validated, the component
fits are decent, most biological laws pass, and the engineering is solid. The model is a credible
foundation. It simply has not yet been measured against the bar that defines the project, and one
biological defect (mortality senescence) must be fixed before any long-horizon claim holds. Fix
those two things and the refinements already underway become the difference between a model that
ties the regional variants and one that clearly beats them.
