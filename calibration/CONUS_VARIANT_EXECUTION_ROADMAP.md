# CONUS variant execution roadmap: from diagnosis to a variant that beats the regionals

**Date:** 2026-06-11
**Purpose:** turn the full stress test and benchmark investigation into a prioritized, bounded
execution plan. Every item below is a defined task with a known cause and a verification gate (the
NE/ACD benchmark, now built). This is the actionable synthesis of the session.

## Where things actually stand

The benchmark to measure success (unified model vs the real regional variants on held-out FIA
remeasurement, stand level) is built and runs. Using it, the diagnostic chain is complete:

- The species-free CONUS equations are **not injected into the runnable engine**. The engine
  applies all calibration as per-species multipliers (MORTMULT, BAIMULT, SDIMAX) on the native
  regional equations. The "calibrated" and "species-free preview" configs reduce to the same
  multipliers. So today's runnable model is per-variant multiplier tuning, not the unified model.
- The current multiplier config **over-thins**, and the cause is isolated by controlled
  experiment: a too-low SDIMAX (mortality multipliers had zero effect; removing SDIMAX made the
  unified run match the native variant). Fixable by recalibrating SDImax.
- Once SDIMAX is corrected, the multiplier config **reproduces the native variant**, it does not
  beat it. Beating the regionals requires the trait-driven equations to drive the engine.
- The benchmark's **absolute** basal-area level is biased low for all models because the FIA
  variable-radius design (microplot vs subplot vs macroplot expansion) is not reconciled by FVS's
  single plot-size factor. This hits all models identically, so the relative comparison is valid;
  the absolute level needs the expansion fix.
- The species-free mortality has a **senescence-shape** error (large-tree hazard collapses to near
  zero; the validated fix is a relative-size term, DBH over species-maximum DBH).

## The roadmap (ordered by leverage and dependency)

### Tier 0 — make the benchmark give a clean, trustworthy absolute verdict
1. **Fix the FIA-to-FVS plot expansion** in the harness so year-0 reproduces observed t1 in both
   density and basal area. The single inv_plot_size cannot represent FIA's nested fixed-radius
   plots; the fix is to set the FVS plot design per FIA's microplot/subplot/macroplot radii (or
   pre-expand per tree to true per-acre and tell FVS not to re-expand). Bounded FVS-internals task.
   Gate: FVS year-0 BA and TPH within a few percent of FIA observed t1.
2. **Recalibrate SDIMAX.** The current values are too low and drive excess self-thinning. Either
   refit the maximum SDI or rescale; verify against the size-density boundary. Gate: the multiplier
   config's TPH bias matches the native variant (about -12 percent on the current NE sample).

### Tier 1 — the actual differentiator (the make-or-break)
3. **Inject the species-free equations into the engine.** This is the central, largest task and the
   only path to beating the regionals. The trait-driven DG (Kuehne), HG (ORGANON v8rd), HT-DBH,
   HCB, survival, CR, and ingrowth must drive FVS as functions, not as multipliers on the native
   equations. Two routes: extend the FVS Fortran to evaluate the species-free forms (most faithful,
   largest effort), or run the species-free projection in the existing R/Python engine
   (17_stand_projection_engine.R) and benchmark that directly against the native FVS variants. The
   second route is faster to a verdict and should come first. Gate: species-free projection scored
   against native NE/ACD and observed on the clean benchmark.
4. **Fix the mortality.** Apply the validated relative-size senescence term (raises large-tree
   hazard to the observed ~1 to 2 percent floor and recovers the U-shape), and confirm the overall
   level is not over-aggressive once SDIMAX is right. Gate: Bakuzis U-shape passes and stand-level
   mortality matches observed.

### Tier 2 — extend and harden
5. **Scale the benchmark** beyond NE/ACD to all variants and add spatial blocking, so the verdict
   is national and honest, not a 20-plot NE pilot.
6. **Confirm trait substitution and the species blend** on the clean benchmark (the held-out-species
   test already supports the trait premise; verify it survives at the stand level).
7. **Wire uncertainty** (posterior draws plus CSPI v7 QRF site draws) once the point verdict is
   favorable, for credible bands on the projections.

## The one-line critical path

Fix the harness expansion and SDIMAX so the benchmark is clean (Tier 0), then get the species-free
equations driving a projection and score them against the native variants (Tier 1, item 3) with the
senescence fix in place (item 4). That sequence is the difference between a model that ties the
regional variants and one that beats them, and the scoreboard to prove it is built.

## What this session delivered toward it

The benchmark harness (real variants via fvs2py, FIA remeasurement linkage), the validated
senescence fix, the isolated SDIMAX cause, the capstone finding that the engine integration is
multiplier-only, the García-grounded equation refinements (annualized CR2, SDImax-driven basal
area, density self-thinning), and this roadmap. Every remaining task is bounded and has a gate.
