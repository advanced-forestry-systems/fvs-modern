# Why AcadianGY.R standalone output does not match FVS-ACD

Aaron Weiskittel | 2026-05-21 | branch context: holoros/fvs-modern, acd-bridge-followup

## Bottom line

The Acadian growth equations themselves are not the problem. The fork's
`fvsOL/inst/extdata/AcadianGY.R` is byte for byte identical to the
`AcadianGY_12.3.5.r` you uploaded (verified: 0 differing lines, both 3,968
lines, same `AcadianVersionTag = "AcadianV12.3.5"`). So a divergence between
"the R code" and "FVS-ACD" cannot come from different equations or coefficients.

The gap is created entirely in the bridge layer, `customRun_fvsRunAcadian.R`,
which is the code FVS uses to drive the Acadian R model one cycle at a time.
That wrapper does three things that change the answer relative to a clean
standalone `AcadianGYOneStand()` call: it injects FVS species calibration
multipliers into the Acadian increments, it changes how height, crown ratio and
climate site index are initialized, and it overrides several run options. Each
of these is examined below, ranked by expected magnitude.

Your two observations are both confirmed by reading the code. Annual versus
periodic stepping is not the primary driver, because the wrapper already loops
annually inside each FVS cycle (`for (year in stdInfo["year"]:stdInfo["cendyear"])`,
customRun line 226, with `cyclen = 1` hardwired at AcadianGY line 1874).
Ingrowth is real but secondary, for the reasons in Driver E.

## Ranked drivers

### A. FVS calibration multipliers are injected, and on the wrong scale (highest)

In a clean standalone run, `dDBH.mult`, `dHt.mult` and `mort.mult` all default
to 1 (AcadianGY lines 1390 to 1392), so the raw Kuehne increments are used.

Under FVS, the wrapper calls `make_fvs_calib()` (AcadianGY lines 1327 to 1397),
which pulls the live species attributes out of the loaded FVS variant:

```
calib.fvs = fvsGetSpeciesAttrs(c("baimult","htgmult","mortmult","mortdia1",
                                 "mortdia2","maxdbh","maxht","minmort","maxdbhcd"))
```

and then renames them (lines 1367 to 1369):

```
dplyr::rename(dDBH.mult = baimult,
              dHt.mult  = htgmult,
              mort.mult = mortmult, ...)
```

These get multiplied straight onto the annual increments inside
`AcadianGYOneStand()`:

```
dDBH = dDBH * dDBH.thin.mod * dDBH.SBW.mod * dDBH.form.risk.mod * dDBH.mult   # line 2232
dHT  = dHT  * dHT.SBW.mod  * dHT.thin.mod  * dHt.mult                          # line 2268
```

There are three issues here, in order of size:

1. The multiplier is non unity in the bridge and unity standalone (largest).
   Whenever the FVS-ACD variant has been calibrated (and after our bridge work
   the NE based ACD does carry calibration), `baimult`, `htgmult` and `mortmult`
   are not 1, while the standalone run uses 1, so the two runs scale the same
   Acadian increments by different factors. The real model attribution on
   Cardinal (RESULTS_attribution.md) shows multipliers in the 0.88 to 1.18 range
   moving total basal area about 3 percent at 25 years but individual species up
   to 9 percent (red spruce) and down 5 percent (yellow birch). So the level
   effect on the stand is modest but the species composition effect is large.

2. The multipliers come from the wrong model (conceptual). `baimult`, `htgmult`
   and `mortmult` were fit for the FVS-NE Fortran growth model, not the Acadian
   R equations. Injecting them into the Acadian increments mixes two calibration
   systems. The Acadian model should carry its own calibration, fit on its own
   annual per species growth, which is what `annualized_calibration.R` produces.

3. The scale of `baimult` is wrong, but this is second order (smallest).
   `baimult` is a basal area increment multiplier (it scales DDS, the change in
   squared inside bark diameter), yet it is applied directly to a linear
   diameter increment `dDBH`. The correct conversion makes the realized basal
   area increment equal to m times the unscaled increment:

       dDBH_new = sqrt(DBH^2 + m * [(DBH + dDBH)^2 - DBH^2]) - DBH

   `bridge_patch_verify.R` quantifies the legacy error: applying m to the linear
   increment instead drifts the realized basal area factor from m by at most
   about 1.5 percent, even for fast growing small trees, because in a single
   year dDBH/DBH is small. So this is a correctness refinement worth roughly 1.5
   percent, not the main driver. It still matters because the error always has
   the same sign (it overshoots for m above 1 and undershoots below 1) and
   accumulates across years.

The same scale concern applies to `mortmult`, which in FVS modifies a survival
or mortality rate, but is here applied multiplicatively to a per year Acadian
mortality that is itself derived from a stand basal area mortality allocation
(AcadianGY lines 2300 to 2392).

### B. Height and crown ratio are initialized differently

The Acadian dubbing of missing heights and crown ratios is commented out in the
wrapper (customRun lines 95 to 140), with the note "commented out code for
height and crown ratio imputation until errors are resolved" and the matching
changelog line "temporarily disabled height and crown ratio imputation due to
problems with FVS stop point 7."

Consequence: under FVS, any tree that arrives without a measured height or crown
ratio gets a height from the FVS NE base dubbing and a crown ratio from FVS,
not from the Acadian `HTPred()` and `HCBPred()` functions (AcadianGY lines 2104
to 2110 and 2160 to 2166). A standalone run that lets the Acadian code fill
those gaps therefore starts from a different state. This matters more than it
looks, because crown ratio enters diameter growth as `log(CR)` with a large
coefficient (b3 = 0.713) and height enters height growth directly. Different
starting CR and HT propagate and compound every year.

If your test trees all carry measured HT and CR, this driver is muted. If any
are dubbed, it is a leading contributor.

### C. Climate site index is sourced differently

CSI drives diameter growth (b8*log(CSI)), height growth (b7*CSI^2), the dynamic
SDImax (the 222.78/CSI term, line 2048), the stand basal area constraint, and
ingrowth. The wrapper sets CSI from the FVS Event Monitor, and if that is
missing it derives CSI from the FVS site index in feet through a piecewise map
(customRun lines 169 to 171):

```
CSI = fvsGetEventMonitorVariables("csi")
if (is.na(CSI)) { CSI = fvsGetEventMonitorVariables("site")*FTtoM
                  CSI = approxfun(c(0,8,14,20), c(0,8,12,14), rule=2)(CSI) }
```

A standalone run uses `stand$CSI` directly, or the default of 12 (AcadianGY line
1889). If the CSI you pass standalone is not the same number the site index map
produces inside FVS, every CSI dependent term shifts, in the same direction for
all trees, which reads as a level bias. This is a common and easily overlooked
source of a clean offset between the two runs.

### D. Run options the wrapper overrides

The wrapper does not use the standalone defaults. Three options differ:

- INGROWTH. Standalone default is "Y" (AcadianGY line 1871). The wrapper default
  is "N" (customRun line 29) and the UI default is "No" (line 347). A standalone
  run grows recruits; a default FVS run may not.
- CutPoint. Standalone default 0.95 (line 1875). The wrapper sets 0.5
  (customRun line 213), which roughly doubles the rate at which the ingrowth
  threshold is crossed when ingrowth is on.
- mortType. The wrapper passes "continuous" (customRun line 214), but
  `AcadianGYOneStand()` hardwires "discrete" in 12.3.5 (line 1876), so the
  wrapper's intent is silently ignored. Both end up discrete, but this is a
  latent inconsistency worth removing so the two code paths cannot drift.

### E. Ingrowth is added at a different time and through a different path

Even with ingrowth on in both, the placement differs. Standalone with
INGROWTH = "Y" adds recruits inside `AcadianGYOneStand()` every simulated year.
The wrapper instead suppresses per year R ingrowth and adds regeneration once
per FVS cycle through `make_fvs_regen()` and `fvsAddTrees()` at stop point 6
(customRun lines 279 to 297), after the annual loop has finished. Same nominal
model, different timing and a different competitive context for the new trees.
Over several cycles this changes TPH and BA. This is why you correctly sensed
ingrowth contributes but is not the main story: it is a timing and path
difference layered on top of A, B and C, not an equation difference.

### F. Maximum size caps come from a different table

Standalone uses the champion tree caps in `tree.size.cap` (lines 1272 to 1324).
The wrapper coalesces FVS `maxdbh` and `maxht` from the TreeSZCp keyword over
those caps (lines 1384 to 1389). Growth is hard clipped to zero at the cap
(`dDBH = 0` when `DBH + dDBH > max.dbh`, line 2239, and the height analog at
2274). If FVS supplies different caps, large trees stop growing at different
diameters. This is a tail effect, important only in older or larger stands.

### G. Annual versus periodic is not a driver here (confirms your read)

The wrapper iterates the Acadian model annually within each FVS cycle and only
hands the accumulated diameter, height and mortality back to FVS once per cycle
(`dg = (DBH - dbh.fvs)*CMtoIN`, etc., in `make_fvs_tree`, lines 1780 to 1787).
So the internal time step is annual on both sides. Annual versus periodic would
only matter if you were comparing the Acadian R model against the native Fortran
NE/ACD variant, which integrates on the FVS cycle. If that is in fact your
comparison target, see the note below, because then the answer is different and
larger.

### H. Mortality model and handoff

`mortModel` defaults to "Acadian" (customRun line 33), so R mortality is passed
to FVS (lines 1786 to 1803). If a run is set to "Base Model" mortality, FVS
applies its own NE mortality and the Acadian mortality is dropped, which the
standalone never does.

## One thing to confirm: which "FVS-ACD" are you comparing against?

There are two distinct things that get called FVS-ACD in this fork:

1. The customRun path: FVS framework plus `AcadianGY.R` driven by
   `customRun_fvsRunAcadian.R`. This runs the real Acadian equations. Everything
   above (A through H) explains the gap for this case, and the gap is reconcilable.

2. The native compiled variant. In this fork, ACD is carried as an NE subvariant
   (the relabel work from the bridge branch). That compiled variant runs the
   Northeastern Fortran growth equations, not the Acadian Kuehne equations at
   all. If your standalone R is being compared against the compiled FVSacd
   binary, the two are simply different models and will never match by tuning;
   they would only be reconciled by routing ACD through the customRun path or by
   porting the Acadian equations into Fortran.

The fastest way to settle this is the decomposition harness
(`acd_divergence_decomposition.R`): run the customRun path on one stand, run
`AcadianGYOneStand()` standalone on the same tree list, then flip drivers A
through F one at a time and watch which one closes the gap. If neither closes it,
you are in case 2.

## Recommended fixes, in order

1. Calibrate the Acadian model on its own annual per species basis and set
   `dDBH.mult`, `dHt.mult` and `mort.mult` from that table, decoupled from the
   FVS-NE `baimult`, `htgmult` and `mortmult`. The optimizer in
   `annualized_calibration.R` fits those annual factors against FIA
   remeasurement. This is the direct answer to your annualized calibration
   request and removes the conceptual mixing of two calibration systems. If you
   must still honor a user supplied `baimult`, convert it through basal area
   (`bridge_patch_verify.R` shows the formula and that the legacy linear
   application is off by about 1.5 percent), but treat that as a secondary
   correctness fix, not the main lever.

2. Decide CSI once. Make the standalone CSI and the FVS derived CSI identical for
   any validation comparison, then keep them identical in production. If you trust
   the Acadian CSI, push it into the Event Monitor so the wrapper reads it rather
   than re-deriving from site index.

3. Re-enable Acadian height and crown dubbing, or guarantee every input tree
   carries measured HT and CR, so both paths start from the same state.

4. Align the option defaults (INGROWTH, CutPoint, mortType) between the wrapper
   and the standalone so the two cannot silently diverge, and remove the dead
   "continuous" mortType that is ignored.

5. Standardize the maximum size cap source.
