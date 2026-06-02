# Gompit mortality projection integration: wiring + validation findings

Status: wiring complete, harness validated, **naive substitution rejected by
validation**. A modeling decision is required before any CONUS-scale gompit run.

## What was built

Greg Johnson's CONUS gompit survival is now wired into the perseus projection
path so a stand can be projected with gompit mortality substituted for FVS's
native mortality, and compared against the default and old-calibrated arms.

| file | role |
|------|------|
| `greg_mortality.py` | gompit hazard `H = exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4)`, `period_survival(spcd,cr,cch,T)`. 133 fitted species. |
| `cch_organon.py` | ORGANON crown-closure-at-tip (`cch`) port of `CAL_CCH.for`. |
| `project_mortality.py` | per-cycle glue: crown profile -> tip cch -> validated affine map (`CCH_A=0.062, CCH_B=0.0036`, Spearman 0.93) -> survival. **Fix:** `spp_group` now returns an int group (was a str, which `KeyError`-ed against `cch_organon`'s int-keyed parameter tables — the `cycle_survival` path had never run end to end). |
| `run_gompit_projection.py` | the runner: `default`, `calibrated`, `posthoc`, `iterative` arms; closed-cohort (`--regen off`) toggle; FVS `TreeId`-tracked TPA. |

FVS treelist facts confirmed during this work: the live-tree table carries a
stable `TreeId`/`TreeIndex` per tree across all cycles (trees are trackable for
proper iterative mortality), and `PctCr` (percent) not a 0-1 crown ratio, so the
runner maps `CrRatio = PctCr/100`.

## Validation result (NE, 6 stands, 100 yr, 20x5 yr)

Mean AGB (short tons/ac):

| proj_yr | default | old-calibrated | gompit posthoc | gompit iterative |
|--------:|--------:|---------------:|---------------:|-----------------:|
| 0   | 8.2  | 8.2  | 8.2   | 8.2   |
| 25  | 69.2 | 36.4 | 775   | 142.6 |
| 50  | 111.2| 54.2 | 1481  | 353.1 |
| 100 | 158.3| 76.7 | 923   | crash |

Both gompit arms run away to physically impossible biomass; the iterative arm
also crashes FVS after ~yr50 on over-dense stands.

## Why naive substitution fails (the finding)

To let gompit own mortality, both gompit arms set `MORTMULT 0` (FVS mortality
off). But **FVS growth is calibrated to occur alongside FVS mortality**. With
FVS mortality disabled:

1. the initial cohort grows for 100 yr with no competition mortality, so trees
   get far too large (the single-shot `posthoc` arm: 6x default biomass even
   after the gompit TPA reduction is applied);
2. automatic establishment keeps adding stems and `NOAUTOES` did **not** suppress
   it for NE (byte-identical results with/without), so density compounds; and
3. gompit's fitted hazards are too low to regulate that unchecked density, so
   even the correct cycle-by-cycle `iterative` re-entry (which feeds the
   gompit-thinned stand into the next cycle, inventory year advanced) still
   balloons and eventually drives FVS to fail on impossibly dense stands.

The TPA bookkeeping is correct (TreeId-tracked, cumulative survival); the problem
is upstream — you cannot remove FVS mortality and keep FVS growth realistic.

## Options (Aaron's call — this is the mortality science)

1. **Fortran integration.** Replace native mortality inside the FVS growth loop
   so growth, SDImax density mortality, and gompit all interact each cycle.
   Most faithful; largest effort.
2. **Closed-cohort re-entry done right.** Keep the iterative harness but (a)
   find the correct keyword to actually disable NE establishment (NOAUTOES is a
   no-op here) and (b) confirm gompit hazards regulate density, or add the
   SDImax density-mortality term back so density can't run away.
3. **Static pure-demographic comparison (already validated).** The `+41.5% BA`
   result came from holding DBH/Ht fixed and iterating only TPA under each
   mortality model (`project_compare.py`). This isolates the mortality model
   cleanly and is CONUS-cheap, but is demographic, not a growth projection.
4. **Re-examine the gompit fit** if the intended design is mortality-regulated
   density without FVS's density term.

## Recommendation

Do not launch a CONUS gompit campaign on the current wiring — it is invalid.
The harness is correct and reusable; pick a design (1-4) first. Option 3
reproduces the validated manuscript number immediately; option 1 or 2 is needed
for a true gompit growth-and-yield projection.
