# Gompit mortality, Fortran in-engine integration: NE validation

Greg Johnson's gompit survival is now substituted for FVS native mortality
**inside the FVS growth loop** (TREGRO -> MORTS), so growth, density, and the
substituted mortality interact each cycle. This replaces the post-hoc TPA
overlay, which was rejected because disabling FVS mortality made FVS growth
unrealistic (runaway density). See `GOMPIT_INTEGRATION_FINDINGS.md`.

## Validation (NE, 8 stands, 100 yr, 20 x 5 yr, same engine, env-toggled)

Mean AGB (short tons/ac), native FVS mortality vs gompit-in-FVS:

| proj yr | native FVS | gompit-in-FVS |
|--------:|-----------:|--------------:|
| 0   | 23.0  | 23.0  |
| 25  | 81.7  | 86.2  |
| 50  | 124.5 | 117.6 |
| 100 | 173.2 | 147.0 |

## Why this is the right result

* **Bounded and realistic.** 147 t/ac at 100 yr, no runaway. Contrast the
  post-hoc overlay (923 t/ac) and the cycle-by-cycle re-entry overlay (crashed
  on over-dense stands). Letting FVS keep doing growth while gompit owns
  mortality each cycle is what fixes it.
* **The crown-closure signal shows up dynamically.** Gompit kills slightly
  *less* than native early (86 vs 82 at yr 25, young/open stands, low crown
  closure) and *more* late (147 vs 173 at yr 100, older/crowded stands, high
  crown closure). That is exactly the cch behaviour the model-level validation
  found: gompit captures crowding mortality the base rate misses.
* **Identical start (23.0).** Same initial cohort; the arms diverge only through
  the mortality model, as intended.

## How it works

`src-converted/base/gompmort.f90` + `base/common/GOMPMC.f90`, hooked into
`vls/morts.f90` (shared by NE/CS/LS):

1. `GOMPLOAD` (called from MORCON) reads `FVS_GOMPIT` / `FVS_GOMPIT_COEF`,
   loads the 133-species coefficients, and resolves them onto the variant's
   species via `FIAJSP`.
2. each cycle, `GOMPCCH` rebuilds per-tree crown closure at tip (ORGANON crown
   port, affine-mapped to the gompit scale) from the live treelist.
3. in the MORTS tree loop, fitted species get `WK2 = PROB*(1 - exp(-H*FINT))`
   with `H = exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4)`; unfit species keep native
   background; the native SDI redistribution (VARMRT) is bypassed (gompit is
   density-aware via cch).

Activation is env-gated (`FVS_GOMPIT=1`, `FVS_GOMPIT_COEF=<csv>`), so the same
binary does native or gompit with no recompile, giving a clean A/B.

## Aaron's validated assumptions (carried into Fortran unchanged)

Gompit coefficients; ORGANON SWO crown geometry as the cch proxy; coarse
FIA->ORGANON group map (softwood->1 DF, hardwood->16 RA); affine map
CCH=0.062+0.0036*cch_hat (Spearman 0.93). Refine the group map for a tighter fit.

## Multi-variant: NE, CS, LS (all share vls/morts.f90)

Adding `base/gompmort.f90` to the CS and LS source lists was sufficient (they
already use `vls/morts.f90`); both built and validated with no further code
changes. Mean AGB at projection year 100 (10 stands each, same binary,
env-toggled):

| variant | native FVS | gompit-in-FVS | change |
|---------|-----------:|--------------:|-------:|
| NE      | 175.7 | 160.8 | -8%  |
| CS      | 166.5 | 124.0 | -26% |
| LS      | 155.6 |  93.6 | -40% |

All bounded, realistic, no runaway/crash. Gompit consistently trims
late-rotation stocking, variant-specific in magnitude (LS most, NE least),
coherent with the differing species mixes and crowding regimes. See
`gompit_fvs_inengine.png` (`calibration/R/41_gompit_fvs_inengine_figure.R`).

![gompit in-engine, NE/CS/LS](gompit_fvs_inengine.png)

## Next steps

* Other variant families use their own `morts.f90` — same hook pattern, separate
  work (the gompit module is variant-agnostic).
* Replace the env switch with a `GOMPMORT` keyword for user/run reproducibility.
* Scale the NE/CS/LS in-engine A/B across the full CONUS stress sample.
