# Offline constrained projector stress test -- findings (2026-07-04)

Branch reconcile/three-arm-onto-main (detached wt-stress1). 40 project() runs +
25 direct solver probes + mechanism probes. Files in this dir.

## BROKE
1. All-zero-TPA stand (01d): TOPHT = NaN at every year, both modes. BA/TPH/QMD
   correctly 0. Cause: stand_top_height denom==0 -> nan, propagates to quantiles
   and topht_ratchet. Severity: SILENT-WRONG (no crash/warning).
2. Top-height per-draw NON-MONOTONICITY (05a/05b/06a + jitter in 04a/04c).
   At n_draws=1 (q50==realized draw) top height declines on ~28-30 of 60 steps
   over 300 yr; worst interior step -0.67 m (normal), -0.48 m (high site); 05b
   ENDS below start in that draw. Under fixed_sdimax=100 the year-5 drop is
   -1.17 m. Cause: the Track B/D "monotone ratchet" floors only the TARGET fed
   to stand_constrain_topheight; the realized snapshot top height is recomputed
   after the mortality/stems reconciler + density cap RECOMPOSE the top-100/ha
   dominant cohort, and nothing guards the realized value. Confirmed: declines
   persist with the density cap effectively off (fixed_sdimax=1e6, 20 declines,
   worst -0.54 m) -> primary driver is the stems/mortality cohort reshuffle, not
   the SDI cap. Severity: SILENT-WRONG (violates the stated dominant-height
   monotonicity invariant; magnitude ~0.5-1.2 m, recovers, never crashes).

## HELD (finite + bounded, no crash)
- empty list, single tree, two identical, all-zero-CR, only-largest-class
- ultra-dense (5000 sph, SDI>>SDIMAX), ultra-sparse (5 sph)
- all-saplings, old-growth ABOVE the GADA asymptote (H1>b1 held: envelope returns
  H1, no inversion blow-up), init top height == asymptote
- bgi 0.5 / 20 / -3 (beyond fitted range), trt just-happened / very-large
- 300 yr normal + high site: BA bounded, TPH self-thins to a positive floor
  (no zero-division, no negative), top height asymptotes near site b1 in the
  cross-draw median (23.5 m high site) -- aside from the per-draw wobble above
- fixed_sdimax 100 (low) and 2000 (high): finite, bounded (100 forces heavy
  early thinning + the -1.17 m top-height dip noted above; not a crash)
- kappa solve at M->0, M->1, M=1 exact, M=5, zero-haz, tiny-haz-big-M, neg-M:
  all converge or clamp gracefully (max kappa ~448 or 1.4e6, mortality < 1)
- stems N->0, N>start, zero-start; topheight unreachable/below-floor; GADA H1>b1,
  H1==b1, H1~b1-eps, tiny H1, zero years; sf_site_asymptote monotone+clamped;
  density cap SDI=0 / QMD=0 / huge SDI / SDIMAX=0 -- all finite.

## MINOR
- gada_tiny_H1: H1=0.01 m -> H2=3.05 m in one 10-yr step (300x). Finite but a
  large relative jump near the origin; only reachable with sub-cm top height.

## ADDITIONAL (from the results CSV)
3. EMPTY tree list (01a): hard CRASH, both modes. IndexError: index -1 out of
   bounds for axis 0 size 0, in stand_top_height (cum[-1] on an empty cohort),
   reached via _stand_state at year 0. Severity: CRASH (unguarded n==0).
4. Constrained BA sits ~90-130 m2/ha at horizon end across MANY cases (baseline
   93, ultra-sparse 129, saplings 91, old-growth 127). The BA-growth target
   clamp (1.2 m2/ha/yr) + SDIMAX cap admit BA far above realistic NE closed-
   canopy maxima (~50-60 m2/ha). Not NaN/negative, but the constrained band is
   biologically high -- worth a look. Ultra-sparse constrained QMD hits 664 cm
   (single surviving huge tree); finite but nonsensical. Severity:
   SILENT-WRONG / plausibility, not a numerical break.
5. Unconstrained top height blows up with no asymptote (300yr normal 174 m,
   high-site 202 m, saplings 116 m). EXPECTED -- this is the pathology the
   constraints exist to fix; recorded as context, not a defect of the projector.
