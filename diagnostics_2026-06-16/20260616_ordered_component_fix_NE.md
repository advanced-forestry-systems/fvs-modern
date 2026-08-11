# Can we fix FVS via ordered, signed component adjustments? NE test (2026-06-16)

Tested Aaron's plan on NE: updated (localized) max SDI plus signed region-level adjustments to the
component equations, in dependency order allometry -> growth -> mortality. Each lever measured against
observed FIA remeasurement; then a combined config. n = 82 NE remeasurement pairs.

| arm | BA bias | TPH bias | QMD bias |
|---|---|---|---|
| default | +16.8% | -9.3% | +13.1% |
| + allometry (CRNMULT + REGHMULT) | +17.1% | -9.4% | +13.1% |
| + growth (BAIMULT 0.7) | +14.3% | -8.7% | +11.4% |
| + mortality (MORTMULT 0.6) | +16.9% | -9.2% | +13.0% |
| + combined (+ localized max SDI x1.4) | +16.8% | -6.0% | +11.6% |

## What each lever does (measured, not assumed)

- Allometry (crown, height-diameter): no effect on stand BA/TPH/QMD. As shown before, crown and
  height-diameter drive volume and biomass, not basal area. They remain the right foundation (and
  matter for the carbon application), but they do not move these stand metrics.
- Growth (signed BAIMULT): modest. A 30 percent diameter-growth slowdown moves BA 16.8 -> 14.3 and
  QMD 13.1 -> 11.4. The leverage is weak because increment over one interval is a small fraction of
  standing diameter; a multiplier cannot fully correct the over-prediction (the growth EQUATION, i.e.
  injection, is the real lever for BA).
- Mortality (signed MORTMULT): essentially no effect on TPH. Reducing mortality 40 percent did not
  raise density (TPH -9.2 vs -9.3). This is the surprise and it relocates the problem: the density
  under-prediction is NOT excess mortality, it is missing INGROWTH (recruitment). FVS in this
  projection does not add enough new small trees. Mortality multipliers cannot fix a recruitment gap.
- Max SDI (localized, raised level): the density lever that actually moved TPH in the combined run
  (-9.3 -> -6.0), by relaxing self-thinning so more stems are retained.

## The corrected fix plan

The ordered, signed-adjustment framework is sound and now implemented (the keyword field layouts are
all verified: CRNMULT/REGHMULT for allometry, BAIMULT for growth, MORTMULT for mortality, SDIMAX for
the density limit). But the NE data redirect the priorities:

1. Allometry first (crown, height-diameter): keep, for volume/biomass and as the growth base, even
   though it does not move stand BA.
2. Growth: a signed BAIMULT helps modestly, but the basal-area over-prediction needs the refitted
   growth EQUATION injected, not a scalar (weak leverage). This is the fvs2py injection blocker.
3. Ingrowth, not mortality, is the density lever. Add the FIA-derived recruitment (the ingrowth model
   is already fit in fvs-conus: output/conus/ingrowth) so the projection recruits the small trees that
   are currently missing. Mortality multipliers are near-neutral here.
4. Max SDI: localized, co-calibrated, to govern self-thinning. It is the working density knob.

So the answer to "can we fix FVS this way" is: partly, and the framework is built and tested. The
combined signed adjustments move TPH (-9.3 to -6.0) and QMD (13.1 to 11.6) in the right direction, but
basal area stays high because the growth correction needs the equation (injection), and density needs
ingrowth (recruitment), not mortality. The next concrete builds are the ingrowth recruitment wiring and
the growth-equation injection; the allometry, max-SDI, and signed-multiplier machinery is ready.

Script: `fix_ne.py`.
