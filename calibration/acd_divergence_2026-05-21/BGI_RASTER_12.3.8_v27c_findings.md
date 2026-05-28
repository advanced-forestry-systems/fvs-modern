# BGI raster as plot-level diameter multiplier under AcadianGY 12.3.8

2026-05-28. Extracts BGI (Biological Growth Index) from
`ME_BGI_V1.tif` at each of the 100 ME FIA sample plot's LAT/LON (98.3 percent
match rate over the full 21,638-plot inventory), then layers BGI into 12.3.8
production posture two ways. Compares against v27b's three site-index variants.

## Result

| config            | input mean   | BA bias % | R^2     | TPA   | QMD cm |
|-------------------|--------------|-----------|---------|-------|--------|
| csi_baseline      | 14.31 m      | +11.05    | 0.4232  | 1043  | 4.92   |
| sicond_replace    | 16.24 m      | +11.39    | 0.4190  | 1045  | 4.91   |
| fvssi_replace     | 15.38 m      | +11.27    | 0.4190  | 1044  | 4.92   |
| **bgi_as_dmult**  | 4133.9       | **+0.95** | **-10.44** | **822** | 4.82 |
| bgi_recenter      | 1.005        | +11.09    | 0.4168  | 1042  | 4.92  |

Observed: BA 94.72 ft^2/ac, TPA 1029, QMD 4.97 in.

## Interpretation

**bgi_as_dmult is pathological**, not a closure. The raw BGI raster values
(mean 4134) used directly as `dDBH.mult` explode diameter growth, trees hit
the species `max.dbh` cap of 200 cm, mortality cleans them out, TPA crashes
200 trees below observed. BA bias drops to +0.95 percent by accident, but
R^2 collapses to -10.44 (worse than predicting the mean) and the dynamics
are broken. This config is a sanity check for the magnitude of dDBH.mult's
leverage on the model, not a real candidate.

**bgi_recenter is the legitimate test.** BGI normalized to mean=1.0 preserves
mean growth and lets BGI carry plot-level productivity variance. Result:
+11.09 percent BA bias, R^2 0.4168 - statistically indistinguishable from the
ClimateSI baseline (+11.05, 0.4232). The 0.006 drop in R^2 is noise.

**BGI is not a lever on this residual** when layered inside 12.3.8 production
posture as a plot-level dDBH multiplier. The ME BGI raster encodes real
plot-level productivity variation, but the variation is not informative about
the model's overshoot of stand BA on Maine FIA conditions.

## What this rules out

After v25 (CSI scaling, partial), v26 (mort.mult, null), v27b (site-index
metric, no improvement), and v27c (BGI from raster, no improvement), the
external productivity inputs to AcadianGY have been exhausted as candidates
for closing the +11 percent FIA BA residual. The residual is not from picking
the wrong proxy for site quality or from underweighting climate; it is
internal to the model.

## What remains

The +11 percent BA residual under 12.3.8 production posture sits inside one
of:

  1. The intrinsic Kuehne et al. 2020 dDBH base-rate coefficients
     (refit against ME FIA, paper-sized).
  2. The Glover/Hool mortality functional form (Weibull or piecewise
     replacement, paper-sized).
  3. A genuine Acadian-variant structural ceiling on Maine FIA conditions
     (accept and document).

None of these is accessible without coefficient-level refitting work that
sits outside the autopilot scope. The sensitivity sweep work is complete.

## Files

`cardinal_acadgy_bgi_v27c.R`, `acadgy_bgi_v27c_results.csv`,
`me_bgi_by_pltcn.csv` (the gdallocationinfo extract on all 21,638 ME plots),
and the Python helper `extract_bgi.py`. Cardinal run as SLURM job 10876272,
about 23 minutes on c0317.
