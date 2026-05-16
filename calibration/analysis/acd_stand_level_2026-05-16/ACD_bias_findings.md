# ACD BA bias diagnosis — 2026-05-15

**Source:** `validation_data_acd_post.csv` from job 9610424 (30,146 ACD-relabeled NE plots)

## Headline reconciliation

The fia_benchmark headline `BA_bias_pct_calib = +19.73%` is misleading. The underlying calibrated arm is actually **under-predicting** BA, not over-predicting it.

| Metric | Calibrated | Default |
|---|---|---|
| mean obs BA | 98.8 | 98.8 |
| mean pred BA | 97.2 | 88.2 |
| mean raw residual (pred−obs) | **−1.63** sq ft/ac | **−10.64** sq ft/ac |
| median raw residual | −4.89 | (n/a) |
| pred / obs ratio | 0.984 | 0.893 |
| mean of per-record pct err | +19.73% | +9.47% |
| median of per-record pct err | −5.88% | (n/a) |

`calc_bias_pct` (line 1773 of the engine) is defined as `100 * mean((pred-obs)/obs)`. That is mean-of-ratios, which is unstable when `obs` is small. The signed mean of raw residuals and the signed median per-record pct both tell the actual story: calibrated ACD is slightly conservative (under-predicts by ~1.6 sq ft/ac, ~1.6% of mean obs).

## Decomposition by BA stratum (the +19.7% culprit)

| Stratum (sq ft/ac) | n | mean obs | mean pred | raw bias | bias_pct | rmse_pct |
|---|---|---|---|---|---|---|
| [0, 25) | 2860 | 13.7 | 22.4 | **+8.68** | **+210%** | 224% |
| [25, 50) | 3685 | 37.8 | 42.4 | +4.64 | +13.1% | 73.8% |
| [50, 75) | 4101 | 62.7 | 64.7 | +1.93 | +3.2% | 42.3% |
| [75, 100) | 4679 | 87.9 | 86.6 | −1.37 | −1.4% | 26.0% |
| [100, 150) | 9489 | 123.4 | 119.1 | −4.34 | −3.5% | 14.9% |
| [150, 200) | 4326 | 170.1 | 161.4 | −8.62 | −5.0% | 10.8% |
| [200, ∞) | 1006 | 223.9 | 209.9 | −13.93 | −6.2% | 10.0% |

The +210% in the lowest stratum (n=2,860 plots, ~9.5% of the sample) pulls the mean-of-ratios way up. Crossover from over-prediction to under-prediction happens around 75 sq ft/ac.

Coarser split: BA<50 (n=6,545) shows bias_pct = +99.1% from raw bias of +6.4, while BA≥100 (n=14,821) shows bias_pct = −4.1% from raw bias of −6.2.

## Decomposition by remeasurement interval

| interval_years | n | raw bias | bias_pct |
|---|---|---|---|
| 4 | 737 | −4.47 | −2.8% |
| 5 | 19,905 | −1.35 | +23.3% |
| 6 | 4,826 | −2.26 | +10.4% |
| 7 | 4,678 | −1.70 | +17.7% |

Raw bias is stable across intervals (−1 to −4 sq ft/ac). The bias_pct swings come from the low-BA records being concentrated in 5-year cycles.

## Recommendations

1. **Don't fix the model — fix the metric.** The +19.7% headline is a metric artifact. Either:
   - Switch the headline to a more robust dimensionless bias (e.g., mean-relative-to-mean: `100 * mean(pred-obs) / mean(obs)` = **−1.65%**), or
   - Report bias_pct only above a minimum BA threshold (say BA>50), and separately report a low-BA cell.
2. **Conservative ~1.6% under-prediction is small and probably real.** The calibrated ACD model is genuinely tracking observed BA at population level. The default arm under-predicts by 11% — calibration is doing useful work.
3. **Watch the low-BA stratum.** A +8.7 sq ft/ac over-prediction on stands averaging 13.7 sq ft/ac (n=2,860) is the only place the model is materially wrong. Likely candidates:
   - DG (ln DDS) over-predicting growth on small/sparse stands where neighbors and BAL are atypical
   - HD curve over-shooting on tall-but-spindly small stands
   - Mortality model not firing on truly stagnant low-BA stands
4. **If the journal-target headline number must be a single bias_pct, recompute it without ingrowth (post-pass setting matched the engine) vs. with ingrowth from the full engine run, and report both.** The post-pass currently sets `INGROWTH_ENABLED = FALSE` for a documented reason; full-fidelity comparison requires the engine run with `FVS_ACD_RELABEL=TRUE`.
