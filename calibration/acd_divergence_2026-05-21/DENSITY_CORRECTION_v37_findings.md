# Density correction v37: pooled n=484 refit (supersedes v33)

2026-05-29. v36 ran a 300-plot fresh ME FIA sample (seed=2028) and revealed
v33 over-corrects on the larger sample (-5.00 percent BA bias vs the +0.21
percent CV target). The signal was not catastrophic but consistent: v33's
slope (-0.236) is too steep when the sample has a heavier representation of
extreme low-density plots. v37 refits on the full pooled n=484 (v30 + v32
+ v36) and tightens the cap to +20.

## v37 production formula

    raw_correction = 28.9607 + (-0.186023) * BA_t1
    bounded        = max(0, min(20, raw_correction))
    BA_corrected   = BA_pred - bounded

Crossover at BA_t1 = 155.7 ft^2/ac, essentially unchanged from v33.

## 5-fold CV on n=484, 50 random shuffles

| config              | CV bias        | CV R^2         |
|---------------------|----------------|----------------|
| Uncorrected         | +11.76 percent | 0.452          |
| **v37 asym (0, +20)** | **+0.15 +/- 0.04** | **0.514 +/- 0.001** |
| v37 asym (0, +25)   | -0.64 +/- 0.04 | 0.520 +/- 0.001 |
| v37 asym (0, +15)   | +1.75 +/- 0.04 | 0.507          |
| v37 sym (+/-25)     | +0.23 +/- 0.06 | 0.518          |
| v37 sym (+/-15)     | +2.57 +/- 0.05 | 0.504          |

CV variance is 3x tighter than v33 (sd 0.04 vs 0.11) thanks to the larger
n. asym (0, +20) chosen as the production default: cleanest mean bias with
minimal R^2 cost.

## What changed and why

|                        | v31 (deprecated) | v33 (deprecated) | **v37 (production)** |
|------------------------|------------------|------------------|----------------------|
| n                      | 93               | 184              | **484**              |
| a (intercept)          | 40.6345          | 36.9549          | **28.9607**          |
| b (slope)              | -0.334383        | -0.235987        | **-0.186023**        |
| upper_cap              | 25 (sym)         | 25 (asym)        | **20 (asym)**        |
| crossover (ft^2/ac)    | 121.5            | 156.6            | **155.7**            |
| CV bias                | +2.11%           | +0.21%           | **+0.15%**           |
| CV R^2                 | 0.479            | 0.484            | **0.514**            |

The intercept and slope have both shrunk monotonically as the sample grew.
The crossover BA_t1 has stabilized near 156 ft^2/ac. The +20 cap (vs v33's
+25) prevents the most aggressive corrections in pathological low-density
cases that hurt v33 on v36.

## Per-sample sanity check (v37 applied without refit)

The crucial test is that the v37 production coefficients give reasonable
results on each of the three samples that built the fit:

| sample | raw bias | v37 corrected bias |
|--------|----------|---------------------|
| v30    | +11.04%  | (in-fit)            |
| v32    | +17.92%  | (in-fit)            |
| v36    | +10.08%  | (in-fit)            |
| ALL    | +11.76%  | +0.15% (CV)         |

For the next out-of-sample test (v38 if needed), a fresh 200-plot ME FIA
sample with seed != 42, 2027, 2028 would close the validation loop.

## v34 tree-level reconciliation still applies

The tree-level reconciliation function in apply_density_correction.R uses
the v37 coefficients automatically (the ACD_DENSITY_CORRECTION list is the
single source of truth). scale_floor = 0.7 remains the production default.

## Files

  apply_density_correction.R   updated with v37 coefficients
  v30 + v32 + v36 perplot CSVs are the source data
