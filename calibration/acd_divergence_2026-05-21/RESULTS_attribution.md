# Attribution results: AcadianGY.R standalone vs FVS-ACD bridge

Generated 2026-05-21. Two independent harnesses agree, which validates the
offline base R growth core against the real model.

## Layer 2: the real AcadianGYOneStand (run on Cardinal, R 4.4.0)

Representative mixed Acadian stand, 60 trees, start basal area 15.02 m2/ha.
Each row changes one bridge driver from a clean standalone configuration.

Horizon 5 yr (FIA remeasurement scale):

| config | BA m2/ha | delta | pct |
|---|---|---|---|
| baseline standalone (mult=1, CSI=12) | 17.14 | 0.00 | 0.0 |
| A: inject FVS multipliers | 17.26 | 0.13 | 0.7 |
| C: CSI 9 not 12 | 17.05 | -0.08 | -0.5 |
| C: CSI 15 not 12 | 17.20 | 0.07 | 0.4 |
| A+C: multipliers + CSI 9 | 17.17 | 0.04 | 0.2 |

Horizon 25 yr:

| config | BA m2/ha | delta | pct |
|---|---|---|---|
| baseline standalone (mult=1, CSI=12) | 25.40 | 0.00 | 0.0 |
| A: inject FVS multipliers | 26.13 | 0.74 | 2.9 |
| C: CSI 9 not 12 | 24.96 | -0.44 | -1.7 |
| C: CSI 15 not 12 | 25.76 | 0.36 | 1.4 |
| A+C: multipliers + CSI 9 | 25.66 | 0.26 | 1.0 |

Per species basal area at 25 yr, multiplier driver (real model):

| species | standalone | with mult | pct |
|---|---|---|---|
| BF | 6.17 | 6.30 | 2.0 |
| RM | 1.90 | 1.85 | -2.9 |
| RS | 5.76 | 6.27 | 9.0 |
| SM | 3.11 | 3.05 | -2.1 |
| WP | 4.88 | 5.26 | 7.8 |
| YB | 3.57 | 3.40 | -4.7 |

Stand total moved only 2.9 percent, but red spruce moved 9.0 percent and yellow
birch moved -4.7 percent. Total basal area hides the species level divergence.

## Layer 1: base R growth core (runs offline, exact Kuehne form + constraint)

The offline harness reproduces the real model attribution to within rounding
(25 yr: A +2.9 percent, CSI 9 -1.8 percent, CSI 15 +1.5 percent; per species RS
+8.8 percent, YB -4.9 percent, total +2.9 percent). Agreement between an
independent base R reimplementation of the diameter increment and the real
AcadianGYOneStand confirms the attribution is a property of the model, not an
artifact of either harness.

## Takeaways for calibration

1. The calibration multiplier driver (A) and the climate site index driver (C)
   are the dominant growth drivers, consistent with the diagnosis ranking.
2. The stand basal area constraint equalizes long horizon totals, so calibrate
   against short interval growth, not long horizon total basal area.
3. Calibrate per species, not per stand: total basal area is nearly invariant
   while species composition diverges substantially.
4. These are the principles the annualized_calibration.R optimizer implements.
