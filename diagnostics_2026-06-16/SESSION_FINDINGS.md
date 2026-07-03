# fvs-modern diagnostics, 2026-06-16

A verification pass over the three fvs-modern refinement claims (calibration, density, uncertainty),
the SDIMAX over-thinning fix, and the engine-injection blocker. Scripts and the stress-test figure
are in this directory; the narrative handoffs live in the project notes.

## 1. SDIMAX over-thinning bug: FIXED and shipped

The calibrated config over-thinned density 25 to 35 percent because the per-species SDIMAX keyword was
emitted with its fixed-column fields misaligned (`"SDIMAX" + 10 spaces + species + value`), so FVS read
the species index as the max-SDI value and set MAX SDI ~= 1 for all species. Fixed in
`config/config_loader.py`: SDIMAX emission disabled by default (`_emit_sdimax`, commit a8e03cd) and the
column alignment corrected for when it is re-enabled (commit 6c0a2c7). Verified at the engine level and
across all variants (`config_stress.py`: SDIMAX now 0 everywhere).

## 2. Calibration delivery: the multipliers do not move the engine

GMUG-style stress test (`gmug_stresstest_2026-06-16.png` / `.csv`), default vs calibrated through the
FVS engine on FIA remeasurement, 7 variants: with the over-thinning removed, calibrated is at PARITY
with default (basal-area bias within about a point everywhere). The earlier GMUG "BA win" was the
over-thinning artifact, not the multipliers. Cause: every per-species multiplier set (diameter growth,
crown, height-diameter, mortality median) is centered on 1.0; they redistribute among species and net
to zero at the stand level. Crown and height-diameter multipliers are not even emitted to the engine
(no CRNMULT/REGHMULT in generate_keywords, not in set_species_attr). Validation note: the
`19_fia_benchmark_engine.R` table that shows calibrated beating default evaluates the refitted
EQUATIONS in R, not the FVS engine; the equations are good, the engine never receives them.

## 3. Uncertainty: incorporated but NOT calibrated (fix verified)

Interval coverage (PICP) against 433,291 observed conditions: the nominal-95 percent intervals cover
about 13 percent (BA 13.4, VOL 12.8, HT 11.4). They carry parameter uncertainty only and omit the
residual variance. Fix verified: a predictive interval (point +/- residual SD) restores coverage to
93 percent overall, 92 to 95 percent per variant. Production change is bounded: draw a residual per
tree from each component sigma in `21_uncertainty_propagation.R`, propagate, take percentiles.

## 4. Engine injection (refit equations into the engine): blocked at the maintainer level

To actually beat default the refitted equations must run inside the engine, not as multipliers. The
fvs2py shared-library path does not load the stand: `FVS(lib="FVSne.so").load_keyfile().run()` produces
zero trees and an empty output while the subprocess executable loads the same stand fully. Root cause
located: the executable's main reads the I/O filenames and the FVS routine initializes via the program
startup, whereas calling `fvs_` via ctypes runs the input phase not at all (empty .out, rc=2),
regardless of cmdline `--keywordfile` or a stdin-redirect that replicates the subprocess. Both fix
routes exist in the library (`filopn_`, `fvsaddtrees_`); wiring them is a focused ctypes task that
needs the fvs2py maintainer's knowledge of the in-process init protocol. Reproductions: `clean_repro.py`,
`fullrun.py`, `ctypes_fix.py`.

## Bottom line

The refits and the uncertainty draws are real and sound in their own R framework. The three gaps are
all in the delivery and validation layer: get the equations into the engine (injection), report
predictive (not parameter) intervals, and keep the SDIMAX emission corrected. The SDIMAX fix is done;
the other two are bounded, well-characterized tasks.

## Files here

- `gmug_stresstest_2026-06-16.png`, `gmug_stresstest_table.csv`: corrected GMUG-style comparison.
- `overthin_diag.py`, `combined_bench.py`, `growth_lever.py`, `config_stress.py`: engine benchmarks.
- `var_scale_diag.py`: localized max-SDI density sweep.
- `picp.py`, `predint.py`: uncertainty coverage and the predictive-interval fix.
- `clean_repro.py`, `fullrun.py`, `ctypes_fix.py`: engine-injection reproductions.
