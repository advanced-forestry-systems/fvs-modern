# FVS-ACD 12.3.6 comprehensive test report

2026-05-22. Full pass over the FVS-ACD work from this session, run on Cardinal
(R 4.4.0, Python 3.9, FVSne/FVSacd standalone binaries). Result: 8 of 8 PASS.

## Deployment under test

The tag-fixed 12.3.6 (`AcadianVersionTag = "AcadianV12.3.6"`) was deployed as the
working `fvsOL/inst/extdata/AcadianGY.R` in the Interface checkout; the prior
12.3.5 is kept as `AcadianGY_12.3.5_backup.r`. MORTCAL is off by default.

## Results

| # | test | result | key figure |
|---|---|---|---|
| T1 | version tag is AcadianV12.3.6 | PASS | tag corrected |
| T2 | 12.3.6 MORTCAL off == canonical 12.3.5 | PASS | BA 25.397 == 25.397 |
| T3 | 12.3.6 MORTCAL=TRUE lowers BA & TPH | PASS | BA -21.2% (25.40 to 20.02) |
| T4 | Maine FIA cross-validation (in-source v16) | PASS | off +15.4%, on +8.6% (R2 0.42 to 0.48) |
| T5 | Canadian NB MAGPlot cross-validation | PASS | off -0.04%, on -6.51% (#128 over-correction) |
| T6 | annual calibration table + bridge helper | PASS | calib.spp valid; unknown species to 1; FVS caps kept |
| T7 | make_fvs_calib patch applies + parses | PASS | patch -p1 clean, CRLF preserved, parses |
| T8 | MAGPlot to Fortran FVS (NE + ACD standalone) | PASS | both ingest+project; NE != ACD (e.g. 64.3 vs 71.7) |

T1, T2, T3, T6 run by `test_acd_comprehensive.R`; T7 by `patch`/`Rscript -e parse`;
T8 by `magplot_fvs_runner.py`. T4 and T5 are the committed in-source
cross-validation runs (`acadgy_insource_v16` and `magplot_insource_v16`).

## Headline conclusions

- 12.3.6 is a safe drop-in for 12.3.5: off by default, byte-for-byte equivalent
  behaviour when off, verified on both synthetic and real (NB) stands.
- The MORTCAL correction is FIA-specific: it halves a real Maine FIA
  over-projection but over-corrects unbiased Canadian CFI. Keep it off for
  Canada/MAGPlot, on only for FIA-like Maine (this matches the 12.3.6 header).
- The full Acadian toolchain validated together: R model, annual calibration
  table + bridge helper + patch, and the standalone Fortran FVS-NE/FVS-ACD path
  for Canadian MAGPlot tree lists.

## One polish item observed

The standalone FVS_Summary2 horizon/cycle count varied run to run for the
MAGPlot runs (a stand reported 2 vs 3 summary cycles on identical NUMCYCLE=2
input). This does not affect the NE-vs-ACD divergence conclusion but should be
pinned (fix the inv_year / NUMCYCLE / TIMEINT alignment so the summary horizon is
deterministic) before using the standalone BA values quantitatively.
