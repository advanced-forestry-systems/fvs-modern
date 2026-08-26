# calibration/drivers

Promoted copies of `run_driver_final.py` and `arm_driver.py`, previously unversioned
scratch-only scripts at `/fs/scratch/PUOM0008/crsfaaron/mortval_20260825/`. Promoted
here on 2026-08-26 as part of the GREGDG deceleration remediation (see
`fix/gregdg-deceleration-20260826`) solely to fix their hardcoded
`FVS_GREGDG=1` default (now `0`) so future runs through this driver family do not
silently enable the unfixed native Greg DG hook. The scratch originals are untouched
(they remain audit evidence for the confirmed defect). No other behavior was changed.
