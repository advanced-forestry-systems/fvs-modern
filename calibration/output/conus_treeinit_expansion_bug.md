# CONUS campaign treeinit expansion bug (found + fixed 2026-06-03)

## Symptom

Merging the full-FIADB CONUS FVS campaign (default/calibrated arms) onto the
PERSEUS dashboard produced implausible trajectories: Maine default AGC reached
2012 Tg by 2125 from a 221 Tg anchor, and the campaign's inventory (PROJ_YEAR 0)
biomass was ~6x too light (ME mean AGB 8.8 tons/ac, median 0.7, vs ~31 expected).

## Diagnosis (not the join)

The STAND_CN join is fine: 89.7% of Maine standinit stands match the treeinit,
and 100% of campaign-output stands are present in the treeinit. The "40%
with-trees" CONUS figure is dominated by other states/variants, not a key bug.

The real defect is the **per-acre expansion factor** in the DataMart
`<ST>_FVS_TREEINIT_PLOT.csv` files. `TREE_COUNT` for overstory trees (DBH >= 5")
has a median of 1.0 and mean ~2.0, versus FIA's authoritative `TPA_UNADJ` of
~5.85. Confirmed against the raw FIA TREE table for the same Maine plots:

| source | per-plot BA (ft2/ac) |
|--------|---:|
| raw FIA TREE (TPA_UNADJ, live, DIA>=1) | 103.5 (median 102.9) |
| FVS_TREEINIT_PLOT (TREE_COUNT) | 34.8 (median 16.3) |
| **median ratio raw / treeinit** | **6.5x** |

The expansion is under-counted ~6.5x, concentrated in the **eastern** variants
(NE/CS/LS/SN). Western treeinit (WC/PN/etc.) was already correct, so its biomass
was unaffected.

## Fix

`treeinit_fix_tpa.py`: each treeinit row carries `TREE_CN` (the FIA tree CN), so
join it to the raw `<ST>_TREE.csv` and overwrite `TREE_COUNT` with the
authoritative `TPA_UNADJ` for matched live trees. Everything else (DIAMETER, HT,
CRRATIO, SPECIES, STAND_CN) is preserved, so the FVS run engine is unchanged --
just point `--treeinit-dir` at the corrected output.

Result across all 49 states: CONUS mean per-plot BA 48.8 -> 91.3 ft2/ac, all
realistic. Idempotent where expansion was already correct (WA 165->165, OR
142->142, CA 132->132); eastern states corrected (ME 35->113, NH 55->123, VT
53->115, MA 50->110). TPA match rate ~80% (older inventory panels not in the
current TREE.csv snapshot keep their original count; conservative).

Fixed files: `FIA_fresh/treeinit_fixed/`.

## Rerun

The original `out_fvs` / `out_gompit` campaign output is invalid (built on the
broken treeinit) and was not published. Relaunched against `treeinit_fixed` into
fresh dirs:

* `submit_conus_fvs_v2.slurm`    -> `out_fvs_v2`    (default + calibrated)
* `submit_conus_gompit_v2.slurm` -> `out_gompit_v2` (national gompit)

Both `sbatch --array=0-380%40`. Median batch ~45 min (NE ~72 min), so ~overnight.

## Residual

A smaller secondary loss remains: ~19% of trees-bearing stands produced near-zero
AGB despite real basal area (corr BA-vs-AGB 0.69), pointing at NSBE species
crosswalk gaps in `compute_plot_agb`. The expansion fix is the dominant (~6.5x)
correction; the NSBE coverage gap is a separate, smaller follow-up.

## Downstream

Once `out_fvs_v2` / `out_gompit_v2` complete: re-run `fvs_perseus_aggregate.py`
then `fvs_perseus_merge.py` to put corrected default/calibrated/gompit engines on
the dashboard (reserve no-harvest). Harvest/disturbance coupling and the TreeMap
spatially-explicit CONUS layer then build on the corrected per-plot biomass.
