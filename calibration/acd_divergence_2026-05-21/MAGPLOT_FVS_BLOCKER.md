# Fortran FVS-NE/ACD on Canadian MAGPlot: blocker diagnosis

2026-05-21. Goal: run the Fortran FVS engine (FVS-NE and FVS-ACD, default and
calibrated) on Canadian (New Brunswick) MAGPlot tree lists. The user reported
this as "blocked by the fvs2py inventory issue." Here is the full diagnosis.

## Context: the R path already works

The R AcadianGY path on MAGPlot succeeds. SLURM job 10129046 (2026-05-20)
projected 262 NB plots ten years and validated against observed remeasurement:
BA bias about 0 percent, R squared 0.88, QMD +9 percent, TPH -7 percent. So the
Acadian model itself ingests and runs MAGPlot fine. The blocker is specific to
the Fortran engine path through fvs2py.

## The blocker is a chain of three issues

### 1. fvs2py does not import on the cluster default Python (root cause)

The cluster default `python3` is 3.9.21. fvs2py uses `enum.StrEnum` (added in
3.11) and `typing.ParamSpec` (added in 3.10), so the import fails before any
inventory can be read:

    ImportError: cannot import name 'StrEnum' from 'enum'
    ImportError: cannot import name 'ParamSpec' from 'typing'

This is the surface "fvs2py inventory issue": fvs2py never loads, so nothing can
be ingested. FIX (confirmed): run under `module load python/3.12`. Under 3.12,
`from fvs2py import FVS` succeeds and `FVSne.so` loads. A version guard was added
to `magplot_fvs_runner.py` so this fails fast with guidance instead of a cryptic
import error. (A StrEnum-only shim was tried and reverted, because ParamSpec also
fails on 3.9; backporting fvs2py to 3.9 is not worth it.)

### 2. The DATABASE keyfile run produces no output (the deeper issue)

Under python/3.12, fvs2py loads FVSne.so and reads the keyword file, but a
DATABASE/STANDSQL keyfile run (the proven bakuzis/silc pattern) does not execute
the projection:

- `restart_code` stays 0 across repeated `run()` calls (never reaches 100, the
  stand-complete state).
- The output SQLite DB contains only the two input tables `fvs_standinit` and
  `fvs_treeinit`; no `FVS_Summary`, `FVS_TreeList`, or `FVS_Cases` are written.
- `fvs.summary` is empty.

This is NOT MAGPlot-specific. A clean synthetic 3-tree stand with valid FIA SPCD
(12, 97, 316) through the same `run_one` DATABASE path produces the same empty
result. So the SQLite DATABASE inventory ingestion / run loop is non-functional
in the current FVSne.so + fvs2py stack, for any input. This is the real engine
blocker and needs engine-level work (see next steps).

### 3. Multiple runs in one process hard-crash on the keyword unit

Running more than one stand in a single Python process aborts with:

    At line 47 of file ./src-converted/bin/../base/keyrdr.f90 (unit = 15, file = 'fort.15')
    Fortran runtime error: Sequential READ or WRITE not allowed after EOF marker

FVS keeps Fortran global state and unit 15 (the keyword scratch file) open across
runs; the dynamic loader caches the `.so`, so the second run reads past the prior
run's EOF marker. A single run per fresh process does not crash. So any batch
runner must isolate each stand in its own process (subprocess, or multiprocessing
with maxtasksperchild=1), or keyrdr.f90 must REWIND/handle EOF on re-entry.

## What was built and is ready

`magplot_fvs_runner.py` converts MAGPlot NB tree lists into the FVS inventory
schema and runs them through the engine, reusing the proven `run_one` path:

- species crosswalk: MAGPlot botanical genus.species (e.g. ABIE.BAL, PICE.RUB)
  to FIA SPCD, covering the NB species (balsam fir, black/red/white spruce, red
  maple, paper/yellow birch, aspen, cedar, beech, pines, etc.); generics and
  unknowns map to FVS "other" buckets.
- units: dbh cm to inches, height m to feet, `stem_ha` to trees per acre.
- initial measurement (meas_num 0), live trees only.

The converter is correct; it will produce results as soon as issue 2 is fixed.
It currently surfaces issue 2 (empty summary) and, when looping, issue 3.

## Recommended next steps, in priority

1. Confirm the DBS (SQLite DATABASE) extension status in FVSne.so. The clean
   synthetic stand producing no output indicates the DATABASE input or the run
   loop is broken in this build. Check whether the DBS extension is compiled in
   and whether `fvs.run()` needs to be driven to completion differently (the
   docstring says run() returns per stand and sets restart_code 100 when done; it
   never advanced here). This likely needs a focused look at fvs2py's run loop
   and/or the FVS DBS/keyrdr Fortran with a rebuild.
2. If DBS is the broken piece, fall back to FVS plain-text inventory: write a
   classic keyword file plus a fixed-format `.tre` tree file (TREEDATA) and run
   that, bypassing SQLite entirely.
3. For batch, run one stand per fresh process (subprocess or multiprocessing
   maxtasksperchild=1), or fix keyrdr.f90 to handle the unit-15 reuse.

## Reproduction (on Cardinal)

    module load python/3.12
    export FVS_PROJECT_ROOT=/users/PUOM0008/crsfaaron/fvs-modern
    export FVS_LIB_DIR=/users/PUOM0008/crsfaaron/fvs-modern/lib
    python3 magplot_fvs_runner.py --nstands 1 --min-trees 15
    # single stand: loads + runs, empty summary (issue 2)
    # --nstands > 1: keyrdr.f90 EOF crash (issue 3)
