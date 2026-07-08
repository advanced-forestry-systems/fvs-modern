# Master handoff and work order (2026-06-16)

The single document to start a new session from. Supersedes earlier handoffs. Part 1 is how to
reconnect, Part 2 is honest status, Part 3 is the prioritized work order with acceptance criteria,
Part 4 is exactly which files to point to.

---

## Part 1. How to reconnect (environment)

- **Cardinal SSH:** key persisted at `~/Documents/Claude/.ssh-cardinal` (mounted as `.ssh-cardinal`;
  discover via find, copy to `~/.ssh/id_ed25519`, chmod 600; re-do each bash call, no carryover). User
  `crsfaaron`, host `cardinal.osc.edu`.
- **GitHub (holoros):** token at `CRSF-Cowork/_context/.gh-holoros/token`; pipe into
  `gh auth login --with-token`; never print it. Branch `feat/conus-sf-integration` (PR #70, open).
- **Python for fvs2py / the FVS class API (injection):** MUST use Python 3.10+. The Cardinal default
  venv is 3.9 and fvs2py will not import there (union hints + ParamSpec). Use `module load python/3.12`
  and the venv at `~/fvs312` (pandas+numpy installed). This is the single most important environment
  fact discovered this session.
- **Python for the projection harness (benchmarks):** the harness `run_fvs_projection` prefers the FVS
  subprocess executable (`lib/FVS<variant>`), which works under 3.9; that is why the benchmarks ran
  despite fvs2py being broken under 3.9.
- **R on Cardinal:** `module load gcc/12.3.0 R/4.4.0` (4.4.0 has data.table + mgcv; 4.5.2 does not).
- **Operational hard stops:** no production config writes or coefficient changes without Aaron's
  sign-off; do not cancel Aaron's running Cardinal jobs; use login-node nohup, not sbatch, while his
  association limit is hit. Note: heavy repeated git operations on the Cardinal `~/fvs-modern` checkout
  churned the working tree this session and deleted untracked scripts; keep new scripts in `~/` (home),
  not inside the repo, and avoid `git pull --rebase` on Cardinal while scripts are uncommitted.

---

## Part 2. Honest status by workstream

- **Maximum SDI (strong, done).** Species-weighted maximum is biased high and uninformative; a
  localized FIA-derived maximum predicts observed self-thinning better (a relative gain on a weak
  signal). Across 20 variants, a level-calibrated localized maximum is at or below native; the optimal
  level is variant-specific (0.6 to 2.0, magnitudes seed-sensitive, read qualitatively). The level
  must be co-calibrated with mortality. Module + R companion committed.
- **fvs-modern calibration (partial, with a real bug).** Mortality, crown, and height-diameter
  multipliers are populated for all 25 variants; diameter growth only ~7 variants; height growth none.
  The calibrated config improves basal area but OVER-THINS density (NE: BA bias +18 to -4%, but TPH +9
  to -32%). Mortality multipliers are mostly benign (NE median 1.0; 4 of 108 species exceed 3x), so the
  over-thinning is dominated by the density-limit emission (SDIMAX) and/or the calibrated self-thinning
  slope, not broad mortality. This is the blocker to "fully functioning calibrated FVS."
- **Species-free CONUS equations (validated).** All six tree-level bundles banked. Held-out NE
  validation: species-free beats per-species and equals hybrid for HCB (0.140 vs 0.147), HG (0.315 vs
  0.318), HT-DBH (3.55 vs 3.77), PICP 0.93 to 0.95. DG validated vs Greg; survival via mortality AUC.
- **Engine injection (env unblocked; one bug).** fvs2py runs under 3.12; FVS class loads the library
  and `run(stop_point_code=5)` executes. But the fvs2py shared-library path loads zero trees from the
  database, while the same stand loads fully via the subprocess executable (639 TPA, BA 48, SDI 118).
  So data + keyfile are correct; the fvs2py `.so` path does not initialize the DB read. Next: set the
  FVS input/output unit assignments the subprocess provides via stdin, or load trees via `fvsAddTrees`.
- **Deliverables (done).** FVS-team single report + update deck; Greg WMENS slides (4, with held-out
  validation); technical report; stress test; implementation-status and focus-status memos.
- **Zenodo.** v2026.05.4 archived (repo zip). Its description overstates calibration ("growth and
  mortality across all 25 variants"); should be softened to match what is populated.

---

## Part 3. Work order (prioritized, with acceptance criteria)

### WO-1 (HIGH): Make the calibrated config improve density, not over-thin it
- **Isolate the cause.** Run NE (and 2-3 other variants), default vs calibrated vs calibrated with
  SDIMAX neutralized (inject `SDIMAX i 9999` for all species via `extra_keywords`) vs calibrated with
  the native self-thinning slope restored. Compare TPH and BA bias vs observed. The arm that returns
  TPH bias to near the default identifies the cause (expected: the SDIMAX emission, since mortality
  multipliers are benign).
- **Fix.** In `config/config_loader.py` `generate_keywords()` / `_find_sdi_param()`, correct the
  SDIMAX emission (do not emit a too-low value; emit the native or a localized level-calibrated value),
  and/or revert the calibrated self-thinning slope if it is the driver.
- **Acceptance:** on the clean calibrated-vs-default-vs-observed benchmark, the calibrated config has
  TPH bias within roughly +/-10% AND keeps the basal-area improvement, across at least NE, CR, SN, LS.
- **Note:** this is a production-behavior change; get Aaron's sign-off before committing the config
  loader change.

### WO-2 (HIGH): Finish the species-free engine injection
- Under `~/fvs312` (python 3.12), get the fvs2py `.so` path to load a stand. Two candidate fixes:
  (a) provide the FVS input/output unit files the subprocess supplies via stdin (keyfile, tree file,
  output, treelist) which `load_keyfile` does not set; (b) load the initial tree list via the FVS API
  `fvsAddTrees` instead of the DATABASE keyword. Confirm `dims["ntrees"]>0` after `run(stop_point_code=
  5, stop_point_year=-1)`.
- Then shadow mode: read `dbh/ht/cr/species/bal` via `fvsTreeAttr` GET (verify the gfortran ctypes
  signature: explicit args then hidden string lengths last), log species-free vs engine increment, no
  SET. Then single-component injection (diameter growth), project, compare to observed.
- **Acceptance:** a one-stand run that stops at point 5 with ntrees>0 and prints the tree list; then a
  shadow log of species-free vs engine DG per tree. Consider looping in David Diaz (MicroFVS) on the
  fvs2py `.so` initialization.

### WO-3 (MEDIUM, compute-gated): Fill diameter- and height-growth multipliers
- The per-variant component fits feed `calibration/R/06_posterior_to_json.R`. DG is populated for ~7
  variants, HG for none. Run the missing per-variant growth fits when a compute slot frees, serialize
  into the configs. **Acceptance:** dds and htg multipliers non-unity for all 25 variants, with
  provenance.

### WO-4 (MEDIUM): Re-verify the all-variant max-SDI magnitudes
- Re-run the level sweep with multiple seeds per variant (the harness hard-coded a seed; PN varied 49
  to 68% across seeds). Report seed-averaged per-variant levels. **Acceptance:** per-variant optimal
  level with a spread/CI, not a single noisy number.

### WO-5 (LOW, compute-gated): Ingrowth composition fit + uncertainty wiring
- Finalize `36_fit_ingrowth_species_composition.R`; wire posterior draws to projection intervals.

### WO-6 (LOW, non-compute): Communications
- Send Greg the WMENS slides. Soften the Zenodo v2026.05.x description to match actual calibration
  coverage. The FVS-team report + deck are ready when Aaron wants to send them.

---

## Part 4. Files to point to (start here)

**Local memos (in `~/Documents/Claude/fvs-conus/`):**
- `20260616_MASTER_HANDOFF_and_WORK_ORDER.md` (this file).
- `20260615_fvs_conus_focus_status.md` (CONUS species-free + injection detail).
- `20260615_fvs_modern_implementation_status.md` (the calibration-coverage audit).
- `20260615_STRESS_TEST_everything.md` (every claim at true strength).
- `20260615_GMUG_April_what_still_holds.md` (what holds / overstated / over-thinning).
- `20260614_USFS_maxSDI_technical_report.md` (the max-SDI technical report).
- Deliverables: `FVS_team_REPORT_3step_plan.docx`, `FVS_team_UPDATE_3step_plan.pptx`,
  `Weiskittel_slides_for_GregWMENS.pptx`.

**Repo (`github.com/holoros/fvs-modern`, branch `feat/conus-sf-integration`):**
- `config/calibrated/<variant>.json` (the calibration_multipliers + SDIMAX the engine uses).
- `config/config_loader.py` (`generate_keywords`, `_find_sdi_param`) -- the WO-1 fix site.
- `calibration/sdimax/` (localized_sdimax.py/.R, var_scale_diag.py, allvar_calibration.csv, figures,
  the technical report copy, and all the committed memos).
- `calibration/injection/sf_injector.py`, `ENGINE_WIRING_DESIGN.md` (the injection prototype + design).
- `deployment/fvs2py/` (the FVS class; `_base.py` has the 3.9 lazy-annotations fix).

**Cardinal (`~` and `~/fvs-conus`, `~/fvs-modern`):**
- `~/fvs312` (python 3.12 venv for fvs2py/injection).
- `~/fvs-modern/calibration/python/`: `perseus_100yr_projection.py` (run_fvs_projection, the harness),
  `calib_vs_default.py` (calibrated-vs-default benchmark), `var_scale_diag.py` (max-SDI sweep),
  `shadow2.py`/`shadow3.py` (injection diagnostics). Re-create `overthin_diag.py` from WO-1 (it was
  removed by git churn).
- `~/fvs-conus/dev_sf_integration/benchmark_sf_vs_legA.R` (species-free held-out validation).
- `~/fvs-conus/output/conus/sf_integration/*_sf_*.csv` (the banked species-free bundles).
- Data: `~/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds` (8.2M rows; SDImax_brms etc.),
  `brms_SDImax.csv`, `VAR_SDIMAX.csv`; FIA per-state CSVs at `/fs/scratch/PUOM0008/crsfaaron/FIA`;
  TreeMap SDImax surface at Zenodo 10.5281/zenodo.19509367.

**First action in the new session:** WO-1 (re-run the over-thinning isolation cleanly and implement the
SDIMAX/slope fix). It is the single change that turns the calibrated config from "better basal area,
worse density" into "fully functioning," and it is the most consequential open item.
