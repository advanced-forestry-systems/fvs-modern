# validation/

One-time validation, gate, and diagnostic code promoted off Cardinal scratch on
September 6, 2026 so the numbers in Paper 2 and the v1 gate archive can be traced
to versioned scripts. Nothing here ships in the engine or the Python package tree.

- `paper2_gate_balfix_20260805/` the Paper 2 benchmark gate (`gate_balfix.py`) and
  its sanity scripts, previously only at `/fs/scratch/.../balfix_gate_20260805/`.
- `conus_projector/` the three-option CONUS equation projector v4, the greg_par
  variant, and `cch_module.R` (R port of GOMPCCH), previously scratch-only.
- `bakuzis_engine_verdict/` Bakuzis verdict and engine A/B scripts from the
  `wt-engine` worktree (`reconcile/three-arm-onto-main`).
- `htdbh_build/` `30b_build_htdbh_input.R` and its summary from `wt-htdbh`.
