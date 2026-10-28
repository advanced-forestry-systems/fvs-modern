# fia_stand_generator.py: inverse plot size double-expansion (found October 27, 2026)

`build_stand_init()` writes `inv_plot_size: 6.0` while `build_tree_init()` writes `tree_count = TPA_UNADJ`, which is already trees per acre. Through the FVS_StandInit path (the route NE and ACD take in `bakuzis_uncertainty_comparison.py`) FVS multiplies each record by the inverse plot size, so every tree enters at six times its expansion (subplot trees at 36.2 TPA, microplot trees at 451 TPA). A six-cell Northeast smoke run on October 27 with the main-branch library seeded at 450 to 28,788 TPA and 4 to 655 ft² ac⁻¹ against labels of 60 to 180.

Patch applied in the Cardinal worktree `/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main` (not yet committed):

    -        "inv_plot_size": 6.0,
    +        "inv_plot_size": 1.0,   # tree_count is TPA_UNADJ (already per acre)

After the patch the same six cells seed at 0 to 4,798 TPA and 0 to 109 ft² ac⁻¹, so the six-fold inflation is gone but the stands still do not reproduce the condition BALIVE they were selected on (ratios 0.02 to 1.5 of the cell label; one replicate has no trees). The remaining loss is upstream of the keyfile, in which trees `load_fia_trees()` attaches to a sampled condition. The regeneration of the Paper 1 ensemble has to wait on that; the year-0 gate in `bakuzis_stress_v2.R` (T0) and the seed check in `engine_bakuzis_v2.py` will catch it when it is fixed. Outputs: `ens_smoke/` (before patch) and `ens_smoke22/` (after) under `/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/`, with `seedcheck.log` and `seedcheck2.log`.
