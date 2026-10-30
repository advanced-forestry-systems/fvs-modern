# scripts/archive

One-off scripts that used to sit at the repository root, moved here unchanged on October 30, 2026. They are kept for the record and are not part of the build, CI or any workflow; nothing in the repository calls them.

- `patch_*.py`: one-time source patchers used while landing engine and calibration changes. They edit files at hard-coded paths on the Cardinal cluster.
- `build_*.sbatch`, `build_*.slurm`: Slurm jobs that built engine libraries or executables on Cardinal.
- `test_12_3_6.R`, `test_acd_comprehensive.R`: byte-identical copies of the scripts of the same name in `calibration/acd_divergence_2026-05-21/`, which the reports in that folder refer to.

See `docs/archive/` for the dated notes moved from the root at the same time.
