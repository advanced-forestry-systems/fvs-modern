# Archive

This section holds dated working notes, plans, handoffs, assessments and old release notes. They record how fvs-modern reached its current state and why some decisions were made. They are kept unchanged for the record and are not maintained, so statements in them can be out of date. For current behaviour, use the [User guide](../user-guide/running-projections.md) and the [Technical reference](../reference/architecture.md).

## What is here

**Notes moved from the repository root (October 30, 2026).** Until the documentation site was set up, these files sat at the repository root. They were moved here with their content unchanged:

- CONUS variant design set: [`categories_conus` JSON schema](CONUS_BLOCK_SCHEMA.md), [FVS-CONUS development notes](CONUS_DEVELOPMENT_NOTES.md), [CONUS component model specification](CONUS_MODEL_SPECIFICATION.md), [plan for incorporating fvs-conus parameters](FVS_CONUS_INTEGRATION_PLAN.md) and the [CONUS projection handoff](HANDOFF_CONUS_PROJECTION.md);
- [Disturbance-aware validation and the four-lever calibration](DISTURBANCE_AWARE_VALIDATION.md);
- two dated landing notes from August 26, 2026: [CONUS modifier layer re-land](MODIFIER_RELAND_20260826.md) and [landing wt-rescore and wt-htdbh](WT_RESCORE_HTDBH_LANDING_20260826.md).

The one-off patch scripts, test scripts and Slurm build scripts that sat at the root were moved to `scripts/archive/` in the repository; they are not part of this site.

**Earlier archive files.** Integration test reports, ACD bridge pull request notes, a session handoff, the April 2026 release readiness report and the v2026.05 release notes.

**Program notes (June 2026).** Charters, program state reports, results and runbooks under `docs/program/`.

**CONUS native adoption (June 2026).** Specifications and plans for adopting the CONUS equations natively in the engine, under `docs/conus_native_adoption/`.

**Assessments and plans.** The March 22, 2026 comprehensive assessment and verification report, the parameter extraction summary, the NVEL integration recommendations, the history rewrite plan and the original getting started page, which the [Quickstart](../getting-started/quickstart.md) and [Installation](../getting-started/installation.md) pages replace.

## Why keep them

Several of these notes are cited by code comments and by `CALIBRATION.md`, and they explain choices that the code alone does not. Keeping them in one place clears the repository root and makes clear that they are history rather than current documentation.
