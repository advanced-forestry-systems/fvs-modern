# Archive: conus-sf-integration record artifacts

## Provenance

These files were extracted from the head of the now-closed pull request
**#70 "CONUS species-free integration scaffold (both legs, blend, uncertainty)"**
(branch `feat/conus-sf-integration`). That PR was development scaffolding and
was closed without merge because it carried no engine code destined for the
`reconcile/three-arm-onto-main` line. The items below have lasting reference
value and are preserved here so the record survives independently of the
closed branch.

**These are archived reference artifacts, not runtime code.** Nothing in this
directory is on any engine build or execution path. Do not wire it into the
engine, configs, or builds.

## Contents

- `sdimax/` - the max-SDI (maximum stand density index) technical reports and
  decks: the CONUS max-SDI technical report, FVS team briefings and updates
  (Markdown, DOCX, PPTX), status/handoff notes, and the supporting analysis
  scripts, figures, and summary CSVs that accompanied those reports.
- `stan/` - the four species-free diameter-growth Stan model variants from the
  Kuehne (2022) exploration:
  - `dg_kuehne2022_speciesfree_v1_quad.stan`
  - `dg_kuehne2022_speciesfree_v2_l1site.stan`
  - `dg_kuehne2022_speciesfree_v3_traitsite.stan`
  - `dg_kuehne2022_speciesfree_v4_full.stan`

## What was deliberately not archived

Engine code, configuration, build files, and the broader calibration pipeline
outputs (e.g. `calibration/output/`, `calibration/figshare/`, `calibration/R/`,
`calibration/osc/`) from PR #70 were intentionally left behind. The full history
of PR #70 remains available on the `feat/conus-sf-integration` branch in the
remote.
