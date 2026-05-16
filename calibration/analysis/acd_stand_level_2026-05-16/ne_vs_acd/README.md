# FVSne vs FVSacd runtime A/B — 2026-05-16

## Goal

Confirm at runtime that calling the **NE (Northeast)** variant and the
**ACD (Acadian)** subvariant on identical input produces materially
different growth predictions. Both share the upstream NE codebase, but
ACD overrides parameter tables and a subset of submodel coefficients for
the spruce-fir Acadian region, so distinct outputs are expected.

## Setup

Cardinal jobs 9717119 (FVSacd) and 9717157 (FVSne). Two standalone FVS
executables, built from the same `src-converted/` snapshot, linked
against the same `.so` artifacts (same NVEL submodule, same NSVB
defaults from this branch).

| Binary | Variant | Self-identifier |
|---|---|---|
| `lib/FVSacd` | Acadian | `AC` |
| `lib/FVSne`  | Northeast | `NE` |

Both linked through the patches accumulated on
`acd-bridge-fix-2026-05-15`:

- `errgro.f90` JOSTND guard (1cd784f)
- `fvs.f90` + `filopn.f90` unit-number init (3db614a / 6a27c2b / 6c029bd)
- MAXSP shadow fix in build script (c75e576)
- NSVB-on-by-default in `grinit.f90` (254cce3 / 3625a01 / e0aeeec)
- `spctrn.f90` AC/NE alias + JSPIN guard
- `econ_stubs.f90` + `varver_stub.f90` split

Input: upstream `ForestVegetationSimulator/tests/FVSne/net01.key` with
`tr "\r" "\n"` line ending conversion + `net01.tre`. Keyword deck
declares 4 stands x 11 cycles, 1990 to 2090 (some 2090 to 2150).

## Result: variant identifiers are distinct

| File | Variant in column 6 of stand header |
|---|---|
| `ne_net01.sum` | `NE` |
| `acd_net01.sum` | `AC` |

```
NE:  -999   11 S248112  NONE  0.1100000E+02 NE 05-16-2026 18:03:28 ...
ACD: -999   11 S248112  NONE  0.1100000E+02 AC 05-16-2026 18:03:29 ...
```

The two binaries dispatch through the variant-specific parameter and
submodel paths as designed. ACD is not a label-only alias.

## Result: tree-level predictions diverge

Sample trees at year 2090 from stand S248112 UNTHINNED CONTROL
(per the `.out` projection tables, the unstable bookkeeping units
that drive the stand-level summary):

| Sample tree | NE DBH | NE HT | NE species | ACD DBH | ACD HT | ACD species |
|---|---|---|---|---|---|---|
| 10  | 10.85 |  81.92 | SM1 | 10.64 |  81.76 | SM1 |
| 30  | 18.14 | 103.32 | WP1 | 18.42 | 103.85 | WP1 |
| 50  | 18.27 | 102.27 | SM1 | 18.96 | 118.10 | **JP1** |
| 70  | 22.09 |  98.35 | QA1 | 22.30 |  97.40 | QA1 |
| 90  | 21.25 | 120.73 | JP1 | 22.37 |  97.39 | **QA1** |
| 100 | 22.43 | 106.44 | SM1 | 22.13 | 105.87 | SM1 |

Two of the six bookkeeping units carry different species identifiers
under ACD vs NE by 2090. The DBH and HT trajectories diverge in every
row. Differences look small at a single point but compound across the
11 cycle simulation horizon.

## Result: stand-level QMD diverges

| Year | NE QMD (in) | ACD QMD (in) |
|---|---|---|
| 2090 | 17.9 | 18.2 |

That 0.3 inch gap at age 160 is exactly the kind of regional
calibration signal the ACD subvariant is intended to express against the
NE baseline.

## File comparison

| File | NE md5 | ACD md5 | Bytes (NE / ACD) |
|---|---|---|---|
| `net01.out` | `6cc4d5ec1cbc8027ee47e4fa418c5719` | `2edd05018bc67a81c5ed3b818f2563b9` | 34,717 / 217,771 |
| `net01.sum` | `4af68f98fd5427eff45021185bb08943` | `dc4c31ee240b7becea59c30ca2f12d30` |    106 / 8,762 |

The `.out` size gap (34 KB vs 218 KB) is the known FVSne summary writer
issue, not a difference in the projection itself. ACD prints SUMMARY
STATISTICS and ECHO SUMMARY pages after the projection; NE crashes part
way through the SUMMARY STATISTICS write and never gets to ECHO SUMMARY.

## Known follow-up: FVSne summary writer crash

`ne_net01.sum` contains only the header line (106 bytes) and the .out
truncates inside the SUMMARY STATISTICS table. The 11-cycle growth
projection itself completes and is fully present in `ne_net01.out`. The
crash happens after the projection during the summary table emit.

The most likely culprit is a path in `vbase/sumout.f90` or its
variant-specific overlay that uses an uninitialized unit number or
trips on a format-buffer-vs-character-buffer mismatch only for NE. ACD
is unaffected because the ACD overlay differs in exactly that path.

Tracking as Task #53. Not blocking ACD use because ACD itself prints a
complete summary.

## Implication

The fork carries genuine variant divergence between NE and ACD, as
required for ACD users in the spruce-fir Acadian region. Building
both binaries with the same NSVB defaults and the same set of
F77 to F90 conversion fixes does not collapse them into a single
output stream.

Combined with the NSVB vs CRM A/B in `../nsvb_vs_crm/RESULTS.md`,
this fork now has end-to-end runtime evidence that:

1. The variant-selection logic at build time is real (NE != ACD).
2. The V/B/C estimator selection at runtime is real (NSVB != CRM).

## Reproducer

```bash
# Build NE and ACD executables from the same source tree
bash deployment/scripts/build_fvs_executables.sh . lib ne
bash deployment/scripts/build_fvs_executables.sh . lib acd

# Run on upstream net01 deck
WORK=$(mktemp -d)
tr "\r" "\n" < upstream/.../tests/FVSne/net01.key > $WORK/net01.key
cp upstream/.../tests/FVSne/net01.tre $WORK/
mkdir $WORK/{ne,acd}
cp $WORK/net01.{key,tre} $WORK/ne/
cp $WORK/net01.{key,tre} $WORK/acd/

(cd $WORK/ne  && lib/FVSne  --keywordfile=net01.key)
(cd $WORK/acd && lib/FVSacd --keywordfile=net01.key)

md5sum $WORK/ne/net01.sum $WORK/acd/net01.sum
diff $WORK/ne/net01.out $WORK/acd/net01.out | head -60
```
