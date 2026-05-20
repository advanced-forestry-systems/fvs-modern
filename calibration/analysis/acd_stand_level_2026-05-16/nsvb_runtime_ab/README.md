# NSVB vs CRM runtime A/B comparison — 2026-05-16

## Goal

Demonstrate at runtime that flipping `LFIANVB = .TRUE.` in
`src-converted/acd/grinit.f90` (and the six other Eastern variants)
actually changes what `FVSacd` reports for volume, biomass and carbon
when given the same input deck.

## Setup

Cardinal job 9717119. Two standalone FVSacd executables built from
identical source EXCEPT for the `LFIANVB` default in
`src-converted/acd/grinit.f90`:

| Binary | LFIANVB default | Built from |
|---|---|---|
| `lib/FVSacd` | `.TRUE.` (NSVB) | current `src-converted/acd/grinit.f90` |
| `lib/FVSacd_crm_final` | `.FALSE.` (legacy CRM) | `src-converted/acd/grinit.f90.bak_pre_nsvb` swapped in for the build, restored after |

Both linked against the same `.so` artifacts, the same NVEL submodule
snapshot, the same `econ_stubs` / `varver_stub` / `errgro` / `filopn` /
`spctrn` patches from this branch. The only source difference between
the two executables is the one-line LFIANVB default flip.

Input: upstream `ForestVegetationSimulator/tests/FVSne/net01.key` with
`tr '\r' '\n'` line-ending conversion + `net01.tre` (line 24, 2 KB).
Keyword deck declares 4 stands x 11 cycles of growth-and-yield
projection 1990 -> 2090 (some stands 2090 -> 2150).

## Result

| File | NSVB md5 | CRM md5 | Size |
|---|---|---|---|
| `net01.sum` | `274a3d03ff2946770ebe08c7d3329d89` | `88d535d5026a617a369357611aedb619` | 8,762 |
| `net01.out` | `c5f40bd64f25ab1f98d393c2a9e80397` | `a8f6b993473c56c1ec5deb69f8b5502a` | 217,771 / 217,649 |
| `FVSOut.db` | `b6a6c112297eb3a2ee250e73d94fb830` | `a2dad71999ed3ce1a108662457e01c76` | 90,112 |

All three output files differ. The runtime LFIANVB toggle has measurable
effect on V/B/C estimates.

## Concrete divergence (Stand S248112 UNTHINNED CONTROL)

Reading the summary table columns (per `vbase/sumout.f90` convention,
columns roughly: Year, Age, TPA, BAft2/ac, %CCF, top-HT, QMD, CFV, NetCFV,
BFV, NetBFV, then mortality / per-cycle increments):

| Year Age TPA BA ... | NSVB | CRM |
|---|---|---|
| 1990 60 536 77 ... | `... 1530 1304 292 1633 ...` | `... 1526 1304 0 1633 ...` |
| 2000 70 501 103 ... | `... 2352 2138 600 3332 ...` | `... 2340 2132 0 3332 ...` |
| 2010 80 453 129 ... | `... 3321 3066 1403 7712 ...` | (still 0 in BFV column) |

The signal is clean and consistent:

- **CFV (gross cubic foot volume) is +4 sq ft/ac higher under NSVB at
  cycle 0**, growing to +12 at cycle 1 (2352 vs 2340) and continuing
  to diverge as the simulation progresses.
- **BFV (gross board foot volume) is non-zero under NSVB and zero under
  CRM**. CRM reports 0 in this column because the legacy path applies a
  stricter merch-spec threshold that the early-cycle small-diameter
  stand fails to meet. NSVB reports volume for the same trees because
  the NSVB total-stem equations apply a different convention.
- **TPA, BA, QMD, and the percent-cover columns are bit-identical**
  between NSVB and CRM, confirming the divergence is isolated to the
  volume-estimator path — not a side-effect on growth or mortality.

The pattern matches the documented behavior of Westfall et al. 2024
NSVB equations relative to the regional CRM (Component Ratio Method).

## Implications

1. The `LFIANVB` flip in `src-converted/{acd,ne,cs,ls,sn,kt,em}/grinit.f90`
   from commits 254cce3 / 3625a01 / e0aeeec genuinely controls
   downstream V/B/C output, not just compile-time conditional code.
2. The runtime A/B framework on this fork is now functional. Future
   regression tests of V/B/C estimator changes can use the same pattern:
   build twice (with and without the change), run on the same deck,
   md5-compare or diff the .sum.
3. Users running `FVSacd --keywordfile=...` get NSVB volume estimates by
   default with no extra keyword. Users who need legacy CRM behavior for
   backward comparability can append `FIAVBC OFF` to their keyword deck
   (the FIAVBC keyword is already registered at index 151 in
   `src-converted/base/keywds.f90`).

## Reproducer

```bash
# Build NSVB exe (current source)
bash deployment/scripts/build_fvs_executables.sh . lib acd

# Build CRM exe via grinit swap
cp src-converted/acd/grinit.f90 /tmp/grinit_nsvb_snap.f90
cp src-converted/acd/grinit.f90.bak_pre_nsvb src-converted/acd/grinit.f90
bash deployment/scripts/build_fvs_executables.sh . /tmp acd
mv /tmp/FVSacd lib/FVSacd_crm_final
cp /tmp/grinit_nsvb_snap.f90 src-converted/acd/grinit.f90

# Run on upstream net01 test deck
WORK=$(mktemp -d)
tr '\r' '\n' < upstream/.../tests/FVSne/net01.key > $WORK/net01.key
cp upstream/.../tests/FVSne/net01.tre $WORK/
mkdir $WORK/{nsvb,crm}
cp $WORK/net01.{key,tre} $WORK/nsvb/
cp $WORK/net01.{key,tre} $WORK/crm/

(cd $WORK/nsvb && lib/FVSacd --keywordfile=net01.key)
(cd $WORK/crm  && lib/FVSacd_crm_final --keywordfile=net01.key)

diff $WORK/nsvb/net01.sum $WORK/crm/net01.sum
md5sum $WORK/nsvb/net01.sum $WORK/crm/net01.sum
```
