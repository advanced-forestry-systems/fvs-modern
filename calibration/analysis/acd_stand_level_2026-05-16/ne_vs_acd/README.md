# FVSne vs FVSacd runtime A/B — 2026-05-16 (revised)

## Goal

Confirm at runtime that calling the **NE (Northeast)** variant and the
**ACD (Acadian)** subvariant on identical input produces materially
different growth predictions, with summary statistics tables complete
for both binaries.

## Setup

Cardinal interactive runs on login node, 2026-05-16. Two standalone
FVS executables built from the same `src-converted/` snapshot, linked
against the same `.so` artifacts, with all branch fixes applied
including the new `sumout.f90` CASE branch deactivation.

| Binary | Variant | Self-identifier | .sum bytes | .out bytes |
|---|---|---|---|---|
| `lib/FVSacd` | Acadian | `AC` | 8,762 | 217,771 |
| `lib/FVSne`  | Northeast | `NE` | 8,762 | 216,971 |

Input: upstream `tests/FVSne/net01.key` + `net01.tre` with CR to LF
conversion. Keyword deck declares 4 stands x 11 cycles, 1990 to 2090.

## Result: full A/B coverage

Both binaries now emit the complete SUMMARY STATISTICS table for all
four stands. `.sum` files are identical in size (8,762 bytes, 61
lines) and structure but differ in content.

`.out` size delta is 800 bytes (0.4%), driven by per-cycle tree-level
text differences. Both binaries get through `STOP 10` cleanly after
processing all keyword records (post-completion IC code, not a crash).

## Result: variant identifiers are distinct

```
NE:  -999   11 S248112  NONE  0.1100000E+02 NE 05-16-2026 19:27:19 ...
ACD: -999   11 S248112  NONE  0.1100000E+02 AC 05-16-2026 18:03:29 ...
```

## Result: stand-level predictions diverge

Stand S248112 UNTHINNED CONTROL at year 2090 (age 160):

| Metric | NE | ACD | Delta |
|---|---|---|---|
| TPA | 111 | 94 | NE +17% |
| BA (ft2/ac) | 194 | 169 | NE +15% |
| SDI | 279 | 245 | NE +14% |
| QMD (in) | 17.9 | 18.2 | ACD +1.7% |
| CFV (cu ft) | 7,638 | 6,727 | NE +13% |
| BFV (bd ft) | 43,258 | 38,540 | NE +12% |
| Top HT (ft) | 106 | 106 | 0 |

ACD produces a lower-stocked, higher-QMD stand at maturity — the
regional spruce-fir signature expected from the Acadian variant
relative to the broader NE calibration. CRSF intuition holds.

## Tree-level divergence (snapshot from .out)

| Sample tree | NE DBH | NE HT | NE species | ACD DBH | ACD HT | ACD species |
|---|---|---|---|---|---|---|
| 10  | 10.85 |  81.92 | SM1 | 10.64 |  81.76 | SM1 |
| 30  | 18.14 | 103.32 | WP1 | 18.42 | 103.85 | WP1 |
| 50  | 18.27 | 102.27 | SM1 | 18.96 | 118.10 | **JP1** |
| 70  | 22.09 |  98.35 | QA1 | 22.30 |  97.40 | QA1 |
| 90  | 21.25 | 120.73 | JP1 | 22.37 |  97.39 | **QA1** |
| 100 | 22.43 | 106.44 | SM1 | 22.13 | 105.87 | SM1 |

## File md5s

| File | md5 |
|---|---|
| ne_net01.sum  | f2f09a21d89e363893dbeff5ade04b6e |
| acd_net01.sum | dc4c31ee240b7becea59c30ca2f12d30 |

## FVSne summary writer bug (now fixed)

The earlier run on this same input produced a 106-byte ne_net01.sum
(header only) because `src-converted/vbase/sumout.f90` carried an
active `CASE (CS,LS,NE,SN)` branch that the F77 to F90
conversion accidentally re-enabled from a commented-out upstream
block. That branch wrote FORMAT 12 (8-column header) plus a 27-arg
WRITE against FORMAT 20 (29-field). The 2-arg deficit, combined with
the dead-code reactivation, crashed the first per-row write inside
the DO 50 loop. ACD was unaffected because VARACD = "AC" falls
through to CASE DEFAULT.

Fix: comment out both CS/LS/NE/SN CASE branches in
`src-converted/vbase/sumout.f90`, matching upstream
`vbase/sumout.f` exactly. CS, LS, NE, SN (and AC via fallthrough)
now all use FORMAT 14 + the full 29-arg WRITE.

## Implication

The fork delivers genuine NE != ACD output divergence with both
variants producing complete reports. Combined with the NSVB vs CRM
A/B in `../nsvb_runtime_ab/`, this fork has end-to-end runtime
evidence at three resolution levels:

1. Variant dispatch (NE vs ACD).
2. V/B/C estimator (NSVB vs CRM).
3. Calibrated posteriors (NE posterior vs ACD posterior, see
   `../calibrated_ne_vs_acd/`).

## Reproducer

```bash
bash deployment/scripts/build_fvs_executables.sh . lib ne
bash deployment/scripts/build_fvs_executables.sh . lib acd

WORK=$(mktemp -d)
tr "\r" "\n" < upstream/.../tests/FVSne/net01.key > $WORK/net01.key
cp upstream/.../tests/FVSne/net01.tre $WORK/
mkdir $WORK/{ne,acd}
cp $WORK/net01.{key,tre} $WORK/ne/
cp $WORK/net01.{key,tre} $WORK/acd/

(cd $WORK/ne  && lib/FVSne  --keywordfile=net01.key)
(cd $WORK/acd && lib/FVSacd --keywordfile=net01.key)

md5sum $WORK/ne/net01.sum $WORK/acd/net01.sum
```
