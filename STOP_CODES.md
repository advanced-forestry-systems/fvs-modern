# FVS standalone STOP code reference

When `lib/FVS<variant>` exits, the integer exit code conveys the
worst-severity condition encountered during the run. The codes are
emitted by `src-converted/base/main.f90` based on `fvsGetICCode`:

| Exit code | ICCODE | Meaning | Typical cause |
| --- | --- | --- | --- |
| 0   | 0 | Clean run, no notes logged | Pristine deck, all keywords accepted, all species in variant set |
| 10  | 1 | Run completed with informational notes | Species substitutions, regen tree replacement, harvest reselection |
| 20  | 2 | Run completed with non-fatal errors | Keyword parameter missing/incorrect, forest code out of range, ECODIVISION missing |
| 30  | 3 | Serious error encountered (run may be partial) | Major model failure but FVS recovered |
| 40  | 4 | Fatal error during the run | Allocation failure, division by zero, similar |
| 50  | 5 | User-requested abort | STOPRUN keyword or similar |

`errgro.f90` is the single source of these escalations. Each call
into errgro names a severity, and `IF (ICCODE .LT. severity)
ICCODE = severity` ensures the worst condition wins.

## What ICCODE values appear in practice on a clean fork

| Scenario | Typical exit |
| --- | --- |
| Variant's own test deck, all species recognized | 0 |
| Variant's own test deck with species substitutions | 10 |
| Foreign deck (NE deck run on CS, LS, SN, KT, EM) | 20 |
| Genuine variant bug or build problem | 30+ |

## Implication for the integration test rubric

The original v1 rubric treated rc=20 as FAIL, which mislabeled five
variants (CS, LS, SN, KT, EM) that were actually working but had
been fed NE's net01.key. The relevant FVS04 ERROR text in the .out
file makes this explicit:

```
********   FVS04 ERROR:  A REQUIRED PARAMETER IS MISSING OR A
           PARAMETER IS INCORRECT; KEYWORD IGNORED.
********   FVS03 WARNING: FOREST CODE INDICATES THE GEOGRAPHIC
           LOCATION IS OUTSIDE THE RANGE OF THE MODEL.
           DEFAULT CODE IS USED.
********   FVS43 WARNING: FIAVBC REQUESTED BUT ECODIVISION CODE
           NOT RECOGNIZED OR MISSING. USING "0000" BY DEFAULT
```

These warnings/errors all stem from the deck specifying NE forest
location code 922, NE-specific species codes (JP, OT, YB), and an
NE-shaped DESIGN record — none of which are valid in CS/LS/SN/KT/EM.
The variants correctly flag the input as foreign and fall back to
defaults, then complete the projection.

The v2 rubric (in `integration_test_v2.sh`) classifies rc=20 as
WARN rather than FAIL when the .sum file is complete, and runs each
variant against its own upstream test deck where one is available
to drive exit codes down to 0 or 10.

## Quick reference: which severity comes from which errgro call

Search `src-converted/base/errgro.f90` for the integer that follows
`IF (ICCODE .LT. <N>) ICCODE = <N>` to see what each call site
escalates to. The 10 escalation sites in the current source map to:

- Severity 1: most species/substitution/regen notes
- Severity 2: keyword parser errors (FVS04), location errors (FVS03)
- Severity 3+: rare in practice on the Eastern variants

For a fresh diagnostic on any specific exit code, grep the variant's
.out file for `FVS\d\d ERROR` or `FVS\d\d WARNING`.
