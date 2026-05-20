# fvs-modern integration test v2 (variant-specific decks)

Run: SLURM job $SLURM_JOB_ID
Date: 2026-05-17 11:12 UTC
Branch: acd-bridge-fix-2026-05-15
Commit: f0dd4ba Autopilot round 9: line-ending fix + deck picker fix = 11/12 clean runs

## Test scope: 12 variants

**Eastern** (7): ACD, NE, CS, LS, SN, KT, EM
**Western sample** (5): WC (Westside Cascades), OP (Olympic Peninsula), CA (Inland California), BM (Blue Mountains), CR (Central Rockies)

**Test deck:** each variant uses its own upstream tests/FVS<variant>/ deck.
ACD falls back to NE's net01 because no FVSacd test dir exists upstream.

## Refined rubric

- **PASS**: rc in {0, 10} and .sum >= 1KB
- **WARN**: rc=20 (input-variant mismatch or non-fatal error) but .sum complete
- **FAIL**: no .sum, .sum too small, marker missing, or rc >= 30

## Results

| Phase | Variant | Status | Detail |
| --- | --- | --- | --- |
| build | acd | PASS | 8156064B |
| build | ne | PASS | 8169728B |
| build | cs | PASS | 8129504B |
| build | ls | PASS | 8116920B |
| build | sn | PASS | 8157192B |
| build | kt | PASS | 9143240B |
| build | em | PASS | 9270872B |
| build | wc | PASS | 9329560B |
| build | op | PASS | 8949160B |
| build | ca | PASS | 8653744B |
| build | bm | PASS | 9278736B |
| build | cr | PASS | 9313736B |
| run | acd | PASS | rc=10 sum=8762B out=217766B  deck=net01.key |
| run | ne | PASS | rc=10 sum=8762B out=216966B  deck=net01.key |
| run | cs | PASS | rc=10 sum=8762B out=211385B  deck=cst01.key |
| run | ls | PASS | rc=10 sum=8762B out=205447B  deck=lst01.key |
| run | sn | PASS | rc=10 sum=8615B out=212770B  deck=snt01.key |
| run | kt | PASS | rc=10 sum=8762B out=185777B  deck=ktt01.key |
| run | em | PASS | rc=10 sum=8762B out=188821B  deck=emt01.key |
| run | wc | PASS | rc=10 sum=8762B out=193535B  deck=wct01.key |
| run | op | PASS | rc=10 sum=9791B out=231583B  deck=opt01.key |
| run | ca | PASS | rc=10 sum=8762B out=198491B  deck=cat01.key |
| run | bm | PASS | rc=10 sum=8762B out=186866B  deck=bmt01.key |
| run | cr | PASS | rc=10 sum=8762B out=196359B  deck=crt01.key |
| marker | acd | PASS | found AC |
| marker | ne | PASS | found NE |
| marker | cs | PASS | found CS |
| marker | ls | PASS | found LS |
| marker | sn | PASS | found SN |
| marker | kt | PASS | found KT |
| marker | em | PASS | found EM |
| marker | wc | PASS | found WC |
| marker | op | PASS | found OP |
| marker | ca | PASS | found CA |
| marker | bm | PASS | found BM |
| marker | cr | PASS | found CR |
| distinct | all | PASS | 12 variants, 12 unique md5s |
| smoke | postpass | PASS | stratified post-pass OK |

## Pass summary

- Total: 38
- PASS: 38
- WARN: 0
- FAIL: 0
- SKIP: 0

**Overall: PASS** (38 PASS / 0 WARN / 0 SKIP)
