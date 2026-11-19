# Greg 2026 FIADB refit: Bakuzis engine check (9 November 2026)

FVSne built from this branch (Cardinal, job 16458342), engine harness v2, 36 northeastern
scenarios (4 forest types by 3 site classes by 3 densities), 100 years. All arms pass the
seeding and site gates (36 of 36). Jobs 16458389 and 16459961 to 16459965.

| Arm | Description |
|---|---|
| native | FVSne native DG, HG and mortality |
| mortnew | native growth, Greg cr_cch2 gompit mortality refit to the 2026 FIADB |
| gregDrv | Greg DG, HG, crown recession and mortality at the calibration driver climate constants |
| gME_all | as gregDrv at the statewide Maine forest median climate |
| gME_grp | as gregDrv, each forest type at its own Maine median climate |

Flagged trajectories by relation (from 2026-11-09_bakuzis-flags-by-arm_DATA.csv):

| Relation | native | mortnew | gregDrv | gME_all | gME_grp |
|---|---|---|---|---|---|
| T1 monotone dominant height | 0 | 0 | 11 | 0 | 0 |
| T2 site ordering HT / VOL | 0 / 0 | 0 / 0 | 3 / 4 | 4 / 4 | 4 / 4 |
| T3 Reineke slope | 7 | 4 | 5 | 3 | 3 |
| T6 Sukachev | 4 | 4 | 1 | 0 | 0 |
| T8 pinned BA plateau | 0 | 0 | 2 | 3 | 3 |

Climate: 1991 to 2020 normals from the Daymet daily cache of the Maine FIA covariate build
(public coordinates, lapse adjusted to field elevation). DD_0 and DD_18 are degree-days below
0 and 18 C, the ClimateNA convention Greg's fits read; his growth qmd text describes DD_0 as
above 0 C, and using that reading stalled spruce-fir and northern hardwood height growth.

Notes. The mortality refit alone changes stand trajectories little. Site ordering flags are
expected: the harness varies site index and Greg's equations have no site term. At Maine
climate, northern hardwood and oak-pine dominant heights reach 109 and 120 ft at year 100,
25 to 35 percent above native FVSne on the medium site; spruce-fir matches native (86 ft).
