# Dominant height diagnosis, Greg arm at Maine climate (November 10, 2026)

Bakuzis v2 stands re-run with tree lists (same seed and stand order as the harness). Arms: native FVSne,
greg (Greg DG, HG, crown, gompit mortality at each forest type's Maine median climate), and gregbne (greg
with HG B0 set to the northeastern FIA 99.9th percentile height, a sensitivity arm only).

Each NE test species maps to its own row of greg_hg_coefficients.csv through FVS-NE FIAJSP (SM 318, YB 371,
AB 531, WA 541, RO 833, WP 129, RM 316, WO 802), so no maximum height is borrowed from another species.

HG B0 equals the largest measured live, unbroken-top height in the 2026 FIA build for 318 (166 ft), 833 (158),
129 (160) and 97 (124), and sits near it for others; northeastern 99.9th percentiles are 111.6, 120, 121 and
91 ft (fia_obs_height_vs_B0.csv).

Year-100 dominant height (40 largest stems per acre), medium site and medium density (BK14, BK32):
northern hardwood native 81.6, greg 109.2, gregbne 95.7 ft; oak-pine native 87.6, greg 120.3, gregbne 111.2 ft.
Across site classes native northern hardwood spans 67 to 98 ft while greg stays at 111 to 112 ft
(topht_by_arm.csv). In oak-pine, greg red maple (129 ft) and white oak (142 ft) overtop white pine (77 ft).
Numbers were recomputed independently from the raw .trl files and the FIA zip with a separate parser.
