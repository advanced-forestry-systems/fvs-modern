#!/usr/bin/env python3
"""Smoke / unit tests for the conus_greg Python consumer.

Run:  python3 sf_integration_dev/test_conus_greg_projector.py
Checks: loader wiring, kernel sanity (positive, monotone, bounded), species
distinction, fallback path, EMT/TD default flagging, and a hard-coded
regression value that pins the greg_est_dg form against the R reference.
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from config.config_loader import FvsConfigLoader
from sf_integration_dev.conus_greg_projector import GregProjector


def approx(a, b, tol=1e-6):
    return abs(a - b) <= tol


def main():
    L = FvsConfigLoader("ne", version="conus_greg")
    gp = GregProjector(L, emt=-28.8, td=24.8)

    # 1. blocks loaded
    assert len(gp.dg_set) == 84 and len(gp.hg_set) == 95 and len(gp.su_set) == 92, \
        (len(gp.dg_set), len(gp.hg_set), len(gp.su_set))

    # 2. DG positive, finite, monotone increasing in DBH for red spruce (97)
    dgs = [gp.dg_annual(97, d, 0.5, 4.5 + 4.0 * d, 90, 1200) for d in (6, 12, 20)]
    assert all(x > 0 for x in dgs)
    assert dgs[0] < dgs[1] < dgs[2], dgs

    # 3. survival is a bounded probability
    ps = gp.survival_annual(97, 0.5, 0.4)
    assert 0.0 <= ps <= 1.0 and gp.survival_annual(97, 0.0, 0.4) == 0.0

    # 4. period survival compounds below annual
    assert gp.survival_period(97, 0.5, 0.4, years=5) < ps

    # 5. species are distinct (oak grows faster than sugar maple at 12 in)
    assert gp.dg_annual(833, 12, 0.5, 52.5, 90, 1200) > \
        gp.dg_annual(318, 12, 0.5, 52.5, 90, 1200)

    # 6. fallback: an out-of-set SPCD (99999) resolves via SW/HW median, flagged
    assert gp.is_fallback(99999, "diameter_growth")
    fb = gp.dg_annual(99999, 12, 0.5, 52.5, 90, 1200)
    assert fb > 0
    # forcing a modal in-set fallback species reproduces that species exactly
    assert approx(gp.dg_annual(99999, 12, 0.5, 52.5, 90, 1200, fallback_spcd=97),
                  gp.dg_annual(97, 12, 0.5, 52.5, 90, 1200))

    # 7. EMT/TD default flag trips only when no per-stand climate is available
    gp2 = GregProjector(L)
    _ = gp2.climate(stand_cn="does_not_exist")
    assert gp2.used_default_climate is True

    # 8. regression pin vs the R reference projector (greg_est_dg, SPCD 97,
    #    dbh=12, cr=0.5, ht=52.5, bal=90, elev=1200, EMT=-28.8). R gives 0.09835.
    assert approx(gp.dg_annual(97, 12, 0.5, 52.5, 90, 1200, emt=-28.8),
                  0.098349, tol=1e-4), gp.dg_annual(97, 12, 0.5, 52.5, 90, 1200, emt=-28.8)

    print("ALL conus_greg projector tests passed.")


if __name__ == "__main__":
    main()
