"""Regression tests for the 2026-10-30 Bakuzis ensemble fixes.

D1  posterior-draw keyword records use the FVS fixed-field layout of the MAP block
D2  the in-process calibrated overlay never writes 0 for a species without a value,
    and posterior-draw runs do not get the MAP overlay
D3  the FIA stand generator seeds a condition's own trees at the condition acre
D4  draw multipliers are on the MAP serializer basis and centred on MAP

Synthetic fixtures only: the draws are generated from the tracked posterior
summaries plus noise, and the FIA plot fixture uses made-up identifiers.

Run:  python3 -m pytest config/tests/test_draws_overlay_standgen.py -q
"""
from __future__ import annotations

import json
import re
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
for p in (ROOT, ROOT / "calibration" / "python", ROOT / "deployment" / "fvs2py"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from config.calibration_multipliers import (  # noqa: E402
    baimult_from_dds, multipliers_from_params, params_from_summaries,
)
from config.config_loader import FvsConfigLoader  # noqa: E402
from config.uncertainty import UncertaintyEngine  # noqa: E402

CONFIG_DIR = ROOT / "config"
CAL_DIR = ROOT / "calibration"
VARIANTS = ("ne", "acd")


# --------------------------------------------------------------------------- helpers
def fvs_fields(record: str) -> tuple[str, list[str]]:
    """Split one keyword record by the keyrdr.f90 fixed-field rules:
    KEYWRD = cols 1-8 (we read 1-10), field i = cols 11+10*(i-1) .. 20+10*(i-1)."""
    kw = record[0:10].strip()
    fields = [record[10 + 10 * i: 20 + 10 * i].strip() for i in range(7)]
    return kw, fields


def records(block: str) -> list[tuple[str, list[str]]]:
    return [fvs_fields(ln) for ln in block.splitlines()
            if ln.strip() and not ln.startswith(("!!", "*"))]


def cal_config(v):
    return FvsConfigLoader(v, version="calibrated", config_dir=CONFIG_DIR).config


def synthetic_draws(v: str, n: int = 200, sd: float = 0.05, seed: int = 1) -> dict:
    """Draws file shaped like 06_posterior_to_json.R output: medians + N(0, sd) noise."""
    rng = np.random.default_rng(seed)
    med = params_from_summaries(CAL_DIR / "output" / "variants" / v)
    fia = [x for x in cal_config(v)["categories"]["species_definitions"]["FIAJSP"]]
    codes = sorted({int(x) for x in fia if str(x).strip().isdigit()})
    sd_med = {"Intercept": 6.0, "b_Intercept": 5.5}
    sd_med.update({f"r_SPCD[{c},Intercept]": 0.1 * ((c % 7) - 3) / 3 for c in codes})
    base = {"mortality": {k: v_ for k, v_ in med["mortality"].items()
                          if k == "Intercept" or k.startswith("r_SPCD[")},
            "diameter_growth": {k: v_ for k, v_ in med["diameter_growth"].items() if re.match(r"^b0\[\d+\]$", k)},
            "stand_density": sd_med}
    draws = {comp: [{k: float(x + rng.normal(0, sd)) for k, x in pars.items()} for _ in range(n)]
             for comp, pars in base.items()}
    return {"variant": v, "n_draws": n, "components": list(draws), "draws": draws}


@pytest.fixture(scope="module", params=VARIANTS)
def engine(request, tmp_path_factory):
    v = request.param
    f = tmp_path_factory.mktemp(f"draws_{v}") / f"{v}_draws.json"
    f.write_text(json.dumps(synthetic_draws(v)))
    return UncertaintyEngine(v, config_dir=CONFIG_DIR, draws_path=f, seed=42, calibration_dir=CAL_DIR)


# --------------------------------------------------------------------------- D4 basis
@pytest.mark.parametrize("v", VARIANTS)
def test_python_port_reproduces_shipped_map_multipliers(v):
    """multipliers_from_params on the posterior medians == calibration_multipliers in <v>.json."""
    cfg = cal_config(v)
    got = multipliers_from_params(params_from_summaries(CAL_DIR / "output" / "variants" / v), cfg, v, CAL_DIR)
    shipped = cfg["calibration_multipliers"]
    for k in ("mort_multiplier", "dds_multiplier", "htg_multiplier"):
        np.testing.assert_allclose(got[k], np.asarray(shipped[k], float), atol=1e-4, err_msg=k)


def test_draw_multiplier_medians_are_centred_on_map(engine):
    cfg = engine._calibrated_config()
    cm = cfg["calibration_multipliers"]
    map_bai = baimult_from_dds(np.asarray(cm["dds_multiplier"], float))
    map_mort = np.asarray(cm["mort_multiplier"], float)
    B, M, S = [], [], []
    for i in range(engine.n_draws):
        m = engine.compute_multipliers_for_draw(engine.get_draw(i))
        B.append(m["growth_mult"]); M.append(m["mort_mult"])
        if "sdi_values" in m:
            S.append(m["sdi_values"])
    np.testing.assert_allclose(np.median(B, axis=0), map_bai, rtol=0.02)
    np.testing.assert_allclose(np.median(M, axis=0), map_mort, rtol=0.03)
    # species that carry no calibrated multiplier (slot never assigned by the
    # serializer, so exactly 1.0 before rounding) stay exactly at unity in every draw
    port = multipliers_from_params(params_from_summaries(CAL_DIR / "output" / "variants" / engine.variant),
                                   cfg, engine.variant, CAL_DIR)
    assert np.all(np.asarray(B)[:, port["dds_multiplier"] == 1.0] == 1.0)
    assert np.all(np.asarray(M)[:, port["mort_multiplier"] == 1.0] == 1.0)
    sdi_map = FvsConfigLoader._find_sdi_param(cfg["categories"])
    if sdi_map is not None:
        base = np.array([np.nan if isinstance(x, str) or x is None else float(x) for x in sdi_map])
        med = np.nanmedian(np.asarray(S), axis=0)[: base.size]
        ok = np.isfinite(base)
        np.testing.assert_allclose(med[ok], base[ok], rtol=0.01)
        assert np.all(np.isnan(np.asarray(S)[:, ~ok[: np.asarray(S).shape[1]]]))


# --------------------------------------------------------------------------- D1 layout
def test_draw_keyfile_fields_parse_per_species(engine):
    draw = engine.get_draw(3)
    m = engine.compute_multipliers_for_draw(draw)
    blk = engine.generate_keywords_for_draw(draw, draw_idx=3, include_comments=False)
    recs = records(blk)
    sdimax = [f for kw, f in recs if kw == "SDIMAX"]
    sdi = m.get("sdi_values")
    if sdi is None:
        assert not sdimax  # variant without a calibrated SDI vector (ACD): none written, as for MAP
    else:
        want = {i + 1: round(float(x), 1) for i, x in enumerate(sdi) if np.isfinite(x) and x > 0}
        got = {int(f[0]): float(f[1]) for f in sdimax}   # option 89: species field 1, value field 2
        assert got == pytest.approx(want, abs=0.051)
        assert len(set(got.values())) > 1               # not one value for every species
        assert all(f[2] == "" for f in sdimax)
    for kw, key in (("MORTMULT", "mort_mult"), ("BAIMULT", "growth_mult")):
        rows = [f for k, f in recs if k == kw]
        assert rows, kw
        for f in rows:
            assert f[0] == ""                            # field 1 = date, blank -> cycle 1
            sp = int(f[1]); val = float(f[2])            # field 2 species, field 3 multiplier
            assert val == pytest.approx(m[key][sp - 1], abs=6e-5)
        assert {int(f[1]) for f in rows} == {i + 1 for i, x in enumerate(m[key]) if abs(x - 1) > 0.01}


def test_draw_and_map_blocks_share_the_formatter(engine):
    """Same values through generate_keywords and generate_keywords_for_draw give identical records."""
    v = engine.variant
    loader = FvsConfigLoader(v, version="calibrated", config_dir=CONFIG_DIR)
    map_recs = {kw for kw, _ in records(loader.generate_keywords(include_comments=False))}
    blk = engine.generate_keywords_for_draw(engine.get_draw(0), include_comments=False)
    for ln in blk.splitlines():
        kw = ln[:10].strip()
        assert kw in map_recs, kw
        if kw == "SDIMAX":
            assert re.fullmatch(r"SDIMAX {4}[ \d]{10}[ \d.]{10}", ln)
        else:
            assert re.fullmatch(rf"{kw:<10} {{10}}[ \d]{{10}}[ \d.]{{10}}([ \d.]{{10}}){{0,2}}", ln)


# --------------------------------------------------------------------------- D2 overlay
class FakeFVS:
    def __init__(self, n=108, fill=777.0):
        self.dims = {"maxspecies": n}
        self.attrs = {"spsdi": np.full(n, fill), "baimult": np.ones(n), "mortmult": np.ones(n), "htgmult": np.ones(n)}

    def get_species_attr(self, a):
        return self.attrs[a].copy()

    def set_species_attr(self, a, arr):
        self.attrs[a] = np.asarray(arr, float).copy()


def test_overlay_skips_na_species():
    loader = FvsConfigLoader("ne", version="calibrated", config_dir=CONFIG_DIR)
    sdi = FvsConfigLoader._find_sdi_param(loader.config["categories"])
    na = [i for i, x in enumerate(sdi) if isinstance(x, str) or x is None]
    assert na, "fixture expects NE SDICON to carry NA entries"
    fake = FakeFVS()
    loader.apply_to_fvs(fake)
    out = fake.attrs["spsdi"]
    assert not np.any(out == 0)
    assert np.all(out[na] == 777.0)
    for i, x in enumerate(sdi):
        if not (isinstance(x, str) or x is None) and x > 0:
            assert out[i] == pytest.approx(float(x))


def test_apply_draw_skips_missing_sdi(engine):
    fake = FakeFVS()
    engine.apply_draw_to_fvs(fake, engine.get_draw(1), {})
    assert not np.any(fake.attrs["spsdi"] == 0)


def test_posterior_runs_do_not_get_the_map_overlay(monkeypatch, tmp_path):
    import bakuzis_uncertainty_comparison as buc
    seen = []

    class StubFVS:
        def __init__(self, lib_path, config_version=None, config_dir=None):
            seen.append(config_version)

        def load_keyfile(self, p):
            pass

        def run(self):
            pass

        summary = pd.DataFrame({"year": [2000], "atba": [1.0]})

    monkeypatch.setitem(sys.modules, "fvs2py", types.SimpleNamespace(FVS=StubFVS))
    (tmp_path / "FVSne.so").write_bytes(b"")
    monkeypatch.setattr(buc, "FVS_LIB_DIR", str(tmp_path))
    sdf = pd.DataFrame([{"stand_id": "S1"}]); tdf = pd.DataFrame([{"stand_id": "S1"}])
    for cv in ("default", "calibrated", "posterior"):
        buc.run_scenario(sdf, tdf, "S1", "ne", cv, draw_keywords="!! x")
    assert seen == ["default", "calibrated", "default"]


# --------------------------------------------------------------------------- D3 stand generator
def _fixture_fia(tmp_path):
    """One forested plot with two conditions (0.6 / 0.4) plus a single-condition plot."""
    cond = pd.DataFrame([
        dict(STATECD=23, COUNTYCD=1, PLOT=1, INVYR=2015, PLT_CN=1001, CONDID=1, CONDPROP_UNADJ=0.6, COND_STATUS_CD=1, SICOND=60, SIBASE=50, BALIVE=np.nan, FORTYPCD=121, STDAGE=50, ASPECT=0, SLOPE=5, ELEV=900),
        dict(STATECD=23, COUNTYCD=1, PLOT=1, INVYR=2015, PLT_CN=1001, CONDID=2, CONDPROP_UNADJ=0.4, COND_STATUS_CD=1, SICOND=62, SIBASE=50, BALIVE=np.nan, FORTYPCD=801, STDAGE=70, ASPECT=0, SLOPE=5, ELEV=900),
        dict(STATECD=23, COUNTYCD=1, PLOT=2, INVYR=2015, PLT_CN=1002, CONDID=1, CONDPROP_UNADJ=1.0, COND_STATUS_CD=1, SICOND=58, SIBASE=50, BALIVE=np.nan, FORTYPCD=121, STDAGE=40, ASPECT=0, SLOPE=5, ELEV=900),
    ])
    rng = np.random.default_rng(7)
    trees = []
    for cn, cid, n in ((1001, 1, 18), (1001, 2, 14), (1002, 1, 22)):
        for k in range(n):
            trees.append(dict(PLT_CN=cn, CONDID=cid, STATECD=23, COUNTYCD=1, PLOT=1, SUBP=1 + k % 4, TREE=k + 1, STATUSCD=1,
                              SPCD=(12 if cid == 1 else 316), DIA=round(float(rng.uniform(5, 16)), 1), HT=50.0, CR=40.0, TPA_UNADJ=6.018046))
    # one P2A_GRM_FLG style live record with a null TPA_UNADJ (PR 129 drop must survive)
    trees.append(dict(PLT_CN=1001, CONDID=1, STATECD=23, COUNTYCD=1, PLOT=1, SUBP=7, TREE=99, STATUSCD=1, SPCD=12, DIA=12.0, HT=50.0, CR=40.0, TPA_UNADJ=np.nan))
    tr = pd.DataFrame(trees)
    live = tr[tr["TPA_UNADJ"].notna()].copy()
    live["ba"] = 0.005454154 * live["DIA"] ** 2 * live["TPA_UNADJ"]
    bal = live.groupby(["PLT_CN", "CONDID"])["ba"].sum()
    cond["BALIVE"] = [bal[(r.PLT_CN, r.CONDID)] / r.CONDPROP_UNADJ for r in cond.itertuples()]
    cond.to_csv(tmp_path / "ME_COND.csv", index=False)
    tr.to_csv(tmp_path / "ME_TREE.csv", index=False)
    return cond


def test_stand_generator_reproduces_balive_on_multicondition_plot(tmp_path):
    import fia_stand_generator as g
    cond = _fixture_fia(tmp_path)
    trees = g.load_fia_trees("acd", cond["PLT_CN"].tolist(), tmp_path)
    loaded = g.load_fia_conditions("acd", tmp_path)
    assert {"CONDID", "CONDPROP_UNADJ"}.issubset(loaded.columns)
    for _, r in loaded.iterrows():
        ti = g.build_tree_init(g.condition_trees(trees, r), "S")
        ba = float((0.005454154 * ti["diameter"] ** 2 * ti["tree_count"]).sum())
        assert ba == pytest.approx(float(r["BALIVE"]), rel=0.01), (r["PLT_CN"], r["CONDID"])
        assert ti["tree_count"].notna().all()
    # end to end through generate_real_stand: every returned stand matches its condition's BALIVE
    n_checked = 0
    for dens in ("Low", "Medium", "High"):
        for stand_df, tree_df, row in g.generate_real_stand("acd", "Medium", dens, fia_dir=tmp_path, n_plots=5, seed=1):
            ba = float((0.005454154 * tree_df["diameter"] ** 2 * tree_df["tree_count"]).sum())
            assert ba == pytest.approx(float(row["BALIVE"]), rel=0.01)
            n_checked += 1
    assert n_checked == len(loaded), (n_checked, cond["BALIVE"].round(1).tolist())


# --------------------------------------------------------------------------- per-stand band output
def test_aggregator_keeps_within_stand_band_separate_from_stand_spread():
    import bakuzis_uncertainty_aggregate as agg
    rng = np.random.default_rng(3)
    rows = []
    for rep, level in ((1, 100.0), (2, 200.0)):          # two stands, very different levels
        for cfg, draws in (("default", [np.nan]), ("calibrated_map", [np.nan]), ("posterior", range(100))):
            for di in draws:
                val = level * (1 + (rng.normal(0, 0.01) if cfg == "posterior" else 0))
                for k, yr in enumerate((2000, 2005)):
                    rows.append(dict(year=yr, scenario=1, species_group="Pine", site_class="Low", density_class="Low",
                                     variant="NE", config=cfg, draw_id=di, replicate=rep, atba=val + k, tpa=500.0))
    df = pd.DataFrame(rows)
    stand, cell = agg.summarize_stand_bands(df)
    c = cell[(cell["variable"] == "atba") & (cell["year"] == 2005)].iloc[0]
    assert c["within_band_pct_med"] < 6            # about 4 percent: 2 x 1.96 x 1 percent noise
    assert c["map_range_pct"] > 60                 # (201 - 101) / 151
    pooled = agg.summarize_ensemble(df, [0.025, 0.975])
    pb = pooled[(pooled["config"] == "posterior_band") & (pooled["variable"] == "atba") & (pooled["year"] == 2005)].iloc[0]
    assert pb["band_basis"] == "pooled_draws_x_stands"
    assert 100 * (pb["q975"] - pb["q025"]) / pb["median"] > 60   # the legacy band mixes in the stand spread


def test_condition_expansion_uses_the_plot_size_proportions():
    """D5: microplot, subplot and macroplot trees are expanded by MICRPROP_UNADJ, SUBPPROP_UNADJ and
    MACRPROP_UNADJ; a missing proportion falls back to CONDPROP_UNADJ."""
    import fia_stand_generator as g
    cond = pd.Series(dict(PLT_CN=7, CONDID=2, CONDPROP_UNADJ=0.5, MICRPROP_UNADJ=0.25, SUBPPROP_UNADJ=0.5, MACRPROP_UNADJ=np.nan))
    trees = pd.DataFrame(dict(PLT_CN=[7, 7, 7, 7], CONDID=[2, 2, 2, 1], DIA=[2.0, 8.0, 30.0, 8.0], TPA_UNADJ=[74.965, 6.018, 0.999, 6.018]))
    t = g.condition_trees(trees, cond)
    assert len(t) == 3
    np.testing.assert_allclose(t["TPA_UNADJ"].to_numpy(), [74.965 / 0.25, 6.018 / 0.5, 0.999 / 0.5])
    # BALIVE as FIA defines it (per-size proportions) is reproduced exactly
    ba = (0.005454154 * t["DIA"] ** 2 * t["TPA_UNADJ"]).sum()
    ref = 0.005454154 * (2.0 ** 2 * 74.965 / 0.25 + 8.0 ** 2 * 6.018 / 0.5 + 30.0 ** 2 * 0.999 / 0.5)
    assert abs(ba / ref - 1) < 1e-12
