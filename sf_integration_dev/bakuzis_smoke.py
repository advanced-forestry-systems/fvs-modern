#!/usr/bin/env python3
"""
bakuzis_smoke.py -- one-stand projection smoke test + Bakuzis realism pass
per variant under conus_hybrid (out-of-process FVS executable).

For each variant it projects a synthetic stand across a 3-level site gradient
(low/med/high site index), extracts the FVS summary trajectory, and runs the
law-like (Bakuzis / Leary 1997) checks:
  1. site ordering   : higher site -> higher dominant height & volume at equal age
  2. monotonic dev   : dominant height, BA, volume non-decreasing over time;
                       TPA non-increasing (mortality); QMD increasing
  3. self-thinning   : Reineke SDI stays below a plausible maximum (Reineke rule)
  4. Eichhorn        : volume is ~a function of dominant height across site classes
Flags any biologically implausible trajectory. Exit is non-fatal: a variant that
fails to project is reported, not crashed on.
"""
import os, sys, json, math, argparse, traceback
import numpy as np, pandas as pd

PROJECT_ROOT = os.path.expanduser("~/fvs-modern")
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
sys.path.insert(0, os.path.expanduser("~/fvs-conus/python"))
import perseus_100yr_projection as P
# point perseus at the freshly built executables (subprocess path, no in-proc segfault)
P.FVS_LIB_DIR = os.environ.get("FVS_BIN_DIR", "/fs/scratch/PUOM0008/crsfaaron/bin-sf85")
CONFIG_DIR = os.environ.get("FVS_CONFIG_DIR", "/fs/scratch/PUOM0008/crsfaaron/wt-gompit/config")

# per-variant species (FIA SPCD) + max-SDI ceiling, derived from the configs
# (categories_conus.diameter_growth.species_intercepts + site_index.SDICON/FMSDI).
# Loaded from variant_params.json; falls back to a generic set if absent.
import json as _json
_PARAMS_PATH = os.environ.get("VARIANT_PARAMS",
                              "/users/PUOM0008/crsfaaron/variant_params.json")
try:
    VARIANT_PARAMS = _json.load(open(_PARAMS_PATH))
except Exception:
    VARIANT_PARAMS = {}
# Sanity ceiling for the self-thinning check when the config carries no SDI max.
# Western conifer variants tolerate higher SDI than eastern; used only as an
# upper bound on runaway density, not as a hard Reineke max.
_WEST = {"pn", "wc", "ie", "ca", "nc", "so", "op", "ec", "wsc", "ws", "ut", "tt",
         "cr", "em", "bm", "ci", "kt", "oc", "ak", "bc", "cs"}
# Northern Rockies habitat-type variants: height growth keyed to habitat type,
# not the site-index field (verified: topht invariant to SICOND and site species).
HABITAT_TYPE_VARIANTS = {"ie", "ci", "kt"}
def variant_max_sdi(variant):
    p = VARIANT_PARAMS.get(variant.lower(), {})
    m = p.get("max_sdi")
    if m and m > 0:
        return float(m)
    return 1100.0 if variant.lower() in _WEST else 800.0


def synth_stand(variant, sicond):
    spp = (VARIANT_PARAMS.get(variant.lower(), {}).get("spp")
           or [202, 93, 122])
    plot = {"INVYR": 2010, "LAT": 45.0, "LON": -110.0 if variant not in ("ne","ls","sn") else -69.0,
            "ELEV": 500, "SLOPE": 15, "ASPECT": 180, "STDAGE": 30,
            "SICOND": sicond, "SISP": spp[0]}
    rows = []
    dbhs = [4.0, 6.0, 8.0, 10.0]; tpas = [80, 60, 40, 25]; hts = [24, 34, 45, 55]
    tid = 0
    for sp in spp:
        for dbh, tpa, ht in zip(dbhs, tpas, hts):
            tid += 1
            rows.append({"DIA": dbh, "HT": ht, "CR": 45, "SPCD": sp,
                         "TPA_UNADJ": tpa / len(spp), "SUBP": 1, "TREE": tid})
    return plot, pd.DataFrame(rows)


def hybrid_keywords(variant):
    try:
        from config.config_loader import FvsConfigLoader
        return FvsConfigLoader(variant.lower(), version="conus_hybrid",
                               config_dir=CONFIG_DIR).generate_keywords(include_comments=False)
    except Exception as e:
        print(f"  [warn] hybrid keywords unavailable for {variant}: {e}")
        return None


def qmd(ba, tpa):
    return math.sqrt(ba / (0.005454 * tpa)) if (tpa and ba and tpa > 0) else float("nan")


def project(variant, sicond, num_cycles=20, cycle_length=5):
    plot, trees = synth_stand(variant, sicond)
    sdf = P.build_fvs_standinit(plot, f"S_{variant}_{sicond}", variant.lower())
    tdf = P.build_fvs_treeinit(trees, f"S_{variant}_{sicond}")
    kw = hybrid_keywords(variant) or "** DEFAULT PARAMETERS"
    res = P.run_fvs_projection(sdf, tdf, f"S_{variant}_{sicond}", variant.lower(),
                               config_version=None, num_cycles=num_cycles,
                               cycle_length=cycle_length,
                               extra_keywords=(kw if kw and kw != "** DEFAULT PARAMETERS" else ""))
    s = res["summary"]
    if s is None or not len(s):
        return None
    s = s.copy()
    cols = {c.lower(): c for c in s.columns}
    def col(*names):
        for n in names:
            if n in cols: return s[cols[n]].astype(float)
        return pd.Series([float("nan")] * len(s))
    out = pd.DataFrame({
        "year": col("year"), "age": col("age"),
        "tpa": col("tpa"), "ba": col("atba", "ba"),
        "topht": col("attopht", "topht"), "tcuft": col("tcuft"),
    })
    out["qmd"] = [qmd(b, t) for b, t in zip(out["ba"], out["tpa"])]
    out["sdi"] = [t * (q / 10.0) ** 1.605 if (q == q and t) else float("nan")
                  for t, q in zip(out["tpa"], out["qmd"])]
    # engine-reported Reineke SDI and the variant's own SDIMax (authoritative;
    # avoids the hand-rolled SDI + guessed-ceiling unit issues)
    out["reineke_sdi"] = col("reinekesdi", "sdi")
    out["sdimax"] = col("sdimax")
    out["site"] = sicond
    return out


def bakuzis(traj, variant, max_sdi=1000.0):
    """Run the law-like checks across the site gradient. Returns (flags, notes, eich)."""
    flags = []
    notes = []
    sites = sorted(traj["site"].unique())
    # 0. site RESPONSE: dominant height should differ across the site gradient.
    #    If it is flat, the site index is not reaching the height-growth model.
    #    EXCEPTION: the Northern Rockies habitat-type variants (IE, CI, KT) drive
    #    height growth on a habitat-type code, not the site-index field, so a
    #    SICOND gradient legitimately does not move them. Verified empirically:
    #    topht is invariant to both SICOND (30/55/80) and site species. For these,
    #    a flat response is expected, not a flag — vary habitat type to test them.
    late0 = traj[traj["age"] >= traj["age"].max() - 5]
    ht_by_site = late0.groupby("site")["topht"].mean()
    if len(ht_by_site) >= 2:
        spread = float(ht_by_site.max() - ht_by_site.min())
        if spread < 2.0:
            if variant.lower() in HABITAT_TYPE_VARIANTS:
                notes.append(f"[{variant}] site index does not drive height growth (expected: "
                             f"habitat-type variant; productivity keyed to habitat type, not SICOND). "
                             f"Vary habitat type to exercise a productivity gradient.")
            else:
                flags.append(f"[{variant}] site index not driving height growth: dominant height "
                             f"flat across sites ({dict(ht_by_site.round(1))}) — check site species/site input")
    # 2. monotonic development within each site. Dominant height, BA, and volume
    #    should rise. TPA and QMD are handled JOINTLY: TPA rising while QMD falls
    #    is the signature of natural regeneration / ingrowth (FVS ESTAB), which is
    #    biologically realistic, so it is a note, not a flag. Only genuinely
    #    implausible motion is flagged: TPA rising without QMD falling (trees
    #    appearing with no size dilution) or QMD falling without ingrowth.
    ingrowth_seen = False
    for si in sites:
        d = traj[traj["site"] == si].sort_values("age")
        for v, lab in [("topht", "dominant height"), ("ba", "basal area"), ("tcuft", "volume")]:
            x = d[v].to_numpy(); x = x[np.isfinite(x)]
            if len(x) < 3: continue
            tol = 0.02 * np.nanmax(np.abs(x)) + 1e-6
            if np.any(np.diff(x) < -tol):
                flags.append(f"[{variant} site {si}] {lab} decreases over time (non-monotonic)")
        tpa = d["tpa"].to_numpy(); qmd = d["qmd"].to_numpy()
        ok = np.isfinite(tpa) & np.isfinite(qmd)
        tpa, qmd = tpa[ok], qmd[ok]
        if len(tpa) < 3: continue
        dt = np.diff(tpa); dq = np.diff(qmd)
        ttol = 0.02 * np.nanmax(np.abs(tpa)) + 1e-6
        qtol = 0.02 * np.nanmax(np.abs(qmd)) + 1e-6
        # implausible: TPA rises while QMD does not fall (no ingrowth to explain it)
        if np.any((dt > ttol) & (dq >= -qtol)):
            flags.append(f"[{variant} site {si}] TPA rises without QMD falling (unexplained tree gain)")
        # implausible: QMD falls while TPA does not rise (size decline without ingrowth)
        if np.any((dq < -qtol) & (dt <= ttol)):
            flags.append(f"[{variant} site {si}] QMD falls without ingrowth (unexplained size loss)")
        # benign: concurrent TPA-up / QMD-down = regeneration establishing
        if np.any((dt > ttol) & (dq < -qtol)):
            ingrowth_seen = True
    if ingrowth_seen:
        notes.append(f"[{variant}] TPA rises with QMD falling in places — natural regeneration/ingrowth active (expected)")
    # 1. site ordering at a common late age
    late = traj[traj["age"] >= traj["age"].max() - 5]
    piv = late.groupby("site")[["topht", "tcuft"]].mean()
    for v, lab in [("topht", "dominant height"), ("tcuft", "volume")]:
        vals = piv[v].to_numpy()
        if len(vals) >= 2 and not np.all(np.diff(vals) >= -1e-6):
            flags.append(f"[{variant}] site ordering violated: {lab} not increasing with site index "
                         f"({dict(zip(piv.index.tolist(), np.round(vals,1).tolist()))})")
    # 3. Reineke self-thinning: a stand that reaches its max density should
    #    PLATEAU (mortality balances growth), not grow density without bound. So
    #    flag only runaway density: SDI above the variant ceiling AND still
    #    climbing late in the projection. A stand sitting at a high but stable
    #    self-thinning limit (flat max-BA across sites) is correct, not a flag.
    for si in sites:
        d = traj[traj["site"] == si].sort_values("age")
        # Prefer the engine's own ReinekeSDI vs SDIMax (authoritative, unit-exact);
        # fall back to the hand-rolled SDI + config/default ceiling only if absent.
        rs = d["reineke_sdi"].to_numpy(); sm = d["sdimax"].to_numpy()
        if np.isfinite(rs).any() and np.isfinite(sm).any() and np.nansum(sm) > 0:
            series = rs[np.isfinite(rs)]; ceiling = float(np.nanmax(sm)); src = "engine ReinekeSDI/SDIMax"
        else:
            series = d["sdi"].to_numpy(); series = series[np.isfinite(series)]
            ceiling = max_sdi; src = "hand-rolled SDI/ceiling"
        if len(series) < 4:
            continue
        mx = float(np.nanmax(series))
        tail = series[-3:]
        still_climbing = (tail[-1] - tail[0]) / max(tail[0], 1.0) > 0.05
        if mx > ceiling * 1.02 and still_climbing:
            flags.append(f"[{variant} site {si}] runaway density: SDI {mx:.0f} > SDIMax {ceiling:.0f} "
                         f"({src}) and still climbing late (not reaching a self-thinning limit)")
    # 4. Eichhorn: volume vs dominant height similar across sites (compare volume at
    #    matched top height); large spread => flag
    eich_note = ""
    try:
        thts = np.linspace(traj["topht"].quantile(0.4), traj["topht"].quantile(0.8), 4)
        spreads = []
        for th in thts:
            vv = []
            for si in sites:
                d = traj[traj["site"] == si]
                i = (d["topht"] - th).abs().idxmin()
                vv.append(d.loc[i, "tcuft"])
            vv = np.array(vv, float)
            if np.nanmean(vv) > 0:
                spreads.append(np.nanstd(vv) / np.nanmean(vv))
        cv = float(np.nanmean(spreads)) if spreads else float("nan")
        eich_note = f"Eichhorn vol~height CV across sites = {cv:.2f}"
        if cv > 0.35:
            flags.append(f"[{variant}] Eichhorn weak: volume at matched height varies {cv*100:.0f}% across sites")
    except Exception:
        pass
    return flags, notes, eich_note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants",
                    default="ne,pn,wc,ls,sn,ie,ca,em,nc,oc,op,tt,ut,ws,bm,cs,ec,ci,acd,kt,so,cr")
    ap.add_argument("--sites", default="40,55,70")
    ap.add_argument("--out", default=os.path.expanduser("~/bakuzis_smoke_out"))
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    sites = [int(x) for x in a.sites.split(",")]
    report = {}
    for v in a.variants.split(","):
        v = v.strip()
        print(f"\n===== {v.upper()} =====")
        exe = os.path.join(P.FVS_LIB_DIR, f"FVS{v.lower()}")
        if not (os.path.exists(exe) and os.access(exe, os.X_OK)):
            print(f"  [skip] no executable at {exe}")
            report[v] = {"status": "no_executable"}
            continue
        trajs = []
        for si in sites:
            try:
                t = project(v, si)
                if t is None or not len(t):
                    print(f"  [site {si}] projection returned no summary")
                    continue
                nyr = len(t)
                print(f"  [site {si}] projected {nyr} periods; "
                      f"age {t['age'].min():.0f}->{t['age'].max():.0f}, "
                      f"topht {t['topht'].iloc[0]:.0f}->{t['topht'].iloc[-1]:.0f}ft, "
                      f"BA {t['ba'].iloc[0]:.0f}->{t['ba'].iloc[-1]:.0f}, "
                      f"vol {t['tcuft'].iloc[-1]:.0f}cuft, maxSDI {t['sdi'].max():.0f}")
                trajs.append(t)
            except Exception as e:
                print(f"  [site {si}] PROJECTION FAILED: {e}")
        if not trajs:
            report[v] = {"status": "projection_failed"}
            continue
        traj = pd.concat(trajs, ignore_index=True)
        traj.to_csv(os.path.join(a.out, f"traj_{v}.csv"), index=False)
        flags, notes, eich = bakuzis(traj, v, max_sdi=variant_max_sdi(v))
        print(f"  {eich}")
        for n in notes: print("  note:", n)
        if flags:
            print(f"  *** {len(flags)} REALISM FLAG(S):")
            for f in flags: print("     -", f)
        else:
            print("  OK: all Bakuzis law-like checks pass"
                  + ("  (with expected-behavior note)" if notes else ""))
        report[v] = {"status": "ok", "n_sites": len(trajs),
                     "flags": flags, "notes": notes, "eichhorn": eich,
                     "max_sdi_ceiling": variant_max_sdi(v)}
    with open(os.path.join(a.out, "bakuzis_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[report] {a.out}/bakuzis_report.json")


if __name__ == "__main__":
    main()
