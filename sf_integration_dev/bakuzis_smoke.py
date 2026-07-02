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

# a couple of region-appropriate species per variant (FIA SPCD) for the stand
VARIANT_SPP = {
    "ne": [12, 97, 316], "pn": [202, 17, 242], "wc": [202, 17, 263],
    "ls": [12, 95, 371], "sn": [131, 110, 611], "ca": [15, 122, 81],
    "em": [93, 19, 202], "ie": [202, 93, 19],
}


def synth_stand(variant, sicond):
    spp = VARIANT_SPP.get(variant.lower(), [202, 93, 122])
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
    out["site"] = sicond
    return out


def bakuzis(traj, variant, max_sdi=450.0):
    """Run the law-like checks across the site gradient. Returns (flags, summary)."""
    flags = []
    sites = sorted(traj["site"].unique())
    # 2. monotonic development within each site
    for si in sites:
        d = traj[traj["site"] == si].sort_values("age")
        for v, lab, want_up in [("topht", "dominant height", True),
                                ("ba", "basal area", True),
                                ("tcuft", "volume", True),
                                ("qmd", "QMD", True),
                                ("tpa", "TPA", False)]:
            x = d[v].to_numpy()
            x = x[np.isfinite(x)]
            if len(x) < 3: continue
            dif = np.diff(x)
            tol = 0.02 * np.nanmax(np.abs(x)) + 1e-6
            if want_up and np.any(dif < -tol):
                flags.append(f"[{variant} site {si}] {lab} decreases over time (non-monotonic)")
            if not want_up and np.any(dif > tol):
                flags.append(f"[{variant} site {si}] {lab} increases over time (TPA should not rise)")
    # 1. site ordering at a common late age
    late = traj[traj["age"] >= traj["age"].max() - 5]
    piv = late.groupby("site")[["topht", "tcuft"]].mean()
    for v, lab in [("topht", "dominant height"), ("tcuft", "volume")]:
        vals = piv[v].to_numpy()
        if len(vals) >= 2 and not np.all(np.diff(vals) >= -1e-6):
            flags.append(f"[{variant}] site ordering violated: {lab} not increasing with site index "
                         f"({dict(zip(piv.index.tolist(), np.round(vals,1).tolist()))})")
    # 3. Reineke self-thinning ceiling
    mx = traj["sdi"].max()
    if mx > max_sdi * 1.15:
        flags.append(f"[{variant}] max SDI {mx:.0f} exceeds plausible Reineke max ~{max_sdi:.0f}")
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
    return flags, eich_note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--variants", default="ne")
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
        flags, eich = bakuzis(traj, v)
        print(f"  {eich}")
        if flags:
            print(f"  *** {len(flags)} REALISM FLAG(S):")
            for f in flags: print("     -", f)
        else:
            print("  OK: all Bakuzis law-like checks pass")
        report[v] = {"status": "ok", "n_sites": len(trajs),
                     "flags": flags, "eichhorn": eich}
    with open(os.path.join(a.out, "bakuzis_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(f"\n[report] {a.out}/bakuzis_report.json")


if __name__ == "__main__":
    main()
