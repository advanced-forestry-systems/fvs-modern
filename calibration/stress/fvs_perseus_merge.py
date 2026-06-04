#!/usr/bin/env python3
"""fvs_perseus_merge.py -- inject the national FVS engines into PERSEUS api/series.

Mirrors ycx_merge_perseus.py but for the full-FIADB CONUS FVS campaign. Three
national models share cls "FVS":

  fvs_national_default_v1     default growth + native (Dixon/VARMRT) mortality
  fvs_national_calibrated_v1  Bayesian-calibrated growth + native mortality
  fvs_national_gompit_v1      default growth + Johnson national gompit mortality

Each model is one campaign arm aggregated to per-state densities by
fvs_perseus_aggregate.py into
  <series_root>/perseus_series_<config>/ycx_<ST>_state_series.csv   (density)
  <series_root>/perseus_series_<config>/fvs_<config>_state_counts.csv (n_plots)

Stock metrics, native densities:
  agc_live_total  Mg C/ha -> Tg C
  agb_dry         Mg/ha   -> Tg

Both are physical (no unit calibration). State totals use a single fixed area
model per state, anchored so the 2025 inventory carbon reproduces fia.json
tg_agc. Because all three arms share the same FIADB inventory, the same area is
applied to every arm; they all start at the FIA anchor in 2025 and diverge over
the projection. Scenario bucket: reserve (no harvest).

Usage:
  python3 fvs_perseus_merge.py <repo_dir> <series_root> [config,config,...]
  (default configs: default,calibrated,gompit -- missing dirs are skipped)
"""
from __future__ import annotations
import csv, json, sys, os, glob, statistics
from collections import defaultdict

repo = sys.argv[1]
series_root = sys.argv[2]
CONFIGS_ARG = (sys.argv[3].split(",") if len(sys.argv) > 3
               else ["default", "calibrated", "gompit"])

api = os.path.join(repo, "public", "api")
fia = json.load(open(os.path.join(api, "fia.json")))

CLS = "FVS"
# FVS reports the inventory (PROJ_YEAR 0 = 2025) treelist before it imputes
# missing tree heights, and NSBE biomass needs height, so the 2025 point is
# systematically understated (50-90% low, in proportion to each state's
# missing-height fraction) and jumps by 2030. We therefore anchor and report
# the FVS series from 2030, the first year with complete (FVS-filled) heights.
START = 2030
# config -> (model id, label)
MODELS = {
    "default":    ("fvs_national_default_v1",
                   "FVS national default (FIADB-initialized, native mortality, no harvest)"),
    "calibrated": ("fvs_national_calibrated_v1",
                   "FVS national calibrated (Bayesian posterior, native mortality, no harvest)"),
    "gompit":     ("fvs_national_gompit_v1",
                   "FVS national gompit mortality (Johnson national hazard, no harvest)"),
}
# metric -> kind. both physical Tg, no unit calibration.
MET = {"agc_live_total": "tg", "agb_dry": "tg"}
BUCKET = "reserve (no harvest)"
MANAGED_BUCKET = "managed (harvest)"
# optional: dir containing managed_<cfg>/managed_<ST>.csv (harvest+disturbance
# scenario from fvs_managed_scenario.py). When set, a managed (harvest) bucket
# and harvest_c_yr flux are injected alongside the reserve trajectory.
MANAGED_ROOT = os.environ.get("FVS_MANAGED_ROOT")


def load_native(path):
    """metric -> mgmt -> year -> (v, lo, hi)."""
    d = defaultdict(lambda: defaultdict(dict))
    for r in csv.DictReader(open(path)):
        d[r["metric"]][r["mgmt"]][int(r["year"])] = (
            float(r["value"]), float(r["value_lo"]), float(r["value_hi"]))
    return d


def load_counts(path):
    out = {}
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            out[r["state"]] = int(r["n_plots"])
    return out


def load_managed(path):
    """metric -> {year: density} for the managed scenario (one state)."""
    d = defaultdict(dict)
    if os.path.exists(path):
        for r in csv.DictReader(open(path)):
            d[r["metric"]][int(r["year"])] = float(r["value"])
    return d


# ---- load every available config arm ----
arms = {}          # config -> {st -> native}
counts = {}        # config -> {st -> n_plots}
for cfg in CONFIGS_ARG:
    d = os.path.join(series_root, f"perseus_series_{cfg}")
    if not os.path.isdir(d):
        print(f"  skip {cfg}: no dir {d}")
        continue
    files = sorted(glob.glob(os.path.join(d, "ycx_*_state_series.csv")))
    if not files:
        print(f"  skip {cfg}: no series files in {d}")
        continue
    arms[cfg] = {os.path.basename(f).split("_")[1]: load_native(f) for f in files}
    counts[cfg] = load_counts(os.path.join(d, f"fvs_{cfg}_state_counts.csv"))
    print(f"  {cfg}: {len(arms[cfg])} states, "
          f"{len(counts[cfg])} count rows")

if not arms:
    sys.exit("no FVS config arms found; nothing to merge")

# reference arm for the shared area model (prefer default)
REF = "default" if "default" in arms else next(iter(arms))
ref_native, ref_counts = arms[REF], counts[REF]
all_states = sorted(set().union(*[set(a) for a in arms.values()]))

# ---- fixed state area model anchored to FIA carbon ----
# total_area_ha[st] = fia_tg[st]*1e6 / agc_density_2025(REF)
# (= the area whose REF 2025 carbon density integrates to the FIA total).
area_ha = {}
ratios = []  # ha per plot, for fallback
for st in all_states:
    tg = fia.get(st, {}).get("tg_agc")
    nat = ref_native.get(st)
    if tg and nat and "agc_live_total" in nat:
        d25 = nat["agc_live_total"][BUCKET][START][0]
        if d25 > 0:
            area_ha[st] = tg * 1e6 / d25
            npl = ref_counts.get(st, 0)
            if npl > 0:
                ratios.append(area_ha[st] / npl)
A0_med = statistics.median(ratios) if ratios else 6000.0
for st in all_states:                     # fallback: median ha/plot * n_plots
    if st not in area_ha:
        npl = ref_counts.get(st) or max(
            (counts[c].get(st, 0) for c in counts), default=0)
        area_ha[st] = npl * A0_med
print(f"area model: {len(ratios)} FIA-anchored states, "
      f"median {A0_med:.0f} ha/plot")


def phys_total(density, st):
    return density * area_ha[st] / 1e6     # Mg/ha -> Tg


OUR_MODELS = {MODELS[c][0] for c in arms}


def fvs_rows(ser):
    """(distinct OUR_MODELS present, total pts of OUR_MODELS) in a series dict."""
    present, pts = set(), 0
    for met in ser:
        for bk in ser[met]:
            for s in ser[met][bk]:
                if s.get("model") in OUR_MODELS:
                    present.add(s["model"])
                    pts += len(s.get("pts", []))
    return present, pts


# ---- inject (incremental, non-destructive bookkeeping like ycx_merge) ----
stmeta = json.load(open(os.path.join(api, "states.json")))
META = json.load(open(os.path.join(api, "meta.json")))
touched = 0
delta_rows_total = 0
globally_new = set()               # OUR_MODELS not present anywhere before
for st in all_states:
    spath = os.path.join(api, "series", f"{st}.json")
    ser = json.load(open(spath)) if os.path.exists(spath) else {}
    before_models, before_pts = fvs_rows(ser)
    metrics_here = set()
    for cfg in arms:
        nat = arms[cfg].get(st)
        if not nat:
            continue
        model, label = MODELS[cfg]
        for metric in MET:
            if metric not in nat or BUCKET not in nat[metric]:
                continue
            nb = nat[metric][BUCKET]
            pts = []
            for y in sorted(nb):
                if y < START:           # drop pre-2030 height-fill artifact
                    continue
                v, lo, hi = nb[y]
                pts.append([y, round(phys_total(v, st), 3),
                               round(phys_total(lo, st), 3),
                               round(phys_total(hi, st), 3)])
            node = ser.setdefault(metric, {}).setdefault(BUCKET, [])
            node[:] = [s for s in node if s.get("model") != model]   # idempotent
            node.append({"model": model, "cls": CLS, "label": label, "pts": pts})
            metrics_here.add(metric)

        # ---- managed (harvest) scenario, if available ----
        if MANAGED_ROOT:
            mpath = os.path.join(MANAGED_ROOT, f"managed_{cfg}",
                                 f"managed_{st}.csv")
            mg = load_managed(mpath)
            mlabel = label.replace("no harvest",
                                   "harvest+disturbance, conus_hcs")
            for metric in ("agc_live_total", "agb_dry", "harvest_c_yr"):
                if metric not in mg:
                    continue
                nb = mg[metric]
                pts = [[y, round(phys_total(v, st), 3)]
                       for y, v in sorted(nb.items()) if y >= START]
                node = ser.setdefault(metric, {}).setdefault(MANAGED_BUCKET, [])
                node[:] = [s for s in node if s.get("model") != model]
                node.append({"model": model, "cls": CLS, "label": mlabel,
                             "pts": pts})
                metrics_here.add(metric)
    json.dump(ser, open(spath, "w"), separators=(",", ":"))
    touched += 1

    after_models, after_pts = fvs_rows(ser)
    delta_eng = len(after_models) - len(before_models)
    delta_rows = after_pts - before_pts
    delta_rows_total += delta_rows
    globally_new |= (after_models - before_models)

    sm = stmeta.get(st)
    if sm is not None:
        if delta_eng > 0:
            sm["engines"] = sm.get("engines", 0) + delta_eng
        sm["rows"] = sm.get("rows", 0) + delta_rows
        sm["has_series"] = True
        new_metrics = set(sm.get("series_metrics", [])) | metrics_here
        if len(new_metrics) > len(sm.get("series_metrics", [])):
            sm["metrics"] = sm.get("metrics", 0) + (
                len(new_metrics) - len(set(sm.get("series_metrics", []))))
        sm["series_metrics"] = sorted(new_metrics)

print(f"injected {len(arms)} FVS arm(s) {list(arms)} into {touched} states")

json.dump(stmeta, open(os.path.join(api, "states.json"), "w"),
          indent=1, ensure_ascii=False)
open(os.path.join(api, "states.json"), "a").write("\n")

# meta engines is a global distinct-model counter; bump only by globally-new
# models, and add the row delta. Leave states/metrics curated counters intact.
META.setdefault("stats", {})
if globally_new:
    META["stats"]["engines"] = META["stats"].get("engines", 0) + len(globally_new)
META["stats"]["rows"] = META["stats"].get("rows", 0) + delta_rows_total
json.dump(META, open(os.path.join(api, "meta.json"), "w"),
          indent=1, ensure_ascii=False)
open(os.path.join(api, "meta.json"), "a").write("\n")
print(f"bookkeeping: +{len(globally_new)} new models {sorted(globally_new)}, "
      f"+{delta_rows_total} rows; meta now engines={META['stats']['engines']}, "
      f"rows={META['stats']['rows']}")

# ---- sanity print ----
for st in ["ME", "GA", "CA", "TX", "OR"]:
    if st not in area_ha:
        continue
    bits = [f"area {area_ha[st]/1e6:.2f} Mha"]
    for cfg in arms:
        nat = arms[cfg].get(st)
        if not nat or "agc_live_total" not in nat:
            continue
        r = nat["agc_live_total"][BUCKET]
        y0, ymid = phys_total(r[2025][0], st), phys_total(r.get(2075, r[max(r)])[0], st)
        bits.append(f"{cfg} agc {y0:.0f}->{ymid:.0f}")
    print(f"  {st}: " + "  ".join(bits))
