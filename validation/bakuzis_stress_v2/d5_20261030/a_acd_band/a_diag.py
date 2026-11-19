"""diagA: why 81 of 180 ACD ensemble stands give identical year-100 BA across 100 posterior draws.
Three of the 81 stands (one per species group), MAP + 3 draws + 4 probes (SDIMAX 150 / 600 all species,
BAIMULT 1.5, MORTMULT 2.0 all species, on top of MAP), keyfile echo, FVS_TreeList at cycles 1-2.
No identifiers written (stands named by scenario/replicate)."""
import os, sys, glob, sqlite3, contextlib, logging, warnings
import numpy as np, pandas as pd
warnings.filterwarnings("ignore")
W = "/fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030"; B = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027"; OUT = B + "/diagA/runs"
for p in (W, W + "/deployment/fvs2py", W + "/calibration/python"): sys.path.insert(0, p)
os.makedirs(OUT, exist_ok=True)
import bakuzis_uncertainty_comparison as buc; logging.disable(logging.CRITICAL)
from fia_stand_generator import generate_real_stand
from config.uncertainty import UncertaintyEngine
from config.config_loader import FvsConfigLoader
state = {"dir": None}
class Keep:
    def __init__(self, *a, **k): pass
    def __enter__(self): os.makedirs(state["dir"], exist_ok=True); return state["dir"]
    def __exit__(self, *a): return False
buc.tempfile.TemporaryDirectory = Keep
s = pd.read_csv(B + "/diagA/yr100_draw_spread_by_stand.csv"); s = s[(s.variant == "ACD") & (s.uBA == 1)]
scen = {num: (sp, site, dens) for (sp, site, dens, num) in buc.enumerate_scenarios(buc.get_species_groups("acd"))}
s["sp"] = s.scenario.map(lambda n: scen[n][0]); picks = s.groupby("sp").head(1).head(3)
eng = UncertaintyEngine("acd", config_dir=buc.CONFIG_DIR, seed=42)
draws = list(np.random.default_rng(42).choice(eng.n_draws, size=100, replace=False))[:3]
mapkw = FvsConfigLoader("acd", version="calibrated", config_dir=buc.CONFIG_DIR).generate_keywords(include_comments=False)
probes = {"sdimax150": f"{'SDIMAX':<10}{0:10d}{150.0:10.1f}", "sdimax600": f"{'SDIMAX':<10}{0:10d}{600.0:10.1f}",
          "baimult1.5": f"{'BAIMULT':<10}{'':>10}{0:10d}{1.5:10.2f}", "mortmult2": f"{'MORTMULT':<10}{'':>10}{0:10d}{2.0:10.2f}"}
rows, trees, echo = [], [], []
for _, r in picks.iterrows():
    sp, site, dens = scen[int(r.scenario)]
    st = generate_real_stand("acd", site, dens, n_plots=5, seed=42 + int(r.scenario))
    sdf, tdf, _ = st[int(r.replicate) - 1]; sid = f"BK{int(r.scenario):02d}R{int(r.replicate)}"
    sdf = sdf.copy(); tdf = tdf.copy(); sdf["stand_id"] = sid; tdf["stand_id"] = sid
    runs = [("map", "calibrated", None)] + [(f"draw{int(d)}", "posterior", eng.generate_keywords_for_draw(eng.get_draw(int(d)), None, draw_idx=int(d))) for d in draws] \
         + [(f"probe_{k}", "posterior", mapkw + "\n" + v) for k, v in probes.items()]
    for lab, cv, kw in runs:
        state["dir"] = f"{OUT}/{sid}_{lab}"
        summ = buc.run_scenario(sdf, tdf, sid, "acd", cv, draw_keywords=kw)
        summ = summ.assign(stand=sid, run=lab, species_group=sp); rows.append(summ[["stand", "species_group", "run", "year", "atba", "tpa", "attopht"]])
        try:
            with contextlib.closing(sqlite3.connect(f"{state['dir']}/FVS_Data.db")) as con:
                tl = pd.read_sql("select * from FVS_TreeList", con)
            yrs = sorted(tl["Year"].unique())[:3]; tl = tl[tl["Year"].isin(yrs)]
            trees.append(tl.groupby(["Year", "SpeciesFIA"]).apply(lambda g: pd.Series(dict(tpa=g.TPA.sum(), dg_mean=np.average(g.DG, weights=g.TPA) if g.TPA.sum() > 0 else np.nan))).reset_index().assign(stand=sid, run=lab))
        except Exception as e:
            print("treelist read failed", sid, lab, e)
        outs = glob.glob(state["dir"] + "/*.out")
        txt = open(outs[0], errors="ignore").read() if outs else ""
        echo.append(dict(stand=sid, run=lab, out_found=bool(outs), n_baimult=txt.count("BAIMULT"), n_mortmult=txt.count("MORTMULT"), n_sdimax=txt.count("SDIMAX"),
                         n_bamax=txt.count("BAMAX"), kw_records=0 if kw is None else sum(1 for l in kw.splitlines() if l[:8].strip() in ("SDIMAX", "BAIMULT", "MORTMULT", "BAMAX"))))
R = pd.concat(rows); E = pd.DataFrame(echo)
R.to_csv(B + "/diagA/a_runs_summary.csv", index=False); E.to_csv(B + "/diagA/a_keyword_echo.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500); pd.set_option("display.max_columns", 30)
print(E.to_string())
w = R.pivot_table(index=["stand", "year"], columns="run", values="atba").reset_index(); print(w.to_string())
t = R.pivot_table(index=["stand", "year"], columns="run", values="tpa").reset_index(); print(t[t.year.isin([2000, 2005, 2010, 2050, 2100])].to_string())
if trees:
    T = pd.concat(trees); T.to_csv(B + "/diagA/a_treelists_cycles.csv", index=False)
    g = T.groupby(["stand", "run", "Year"]).agg(tpa=("tpa", "sum")).reset_index().pivot_table(index=["stand", "Year"], columns="run", values="tpa"); print(g.to_string())
print("A_DIAG_DONE")
