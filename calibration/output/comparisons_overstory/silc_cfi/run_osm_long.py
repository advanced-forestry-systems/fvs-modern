#!/usr/bin/env python3
"""run_osm_long.py - OSM-ACD long-horizon driver (YEARS 20)"""
import logging, sqlite3, subprocess
from io import StringIO
from pathlib import Path
import numpy as np, pandas as pd
logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("osm_long")

OSM_BINARY = Path.home() / "OSM/v2.26.1/OSMv2.26.1_Linux64/OSM.ConsoleApp"
BGI_DEFAULT = 3000
BA_PER_TREE_CM = 0.00007854
FT2_AC_PER_M2_HA = 4.35
ACRES_PER_HA = 2.4710538147
N_YEARS = 20  # max horizon ~ 14 yr

OSM_REMAP = {"PB": "WB", "TA": "TL", "PC": "PR"}
acadgy_species = {"AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
                  "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM",
                  "RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP",
                  "WS","YB"}
cfi_to_agy = {"BF":"BF","RS":"RS","PB":"PB","YB":"YB","RM":"RM","SM":"SM",
              "WA":"WA","RO":"RO","BS":"BS","WS":"WS","JP":"JP","RP":"RP",
              "WP":"WP","EH":"EH","HM":"EH","CE":"WC","WC":"WC","BE":"BC",
              "QA":"QA","GB":"GB","BC":"BC","BA":"BA","AB":"AB","ST":"ST",
              "TA":"TA"}
spcd_to_agy = {"12":"BF","97":"RS","375":"PB","371":"YB","316":"RM","241":"WC",
               "261":"EH","95":"BS","91":"WS","105":"JP","129":"WP","318":"SM",
               "531":"BC","746":"QA","833":"RO","541":"WA","934":"GB"}

def cfi_to_osm(sp_raw, spcd):
    sp = cfi_to_agy.get(str(sp_raw), None)
    if sp is None:
        sp = spcd_to_agy.get(str(int(spcd)) if pd.notna(spcd) else "", None)
    if sp is None or sp not in acadgy_species:
        sp = "OH" if (pd.notna(spcd) and float(spcd) >= 300) else "OS"
    return OSM_REMAP.get(sp, sp)

def build_inputs(manifest, pair_input_dir, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = out_dir / "osm_long.sqlite"
    osmc_path = out_dir / "osm_long.osmc"
    stand_csv = out_dir / "osm_long_StandList.csv"
    tree_csv = out_dir / "osm_long_TreeList.csv"
    if sqlite_path.exists(): sqlite_path.unlink()

    sid_map = {}; stand_rows = []; tree_rows = []
    for i, p in manifest.iterrows():
        sid = i + 1
        key = (int(p["PLOT"]), int(p["YEAR_PREV"]))
        sid_map[key] = sid
        stand_rows.append((sid, 1, "ME", "None", BGI_DEFAULT, 0))
        tf = pair_input_dir / Path(p["tree_list_file"]).name
        if not tf.exists(): continue
        td = pd.read_csv(tf)
        for _, r in td.iterrows():
            d = float(r["DBH"]) if pd.notna(r["DBH"]) else 0
            if d <= 0: continue
            dbh_cm = d * 2.54
            ht = float(r["HT"]) if pd.notna(r["HT"]) and r["HT"] > 0 else 0
            if ht <= 0:
                ht_ft = max(6, 4.27 + 82 * (1 - np.exp(-0.04 * dbh_cm)))
                ht_m = ht_ft * 0.3048
            else:
                ht_m = ht * 0.3048
            expf_ha = float(r["EXPF"]) * ACRES_PER_HA
            sp = cfi_to_osm(r.get("SP", ""), r.get("SPCD"))
            tree_rows.append((sid, sp, dbh_cm, ht_m, expf_ha))

    con = sqlite3.connect(str(sqlite_path)); cur = con.cursor()
    cur.execute("CREATE TABLE OSM_StandList (SurveyID INT, Plots INT, Zone TEXT, "
                "Management TEXT, BGI INT, PoorSite SMALLINT)")
    cur.execute("CREATE TABLE OSM_TreeList (SurveyID INT, Species TEXT, DBH REAL, "
                "HT REAL, Stems REAL)")
    cur.execute("CREATE INDEX idx_tl_sid ON OSM_TreeList(SurveyID)")
    cur.executemany("INSERT INTO OSM_StandList VALUES (?, ?, ?, ?, ?, ?)", stand_rows)
    cur.executemany("INSERT INTO OSM_TreeList (SurveyID, Species, DBH, HT, Stems) VALUES (?, ?, ?, ?, ?)", tree_rows)
    con.commit(); con.close()
    lines = ["SIMULATION", f" YEARS {N_YEARS}", " YPC 1", "",
             f'INPUTS.SOURCE "{sqlite_path}"', "", "OUTPUTS", " Messages FALSE",
             " StandSummary.ConsoleOn FALSE",
             f' StandSummary.FilePath "{stand_csv}"',
             f' TreeList.FilePath "{tree_csv}"', "", "SIMULATION.Scenario BASE"]
    for sid in sorted(sid_map.values()): lines.append(f"SIMULATE {sid}")
    osmc_path.write_text("\n".join(lines) + "\n")
    return sqlite_path, osmc_path, stand_csv, tree_csv, sid_map

def parse_stand(stand_csv):
    raw = Path(stand_csv).read_text(errors="replace").replace("\x00", "")
    df = pd.read_csv(StringIO(raw))
    def col(*names):
        for n in names:
            if n in df.columns: return n
        return None
    c_sid = col("SurveyID","Survey_ID","Id")
    c_yr  = col("Year","Yr","Cycle","Age")
    c_ba  = col("BA","BasalArea","BasalArea_m2ha","GBA")
    c_tph = col("Trees","TPH","Stems","Stems_per_ha","Density")
    c_qmd = col("QMD","Qmd","Qmd_cm")
    df = df[df[c_sid].notna()].copy()
    df["SurveyID"] = df[c_sid].astype(float).astype(int)
    df["yr_off"] = df.groupby("SurveyID")[c_yr].transform(
        lambda s: (s.astype(float) - s.astype(float).min()).round().astype(int))
    df["BA_m2ha"] = df[c_ba].astype(float)
    df["TPH_ha"]  = df[c_tph].astype(float)
    df["QMD_cm"]  = (df[c_qmd].astype(float) if c_qmd
                     else np.sqrt((df["BA_m2ha"] / df["TPH_ha"].replace(0, np.nan))
                                  / BA_PER_TREE_CM))
    return df[["SurveyID","yr_off","BA_m2ha","TPH_ha","QMD_cm"]]

def main():
    wd = Path.cwd()
    manifest = pd.read_csv(wd / "silc_cfi_longhorizon_pairs.csv")
    out_dir = wd / "osm_long_work"
    sqlite_path, osmc_path, stand_csv, tree_csv, sid_map = build_inputs(
        manifest, wd / "pair_input_long", out_dir)
    log.info(f"Built inputs for {len(sid_map)} pairs")
    res = subprocess.run([str(OSM_BINARY), "Acadian", str(osmc_path)],
                         capture_output=True, text=True, timeout=600)
    if res.returncode != 0:
        log.error(f"OSM failed: {res.stderr[-500:]}")
        raise SystemExit(1)
    sp = parse_stand(stand_csv)
    log.info(f"Parsed {len(sp)} stand-year rows")
    sid_df = pd.DataFrame([(s, k[0], k[1]) for k, s in sid_map.items()],
                          columns=["SurveyID","PLOT","YEAR_PREV"])
    m = manifest[["PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR",
                  "BA_PREV_FT2AC","BA_CURR_FT2AC"]].copy()
    m["PERIOD_YR"] = m["PERIOD_YR"].clip(1, 20).astype(int)
    m = m.merge(sid_df, on=["PLOT","YEAR_PREV"], how="inner")
    sp_at = m.merge(sp, on="SurveyID", how="left")
    sp_at = sp_at[sp_at["yr_off"] == sp_at["PERIOD_YR"]].copy()
    sp_at["BA_PRED_ft2ac"] = sp_at["BA_m2ha"] * FT2_AC_PER_M2_HA
    sp_at["TPA_PRED"]      = sp_at["TPH_ha"]  / ACRES_PER_HA
    sp_at["QMD_PRED_in"]   = sp_at["QMD_cm"]  / 2.54
    sp_at["variant"] = "ACD"; sp_at["config"] = "osm"
    sp_at = sp_at.rename(columns={"BA_PREV_FT2AC":"BA_OBS_PREV",
                                   "BA_CURR_FT2AC":"BA_OBS_CURR"})
    sp_at[["PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","variant","config",
           "BA_PRED_ft2ac","TPA_PRED","QMD_PRED_in",
           "BA_OBS_PREV","BA_OBS_CURR"]].to_csv(
        wd / "silc_cfi_long_osm_results.csv", index=False)
    log.info(f"Wrote silc_cfi_long_osm_results.csv: {len(sp_at)} rows")

if __name__ == "__main__":
    main()
