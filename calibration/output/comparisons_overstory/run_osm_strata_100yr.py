#!/usr/bin/env python3
"""run_osm_strata_100yr.py - OSM-ACD 100-yr driver for SILC byStrata
stands. Uses GrownDB year-2023 snapshot for the starting tree list
(matches AGM/FVS-NE), with refined per-stand BGI=3902 from ME_BGI_V1
raster at Davistown lat/long and CSI from StandInit."""
import sqlite3, subprocess
from io import StringIO
from pathlib import Path
import numpy as np, pandas as pd

OSM_BINARY = Path.home() / "OSM/v2.26.1/OSMv2.26.1_Linux64/OSM.ConsoleApp"
BGI_REFINED = 3902  # from ME_BGI_V1.tif sampled at lat 46.4628, lon -68.4253
BA_PER_TREE_CM = 0.00007854
FT2_AC_PER_M2_HA = 4.35
ACRES_PER_HA = 2.4710538147
N_YEARS = 100

OSM_REMAP = {"PB": "WB", "TA": "TL", "PC": "PR"}
acadgy_species = {"AB","AS","BA","BC","BF","BP","BS","BT","EC","EH","GA","GB",
                  "HH","HW","JP","NS","OH","OS","PB","PC","PR","QA","RB","RM",
                  "RN","RO","RP","RS","SB","SM","ST","SW","TA","WA","WC","WP",
                  "WS","YB"}

def to_osm_sp(sp):
    s = str(sp)
    if s in acadgy_species:
        return OSM_REMAP.get(s, s)
    return "OS"

def build_inputs(stand_ids, gr2023, si_lookup, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = out_dir / "osm_strata.sqlite"
    osmc_path   = out_dir / "osm_strata.osmc"
    stand_csv   = out_dir / "osm_strata_StandList.csv"
    tree_csv    = out_dir / "osm_strata_TreeList.csv"
    if sqlite_path.exists(): sqlite_path.unlink()

    sid_map = {}
    stand_rows = []
    tree_rows  = []
    for i, sid in enumerate(stand_ids):
        n = i + 1
        sid_map[sid] = n
        stand_rows.append((n, 1, "ME", "None", BGI_REFINED, 0))
        td = gr2023[gr2023["StandID"] == sid]
        for _, r in td.iterrows():
            dbh = float(r["DBH"]) if pd.notna(r["DBH"]) else 0
            # OSM rejects trees below the ingrowth threshold; drop seedlings (<0.5 in)
            # and unrealistic upper outliers (>40 in is likely data entry error)
            if dbh < 0.5 or dbh > 40: continue
            dbh_cm = dbh * 2.54
            ht = float(r["Ht"]) if pd.notna(r["Ht"]) and r["Ht"] > 0 else 0
            if ht <= 0:
                ht_ft = max(6, 4.27 + 82 * (1 - np.exp(-0.04 * dbh_cm)))
                ht_m = ht_ft * 0.3048
            else:
                ht_m = ht * 0.3048
            expf_ha = float(r["TPA"]) * ACRES_PER_HA
            sp = to_osm_sp(str(r["Species"]))
            tree_rows.append((n, sp, dbh_cm, ht_m, expf_ha))

    con = sqlite3.connect(str(sqlite_path)); cur = con.cursor()
    cur.execute("CREATE TABLE OSM_StandList (SurveyID INT, Plots INT, Zone TEXT, "
                "Management TEXT, BGI INT, PoorSite SMALLINT)")
    cur.execute("CREATE TABLE OSM_TreeList (SurveyID INT, Species TEXT, DBH REAL, "
                "HT REAL, Stems REAL)")
    cur.execute("CREATE INDEX idx_tl_sid ON OSM_TreeList(SurveyID)")
    cur.executemany("INSERT INTO OSM_StandList VALUES (?,?,?,?,?,?)", stand_rows)
    cur.executemany("INSERT INTO OSM_TreeList (SurveyID, Species, DBH, HT, Stems) VALUES (?,?,?,?,?)", tree_rows)
    con.commit(); con.close()

    lines = ["SIMULATION", f" YEARS {N_YEARS}", " YPC 5", "",
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
    c_sid = col("SurveyID","Survey_ID","Id")
    c_yr  = col("Year","Yr","Cycle","Age")
    c_ba  = col("BA","BasalArea","BasalArea_m2ha","GBA")
    c_tph = col("Trees","TPH","Stems","Stems_per_ha","Density")
    c_qmd = col("QMD","Qmd","Qmd_cm")
    df = df[df[c_sid].notna()].copy()
    df["SurveyID"] = df[c_sid].astype(float).astype(int)
    df["yr_off"]   = df.groupby("SurveyID")[c_yr].transform(
        lambda s: (s.astype(float) - s.astype(float).min()).round().astype(int))
    df["BA_m2ha"]  = df[c_ba].astype(float)
    df["TPH_ha"]   = df[c_tph].astype(float)
    df["QMD_cm"]   = (df[c_qmd].astype(float) if c_qmd
                     else np.sqrt((df["BA_m2ha"] / df["TPH_ha"].replace(0, np.nan)) / BA_PER_TREE_CM))
    return df[["SurveyID","yr_off","BA_m2ha","TPH_ha","QMD_cm"]]

def main():
    wd = Path.cwd()
    gr = pd.read_csv(wd / "GrownDB_byStrata_ALL.csv")
    gr2023 = gr[(gr["Year"] == 2023) & gr["DBH"].notna() & gr["TPA"].notna() & (gr["TPA"] > 0)]
    si = pd.read_csv(wd / "Acadian_Matrix_StandInit_2023.csv")
    si_lookup = dict(zip(si["STAND_ID"], si["ClimateSiteIndexMeters"]))
    stand_ids = sorted(gr2023["StandID"].unique())
    print(f"Running OSM-ACD 100-yr on {len(stand_ids)} byStrata stands, BGI={BGI_REFINED}")

    out_dir = wd / "osm_strata_work"
    sqlite_path, osmc_path, stand_csv, tree_csv, sid_map = build_inputs(
        stand_ids, gr2023, si_lookup, out_dir)
    print(f"Built inputs for {len(sid_map)} stands")
    res = subprocess.run([str(OSM_BINARY), "Acadian", str(osmc_path)],
                         capture_output=True, text=True, timeout=1800)
    if res.returncode != 0:
        print(f"OSM failed: {res.stderr[-500:]}")
        raise SystemExit(1)
    sp = parse_stand(stand_csv)
    print(f"Parsed {len(sp)} stand-year rows")

    # Convert to imperial and map back to StandID
    sid_df = pd.DataFrame([(s, sid) for sid, s in sid_map.items()],
                           columns=["StandID","SurveyID"])
    sp = sp.merge(sid_df, on="SurveyID", how="left")
    sp["Year"]     = 2023 + sp["yr_off"]
    sp["BA"]       = sp["BA_m2ha"]   * FT2_AC_PER_M2_HA
    sp["Tpa"]      = sp["TPH_ha"]    / ACRES_PER_HA
    sp["QMD"]      = sp["QMD_cm"]    / 2.54
    sp = sp[["StandID","Year","BA","Tpa","QMD"]]
    sp.to_csv(wd / "silc_strata_100yr_osmacd_results.csv", index=False)
    print(f"Wrote silc_strata_100yr_osmacd_results.csv: {len(sp)} rows")

if __name__ == "__main__":
    main()
