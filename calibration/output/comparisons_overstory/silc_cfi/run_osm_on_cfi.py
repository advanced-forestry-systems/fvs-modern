#!/usr/bin/env python3
"""
run_osm_on_cfi.py
=================
Cardinal-side OSM-ACD driver for the 24 SILC CFI prediction pairs.

Pattern follows cardinal_magplot_osm_verify.py exactly: one SurveyID per
CFI plot+year_prev pair, annual snapshots (YPC 1) over the maximum
PERIOD_YR window, then each pair sliced to its own interval. BGI defaults
to 3000 (the OSM demo / median placeholder) since SILC CFI plots have
no Brunswick Growth Index attached.

Outputs:
  silc_cfi_osm_results.csv  one row per pair with OSM-ACD predicted year_curr
                            BA / TPH / QMD (metric), converted to imperial
                            for the cross-model scorecard.
"""
import logging, sqlite3, subprocess
from io import StringIO
from pathlib import Path
import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(message)s")
log = logging.getLogger("osm_cfi")

OSM_BINARY = Path.home() / "OSM/v2.26.1/OSMv2.26.1_Linux64/OSM.ConsoleApp"
BGI_DEFAULT = 3000
BA_PER_TREE_CM = 0.00007854  # m2 per (cm^2)
FT2_AC_PER_M2_HA = 4.35
ACRES_PER_HA = 2.4710538147

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


def cfi_species_to_osm(sp_raw, spcd):
    sp = cfi_to_agy.get(str(sp_raw), None)
    if sp is None:
        sp = spcd_to_agy.get(str(int(spcd)) if pd.notna(spcd) else "", None)
    if sp is None or sp not in acadgy_species:
        sp = "OH" if (pd.notna(spcd) and float(spcd) >= 300) else "OS"
    return OSM_REMAP.get(sp, sp)


def build_inputs(manifest, pair_input_dir, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    sqlite_path = out_dir / "silc_cfi_osm.sqlite"
    osmc_path = out_dir / "silc_cfi_osm.osmc"
    stand_csv = out_dir / "silc_cfi_StandListProjections.csv"
    tree_csv = out_dir / "silc_cfi_TreeListProjections.csv"
    if sqlite_path.exists():
        sqlite_path.unlink()

    sid_map = {}  # (PLOT, YEAR_PREV) -> SurveyID
    stand_rows = []
    tree_rows = []
    for i, p in manifest.iterrows():
        sid = i + 1
        key = (int(p["PLOT"]), int(p["YEAR_PREV"]))
        sid_map[key] = sid
        stand_rows.append((sid, 1, "ME", "None", BGI_DEFAULT, 0))
        tf = pair_input_dir / Path(p["tree_list_file"]).name
        if not tf.exists():
            log.warning(f"  missing tree list: {tf}")
            continue
        td = pd.read_csv(tf)
        for _, r in td.iterrows():
            dbh_in = float(r["DBH"]) if pd.notna(r["DBH"]) else 0.0
            if dbh_in <= 0:
                continue
            dbh_cm = dbh_in * 2.54
            ht_ft = float(r["HT"]) if pd.notna(r["HT"]) and r["HT"] > 0 else 0
            # impute missing heights
            if ht_ft <= 0:
                ht_ft = max(6, 4.27 + 82 * (1 - np.exp(-0.04 * dbh_cm))) / 3.28084
                # That gave m; convert to ft
                ht_ft = ht_ft * 3.28084
            ht_m = ht_ft * 0.3048
            expf_ac = float(r["EXPF"])    # trees/ac (5 for CFI)
            expf_ha = expf_ac * ACRES_PER_HA
            sp = cfi_species_to_osm(r.get("SP", ""), r.get("SPCD"))
            tree_rows.append((sid, sp, dbh_cm, ht_m, expf_ha))

    con = sqlite3.connect(str(sqlite_path))
    cur = con.cursor()
    cur.execute("CREATE TABLE OSM_StandList (SurveyID INT, Plots INT, Zone TEXT, "
                "Management TEXT, BGI INT, PoorSite SMALLINT)")
    cur.execute("CREATE TABLE OSM_TreeList (SurveyID INT, Species TEXT, DBH REAL, "
                "HT REAL, Stems REAL)")
    cur.execute("CREATE INDEX idx_tl_sid ON OSM_TreeList(SurveyID)")
    cur.executemany("INSERT INTO OSM_StandList VALUES (?, ?, ?, ?, ?, ?)", stand_rows)
    cur.executemany("INSERT INTO OSM_TreeList (SurveyID, Species, DBH, HT, Stems) "
                    "VALUES (?, ?, ?, ?, ?)", tree_rows)
    con.commit()
    con.close()
    log.info(f"Wrote {sqlite_path}: {len(stand_rows)} stands, {len(tree_rows)} tree records")

    # 10 yr horizon to cover the longest CFI period (9 yr)
    lines = ["SIMULATION", " YEARS 10", " YPC 1", "",
             f'INPUTS.SOURCE "{sqlite_path}"', "", "OUTPUTS", " Messages FALSE",
             " StandSummary.ConsoleOn FALSE",
             f' StandSummary.FilePath "{stand_csv}"',
             f' TreeList.FilePath "{tree_csv}"', "", "SIMULATION.Scenario BASE"]
    for sid in sorted(sid_map.values()):
        lines.append(f"SIMULATE {sid}")
    osmc_path.write_text("\n".join(lines) + "\n")
    return sqlite_path, osmc_path, stand_csv, tree_csv, sid_map


def run_osm(osmc, timeout=1800):
    log.info(f"Running: {OSM_BINARY} Acadian {osmc}")
    res = subprocess.run([str(OSM_BINARY), "Acadian", str(osmc)],
                         capture_output=True, text=True, timeout=timeout)
    if res.returncode != 0:
        log.error(f"OSM returned {res.returncode}")
        log.error(f"stdout tail: {res.stdout[-1500:]}")
        log.error(f"stderr tail: {res.stderr[-1500:]}")
        return False
    log.info("OSM run complete")
    return True


def parse_stand(stand_csv):
    raw = Path(stand_csv).read_text(errors="replace").replace("\x00", "")
    df = pd.read_csv(StringIO(raw))

    def col(*names):
        for n in names:
            if n in df.columns:
                return n
        return None
    c_sid = col("SurveyID", "Survey_ID", "Id")
    c_yr  = col("Year", "Yr", "Cycle", "Age")
    c_ba  = col("BA", "BasalArea", "BasalArea_m2ha", "GBA")
    c_tph = col("Trees", "TPH", "Stems", "Stems_per_ha", "Density")
    c_qmd = col("QMD", "Qmd", "Qmd_cm")
    df = df[df[c_sid].notna()].copy()
    df["SurveyID"] = df[c_sid].astype(float).astype(int)
    df["yr_off"] = df.groupby("SurveyID")[c_yr].transform(
        lambda s: (s.astype(float) - s.astype(float).min()).round().astype(int))
    df["BA_m2ha"]  = df[c_ba].astype(float)
    df["TPH_ha"]   = df[c_tph].astype(float)
    df["QMD_cm"]   = (df[c_qmd].astype(float) if c_qmd
                      else np.sqrt((df["BA_m2ha"]
                                    / df["TPH_ha"].replace(0, np.nan))
                                   / BA_PER_TREE_CM))
    return df[["SurveyID", "yr_off", "BA_m2ha", "TPH_ha", "QMD_cm"]]


def main():
    wd = Path.cwd()
    manifest = pd.read_csv(wd / "silc_cfi_pair_summary.csv")
    out_dir = wd / "osm_work"
    sqlite_path, osmc_path, stand_csv, tree_csv, sid_map = build_inputs(
        manifest, wd / "pair_input", out_dir)
    if not run_osm(osmc_path):
        raise SystemExit("OSM failed")
    sp = parse_stand(stand_csv)
    log.info(f"Parsed {len(sp)} OSM stand-year rows for {sp['SurveyID'].nunique()} stands")

    # Match each pair at its own interval
    sid_df = pd.DataFrame([(sid, k[0], k[1]) for k, sid in sid_map.items()],
                          columns=["SurveyID", "PLOT", "YEAR_PREV"])
    pair_keys = manifest[["PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR",
                          "BA_PREV_FT2AC","BA_CURR_FT2AC"]].copy()
    pair_keys["PERIOD_YR"] = pair_keys["PERIOD_YR"].clip(1, 10).astype(int)
    m = pair_keys.merge(sid_df, on=["PLOT","YEAR_PREV"], how="inner")
    sp_at = m.merge(sp, on="SurveyID", how="left")
    sp_at = sp_at[sp_at["yr_off"] == sp_at["PERIOD_YR"]].copy()
    # Convert to imperial for cross-model
    sp_at["BA_PRED_ft2ac"] = sp_at["BA_m2ha"] * FT2_AC_PER_M2_HA
    sp_at["TPA_PRED"]      = sp_at["TPH_ha"]  / ACRES_PER_HA
    sp_at["QMD_PRED_in"]   = sp_at["QMD_cm"]  / 2.54
    sp_at["variant"] = "ACD"
    sp_at["config"]  = "osm"
    cols = ["PLOT","YEAR_PREV","YEAR_CURR","PERIOD_YR","variant","config",
            "BA_PRED_ft2ac","TPA_PRED","QMD_PRED_in",
            "BA_OBS_PREV","BA_OBS_CURR"]
    sp_at = sp_at.rename(columns={"BA_PREV_FT2AC":"BA_OBS_PREV",
                                   "BA_CURR_FT2AC":"BA_OBS_CURR"})
    sp_at[cols].to_csv(wd / "silc_cfi_osm_results.csv", index=False)
    log.info(f"Wrote silc_cfi_osm_results.csv: {len(sp_at)} pairs")

    ok = sp_at["BA_PRED_ft2ac"].notna() & sp_at["BA_OBS_CURR"].notna()
    if ok.sum() > 0:
        bias_pct = 100*(sp_at.loc[ok,"BA_PRED_ft2ac"].mean() /
                        sp_at.loc[ok,"BA_OBS_CURR"].mean() - 1)
        rmse = float(np.sqrt(((sp_at.loc[ok,"BA_PRED_ft2ac"]
                              - sp_at.loc[ok,"BA_OBS_CURR"])**2).mean()))
        log.info(f"OSM-ACD on CFI: n={ok.sum()}  BA bias {bias_pct:+.2f}%  "
                 f"RMSE {rmse:.2f} ft^2/ac")

if __name__ == "__main__":
    main()
