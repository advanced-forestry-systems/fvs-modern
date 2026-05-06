"""
Convert Marshall-style FIA CSVs to standard FIA DataMart-format
per-state COND and TREE CSVs that fia_stand_generator.py expects.

Marshall CSVs are a curated/cleaned export with derived columns
(BAPA, QMD, SDI, RD) and missing PLT_CN. fia_stand_generator.py
needs the standard FIA DataMart schema:
  COND: STATECD, COUNTYCD, PLOT, INVYR, COND_STATUS_CD, SICOND,
        SIBASE, BALIVE, FORTYPCD, STDAGE, ASPECT, SLOPE, ELEV, PLT_CN
  TREE: PLT_CN, STATECD, COUNTYCD, PLOT, SUBP, TREE, STATUSCD,
        SPCD, DIA, HT, CR, TPA_UNADJ

PLT_CN is synthesized from (STATECD, COUNTYCD, PLOT, INVYR) so the
COND-to-TREE join works.

Usage:
  python3 marshall_to_fia_csv.py --state OR \
      --in-dir ~/fvs-modern/calibration/data/marshall_csvs \
      --out-dir ~/fia_data

Author: A. Weiskittel
Date: 2026-04-25
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)


def synth_plt_cn(df: pd.DataFrame) -> pd.Series:
    """Synthesize a unique PLT_CN from STATECD/COUNTYCD/PLOT/INVYR."""
    return (
        df["STATECD"].astype(int).astype(str).str.zfill(2)
        + df["COUNTYCD"].astype(int).astype(str).str.zfill(3)
        + df["PLOT"].astype(int).astype(str).str.zfill(6)
        + df["INVYR"].astype(int).astype(str)
    )


def convert_state(state: str, in_dir: Path, out_dir: Path) -> None:
    plot_csv = in_dir / f"{state}_PLOTdata.CSV"
    tree_csv = in_dir / f"{state}_TREEdata.CSV"
    if not plot_csv.exists() or not tree_csv.exists():
        raise FileNotFoundError(f"Missing inputs for {state}: {plot_csv}, {tree_csv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    logger.info(f"Reading {plot_csv.name}")
    plot = pd.read_csv(plot_csv, low_memory=False)
    logger.info(f"Read {len(plot)} plot rows")

    plot["PLT_CN"] = synth_plt_cn(plot)
    cond_cols_map = {
        "STATECD": "STATECD",
        "COUNTYCD": "COUNTYCD",
        "PLOT": "PLOT",
        "INVYR": "INVYR",
        "COND_STATUS_CD": "COND_STATUS_CD",
        "SICOND": "SICOND",
        "SIBASE": "SIBASE",
        "BAPA": "BALIVE",      # Marshall BAPA == FIA BALIVE
        "FORTYPCD": "FORTYPCD",
        "STDAGE": "STDAGE",
        "ASPECT": "ASPECT",
        "SLOPE": "SLOPE",
        "ELEV": "ELEV",
        "PLT_CN": "PLT_CN",
    }
    cond_present = {k: v for k, v in cond_cols_map.items() if k in plot.columns}
    cond = plot[list(cond_present.keys())].rename(columns=cond_present)
    if "BALIVE" not in cond.columns:
        logger.warning("BAPA missing; BALIVE column will be NaN-only")
    cond_out = out_dir / f"{state}_COND.csv"
    cond.to_csv(cond_out, index=False)
    logger.info(f"Wrote {cond_out} with {len(cond)} rows")

    logger.info(f"Reading {tree_csv.name}")
    tree = pd.read_csv(tree_csv, low_memory=False)
    logger.info(f"Read {len(tree)} tree rows")

    tree["PLT_CN"] = synth_plt_cn(tree)
    tree_cols_map = {
        "PLT_CN": "PLT_CN",
        "STATECD": "STATECD",
        "COUNTYCD": "COUNTYCD",
        "PLOT": "PLOT",
        "SUBP": "SUBP",
        "TREE": "TREE",
        "STATUSCD": "STATUSCD",
        "SPCD": "SPCD",
        "DIA": "DIA",
        "HT": "HT",
        "CR": "CR",
        "TPA_UNADJ": "TPA_UNADJ",
    }
    tree_present = {k: v for k, v in tree_cols_map.items() if k in tree.columns}
    tree_out_df = tree[list(tree_present.keys())].rename(columns=tree_present)
    tree_out = out_dir / f"{state}_TREE.csv"
    tree_out_df.to_csv(tree_out, index=False)
    logger.info(f"Wrote {tree_out} with {len(tree_out_df)} rows")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--state", required=True, help="USPS state abbrev like OR")
    parser.add_argument("--in-dir", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    convert_state(args.state.upper(), args.in_dir, args.out_dir)


if __name__ == "__main__":
    main()
