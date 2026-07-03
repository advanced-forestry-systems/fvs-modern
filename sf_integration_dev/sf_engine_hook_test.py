#!/usr/bin/env python3
"""Minimal proof that the fvs2py stop-point-5 hook works IN-PROCESS.
Uses the known-good NE test keyfile (inline TREEDATA), prepends a DSNOUT to a
temp DB exactly like the fvs2py test fixture, runs the stop-point-5 loop, and
reads per-tree engine increments. Isolates the injection mechanism from the
SQLite-input plumbing that crashed the perseus keyfile path.
"""
import os, sys, tempfile
import numpy as np
sys.path.insert(0, os.path.expanduser("~/fvs-modern/deployment/fvs2py"))
from fvs2py import FVS

LIB = os.path.expanduser(os.environ.get("FVSNE_SO", "~/fvs-modern/lib/FVSne.so"))
BASE_KEY = os.path.expanduser(
    "~/fvs-modern/deployment/fvs2py/fvs2py/tests/keyfiles/NE.key")

def main():
    with tempfile.TemporaryDirectory() as tmp:
        # NO DATABASE block: shadow logging reads tree state from engine memory
        # via get_tree_attr, so the SQLite DSNOUT (which fails in-process here)
        # is unnecessary. Use the inline-inventory keyfile as-is.
        content = open(BASE_KEY).read()
        kpath = os.path.join(tmp, "ne_hook.key")
        open(kpath, "w").write(content)
        os.chdir(tmp)  # clean CWD so FVS creates a fresh default FVSOut.db

        fvs = FVS(lib_path=LIB, config_version=None,
                  config_dir=os.path.expanduser("~/fvs-modern/config"))
        fvs.load_keyfile(kpath)

        def ga(a):
            try: return np.asarray(fvs.get_tree_attr(a), dtype=float)
            except Exception as e: print("  attr", a, "err", e); return np.array([])

        stop = 0
        fvs.run(stop_point_code=5, stop_point_year=-1)
        while getattr(fvs, "restart_code", 100) == 5:
            stop += 1
            dbh, dg, htg = ga("dbh"), ga("dg"), ga("htg")
            n = len(dbh)
            md = float(np.nanmean(dbh)) if n else float("nan")
            mdg = float(np.nanmean(dg)) if len(dg) else float("nan")
            mhg = float(np.nanmean(htg)) if len(htg) else float("nan")
            print(f"[stop5 #{stop:2d}] ntrees={n:3d} mean_dbh={md:6.2f}in "
                  f"engine_dg={mdg:7.4f}in engine_htg={mhg:6.3f}ft")
            fvs.run(stop_point_code=5, stop_point_year=-1)
        try: fvs.run()
        except Exception: pass
        print(f"RESULT: {stop} stop-point-5 hooks fired; restart_code={getattr(fvs,'restart_code',None)}")
        print("MECHANISM_OK" if stop > 0 else "MECHANISM_FAILED")

if __name__ == "__main__":
    main()
