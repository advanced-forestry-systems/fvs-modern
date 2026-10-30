#!/usr/bin/env python3
"""d2_wrapper_v6.py: the runoff D2 check pattern on the fixed code (fvs-modern main d08b8b96), driven
through the catalog runner. D2_MODE=inproc runs the calibrated arm as the catalog does (FVS
config_version "calibrated": keyword block + the in-process overlay, which now skips species without a
calibrated SDImax). D2_MODE=kwonly constructs FVS with config_version "default", so only the keyword
block applies. After each run the per-species SDI maximum (spsdi) is logged to D2_SPSDI_LOG."""
import os, sys, runpy, csv
import fvs2py
MODE = os.environ.get("D2_MODE", "inproc"); LOG = os.environ.get("D2_SPSDI_LOG")
_init, _run = fvs2py.FVS.__init__, fvs2py.FVS.run
def init(self, *a, **k):
    if MODE == "kwonly": k["config_version"] = "default"
    _init(self, *a, **k)
def run(self, *a, **k):
    out = _run(self, *a, **k)
    try:
        s = self.get_species_attr("spsdi"); nz = [x for x in s if x > 0]
        if LOG:
            new = not os.path.exists(LOG)
            with open(LOG, "a", newline="") as f:
                w = csv.writer(f)
                if new: w.writerow(["mode", "config_version", "n_species", "n_spsdi_zero", "min_nonzero", "max"])
                w.writerow([MODE, self._config_version, len(s), int(sum(1 for x in s if x <= 0)), min(nz) if nz else "", max(s) if len(s) else ""])
    except Exception as e:
        print("spsdi read failed:", e, file=sys.stderr)
    return out
fvs2py.FVS.__init__, fvs2py.FVS.run = init, run
R = "/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/conus_catalog_runs.py"
sys.argv = [R] + sys.argv[1:]
runpy.run_path(R, run_name="__main__")
