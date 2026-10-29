"""fvs2py.catalog  (2026-10-28)

Configuration catalog for fvs-modern: which survival model and which parameter set to run for a
given FVS variant, with the long-run realism evidence behind the choice.

The catalog is data, not code: `config/catalog/catalog_recommendation_by_variant.csv` and
`config/catalog/catalog_gates_by_arm.csv`, produced by validation/bakuzis_stress_v2/catalog/
(conus_catalog_runs.py on Cardinal, bakuzis_stress_v2.R, catalog_gates.R). Every number those
files carry is traceable to the run folder named in the README beside them.

Arms are named <parameter set>_<survival>:
  parameter set:  default      fvs-modern default (variant native coefficients, no calibration keywords)
                  calibrated   fvs-modern Bayesian calibrated MAP (config/calibrated/<variant>.json)
                  conus_b1, conus_b2   CONUS trait-mediated projector arms (fvs-conus-pipeline; not run
                               through this engine; listed for comparison only)
  survival:       native       the variant's own mortality routine
                  gompit       GOMPIT crown-closure survival (FVS_GOMPIT=1, FVS_GOMPIT_COEF=<csv>)

Usage:
    from fvs2py.catalog import select_config, catalog_table
    choice = select_config("ne")                      # the catalog recommendation
    choice = select_config("sn", survival="gompit")   # force the survival model, keep the recommended set
    fvs = FVS(lib_path=..., config_version=choice.config_version, config_dir=...)
    os.environ.update(choice.env)                     # before the first run(); GOMPIT is env gated

`select_config` never silently returns an arm the catalog did not admit: when the variant is
NOT RECOMMENDED for long horizons, `choice.admissible` is False and `choice.note` says why, and
the caller decides.
"""
from __future__ import annotations
import csv, os
from dataclasses import dataclass, field

_HERE = os.path.dirname(os.path.abspath(__file__))
_DEFAULT_DIR = os.path.normpath(os.path.join(_HERE, "..", "..", "..", "config", "catalog"))
GOMPIT_COEF_DEFAULT = "greg_mortality_coefficients.csv"


@dataclass
class CatalogChoice:
    variant: str
    arm: str
    parameter_set: str
    survival: str
    config_version: str
    env: dict = field(default_factory=dict)
    admissible: bool = False
    status: str = ""
    admissible_arms: list = field(default_factory=list)
    nearterm_evidence: str = ""
    evidence: dict = field(default_factory=dict)
    note: str = ""


def _read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def catalog_table(catalog_dir: str | None = None):
    """All arm by variant gate rows (list of dicts)."""
    d = catalog_dir or os.environ.get("FVS_CATALOG_DIR", _DEFAULT_DIR)
    return _read(os.path.join(d, "catalog_gates_by_arm.csv"))


def _split(arm: str):
    surv = "gompit" if arm.endswith("gompit") else "native"
    pset = arm[: -(len(surv) + 1)] if arm.endswith("_" + surv) else arm.replace("_gompit", "")
    return pset, surv


def select_config(variant: str, survival: str | None = None, parameter_set: str | None = None,
                  catalog_dir: str | None = None, config_dir: str | None = None) -> CatalogChoice:
    v = variant.upper()
    d = catalog_dir or os.environ.get("FVS_CATALOG_DIR", _DEFAULT_DIR)
    recs = {r["group"]: r for r in _read(os.path.join(d, "catalog_recommendation_by_variant.csv"))}
    gates = [r for r in _read(os.path.join(d, "catalog_gates_by_arm.csv")) if r["group"] == v]
    if v not in recs:
        raise KeyError(f"variant {v} is not in the catalog (no fitting plots assigned or not run)")
    rec = recs[v]
    if not rec["recommended"]:
        raise KeyError(f"variant {v}: {rec['status']}")
    adm = [a for a in rec["admissible_arms"].split("; ") if a]
    arm = rec["recommended"]
    prefer = "gompit" if "NE cohort" in rec.get("nearterm_evidence", "") else None   # the only per-variant near-term evidence
    # caller constraints: pick the simplest admissible arm that satisfies them, else the matching
    # non-admissible arm with admissible=False
    cands = [a for a in adm if (survival is None or a.endswith(survival)) and (parameter_set is None or a.startswith(parameter_set))]
    if survival is not None or parameter_set is not None:
        if cands:
            pref = [a for a in cands if prefer and a.endswith(prefer)]
            arm = (pref or cands)[0]
        else:
            arm = f"{parameter_set or _split(rec['recommended'])[0]}_{survival or _split(rec['recommended'])[1]}"
    pset, surv = _split(arm)
    if pset.startswith("conus"):
        raise ValueError(f"{arm} is a CONUS projector arm; it is not runnable through fvs2py. Use fvs-conus-pipeline.")
    env = {}
    if surv == "gompit":
        cdir = config_dir or os.environ.get("FVS_CONFIG_DIR", os.path.normpath(os.path.join(_HERE, "..", "..", "..", "config")))
        env = {"FVS_GOMPIT": "1", "FVS_GOMPIT_COEF": os.path.join(cdir, GOMPIT_COEF_DEFAULT)}
    ev = next((g for g in gates if g["arm"] == arm), {})
    admissible = arm in adm
    note = "" if admissible else (f"{arm} is not an admissible arm for {v}: {ev.get('why_not', '')}".strip() if ev else f"{arm} has no gate row for {v}")
    if ev.get("reineke_unknown", "") in ("TRUE", "True"):
        note = (note + " " if note else "") + f"Reineke self-thinning not testable ({ev.get('reineke_eval', '0')} trajectories reached the relative-density gate)."
    if ev.get("calibrated_distinct", "") in ("FALSE", "False"):
        note = (note + " " if note else "") + "The calibrated config did not load for this variant; this arm is identical to the default."
    if rec.get("rule_sensitive", "") in ("TRUE", "True"):
        note = (note + " " if note else "") + "Verdict is rule-sensitive (admissible set differs under the alternative Reineke gate)."
    if not rec["status"].startswith("RECOMMENDED"):
        note = (note + " " if note else "") + rec["status"]
    return CatalogChoice(variant=v, arm=arm, parameter_set=pset, survival=surv,
                         config_version="calibrated" if pset == "calibrated" else "default", env=env,
                         admissible=admissible, status=rec["status"], admissible_arms=adm,
                         nearterm_evidence=rec.get("nearterm_evidence", ""), evidence=ev, note=note)


if __name__ == "__main__":
    import sys
    for v in sys.argv[1:] or ["ne"]:
        c = select_config(v)
        print(v.upper(), c.arm, c.config_version, c.env, "admissible" if c.admissible else "NOT admissible", c.status, c.note)
