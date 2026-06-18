#!/usr/bin/env python3
"""Select a CONUS-consistent MCW per species from the recovered per-variant curves (2026-06-18).

Reads mcw_by_variant_species.csv (per-variant MCW = B0 + B1*D recovered from CCF coefficients), groups by
species, and chooses a consensus curve: the cross-variant median B0 and B1 (linear form), with the median
MCW at 20 in DBH and the spread. This is the unified crown-width equation per species to feed into the
crown / CCF competition term of both the calibrated engine and the fvs-conus crown component. For the
Acadian and Northeast region, Russell and Weiskittel (2010) is the preferred source (those variants do not
carry a parabolic MCW; see the recovery note). Output: mcw_conus_consensus.csv.
"""
import csv, statistics as st
IN="mcw_out/mcw_by_variant_species.csv"; OUT="mcw_out/mcw_conus_consensus.csv"
rows=[r for r in csv.DictReader(open(IN))]
by={}
for r in rows:
    ab=r["sp_abbrev"].strip()
    if not ab: continue
    by.setdefault(ab,[]).append(r)
def med(xs):
    xs=[float(x) for x in xs if x not in ("","nan")]
    return round(st.median(xs),4) if xs else float("nan")
out=[]
for ab,rs in sorted(by.items()):
    lin=[r for r in rs if r.get("form")=="linear"]
    mcw20=[float(r["MCW_at_20in_ft"]) for r in rs if r["MCW_at_20in_ft"] not in ("","nan")]
    src=lin if lin else rs
    rec={
        "sp_abbrev":ab,
        "n_variants":len({r["variant"] for r in rs}),
        "n_curves":len(rs),
        "form":"linear" if lin else (rs[0].get("form") or ""),
        "B0_consensus_ft":med([r["B0_or_A1"] for r in src]),
        "B1_consensus_ft_per_in":med([r["B1_or_A2"] for r in src]),
        "MCW20_consensus_ft":round(st.median(mcw20),2) if mcw20 else float("nan"),
        "MCW20_min":round(min(mcw20),2) if mcw20 else float("nan"),
        "MCW20_max":round(max(mcw20),2) if mcw20 else float("nan"),
        "MCW20_spread":round(max(mcw20)-min(mcw20),2) if mcw20 else float("nan"),
        "variants":";".join(sorted({r["variant"] for r in rs})),
    }
    out.append(rec)
out.sort(key=lambda x:-x["n_variants"])
with open(OUT,"w",newline="") as fh:
    w=csv.DictWriter(fh,fieldnames=list(out[0].keys())); w.writeheader(); w.writerows(out)
multi=[r for r in out if r["n_variants"]>1]
print("species with a consensus curve:",len(out),"; multi-variant species:",len(multi))
print("consensus MCW = B0 + B1*D (ft, D in inches); selected as cross-variant median")
print("%-6s %4s %8s %9s %9s %8s"%("spp","nvar","B0","B1","MCW@20","spread"))
for r in sorted(multi,key=lambda x:-x["n_variants"])[:18]:
    print("%-6s %4d %8.2f %9.3f %9.2f %8.2f"%(r["sp_abbrev"],r["n_variants"],r["B0_consensus_ft"],r["B1_consensus_ft_per_in"],r["MCW20_consensus_ft"],r["MCW20_spread"]))
print("wrote",OUT)
