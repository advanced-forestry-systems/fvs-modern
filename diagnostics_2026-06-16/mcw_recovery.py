#!/usr/bin/env python3
"""Recover Maximum Crown Width (MCW) from each variant's FVS CCF coefficients (Marshall, 2026-06-17).

FVS stores per-tree CCF as CCFT = RD1 + RD2*D + RD3*D^2 (D in inches, D>=1), which under the open-grown
identity CCF = k*MCW^2 (k = 0.001803026, MCW in feet) implies a linear MCW = B0 + B1*D with
  B0 = sqrt(RD1/k),  B1 = sqrt(RD3/k).
(Power-form variants: A1 = sqrt(R4/k), A2 = R5/2; not present in the quadratic-form variant files.)

This script parses every <variant>/ccfcal.f90, recovers B0/B1 and MCW at a reference DBH per species,
and tabulates the cross-variant spread per species abbreviation (e.g. Douglas-fir across variants).
Outputs: mcw_by_variant_species.csv (per row) and mcw_cross_variant_spread.csv (per species).

Usage: python3 mcw_recovery.py <SRC_DIR> <OUT_DIR>   (SRC_DIR = src-converted/)
"""
import os, re, sys, math, csv, glob

K = 0.001803026          # CCF-to-crown-area constant (Marshall)
REF_DBH = 20.0           # inches, reference DBH for a single comparable MCW value
SRC = sys.argv[1] if len(sys.argv) > 1 else "src-converted"
OUT = sys.argv[2] if len(sys.argv) > 2 else "."

def parse_data_array(text, name):
    """Return the list of numeric tokens in a Fortran `DATA <name>/ ... /` block, or None."""
    m = re.search(r"DATA\s+" + re.escape(name) + r"\s*/(.*?)/", text, re.S | re.I)
    if not m:
        return None
    body = re.sub(r"&", " ", m.group(1))          # line continuations
    body = re.sub(r"!.*", "", body)               # strip trailing comments
    toks = re.findall(r"[-+]?\d*\.?\d+(?:[EeDd][-+]?\d+)?", body)
    out = []
    for t in toks:
        t = t.replace("D", "E").replace("d", "e")
        try:
            out.append(float(t))
        except ValueError:
            pass
    return out

def parse_species_order(text):
    """Parse the 'VARIANT SPECIES ORDER' header comment into {index: abbrev}. Tolerant of layout."""
    order = {}
    m = re.search(r"VARIANT SPECIES ORDER(.*?)(?:CCF EQUATIONS ORDER|SOURCES OF COEFF|----------)", text, re.S | re.I)
    block = m.group(1) if m else text[:1500]
    for idx, ab in re.findall(r"(\d+)\s*=\s*([A-Za-z][A-Za-z0-9 ]{0,3}?)(?=[,\n])", block):
        ab = ab.strip()
        if ab and ab not in (",",):
            order[int(idx)] = ab
    return order

rows = []
flagged = []
variants = sorted(d for d in os.listdir(SRC) if os.path.isfile(os.path.join(SRC, d, "ccfcal.f90")))
for var in variants:
    text = open(os.path.join(SRC, var, "ccfcal.f90"), errors="ignore").read()
    indccf = parse_data_array(text, "INDCCF")
    rd1 = parse_data_array(text, "RD1")
    rd2 = parse_data_array(text, "RD2")
    rd3 = parse_data_array(text, "RD3")
    rda = parse_data_array(text, "RDA")   # power-form coefficient (CCF = RDA*D^RDB)
    rdb = parse_data_array(text, "RDB")   # power-form exponent
    b1a = parse_data_array(text, "B1")    # direct-MCW power coef (MCW = B1*D^B2), e.g. AK
    b2a = parse_data_array(text, "B2")
    eqmap = parse_data_array(text, "EQMAP")
    spo = parse_species_order(text)
    # layouts: INDCCF-mapped CCF, direct-indexed CCF, or direct-MCW power form (AK: MCW=B1*D^B2 via EQMAP)
    if indccf and rd1 and rd3:
        nsp = len(indccf); mapped = True; direct_mcw = False
    elif rd1 and rd3:
        nsp = len(rd1); mapped = False; direct_mcw = False
    elif rda and rdb:
        nsp = len(rda); mapped = False; direct_mcw = False
    elif eqmap and b1a and b2a:
        nsp = len(eqmap); mapped = False; direct_mcw = True
    else:
        flagged.append((var, "CCF delegated to crown-width routine (R5CRWD/CWCALC); MCW lives there - read separately"))
        continue
    for i in range(nsp):
        if direct_mcw:
            eq = int(eqmap[i])
            if eq < 1 or eq > len(b1a):
                continue
            a1, a2 = b1a[eq - 1], b2a[eq - 1]
            if a1 <= 0:
                continue
            rows.append({
                "variant": var, "sp_index": i + 1, "sp_abbrev": spo.get(i + 1, ""),
                "eq_index": eq, "form": "power_mcw", "RD1": float("nan"), "RD2": float("nan"), "RD3": float("nan"),
                "B0_or_A1": round(a1, 4), "B1_or_A2": round(a2, 5),
                "MCW_at_20in_ft": round(a1 * REF_DBH ** a2, 3), "RD2_consistency": float("nan")})
            continue
        ic = int(indccf[i]) if mapped else i + 1
        if ic < 1:
            continue
        r1 = rd1[ic - 1] if rd1 and ic - 1 < len(rd1) else 0.0
        r2 = rd2[ic - 1] if rd2 and ic - 1 < len(rd2) else float("nan")
        r3 = rd3[ic - 1] if rd3 and ic - 1 < len(rd3) else 0.0
        form = ""; b0 = b1 = mcw20 = float("nan"); cons = float("nan")
        if r1 > 0 and r3 > 0:                                   # quadratic large-tree form -> linear MCW
            b0 = math.sqrt(r1 / K); b1 = math.sqrt(r3 / K); mcw20 = b0 + b1 * REF_DBH; form = "linear"
            ir2 = 2.0 * math.sqrt(r1 * r3)
            cons = r2 / ir2 if ir2 > 0 and r2 == r2 else float("nan")
        elif rda and rdb and i < len(rda) and i < len(rdb) and rda[i] > 0:   # power form -> MCW = A1*D^A2
            a1 = math.sqrt(rda[i] / K); a2 = rdb[i] / 2.0; mcw20 = a1 * REF_DBH ** a2
            b0, b1, form = a1, a2, "power"   # B0/B1 columns reused as A1/A2 for power rows
        else:
            continue
        rows.append({
            "variant": var, "sp_index": i + 1, "sp_abbrev": spo.get(i + 1, ""),
            "eq_index": ic, "form": form, "RD1": round(r1, 6), "RD2": r2, "RD3": round(r3, 6),
            "B0_or_A1": round(b0, 4), "B1_or_A2": round(b1, 5),
            "MCW_at_20in_ft": round(mcw20, 3),
            "RD2_consistency": round(cons, 3) if cons == cons else float("nan"),
        })

# write per-row table
os.makedirs(OUT, exist_ok=True)
f1 = os.path.join(OUT, "mcw_by_variant_species.csv")
with open(f1, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

# cross-variant spread by species abbreviation
by_ab = {}
for r in rows:
    ab = r["sp_abbrev"]
    if not ab:
        continue
    by_ab.setdefault(ab, []).append(r)
spread = []
for ab, rs in sorted(by_ab.items()):
    mcw = [r["MCW_at_20in_ft"] for r in rs]
    vars_ = sorted({r["variant"] for r in rs})
    spread.append({
        "sp_abbrev": ab, "n_variants": len(vars_), "n_curves": len(rs),
        "MCW20_min": round(min(mcw), 2), "MCW20_max": round(max(mcw), 2),
        "MCW20_range": round(max(mcw) - min(mcw), 2),
        "MCW20_mean": round(sum(mcw) / len(mcw), 2),
        "variants": ";".join(vars_),
    })
spread.sort(key=lambda x: -x["MCW20_range"])
f2 = os.path.join(OUT, "mcw_cross_variant_spread.csv")
with open(f2, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(spread[0].keys()))
    w.writeheader(); w.writerows(spread)

print("variants parsed (quadratic ccfcal.f90):", len(variants) - len(flagged), "of", len(variants))
print("rows (variant x species):", len(rows))
print("species abbrevs with >1 curve across variants:",
      sum(1 for s in spread if s["n_curves"] > 1))
print("\nTop 12 widest cross-variant MCW spread (at 20in DBH):")
print("%-6s %4s %4s %8s %8s %8s   %s" % ("spp", "nvar", "ncrv", "min", "max", "range", "variants"))
for s in spread[:12]:
    print("%-6s %4d %4d %8.2f %8.2f %8.2f   %s" %
          (s["sp_abbrev"], s["n_variants"], s["n_curves"], s["MCW20_min"], s["MCW20_max"], s["MCW20_range"], s["variants"]))
if flagged:
    print("\nFLAGGED (need manual handling):")
    for v, why in flagged:
        print("  %-5s %s" % (v, why))
print("\nwrote", f1, "and", f2)
