#!/usr/bin/env python3
"""Crown closure at tree tip (cch), ported from CAL_CCH.for (ORGANON crown model).

Algorithm (faithful to the Fortran):
  CRNCLO   builds a 41-element crown-area profile CCH[0..40]; CCH[40] = tallest
           HT in the stand; CCH[i] = total crown area per acre (CCF% units,
           CA = 0.001803 * CW^2 * EXPAN) present at height stratum i*(maxHT/40),
           summed over all trees whose crown reaches that height.
  tree_cch interpolates that profile to a subject tree's tip height -> the crown
           closure experienced AT the tip (0 for the tallest tree, high for
           understory). This is the `cch` covariate in Greg's mortality model.

Per-tree crown geometry (ORGANON SWO = version 1, 18 species groups):
  MCW  maximum crown width = b0 + b1*DBH + b2*DBH^2 (capped at PKDBH; small-tree
       rule below 4.5 ft).
  LCW  largest crown width = MCW * CR^(b1 + b2*CL + b3*DBH/HT).
  HLCW height to largest crown width = HT - (1-b1)*CL.
  CW   crown width above HLCW = LCW * RP^(b1 + b2*sqrt(RP) + b3*HT/DBH),
       RP = (HT-XL)/(HT-HLCW); below HLCW the crown width is LCW.

NOTE ON SCOPE: CAL_CCH.for carries Pacific-NW ORGANON species-group parameters
(18 SWO groups). CONUS-wide application needs an FIA-species -> group crosswalk
(and ideally national crown-width coefficients); the panel's stored CCH1 was
built upstream with such a mapping. Validate this port's cch against the stored
CCH1 on real tree lists before trusting projection-level numbers; supply
SPP_TO_GROUP for the species present.
"""
from __future__ import annotations
import numpy as np

# ---- ORGANON SWO parameters (version 1), indexed by group 1..18 ----
# groups: 1 DF 2 GW 3 PP 4 SP 5 IC 6 WH 7 RC 8 PY 9 MD 10 GC 11 TA 12 CL
#         13 BL 14 WO 15 BO 16 RA 17 PD 18 WI
MCWPAR = {  # b0, b1, b2, PKDBH
 1:(4.6366,1.6078,-0.009625,88.52), 2:(6.1880,1.0069,0.0,999.99),
 3:(3.4835,1.343,-0.0082544,81.35), 4:(4.6600546,1.0701859,0.0,999.99),
 5:(3.2837,1.2031,-0.0071858,83.71), 6:(4.5652,1.4147,0.0,999.99),
 7:(4.0,1.65,0.0,999.99), 8:(4.5652,1.4147,0.0,999.99),
 9:(3.4298629,1.3532302,0.0,999.99), 10:(2.9793895,1.5512443,-0.01416129,54.77),
 11:(4.4443,1.7040,0.0,999.99), 12:(4.4443,1.7040,0.0,999.99),
 13:(4.0953,2.3849,-0.011630,102.53), 14:(3.0785639,1.9242211,0.0,999.99),
 15:(3.3625,2.0303,-0.0073307,138.93), 16:(8.0,1.53,0.0,999.99),
 17:(2.9793895,1.5512443,-0.01416129,54.77), 18:(2.9793895,1.5512443,-0.01416129,54.77)}
LCWPAR = {  # b1, b2, b3
 1:(0.0,0.00371834,0.808121), 2:(0.0,0.00308402,0.0), 3:(0.355532,0.0,0.0),
 4:(0.0,0.00339675,0.532418), 5:(-0.251389,0.00692512,0.985922), 6:(0.0,0.0,0.0),
 7:(-0.251389,0.00692512,0.985922), 8:(0.0,0.0,0.0), 9:(0.118621,0.00384872,0.0),
 10:(0.0,0.0,1.161440), 11:(0.0,0.0111972,0.0), 12:(0.0,0.0207676,0.0),
 13:(0.0,0.0,1.47018), 14:(0.364811,0.0,0.0), 15:(0.0,0.0,1.27196),
 16:(0.3227140,0.0,0.0), 17:(0.0,0.0,1.161440), 18:(0.0,0.0,1.161440)}
CWAPAR = {  # b1, b2, b3 (crown width above LCW)
 1:(0.929973,-0.135212,-0.0157579), 2:(0.999291,0.0,-0.0314603),
 3:(0.755583,0.0,0.0), 4:(0.755583,0.0,0.0), 5:(0.629785,0.0,0.0),
 6:(0.629785,0.0,0.0), 7:(0.629785,0.0,0.0), 8:(0.629785,0.0,0.0),
 9:(0.5,0.0,0.0), 10:(0.5,0.0,0.0), 11:(0.5,0.0,0.0), 12:(0.5,0.0,0.0),
 13:(0.5,0.0,0.0), 14:(0.5,0.0,0.0), 15:(0.5,0.0,0.0), 16:(0.5,0.0,0.0),
 17:(0.5,0.0,0.0), 18:(0.5,0.0,0.0)}
DACBPAR = {1:0.062,2:0.028454,3:0.05,4:0.05,5:0.20,6:0.209806,7:0.20,8:0.209806,
 9:0.0,10:0.0,11:0.0,12:0.0,13:0.0,14:0.0,15:0.0,16:0.0,17:0.0,18:0.0}


def mcw(g, D, H):
    b0, b1, b2, pk = MCWPAR[g]
    d = min(D, pk)
    return (H / 4.5 * b0) if H < 4.501 else (b0 + b1 * d + b2 * d * d)

def hlcw(g, H, CR):
    return H - (1.0 - DACBPAR[g]) * CR * H

def lcw(g, MCW, CR, DBH, HT):
    b1, b2, b3 = LCWPAR[g]; CL = CR * HT
    return MCW * CR ** (b1 + b2 * CL + b3 * (DBH / HT))

def cw_above(g, HLCW, LCW, HT, DBH, XL):
    b1, b2, b3 = CWAPAR[g]
    rp = (HT - XL) / (HT - HLCW)
    if rp <= 0: return 0.0
    return LCW * rp ** (b1 + b2 * np.sqrt(rp) + b3 * (HT / DBH))


def crown_closure(trees):
    """trees: list of dicts with group, DBH, HT, CR (0-1), EXPAN (TPA).
    Returns the 41-element CCH profile (index 0..40)."""
    cch = np.zeros(41)
    cch[40] = max(t["HT"] for t in trees)
    for t in trees:
        g, DBH, HT, CR, EXPAN = t["group"], t["DBH"], t["HT"], t["CR"], t["EXPAN"]
        CL = CR * HT; HCB = HT - CL
        MCW = mcw(g, DBH, HT); LCW = lcw(g, MCW, CR, DBH, HT); HL = hlcw(g, HT, CR)
        for ii in range(40, 0, -1):          # Fortran II=40..1
            xl = (ii - 1) * (cch[40] / 40.0)
            if xl <= (HCB if HCB > HL else HL):
                # below HLCW: crown width = LCW (with HCB>HLCW edge handled below)
                cw = LCW if HCB <= HL else cw_above(g, HL, LCW, HT, DBH, max(xl, HCB))
            elif xl < HT:
                cw = cw_above(g, HL, LCW, HT, DBH, xl)
            else:
                cw = 0.0
            cch[ii] += (cw ** 2) * (0.001803 * EXPAN)
    return cch


def tree_cch(HT, cch):
    """Interpolate the profile to a subject tree's tip height (TREE_CCH)."""
    top = cch[40]
    if HT >= top: return 0.0
    xi = 40.0 * (HT / top)
    i = int(xi) + 1           # Fortran INT(XI)+2 with 1-based -> +1 in 0-based
    if i == 40:
        return cch[39] * (40.0 - xi)
    xxi = float(i + 1) - 1.0
    return cch[i] + (cch[i - 1] - cch[i]) * (xxi - xi)


if __name__ == "__main__":
    # smoke: a 3-tree stand (DF group 1). Tall tree -> cch~0; understory -> higher.
    trees = [dict(group=1, DBH=20, HT=100, CR=0.5, EXPAN=20),
             dict(group=1, DBH=12, HT=70,  CR=0.5, EXPAN=30),
             dict(group=1, DBH=4,  HT=30,  CR=0.6, EXPAN=50)]
    cch = crown_closure(trees)
    print("profile max HT:", cch[40], " peak crown-area:", round(cch[:40].max(), 3))
    for t in trees:
        print(f"  HT={t['HT']:>3}: cch(tip)={tree_cch(t['HT'], cch):.4f}")
