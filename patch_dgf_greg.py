#!/usr/bin/env python3
# Wire the Greg DG substitution into ne/dgf.f90: includes, decls, lazy GREGLOADDG,
# and the per-tree WK2 override (Greg annual DG over the 10-yr cycle -> inside-bark
# DDS via BRATIO -> WK2=log(DDS), omitting native COR). Guarded by LGREGDG+GHAVE_DG.
p = "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/ne/dgf.f90"
s = open(p).read(); ch = []

a1 = "INCLUDE 'PDEN.f90'"
b1 = "INCLUDE 'PDEN.f90'\n!\nINCLUDE 'VARCOM.f90'\n!\nINCLUDE 'GREGMC.f90'"
if "GREGMC.f90" not in s and a1 in s:
    s = s.replace(a1, b1, 1); ch.append("includes")

a2 = "REAL DIAGR,DDS"
b2 = "REAL DIAGR,DDS\nREAL DGGD,GGINC,GGBARK,GGDIAGR,GGDDS\nINTEGER IGYR"
if "DGGD" not in s and a2 in s:
    s = s.replace(a2, b2, 1); ch.append("decls")

a3 = "CALL BADIST(DEBUG)"
b3 = "CALL BADIST(DEBUG)\n!  Greg DG: lazy-load coefficients/climate once (no-op unless FVS_GREGDG set)\nCALL GREGLOADDG"
if "CALL GREGLOADDG" not in s and a3 in s:
    s = s.replace(a3, b3, 1); ch.append("loadcall")

a4 = "WK2(I)=ALOG(DDS)+COR(ISPC)"
b4 = ("WK2(I)=ALOG(DDS)+COR(ISPC)\n"
      "!  --- Greg DG substitution (option 1): override native WK2 for covered species ---\n"
      "IF (LGREGDG .AND. GHAVE_DG(ISPC)) THEN\n"
      "  DGGD = DIAM(I)\n"
      "  DO IGYR = 1, 10\n"
      "    CALL GREGDGV(ISPC, DGGD, FLOAT(ICR(I))/100.0, HT(I), PTBALT(I), GELEV, GEMT, GGINC)\n"
      "    DGGD = DGGD + GGINC\n"
      "  END DO\n"
      "  GGBARK = BRATIO(ISPC, DGGD, HT(I))\n"
      "  GGDIAGR = (DGGD - DIAM(I)) * GGBARK\n"
      "  GGDDS = GGDIAGR * (2.0*DIAM(I)*GGBARK + GGDIAGR)\n"
      "  IF (GGDDS .LT. 1.0E-6) GGDDS = 1.0E-6\n"
      "  WK2(I) = ALOG(GGDDS)\n"
      "ENDIF")
if "Greg DG substitution" not in s and a4 in s:
    s = s.replace(a4, b4, 1); ch.append("hook")

open(p, "w").write(s)
print("changed:", ch, "| verified:",
      all(x in s for x in ["GREGMC.f90","DGGD","CALL GREGLOADDG","Greg DG substitution"]))
