#!/usr/bin/env python3
src = "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v4.R"
dst = "/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v6.R"
s = open(src).read()

marker = "## v4 CHANGES vs v3 (this file):"
v6note = """## conus_eq_projector_v6.R == v4 + HT-DBH MONOTONICITY GUARD (projection-time only).
##
## v6 CHANGE vs v4 (ONLY this; everything else identical to v4):
##   FIX D  HT-DBH MONOTONICITY GUARD. The fvs-conus Wykoff HT-DBH kernel predicts the
##          largest-DBH stems SHORTER than mid-size stems (saturating b1/(DBH+1) term plus a
##          positive a_ba*sqrt(BA) competition term), so top height (HT_M_DOM, defined on the
##          largest stems) plateaus/recedes late in the rotation while mean height stays monotone.
##          This is a PROJECTION-TIME correction, NOT a refit of the published equation. After
##          predicting per-tree HT each stand-cycle, mono_ht() enforces isotonic non-decreasing
##          HT in DBH: each tree's HT is raised to be >= the HT of any smaller-DBH tree
##          (running max over DBH-ascending order = lower monotone envelope). Two hard guards:
##            (a) NO EXPLOSION: each guarded HT is capped at min(1.2*max(raw HT in stand), HT_ABS_MAX=80 m).
##            (b) The guard only RAISES HT for the large trees; the equation is untouched for the
##                smallest tree and for any tree already on the monotone envelope, so AGB/Eichhorn
##                (which key off BA and mean HT) move negligibly.
##          mono_ht() is applied at the two HT-assignment points: the projected-cycle HT
##          (HT2 = ht_from_dbh*ht_ratio) and the recruit/NA backfill. Seed-year (year-0) HT keeps
##          the measured treeinit HT exactly as in v4 (no guard at seeding).
##
"""
s = s.replace(marker, v6note + marker, 1)
s = s.replace("## conus_eq_projector_v3.R -- annualized",
              "## conus_eq_projector_v6.R -- annualized", 1)

htdbh_close = '  1.37+exp(eta+0.5*hc("sigma")^2)}'
assert s.count(htdbh_close) == 1, "ht_from_dbh close count=%d" % s.count(htdbh_close)
mono_fn = htdbh_close + """
## v6 FIX D: isotonic HT-DBH monotonicity guard (projection-time). Raise each tree's HT to be
## >= the HT of any smaller-DBH tree (running max over DBH-ascending order). Cap at
## min(MONO_CAP_MULT * max(raw HT in stand), HT_ABS_MAX) so the guard cannot make HT explode.
MONO_CAP_MULT <- as.numeric(Sys.getenv("MONO_CAP_MULT","1.2"))
HT_ABS_MAX    <- as.numeric(Sys.getenv("HT_ABS_MAX","80"))
mono_ht <- function(dbh, ht){
  n <- length(ht); if(n <= 1) return(ht)
  ok <- is.finite(dbh) & is.finite(ht)
  if(sum(ok) <= 1) return(ht)
  cap <- min(MONO_CAP_MULT * max(ht[ok]), HT_ABS_MAX)
  o   <- order(dbh[ok])
  hh  <- ht[ok][o]
  hm  <- cummax(hh)
  hm  <- pmin(hm, cap)
  hm  <- cummax(hm)
  ht2 <- ht
  idx <- which(ok)
  ht2[idx[o]] <- hm
  ht2
}"""
s = s.replace(htdbh_close, mono_fn, 1)

old389 = "      HT2<-ht_from_dbh(tl)*tl$ht_ratio; tl$CR<-cr_update(tl,HT2); tl$HT<-HT2"
new389 = "      HT2<-ht_from_dbh(tl)*tl$ht_ratio; HT2<-mono_ht(tl$dbh_in,HT2); tl$CR<-cr_update(tl,HT2); tl$HT<-HT2  # v6 FIX D"
assert s.count(old389) == 1, "old389 count=%d" % s.count(old389)
s = s.replace(old389, new389, 1)

old436 = "      tl$HT[!is.finite(tl$HT)]<-ht_from_dbh(tl)[!is.finite(tl$HT)]"
new436 = "      tl$HT[!is.finite(tl$HT)]<-ht_from_dbh(tl)[!is.finite(tl$HT)]\n      tl$HT<-mono_ht(tl$dbh_in,tl$HT)  # v6 FIX D: keep monotone after recruit/NA backfill"
assert s.count(old436) == 1, "old436 count=%d" % s.count(old436)
s = s.replace(old436, new436, 1)

open(dst,"w").write(s)
print("WROTE", dst, "len", len(s))
print("mono_ht defs =", s.count("mono_ht <- function"), "| applied =", s.count("mono_ht(tl$dbh_in"))
