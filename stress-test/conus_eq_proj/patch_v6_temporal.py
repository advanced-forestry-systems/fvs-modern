f="/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v6.R"
s=open(f).read()
old="      HT2<-ht_from_dbh(tl)*tl$ht_ratio; HT2<-mono_ht(tl$dbh_in,HT2); tl$CR<-cr_update(tl,HT2); tl$HT<-HT2  # v6 FIX D"
assert s.count(old)==1, s.count(old)
new=("      HT2<-ht_from_dbh(tl)*tl$ht_ratio; HT2<-mono_ht(tl$dbh_in,HT2)  # v6 FIX D (isotonic-in-DBH)\n"
     "      ## v6 FIX D2 (temporal floor): a surviving tree cannot get SHORTER. Floor each tree's new HT\n"
     "      ## at its previous-cycle HT (tl$HT still holds last cycle's value here). This removes the\n"
     "      ## between-cycle top-height recession the a_ba*sqrt(BA) competition term induces on the\n"
     "      ## largest-DBH stems. Capped (with the isotonic cap) so HT cannot explode.\n"
     "      if(TEMPORAL_FLOOR){ prevH<-tl$HT; ok<-is.finite(prevH)&is.finite(HT2); HT2[ok]<-pmax(HT2[ok],prevH[ok])\n"
     "        HT2<-pmin(HT2, HT_ABS_MAX) }\n"
     "      tl$CR<-cr_update(tl,HT2); tl$HT<-HT2")
s=s.replace(old,new,1)
# add TEMPORAL_FLOOR env switch next to MONO_CAP_MULT def
anchor='MONO_CAP_MULT <- as.numeric(Sys.getenv("MONO_CAP_MULT","1.2"))'
assert s.count(anchor)==1
s=s.replace(anchor, anchor+'\nTEMPORAL_FLOOR <- as.logical(as.integer(Sys.getenv("TEMPORAL_FLOOR","1")))  # v6 FIX D2 default ON',1)
# doc note
s=s.replace("##          (HT2 = ht_from_dbh*ht_ratio) and the recruit/NA backfill. Seed-year (year-0) HT keeps",
 "##          (HT2 = ht_from_dbh*ht_ratio) and the recruit/NA backfill. v6 ALSO applies a per-tree\n"
 "##          TEMPORAL FLOOR (FIX D2): a surviving tree's HT is never allowed to drop below its previous-\n"
 "##          cycle HT (a living tree cannot shrink), which is what actually removes the between-cycle\n"
 "##          top-height recession; the isotonic-in-DBH guard alone only fixes within-snapshot ordering.\n"
 "##          (HT2 = ht_from_dbh*ht_ratio) and the recruit/NA backfill. Seed-year (year-0) HT keeps",1)
open(f,"w").write(s)
print("temporal floor applied:", s.count("FIX D2 (temporal floor)"), "switch:", s.count("TEMPORAL_FLOOR <- as.logical"))
