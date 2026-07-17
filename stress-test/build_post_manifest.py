#!/usr/bin/env python3
# build_post_manifest.py - extend the FVS posterior-uncertainty manifest to all CONUS
# states. For each state (FIPS) finds its dominant FVS variant (most stands across the
# standinit_by_variant files), then emits manifest rows (idx, variant_lower, abbrev, draw)
# for the states not already run, draws 0..NDRAW-1. Appends to manifest_post.tsv numbering.
import glob, os, csv, collections
SD="/fs/scratch/PUOM0008/crsfaaron/fvs_stress"
SI=os.path.join(SD,"standinit_by_variant")
NDRAW=24
DONE={"GA","ID","IN","ME","MN","OR","WA"}
FIPS={1:"AL",4:"AZ",5:"AR",6:"CA",8:"CO",9:"CT",10:"DE",12:"FL",13:"GA",16:"ID",17:"IL",18:"IN",
 19:"IA",20:"KS",21:"KY",22:"LA",23:"ME",24:"MD",25:"MA",26:"MI",27:"MN",28:"MS",29:"MO",30:"MT",
 31:"NE",32:"NV",33:"NH",34:"NJ",35:"NM",36:"NY",37:"NC",38:"ND",39:"OH",40:"OK",41:"OR",42:"PA",
 44:"RI",45:"SC",46:"SD",47:"TN",48:"TX",49:"UT",50:"VT",51:"VA",53:"WA",54:"WV",55:"WI",56:"WY"}

# count stands per (variant, fips)
cnt=collections.defaultdict(lambda: collections.Counter())   # fips -> Counter(variant)
for f in glob.glob(os.path.join(SI,"standinit_*.csv")):
    var=os.path.basename(f)[len("standinit_"):-4]
    if var.startswith("_"): continue
    with open(f, newline="") as fh:
        rd=csv.reader(fh); hdr=next(rd)
        try: si=hdr.index("STATE")
        except ValueError: continue
        for row in rd:
            try: fp=int(float(row[si]))
            except: continue
            cnt[fp][var]+=1
# dominant variant per state
dom={}
for fp,c in cnt.items():
    if fp in FIPS: dom[FIPS[fp]]=c.most_common(1)[0][0]

# current max idx in manifest_post.tsv
mp=os.path.join(SD,"manifest_post.tsv")
maxidx=-1
if os.path.exists(mp):
    for line in open(mp):
        p=line.split()
        if p and p[0].isdigit(): maxidx=max(maxidx,int(p[0]))

rows=[]; idx=0
for ab in sorted(dom):
    if ab in DONE: continue
    var=dom[ab].lower()
    if not os.path.exists(f"/users/PUOM0008/crsfaaron/fvs-modern/config/calibrated/{var}_draws.json"):
        print(f"  skip {ab}: no {var}_draws.json"); continue
    if not os.path.exists(os.path.join(SD,"standinit_by_variant",f"standinit_{var.upper()}.csv")):
        print(f"  skip {ab}: no standinit_{var.upper()}"); continue
    for d in range(NDRAW):
        rows.append((idx,var,ab,d)); idx+=1

ext=os.path.join(SD,"manifest_post_ext.tsv")
with open(ext,"w") as o:
    for r in rows: o.write("\t".join(map(str,r))+"\n")
print(f"dominant variants: {dict(sorted(dom.items()))}")
print(f"wrote {len(rows)} rows ({len(set(r[2] for r in rows))} states x {NDRAW} draws) -> {ext}")
print(f"array range: {maxidx+1}-{idx-1}")
