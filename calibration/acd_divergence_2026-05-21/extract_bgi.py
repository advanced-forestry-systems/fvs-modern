#!/usr/bin/env python3
# Extract BGI from /users/PUOM0008/crsfaaron/raster_layers/bgi/ME_BGI_V1.tif at
# each ME FIA plot's LAT/LON. Uses gdallocationinfo with -wgs84 flag to handle
# the projection from WGS84 to the raster's native CRS (NAD83 UTM 19N).
# Writes /users/PUOM0008/crsfaaron/acadgy_fia_verify/me_bgi_by_pltcn.csv.
import csv
import subprocess
import sys

RASTER = "/users/PUOM0008/crsfaaron/raster_layers/bgi/ME_BGI_V1.tif"
PLOT_CSV = "/users/PUOM0008/crsfaaron/fia_data/ME_PLOT.csv"
OUT_CSV = "/users/PUOM0008/crsfaaron/acadgy_fia_verify/me_bgi_by_pltcn.csv"

with open(PLOT_CSV) as f:
    reader = csv.DictReader(f)
    rows = []
    for r in reader:
        try:
            lat = float(r["LAT"]); lon = float(r["LON"]); cn = r["CN"]
            if -90 < lat < 90 and -180 < lon < 180:
                rows.append((cn, lat, lon))
        except (KeyError, ValueError):
            continue

print(f"loaded {len(rows)} plots with valid LAT/LON", flush=True)

# Batch gdallocationinfo via stdin "lon lat\n" per query, -wgs84 flag
proc = subprocess.Popen(
    ["gdallocationinfo", "-wgs84", "-valonly", RASTER],
    stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
    text=True, bufsize=1
)
input_text = "\n".join(f"{lon} {lat}" for cn, lat, lon in rows) + "\n"
out, err = proc.communicate(input_text)
vals = out.strip().split("\n")
print(f"gdallocationinfo returned {len(vals)} values (rc={proc.returncode})", flush=True)
if proc.returncode != 0:
    print("STDERR:", err[:500], flush=True)

with open(OUT_CSV, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["CN", "LAT", "LON", "BGI"])
    matched = 0
    for (cn, lat, lon), v in zip(rows, vals):
        try:
            fv = float(v)
            w.writerow([cn, lat, lon, fv])
            if fv > 0:
                matched += 1
        except ValueError:
            w.writerow([cn, lat, lon, ""])
print(f"wrote {OUT_CSV} with {matched} non-zero BGI matches of {len(rows)} plots", flush=True)
