#!/usr/bin/env python3
"""lookup_bgi_csi.py - Extract BGI and CSI from rasters at SILC CFI lat/long
and the 11 byStrata stand locations on Cardinal.
Runs on Cardinal."""
import os, sys
try:
    import rasterio
    from rasterio.warp import transform as rio_transform
except Exception:
    print("rasterio missing; install with: pip install --user rasterio")
    sys.exit(1)
import pandas as pd
from pyproj import Transformer

BGI_TIF = os.path.expanduser("~/raster_layers/bgi/ME_BGI_V1.tif")
# Find CSI raster
import glob
CSI_TIFS = glob.glob(os.path.expanduser("~/raster_layers/csi/**/*.tif"), recursive=True)
print("BGI:", BGI_TIF, "exists:", os.path.exists(BGI_TIF))
print("CSI:", CSI_TIFS)

def sample_raster(path, lat, lon):
    with rasterio.open(path) as src:
        # transform lat/long (EPSG:4326) to raster CRS
        xs, ys = rio_transform("EPSG:4326", src.crs, [lon], [lat])
        try:
            v = next(src.sample([(xs[0], ys[0])]))[0]
            return float(v) if v is not None else None
        except StopIteration:
            return None

# === SILC CFI: single approximate location ===
print("\n=== SILC CFI plots (approximate lat 46.4628, lon -68.4253) ===")
bgi_cfi = sample_raster(BGI_TIF, 46.4628, -68.4253)
print(f"  BGI: {bgi_cfi}")
for c in CSI_TIFS[:3]:
    v = sample_raster(c, 46.4628, -68.4253)
    print(f"  CSI from {os.path.basename(c)}: {v}")

# === byStrata stands: from StandInit lat/long if any ===
print("\n=== byStrata stands ===")
si = pd.read_csv(os.path.expanduser("~/silc_strata/Acadian_Matrix_StandInit_2023.csv"))
# Stand init has Climate Site Index already
print(si[["STAND_ID","INV_PLOT_SIZE","ClimateSiteIndexMeters","ElevationMeters"]].to_string(index=False))
print("\nNote: byStrata stands have CSI already in StandInit (ClimateSiteIndexMeters column)")
