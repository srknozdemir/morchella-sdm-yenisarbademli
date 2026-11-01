"""
ADIM 2b — Topoğrafik nemlilik indeksi (TWI), WhiteboxTools bulunmadığında.
pysheds ile D-sonsuz (D-infinity) akış biriktirme; TWI = ln(a / tan β),
a = (biriken hücre + 1) × hücre boyu (özgül katkı alanı, m).
Girdi : data/raw/terrain_*_utm.tif  (02_dem_topography.R çıktısı)
Çıktı : data/raw/twi_native.tif ; 30 m'ye toplama 02_dem_topography.R'de yapılır.
Kullanım: python3 02b_twi.py <proje_klasörü>
Gerekli: pip install pysheds rasterio
"""
import sys, glob, os
import numpy as np
from pysheds.grid import Grid
import rasterio

proj = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
src = sorted(glob.glob(os.path.join(proj, "data", "raw", "terrain_*_utm.tif")))[0]
grid = Grid.from_raster(src)
dem = grid.read_raster(src)
dem = grid.fill_pits(dem)
dem = grid.fill_depressions(dem)
dem = grid.resolve_flats(dem)
fdir = grid.flowdir(dem, routing="dinf")
acc = grid.accumulation(fdir, routing="dinf")
with rasterio.open(src) as r:
    z = r.read(1).astype("float64"); prof = r.profile; res = r.res[0]
    nod = r.nodata
gy, gx = np.gradient(np.where(np.isfinite(z), z, np.nan), res)
slope = np.arctan(np.hypot(gx, gy))
a = (np.asarray(acc, dtype="float64") + 1.0) * res
twi = np.log(a / np.maximum(np.tan(slope), 0.001))
twi[~np.isfinite(z)] = np.nan
if nod is not None:
    twi[z == nod] = np.nan
prof.update(dtype="float32", nodata=np.nan, compress="deflate")
out = os.path.join(proj, "data", "raw", "twi_native.tif")
with rasterio.open(out, "w", **prof) as w:
    w.write(twi.astype("float32"), 1)
print("TWI yazıldı:", out, "aralık", np.nanmin(twi), np.nanmax(twi))
