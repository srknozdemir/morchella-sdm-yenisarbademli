# =============================================================================
#  ADIM 2 — Arazi modeli ve topoğrafik / mikroiklim / kanopi değişkenleri (30 m)
#  Girdi önceliği:
#   1) data/raw/dtm.tif (çıplak arazi) — varsa türevler TOPO_CALC_RES'te hesaplanır,
#      30 m'ye ortalama ile toplanır; ince ölçek heterojenliği SD katmanlarıyla korunur.
#      data/raw/dsm.tif de varsa kanopi yüksekliği (DSM − DTM) metrikleri üretilir.
#   2) Yoksa Copernicus GLO-30 (AWS açık veri). GLO-30 bir yüzey modelidir; ormanda
#      ağaç taçlarını içerir. Makalede bu sınırlılık belirtilmelidir.
#  Çıktılar: data/topo/<değişken>.tif (UTM 36N, 30 m şablon ızgara)
#
#   elev         yükselti (ölçek küçültme ve tanı; modele girmez)
#   slope        eğim (°)                     slope_sd  hücre içi eğim SD (yalnız DTM ile)
#   northness    cos(bakı)                    eastness  sin(bakı)
#   tpi100/300/1000  topoğrafik konum indeksi, üç ölçek (m) — soğuk hava birikimi, sırt/vadi
#   tri          arazi engebelilik indeksi
#   hli          ısı yükü indeksi (McCune & Keon 2002, eşitlik 3)
#   twi          topoğrafik nemlilik indeksi ln(a / tanβ)
#   insol_spr    Mart–Mayıs açık gök doğrudan ışınım (kWh m-2; gölge düşümü yok)
#   chm_mean, chm_sd, canopy_cover   kanopi yüksekliği ortalaması, SD, > CANOPY_MIN_H örtü oranı (DSM ile)
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))

tmpl     <- local_template(RES_LOCAL)
EDGE     <- 3000                                  # kenar etkisi için tampon (m)
work_ext <- ext(tmpl) + EDGE
has_dtm  <- file.exists(DTM_FILE)
has_dsm  <- file.exists(DSM_FILE)

# ---- 1. Arazi modeli --------------------------------------------------------
if (has_dtm) {
  dtm_raw <- rast(DTM_FILE)
  if (is.na(crs(dtm_raw)) || crs(dtm_raw) == "") stop("DTM dosyasında koordinat sistemi tanımlı değil: ", DTM_FILE)
  calc_res <- max(TOPO_CALC_RES, res(dtm_raw)[1])
  tc  <- rast(work_ext, resolution = calc_res, crs = CRS_LOCAL)
  dem <- project(dtm_raw, tc, method = if (calc_res > res(dtm_raw)[1] * 1.5) "average" else "bilinear")
  cov <- global(!is.na(crop(dem, ext(tmpl))), "mean")[1, 1]
  log_msg("DTM kullanıldı: ", basename(DTM_FILE), "; hesap çözünürlüğü ", calc_res, " m; model alanını kaplama oranı ",
          round(100 * cov, 1), "%")
  if (cov < 0.95) log_msg("UYARI: DTM model alanının tamamını kaplamıyor. Model alanı DTM kapsamıyla sınırlandırılacak.")
  src <- "DTM"
} else {
  bb   <- ext(project(as.polygons(work_ext, crs = CRS_LOCAL), "EPSG:4326"))
  tile_url <- function(la, lo) {
    id <- sprintf("Copernicus_DSM_COG_10_%s%02d_00_%s%03d_00_DEM",
                  ifelse(la >= 0, "N", "S"), abs(la), ifelse(lo >= 0, "E", "W"), abs(lo))
    sprintf("https://copernicus-dem-30m.s3.amazonaws.com/%s/%s.tif", id, id)
  }
  dem_dir <- file.path(RAW_DIR, "copdem30"); dir.create(dem_dir, FALSE, TRUE)
  options(timeout = 3600); tiles <- character()
  for (la in floor(bb$ymin):floor(bb$ymax)) for (lo in floor(bb$xmin):floor(bb$xmax)) {
    u <- tile_url(la, lo); f <- file.path(dem_dir, basename(u))
    if (!file.exists(f)) {
      ok <- tryCatch({ download.file(u, f, mode = "wb", quiet = TRUE); TRUE },
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (!ok) { unlink(f); log_msg("DEM karosu indirilemedi: ", u); next }
    }
    tiles <- c(tiles, f)
  }
  if (!length(tiles)) stop("Hiç DEM karosu yok. Ağ erişimini denetleyin ya da data/raw/dtm.tif ekleyin.")
  dem_ll <- if (length(tiles) > 1) merge(sprc(lapply(tiles, rast))) else rast(tiles)
  calc_res <- RES_LOCAL
  tc  <- rast(work_ext, resolution = calc_res, crs = CRS_LOCAL)
  dem <- project(crop(dem_ll, bb + 0.02), tc, method = "bilinear")
  log_msg("Copernicus GLO-30 kullanıldı (yüzey modeli).")
  src <- "GLO-30"
}
names(dem) <- "elev"
writeRaster(dem, file.path(RAW_DIR, sprintf("terrain_%sm_utm.tif", calc_res)), overwrite = TRUE, wopt = WOPT)

# ---- 2. Temel türevler (hesap çözünürlüğünde) --------------------------------
slp_r <- terrain(dem, "slope",  unit = "radians")
asp_r <- terrain(dem, "aspect", unit = "radians")
slope     <- slp_r * 180 / pi;    names(slope) <- "slope"
northness <- cos(asp_r);          names(northness) <- "northness"
eastness  <- sin(asp_r);          names(eastness)  <- "eastness"
tri       <- terrain(dem, "TRI"); names(tri) <- "tri"

# ---- 3. Isı yükü indeksi (McCune & Keon 2002, eşitlik 3) ---------------------
# Alan enlemde ~0.3° kapsadığından enlem sabit (alan merkezi) alınır
lat0  <- mean(local_bbox()[c("ymin", "ymax")])
lat_r <- lat0 * pi / 180
A   <- abs(pi - abs(asp_r - 5 * pi / 4))
hli <- exp(-1.467 + 1.582 * cos(lat_r) * cos(slp_r) - 1.5 * cos(A) * sin(slp_r) * sin(lat_r) -
           0.262 * sin(lat_r) * sin(slp_r) + 0.607 * sin(A) * sin(slp_r))
names(hli) <- "hli"

# ---- 4. Topoğrafik nemlilik indeksi -----------------------------------------
twi <- NULL
wbt_ok <- requireNamespace("whitebox", quietly = TRUE) &&
  isTRUE(tryCatch(file.exists(whitebox::wbt_exe_path(shell_quote = FALSE)), error = function(e) FALSE))
if (wbt_ok) {
  wd <- file.path(RAW_DIR, "wbt"); dir.create(wd, FALSE)
  f_dem <- file.path(wd, "dem.tif"); writeRaster(dem, f_dem, overwrite = TRUE)
  whitebox::wbt_breach_depressions_least_cost(f_dem, file.path(wd, "dem_br.tif"), dist = 50, fill = TRUE)
  whitebox::wbt_d_inf_flow_accumulation(file.path(wd, "dem_br.tif"), file.path(wd, "sca.tif"),
                                        out_type = "Specific Contributing Area")
  whitebox::wbt_slope(file.path(wd, "dem_br.tif"), file.path(wd, "slope.tif"), units = "degrees")
  whitebox::wbt_wetness_index(file.path(wd, "sca.tif"), file.path(wd, "slope.tif"), file.path(wd, "twi.tif"))
  twi <- rast(file.path(wd, "twi.tif")); crs(twi) <- CRS_LOCAL
  log_msg("TWI: WhiteboxTools D-infinity")
} else if ("flowAccumulation" %in% ls("package:terra")) {
  acc <- flowAccumulation(terrain(dem, "flowdir"))
  twi <- log(((acc + 1) * calc_res) / max(tan(slp_r), 0.001))
  log_msg("TWI: terra D8 (WhiteboxTools bulunamadı)")
} else if (file.exists(file.path(RAW_DIR, "twi_native.tif"))) {
  twi <- rast(file.path(RAW_DIR, "twi_native.tif")); crs(twi) <- CRS_LOCAL
  log_msg("TWI: 02b_twi.py çıktısı (pysheds D-infinity)")
} else log_msg("UYARI: TWI hesaplanamadı. install.packages('whitebox'); whitebox::install_whitebox() ya da 02b_twi.py")
if (!is.null(twi)) names(twi) <- "twi"

log_msg("Işınım hesabı başlıyor")
# ---- 5. Fruktifikasyon dönemi güneş ışınımı (Mart–Mayıs) ---------------------
sun_pos <- function(doy, hour, lat_deg) {
  d  <- 23.45 * sin(2 * pi * (284 + doy) / 365) * pi / 180
  ha <- (hour - 12) * 15 * pi / 180; la <- lat_deg * pi / 180
  alt <- asin(sin(la) * sin(d) + cos(la) * cos(d) * cos(ha))
  az  <- acos(pmin(1, pmax(-1, (sin(d) - sin(alt) * sin(la)) / (cos(alt) * cos(la)))))
  if (ha > 0) az <- 2 * pi - az
  c(alt = alt * 180 / pi, az = az * 180 / pi)
}
#  Hesap satır blokları hâlinde yapılır (bellek kullanımı sabit kalır):
#  cos(i) = cos(Z) cos(S) + sin(Z) sin(S) cos(A_güneş − A_yamaç)
sun_tab <- do.call(rbind, lapply(seq(60, 151, by = 5), function(dy) do.call(rbind, lapply(5:19, function(h) {
  sp <- sun_pos(dy, h, lat0); if (sp["alt"] <= 2) return(NULL)
  data.frame(zr = (90 - sp[["alt"]]) * pi / 180, ar = sp[["az"]] * pi / 180,
             w = 1367 * 0.70^(1 / sin(sp[["alt"]] * pi / 180)) * 5)   # Wh m-2 (1 sa × 5 gün)
}))))
f_ins <- file.path(RAW_DIR, "insol_native.tif")
insol <- rast(dem)
bs <- writeStart(insol, f_ins, overwrite = TRUE, wopt = WOPT)
readStart(slp_r); readStart(asp_r)
for (b in seq_len(bs$n)) {
  sv <- readValues(slp_r, bs$row[b], bs$nrows[b]); av <- readValues(asp_r, bs$row[b], bs$nrows[b])
  cs <- cos(sv); ss <- sin(sv); acc <- numeric(length(sv))
  for (k in seq_len(nrow(sun_tab))) {
    ci <- cos(sun_tab$zr[k]) * cs + sin(sun_tab$zr[k]) * ss * cos(sun_tab$ar[k] - av)
    ci[ci < 0] <- 0; acc <- acc + ci * sun_tab$w[k]
  }
  writeValues(insol, acc, bs$row[b], bs$nrows[b])
}
readStop(slp_r); readStop(asp_r); insol <- writeStop(insol)
insol <- insol / 1000; names(insol) <- "insol_spr"

log_msg("30 m'ye toplama başlıyor")
# ---- 6. 30 m'ye toplama -----------------------------------------------------
to_tmpl <- function(x, fun = "average") {
  if (abs(res(x)[1] - RES_LOCAL) < 1e-6) return(resample(x, tmpl, method = "bilinear"))
  resample(x, tmpl, method = fun)
}
out <- list(elev = dem, slope = slope, northness = northness, eastness = eastness,
            tri = tri, hli = hli, insol_spr = insol)
if (!is.null(twi)) out$twi <- twi
out_l <- lapply(out, to_tmpl)
# Hücre içi SD: E[x²] − E[x]² üzerinden; tam sayı olmayan oranlarda da (ör. 12 m → 30 m) doğru çalışır
agg_sd <- function(x) {
  m  <- resample(x,   tmpl, method = "average")
  m2 <- resample(x^2, tmpl, method = "average")
  sqrt(clamp(m2 - m^2, lower = 0, values = TRUE))
}
if (calc_res < RES_LOCAL) out_l$slope_sd <- agg_sd(slope)

# Çok ölçekli TPI: her yarıçap, yarıçap/10 çözünürlükte hesaplanır
tpi_focal <- function(x, radius_m) {
  w <- focalMat(x, radius_m, type = "circle"); w[w > 0] <- 1; w[w == 0] <- NA
  x - focal(x, w = w, fun = "mean", na.rm = TRUE)    # kenarda yalnız mevcut hücrelerin ortalaması
}
for (r in c(100, 300, 1000)) {
  rr <- max(calc_res, r / 10)
  d_r <- if (rr > calc_res) resample(dem, rast(work_ext, resolution = rr, crs = CRS_LOCAL), method = "average") else dem
  out_l[[paste0("tpi", r)]] <- resample(tpi_focal(d_r, r), tmpl, method = "bilinear")
}

log_msg("Kanopi hesabı başlıyor")
# ---- 7. Kanopi metrikleri (DSM − DTM) ---------------------------------------
#  Fark DSM'nin kendi çözünürlüğünde alınır; DTM bu ızgaraya bilineer aktarılır.
#  (DTM DSM'den kabaysa — ör. 12 m'ye karşı 5 m — dik yamaçlarda enterpolasyon hatası
#   1–3 m'lik sahte "kanopi" üretebilir; örtü eşiği bu yüzden CANOPY_MIN_H ile verilir.)
if (has_dtm && has_dsm) {
  dsm_raw <- rast(DSM_FILE)
  if (is.na(crs(dsm_raw)) || crs(dsm_raw) == "") stop("DSM dosyasında koordinat sistemi tanımlı değil: ", DSM_FILE)
  dsm_n <- if (same.crs(dsm_raw, CRS_LOCAL)) crop(dsm_raw, work_ext) else
             crop(project(dsm_raw, CRS_LOCAL, res = res(dsm_raw)[1], method = "bilinear"), work_ext)
  dtm_n <- resample(project(rast(DTM_FILE), CRS_LOCAL, method = "bilinear"), dsm_n, method = "bilinear")
  dh    <- dsm_n - dtm_n
  # Dikey datum / hizalama denetimi: açık alanlarda fark ~0 olmalıdır
  q <- quantile(spatSample(dh, 200000, na.rm = TRUE)[, 1], c(0.01, 0.05, 0.25, 0.5), na.rm = TRUE)
  log_msg("DSM − DTM yüzdelikleri (m): ", paste(names(q), round(q, 2), collapse = "; "))
  write.csv(data.frame(q = names(q), m = q), file.path(LOG_DIR, "chm_quantiles.csv"), row.names = FALSE)
  if (abs(q[["5%"]]) > 3)
    log_msg("UYARI: DSM − DTM'nin alt yüzdelikleri 0'dan belirgin uzak. DSM ve DTM farklı dikey datuma ",
            "(elipsoit / geoit) ya da farklı hizalamaya sahip olabilir; kanopi metrikleri güvenilir değildir.")
  chm <- clamp(dh, lower = 0, upper = 45, values = TRUE); names(chm) <- "chm"
  writeRaster(chm, file.path(RAW_DIR, sprintf("chm_%sm.tif", round(res(chm)[1]))), overwrite = TRUE, wopt = WOPT)
  out_l$chm_mean     <- resample(chm, tmpl, method = "average")
  out_l$chm_sd       <- agg_sd(chm)
  out_l$canopy_cover <- resample(chm > CANOPY_MIN_H, tmpl, method = "average")
  log_msg("Kanopi metrikleri üretildi (DSM ", round(res(dsm_n)[1], 1), " m − DTM ", round(res(rast(DTM_FILE))[1], 1),
          " m); örtü eşiği ", CANOPY_MIN_H, " m.")
} else log_msg("Not: DSM + DTM birlikte bulunmadığı için kanopi metrikleri üretilmedi.")

# ---- 8. Yazma ve özet -------------------------------------------------------
for (nm in names(out_l)) {
  r <- out_l[[nm]]; names(r) <- nm
  writeRaster(r, file.path(TOPO_DIR, paste0(nm, ".tif")), overwrite = TRUE, wopt = WOPT)
}
s <- rast(file.path(TOPO_DIR, paste0(names(out_l), ".tif")))
summ <- data.frame(var = names(s), t(sapply(names(s), function(v)
  unlist(global(s[[v]], c("min", "mean", "max"), na.rm = TRUE)))), kaynak = src)
write.csv(summ, file.path(LOG_DIR, "topo_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)
log_msg("ADIM 2 tamamlandı: ", nlyr(s), " katman, ", RES_LOCAL, " m (kaynak: ", src, ")")
