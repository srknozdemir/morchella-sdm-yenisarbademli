# =============================================================================
#  ADIM 5 — İklim değişkenlerinin 1 km'den 30 m'ye ölçek küçültülmesi
#  Güncel dönem: yükselti regresyonu + kalıntı enterpolasyonu
#     v_1km = a + b · z_1km + e          (z_1km: arazi modelinin CHELSA hücresi ortalaması)
#     v_30m = a + b · z_30m + bilineer(e)
#   - "lin" değişkenleri doğrudan, "log" değişkenleri log(v + 1) ölçeğinde modellenir
#   - diğerleri (ör. bio15) yalnız bilineer enterpolasyonla aktarılır
#   - regresyon, model alanı (+3 km tampon) içindeki CHELSA hücreleriyle kurulur;
#     b katsayısı sıcaklık değişkenlerinde yerel sıcaklık azalma oranıdır (°C/m)
#  Gelecek dönemler: delta (değişim) yöntemi
#     lin:  v_30m,F = v_30m,G + bilineer(v_1km,F − v_1km,G)
#     log:  v_30m,F = v_30m,G · bilineer(v_1km,F / v_1km,G)
#  Böylece arazinin ince ölçek deseni korunur, iklim sinyali GCM'den gelir.
#  Çıktılar: data/climate_30m/<set>/<var>.tif ; logs/downscale_coefficients.csv
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))

tmpl  <- local_template(RES_LOCAL)
z30   <- rast(file.path(TOPO_DIR, "elev.tif"))
if (!compareGeom(z30, tmpl, stopOnError = FALSE)) stop("topo/elev.tif şablonla uyuşmuyor; 02'yi yeniden çalıştırın.")
fz    <- list.files(RAW_DIR, pattern = "^terrain_.*_utm\\.tif$", full.names = TRUE)[1]
if (is.na(fz)) stop("Ham arazi modeli (data/raw/terrain_*_utm.tif) yok; önce 02'yi çalıştırın.")
zfine <- rast(fz)

cur_dir <- file.path(CLIM_DIR, "current")
vars <- sub("\\.tif$", "", list.files(cur_dir, pattern = "\\.tif$"))
if (!length(vars)) stop("data/climate/current boş; önce 01'i çalıştırın.")

# CHELSA ızgarasında arazi yüksekliği (hücre ortalaması), model alanı + 3 km
g1  <- rast(file.path(cur_dir, paste0(vars[1], ".tif")))
win <- ext(project(as.polygons(ext(tmpl) + 3000, crs = CRS_LOCAL), crs(g1)))
g1  <- crop(g1, win, snap = "out")
zC  <- project(zfine, g1, method = "average")
names(zC) <- "z"

rule_of <- function(v) if (v %in% DOWNSCALE_RULES$lin) "lin" else if (v %in% DOWNSCALE_RULES$log) "log" else "interp"
to30 <- function(x) project(x, tmpl, method = "bilinear")

coef_tab <- list()
cur30 <- list()
od <- file.path(CLIM30_DIR, "current"); dir.create(od, FALSE, TRUE)
for (v in vars) {
  vC <- crop(rast(file.path(cur_dir, paste0(v, ".tif"))), g1)
  rl <- rule_of(v)
  if (rl == "interp") {
    v30 <- to30(vC); a <- b <- r2 <- NA
  } else {
    y  <- if (rl == "log") log1p(clamp(vC, lower = 0, values = TRUE)) else vC
    df <- na.omit(data.frame(y = values(y, mat = FALSE), z = values(zC, mat = FALSE)))
    if (nrow(df) < 20 || sd(df$z) < 1) { v30 <- to30(vC); a <- b <- r2 <- NA; rl <- "interp" } else {
      fit <- lm(y ~ z, data = df); a <- coef(fit)[1]; b <- coef(fit)[2]; r2 <- summary(fit)$r.squared
      res <- y - (a + b * zC)
      v30 <- a + b * z30 + to30(res)
      if (rl == "log") v30 <- clamp(expm1(v30), lower = 0, values = TRUE)
    }
  }
  v30 <- mask(v30, z30); names(v30) <- v
  writeRaster(v30, file.path(od, paste0(v, ".tif")), overwrite = TRUE, wopt = WOPT)
  cur30[[v]] <- v30
  coef_tab[[v]] <- data.frame(var = v, rule = rl, intercept = unname(a), slope_per_m = unname(b),
                              slope_per_100m = unname(b) * 100, r2 = unname(r2), n_cells = ncell(vC))
}
coef_df <- do.call(rbind, coef_tab)
write.csv(coef_df, file.path(LOG_DIR, "downscale_coefficients.csv"), row.names = FALSE)
log_msg("Güncel iklim 30 m'ye indirildi: ", length(vars), " değişken.")
print(coef_df[, c("var", "rule", "slope_per_100m", "r2")], row.names = FALSE, digits = 3)

# ---- Gelecek: delta yöntemi ------------------------------------------------
fsets <- list.dirs(CLIM_DIR, recursive = FALSE, full.names = FALSE)
fsets <- grep("^future_", fsets, value = TRUE)
for (fs in fsets) {
  od <- file.path(CLIM30_DIR, fs); dir.create(od, FALSE, TRUE)
  fv <- intersect(vars, sub("\\.tif$", "", list.files(file.path(CLIM_DIR, fs), pattern = "\\.tif$")))
  for (v in fv) {
    out <- file.path(od, paste0(v, ".tif")); if (file.exists(out)) next
    vC <- crop(rast(file.path(cur_dir, paste0(v, ".tif"))), g1)
    vF <- crop(rast(file.path(CLIM_DIR, fs, paste0(v, ".tif"))), g1)
    r  <- if (rule_of(v) == "log") cur30[[v]] * to30((vF + 1) / (vC + 1)) else cur30[[v]] + to30(vF - vC)
    if (rule_of(v) == "log") r <- clamp(r + (to30((vF + 1) / (vC + 1)) - 1), lower = 0, values = TRUE)   # (v+1) oranı
    r <- mask(r, z30); names(r) <- v
    writeRaster(r, out, overwrite = TRUE, wopt = WOPT)
  }
  log_msg(fs, ": ", length(fv), " değişken")
}
log_msg("ADIM 5 tamamlandı.")
