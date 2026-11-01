# =============================================================================
#  ADIM 3 — Kayıtların hazırlanması
#  Yalnız saha kayıtları (kuzu.csv, yalnız var) kullanılır.
#  temizlik → UTM 36N → 30 m model hücresinde tekleme → (isteğe bağlı) mesafe inceltmesi
#  → arka plan noktaları (rastgele ve yanlılık ağırlıklı)
#  Çıktılar (data/occurrences/):
#     occ_local_clean.csv         temizlenmiş tüm kayıtlar (lon, lat, x, y)
#     occ_local_thinned.csv       biomod2-sdm-runner girdisi: Tur, x, y, presence (UTM)
#     bg_local_random.csv, bg_local_bias.csv   arka plan (x, y)
#     bias_surface_local.tif      örnekleme yanlılığı yüzeyi
#     occ_summary.csv             sayılar (Yöntem bölümü için)
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
set.seed(20261030)
summ <- list(); S <- function(k, v) summ[[k]] <<- v

# ---- A1. Okuma ve temizlik --------------------------------------------------
occ <- read_occ()
if (is.null(occ)) stop("Kayıt dosyası bulunamadı: ", OCC_FILE)
S("ham_kayit", nrow(occ))
occ <- occ[occ$lon > 25 & occ$lon < 45 & occ$lat > 35 & occ$lat < 43, ]   # Türkiye dışı / ters girilmiş
occ <- occ[!duplicated(round(occ[, c("lon", "lat")], 6)), ]
S("tam_tekrar_sonrasi", nrow(occ))

pv  <- project(vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326"), CRS_LOCAL)
occ <- cbind(occ, crds(pv))
write.csv(occ[, c("lon", "lat", "x", "y")], file.path(OCC_DIR, "occ_local_clean.csv"), row.names = FALSE)

# 1 km CHELSA hücresi sayısı — yerel iklim çeşitliliğinin sınırını gösterir
S("farkli_1km_iklim_hucresi", nrow(unique(floor(occ[, c("lon", "lat")] * 120))))

# ---- A2. Model hücresinde tekleme ---------------------------------------------
tmpl <- local_template(RES_LOCAL)
cid  <- cellFromXY(tmpl, as.matrix(occ[, c("x", "y")]))
occ  <- occ[!duplicated(cid), ]
S(paste0("hucre_tekleme_", RES_LOCAL, "m"), nrow(occ))

# ---- A3. Mekânsal inceltme (rastgele açgözlü, THIN_REPS tekrar) -------------
thin_once <- function(xy, d) {
  D <- as.matrix(dist(xy)); keep <- integer(0)
  for (i in sample(nrow(xy))) if (!length(keep) || all(D[i, keep] >= d)) keep <- c(keep, i)
  keep
}
thin_best <- function(xy, d, reps = THIN_REPS) {
  best <- integer(0)
  for (k in seq_len(reps)) { kk <- thin_once(xy, d); if (length(kk) > length(best)) best <- kk }
  sort(best)
}
xy   <- as.matrix(occ[, c("x", "y")])
occ_t <- if (THIN_DIST_M > 0) occ[thin_best(xy, THIN_DIST_M), ] else occ
S("modele_giren_kayit", nrow(occ_t))
for (d in c(50, 100)) S(paste0("duyarlilik_inceltme_", d, "m"), length(thin_best(xy, d, 30)))

nnd <- apply(as.matrix(dist(xy)) + diag(Inf, nrow(xy)), 1, min)
S("en_yakin_komsu_medyan_m", round(median(nnd)))

write.csv(data.frame(Tur = SPECIES_NAME, x = round(occ_t$x, 1), y = round(occ_t$y, 1), presence = 1),
          file.path(OCC_DIR, "occ_local_thinned.csv"), row.names = FALSE)

# ---- A4. Arka plan noktaları ------------------------------------------------
#  Alan: yerel şablon. Göl yüzeyi (düz ve LAKE_MAX_ELEV altı) dışarıda bırakılır.
msk <- tmpl; values(msk) <- 1
fe <- file.path(TOPO_DIR, "elev.tif"); fs <- file.path(TOPO_DIR, "slope.tif")
if (file.exists(fe) && file.exists(fs)) {
  el <- rast(fe); sl <- rast(fs)
  if (compareGeom(el, tmpl, stopOnError = FALSE)) {
    msk <- ifel(el <= LAKE_MAX_ELEV & sl < 0.5, NA, 1)
    S("arka_plan_disi_hucre_gol_veya_veri_yok", global(is.na(msk), "sum")[1, 1])
  } else log_msg("UYARI: topo ızgarası şablonla uyuşmuyor; 02_dem_topography.R yeniden çalıştırılmalı.")
} else log_msg("Not: topo katmanları yok; göl maskesi uygulanmadı. 02'den sonra bu betiği yeniden çalıştırın.")

# (a) rastgele arka plan
bg_r <- spatSample(msk, N_BACKGROUND, method = "random", na.rm = TRUE, xy = TRUE, values = FALSE)
write.csv(data.frame(x = bg_r[, 1], y = bg_r[, 2]), file.path(OCC_DIR, "bg_local_random.csv"), row.names = FALSE)

# (b) yanlılık ağırlıklı arka plan: tüm (inceltilmemiş) kayıtların Gauss çekirdek yoğunluğu.
#     Yerel bilgiyle yönlendirilmiş örneklemenin coğrafi yanlılığını arka plana taşır
#     (Phillips et al. 2009; Fourcade et al. 2014).
#     Yoğunluk 150 m'de hesaplanıp model ızgarasına bilineer aktarılır (hız için).
tc   <- aggregate(tmpl, 5)
cnt  <- rasterize(as.matrix(occ[, c("x", "y")]), tc, fun = "count", background = 0)
w    <- focalMat(cnt, BIAS_SIGMA_M, type = "Gauss")
bias <- resample(focal(cnt, w = w, fun = "sum", na.rm = TRUE), tmpl, method = "bilinear")
bias <- (bias - global(bias, "min", na.rm = TRUE)[1, 1]) + 1e-6
bias <- mask(bias / global(bias, "max", na.rm = TRUE)[1, 1], msk)
names(bias) <- "bias"
writeRaster(bias, file.path(OCC_DIR, "bias_surface_local.tif"), overwrite = TRUE)
p    <- values(bias, mat = FALSE); ok <- which(!is.na(p))
ci   <- sample(ok, N_BACKGROUND, replace = TRUE, prob = p[ok])
bxy  <- xyFromCell(bias, ci) + matrix(runif(2 * length(ci), -RES_LOCAL / 2, RES_LOCAL / 2), ncol = 2)
write.csv(data.frame(x = bxy[, 1], y = bxy[, 2]), file.path(OCC_DIR, "bg_local_bias.csv"), row.names = FALSE)

# ---- Özet -------------------------------------------------------------------
summ_df <- data.frame(adim = names(summ), deger = unlist(summ))
write.csv(summ_df, file.path(OCC_DIR, "occ_summary.csv"), row.names = FALSE)
print(summ_df, row.names = FALSE)
log_msg("ADIM 3 tamamlandı.")
