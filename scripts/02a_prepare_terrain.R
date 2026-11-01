# =============================================================================
#  ADIM 2a — HGM SAM (12 m) ve SYM (5 m pafta) verilerinin hazırlanması
#   1. SAM (ArcGRID, Lambert TC1M) model alanı + 3.5 km için kesilir, UTM 36N'ye aktarılır
#   2. SYM paftalarından model alanıyla kesişenler seçilir ve tek mozaikte birleştirilir
#   3. DTM, DSM'ye yatay olarak hizalanır: düz ve açık alanlarda DSM − DTM farkının
#      yayılımını (MAD) en küçük yapan kayma 1 m adımla aranır
#  Çıktılar: data/raw/dtm.tif, data/raw/dsm.tif, logs/terrain_prep.csv
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
set.seed(31)

tmpl     <- local_template(RES_LOCAL)
work_ext <- ext(tmpl) + 3500
work_pol <- as.polygons(work_ext, crs = CRS_LOCAL)
info <- list(); I <- function(k, v) info[[k]] <<- as.character(v)

# ---- 1. SAM -----------------------------------------------------------------
if (!dir.exists(HGM_SAM_PATH) && !file.exists(HGM_SAM_PATH)) stop("SAM bulunamadı: ", HGM_SAM_PATH)
sam  <- rast(HGM_SAM_PATH)
samc <- crop(sam, project(work_pol, crs(sam)), snap = "out")
r_sam <- res(sam)[1]
dtm  <- project(samc, rast(work_ext, resolution = r_sam, crs = CRS_LOCAL), method = "bilinear")
names(dtm) <- "dtm"
I("sam_kaynak", HGM_SAM_PATH); I("sam_cozunurluk_m", round(r_sam, 3))
I("sam_kapsam_orani", round(global(!is.na(crop(dtm, ext(tmpl))), "mean")[1, 1], 4))

# ---- 2. SYM paftaları --------------------------------------------------------
tifs <- unlist(lapply(HGM_SYM_DIRS[dir.exists(HGM_SYM_DIRS)], list.files,
                      pattern = "\\.tif$", full.names = TRUE, ignore.case = TRUE))
hit <- character()
for (f in tifs) {
  r <- tryCatch(rast(f), error = function(e) NULL); if (is.null(r)) next
  e <- if (same.crs(r, CRS_LOCAL)) ext(r) else ext(project(as.polygons(ext(r), crs = crs(r)), CRS_LOCAL))
  if (relate(ext(e), work_ext, "intersects")) hit <- c(hit, f)
}
# Aynı paftanın yinelenen kopyalarını (ör. N26-d3_5m.tif / N26d3_5.tif) ayıkla
key <- toupper(gsub("[^A-Za-z0-9]", "", sub("_5m?\\.tif$", "", basename(hit), ignore.case = TRUE)))
hit <- hit[!duplicated(key)]
if (!length(hit)) stop("Model alanıyla kesişen SYM paftası bulunamadı.")
log_msg("SYM paftaları (", length(hit), "): ", paste(basename(hit), collapse = ", "))
I("sym_paftalar", paste(basename(hit), collapse = " "))

crs_ok <- vapply(hit, function(f) same.crs(rast(f), CRS_LOCAL), logical(1))
grid5  <- rast(work_ext, resolution = DSM_RES, crs = CRS_LOCAL)
if (all(crs_ok)) {
  v   <- vrt(hit, filename = file.path(RAW_DIR, "sym_mosaic.vrt"), overwrite = TRUE)
  dsm <- resample(crop(v, work_ext, snap = "out"), grid5, method = if (DSM_RES > 5) "average" else "bilinear")
} else {
  dsm <- mosaic(sprc(lapply(hit, function(f) project(rast(f), grid5, method = "bilinear"))), fun = "mean")
}
names(dsm) <- "dsm"
I("sym_kapsam_orani", round(global(!is.na(crop(dsm, ext(tmpl))), "mean")[1, 1], 4))

# ---- 3. Yatay hizalama (DTM → DSM) -----------------------------------------
#  Örnek noktalar: eğimi < 10° olan ve DSM − DTM farkı küçük (açık alan adayı) hücreler
slp <- terrain(dtm, "slope")
pts <- spatSample(dsm, 60000, method = "random", na.rm = TRUE, xy = TRUE)
pts$slope <- extract(slp, as.matrix(pts[, c("x", "y")]))[, 1]
pts <- pts[!is.na(pts$slope) & pts$slope < 10, ]
mad_at <- function(dx, dy) {
  z <- extract(dtm, cbind(pts$x + dx, pts$y + dy), method = "bilinear")[, 1]
  d <- pts$dsm - z; d <- d[is.finite(d)]
  d <- d[d < quantile(d, 0.5)]           # kanopi etkisini dışarıda tutmak için alt yarı
  median(abs(d - median(d)))
}
grid_search <- function(rng, step) {
  g <- expand.grid(dx = seq(-rng, rng, step), dy = seq(-rng, rng, step))
  g$mad <- mapply(mad_at, g$dx, g$dy); g[which.min(g$mad), ]
}
mad0 <- mad_at(0, 0)
b1   <- grid_search(COREG_MAX_SHIFT_M, 3)
g2   <- expand.grid(dx = b1$dx + (-3:3), dy = b1$dy + (-3:3)); g2$mad <- mapply(mad_at, g2$dx, g2$dy)
best <- g2[which.min(g2$mad), ]
log_msg(sprintf("Hizalama: kayma dx = %+d m, dy = %+d m; MAD %.2f m → %.2f m", best$dx, best$dy, mad0, best$mad))
I("hizalama_dx_m", best$dx); I("hizalama_dy_m", best$dy)
I("mad_once_m", round(mad0, 3)); I("mad_sonra_m", round(best$mad, 3))
if (abs(best$dx) >= COREG_MAX_SHIFT_M || abs(best$dy) >= COREG_MAX_SHIFT_M)
  log_msg("UYARI: kayma arama sınırına dayandı; COREG_MAX_SHIFT_M'yi büyütüp yeniden çalıştırın.")
# DSM(x) ≈ DTM(x + d)  →  hizalı DTM = DTM'nin −d kadar ötelenmesi
if (best$dx != 0 || best$dy != 0) {
  dtm <- resample(shift(dtm, dx = -best$dx, dy = -best$dy), rast(work_ext, resolution = r_sam, crs = CRS_LOCAL),
                  method = "bilinear")
}

# Dikey fark özeti (hizalama sonrası, düz alanlar)
z  <- extract(dtm, as.matrix(pts[, c("x", "y")]), method = "bilinear")[, 1]
dq <- quantile(pts$dsm - z, c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), na.rm = TRUE)
log_msg("Düz alanlarda DSM − DTM yüzdelikleri (m): ", paste(names(dq), round(dq, 2), collapse = "; "))
for (k in names(dq)) I(paste0("fark_", k), round(dq[[k]], 2))

writeRaster(dtm, DTM_FILE, overwrite = TRUE, wopt = WOPT)
writeRaster(dsm, DSM_FILE, overwrite = TRUE, wopt = WOPT)
write.csv(data.frame(alan = names(info), deger = unlist(info)), file.path(LOG_DIR, "terrain_prep.csv"), row.names = FALSE)
log_msg("ADIM 2a tamamlandı: ", basename(DTM_FILE), " (", round(r_sam, 2), " m), ", basename(DSM_FILE), " (", DSM_RES, " m)")
