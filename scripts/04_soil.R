# =============================================================================
#  ADIM 4 — Edafik değişkenler: SoilGrids 2.0 (250 m), 0–30 cm ağırlıklı ortalama
#  Katmanlar uzaktan (/vsicurl/) okunur, yalnız çalışma penceresi kesilir.
#  250 m veri 30 m ızgaraya bilineer aktarılır: bu çözünürlük artışı yeni bilgi
#  üretmez; makalede "250 m'den yeniden örneklenmiştir" diye belirtilmelidir.
#  Çıktılar: data/soil/soil_<değişken>.tif (30 m, UTM 36N)
#  Kaynak: Poggio et al. (2021) SOIL 7:217–240
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))

SG_BASE <- "https://files.isric.org/soilgrids/latest/data"
# SoilGrids depolama birimi → raporlama birimi
conv <- c(phh2o = 10, soc = 10, clay = 10, sand = 10, silt = 10, cfvo = 10, bdod = 100, nitrogen = 100, cec = 10)
unit <- c(phh2o = "pH (H2O)", soc = "g/kg", clay = "%", sand = "%", silt = "%", cfvo = "hacim %",
          bdod = "g/cm3", nitrogen = "g/kg", cec = "cmol(c)/kg")
dw   <- c("0-5cm" = 5, "5-15cm" = 10, "15-30cm" = 15)    # katman kalınlıkları (cm)

tmpl <- local_template(RES_LOCAL)
win  <- as.polygons(ext(tmpl) + 1000, crs = CRS_LOCAL)
raw_dir <- file.path(RAW_DIR, "soilgrids"); dir.create(raw_dir, FALSE, TRUE)
reg <- list()

for (v in SOIL_VARS) {
  out <- file.path(SOIL_DIR, paste0("soil_", v, ".tif"))
  if (file.exists(out)) { log_msg("var: ", basename(out)); next }
  lay <- list()
  for (d in SOIL_DEPTHS) {
    fr <- file.path(raw_dir, sprintf("%s_%s_mean.tif", v, d))
    if (!file.exists(fr)) {
      u <- sprintf("%s/%s/%s_%s_mean.vrt", SG_BASE, v, v, d)
      r <- tryCatch({
        x  <- rast(paste0("/vsicurl/", u))
        ex <- ext(project(win, crs(x)))
        writeRaster(crop(x, ex, snap = "out"), fr, overwrite = TRUE, wopt = WOPT)
      }, error = function(e) { log_msg("SoilGrids okunamadı: ", v, " ", d, " — ", conditionMessage(e)); NULL })
      if (is.null(r)) next
    }
    lay[[d]] <- rast(fr)
  }
  if (length(lay) == 0) { reg[[v]] <- "missing"; next }
  w  <- dw[names(lay)]
  st <- rast(lay)
  st <- classify(st, cbind(0, NA))  # SoilGrids'te 0 = veri yok (su, yerleşim)
  wm <- sum(st * w, na.rm = TRUE) / sum((!is.na(st)) * w)
  wm <- wm / conv[[v]]
  r30 <- project(wm, tmpl, method = "bilinear")
  names(r30) <- paste0("soil_", v)
  writeRaster(r30, out, overwrite = TRUE, wopt = WOPT)
  g <- global(r30, c("min", "mean", "max"), na.rm = TRUE)
  log_msg(sprintf("soil_%s (%s): %.2f / %.2f / %.2f  [%d katman]", v, unit[[v]], g[1, 1], g[1, 2], g[1, 3], length(lay)))
  reg[[v]] <- "ok"
}

# Kum + silt + kil toplam denetimi (~%100 olmalı)
f3 <- file.path(SOIL_DIR, paste0("soil_", c("sand", "silt", "clay"), ".tif"))
if (all(file.exists(f3))) {
  tot <- sum(rast(f3))
  log_msg("Doku toplamı (kum+silt+kil) ortalaması: ", round(global(tot, "mean", na.rm = TRUE)[1, 1], 1), " %")
}
write.csv(data.frame(var = names(reg), status = unlist(reg)), file.path(LOG_DIR, "soil_registry.csv"), row.names = FALSE)
log_msg("ADIM 4 tamamlandı.")
