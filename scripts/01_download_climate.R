# =============================================================================
#  ADIM 1 — CHELSA v2.1 iklim değişkenlerinin indirilmesi
#  - Güncel dönem (1981–2010) ve gelecek (CHELSA_PERIODS × CHELSA_SSPS × CHELSA_GCMS)
#  - Pencere: kayıtlar + CLIM_BUFFER_KM (30 m'ye ölçek küçültme regresyonu için)
#  - Dosyalar uzaktan /vsicurl/ ile okunur ve yalnızca gereken pencere kesilir;
#    küresel dosyanın tamamı indirilmez. Kesinti olursa betiği yeniden çalıştırın,
#    var olan dosyalar atlanır.
#  Çıktı yapısı (biomod2-sdm-runner "Automatic" modu ile uyumlu):
#    data/climate/current/<var>.tif
#    data/climate/future_<dönem>_<ssp>_<GCM>/<var>.tif
#  Kaynak: Karger et al. (2017) Sci Data 4:170122; Karger et al. (2023) Ecography e06535
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))

urls_current <- function(v, base) {
  c(file.path(base, CHELSA_BASELINE, "bio",
              sprintf("CHELSA_%s_%s_V.2.1.tif", v, CHELSA_BASELINE)))
}
urls_future <- function(v, period, ssp, gcm, base) {
  file.path(base, period, gcm, ssp, "bio",
            sprintf("CHELSA_%s_%s_%s_%s_V.2.1.tif", v, period, tolower(gcm), ssp))
}
urls_monthly <- function(var, m, period = NULL, ssp = NULL, gcm = NULL, base) {
  mm <- sprintf("%02d", m)
  if (is.null(period)) {
    file.path(base, CHELSA_BASELINE, var,
              sprintf("CHELSA_%s_%s_%s_V.2.1.tif", var, mm, CHELSA_BASELINE))
  } else {
    p2 <- sub("-", "_", period)
    file.path(base, period, gcm, ssp, var,
              sprintf("CHELSA_%s_r1i1p1f1_w5e5_%s_%s_%s_%s_norm.tif",
                      tolower(gcm), ssp, var, mm, p2))
  }
}

# Uzak dosyadan pencere kesme; başarısız olursa tam indirip keser
fetch_crop <- function(url_fun, out, bbox) {
  if (file.exists(out)) return("exists")
  e <- bbox_ext(bbox)
  for (base in CHELSA_BASE) {
    u <- url_fun(base)
    r <- tryCatch({
      x <- crop(rast(paste0("/vsicurl/", u)), e, snap = "out")
      names(x) <- sub("\\.tif$", "", basename(out))
      writeRaster(x, out, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
    }, error = function(err) NULL)
    if (!is.null(r)) return(paste("ok", base))
    # yedek: tam indirme
    tmp <- tempfile(fileext = ".tif")
    ok <- tryCatch({ download.file(u, tmp, mode = "wb", quiet = TRUE); TRUE },
                   error = function(err) FALSE, warning = function(w) FALSE)
    if (ok) {
      x <- crop(rast(tmp), e, snap = "out")
      names(x) <- sub("\\.tif$", "", basename(out))
      writeRaster(x, out, overwrite = TRUE, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
      unlink(tmp)
      return(paste("ok-full", base))
    }
  }
  "missing"
}

registry <- list()
add_reg <- function(...) registry[[length(registry) + 1]] <<- data.frame(..., stringsAsFactors = FALSE)
options(timeout = 3600)

# ---- 1. Güncel dönem --------------------------------------------------------
targets_now <- list(current = clim_bbox())
for (tg in names(targets_now)) {
  od <- file.path(CLIM_DIR, tg); dir.create(od, showWarnings = FALSE)
  for (v in CHELSA_VARS) {
    st <- fetch_crop(function(b) urls_current(v, b), file.path(od, paste0(v, ".tif")), targets_now[[tg]])
    log_msg(tg, " ", v, ": ", st); add_reg(set = tg, var = v, status = st)
  }
  if (CHELSA_MONTHLY) for (var in c("tas", "pr")) for (m in FRUITING_MONTHS) {
    nm <- sprintf("%s_%02d", var, m)
    dir.create(file.path(od, "monthly"), showWarnings = FALSE)
    st <- fetch_crop(function(b) urls_monthly(var, m, base = b),
                     file.path(od, "monthly", paste0(nm, ".tif")), targets_now[[tg]])
    log_msg(tg, " ", nm, ": ", st); add_reg(set = tg, var = nm, status = st)
  }
}

# ---- 2. Gelecek dönemler ----------------------------------------------------
fut_bbox <- clim_bbox()
for (p in CHELSA_PERIODS) for (s in CHELSA_SSPS) for (g in CHELSA_GCMS) {
  tg <- sprintf("future_%s_%s_%s", p, s, g)
  od <- file.path(CLIM_DIR, tg); dir.create(od, showWarnings = FALSE)
  for (v in CHELSA_VARS) {
    st <- fetch_crop(function(b) urls_future(v, p, s, g, b), file.path(od, paste0(v, ".tif")), fut_bbox)
    log_msg(tg, " ", v, ": ", st); add_reg(set = tg, var = v, status = st)
  }
  if (CHELSA_MONTHLY) for (var in c("tas", "pr")) for (m in FRUITING_MONTHS) {
    nm <- sprintf("%s_%02d", var, m)
    dir.create(file.path(od, "monthly"), showWarnings = FALSE)
    st <- fetch_crop(function(b) urls_monthly(var, m, p, s, g, b),
                     file.path(od, "monthly", paste0(nm, ".tif")), fut_bbox)
    log_msg(tg, " ", nm, ": ", st); add_reg(set = tg, var = nm, status = st)
  }
}

# ---- 3. Fruktifikasyon dönemi türev değişkenleri ----------------------------
#  tas_spr : Mart–Mayıs ortalama sıcaklık (°C)
#  pr_spr  : Mart–Mayıs toplam yağış (mm)
#  CHELSA aylık dosyalarında ölçek/ofset meta verisi terra tarafından uygulanır;
#  aşağıdaki birim denetimi yine de yapılır.
make_spring <- function(od) {
  md <- file.path(od, "monthly")
  ft <- file.path(md, sprintf("tas_%02d.tif", FRUITING_MONTHS))
  fp <- file.path(md, sprintf("pr_%02d.tif",  FRUITING_MONTHS))
  if (all(file.exists(ft)) && !file.exists(file.path(od, "tas_spr.tif"))) {
    t <- mean(rast(ft))
    if (global(t, "mean", na.rm = TRUE)[1, 1] > 200) t <- t - 273.15   # Kelvin ise
    writeRaster(t, file.path(od, "tas_spr.tif"), overwrite = TRUE)
  }
  if (all(file.exists(fp)) && !file.exists(file.path(od, "pr_spr.tif")))
    writeRaster(sum(rast(fp)), file.path(od, "pr_spr.tif"), overwrite = TRUE)
}
for (od in list.dirs(CLIM_DIR, recursive = FALSE)) make_spring(od)

# ---- 4. Özet ve denetimler --------------------------------------------------
reg <- do.call(rbind, registry)
write.csv(reg, file.path(LOG_DIR, "climate_download_registry.csv"), row.names = FALSE)
miss <- subset(reg, status == "missing")
if (nrow(miss)) {
  log_msg("UYARI: ", nrow(miss), " katman bulunamadı. Ayrıntı: logs/climate_download_registry.csv")
  print(table(miss$var))
}
# Gelecekte eksik olan değişken, modelde kullanılamaz: listeyi yaz
fut_sets <- grep("^future_", unique(reg$set), value = TRUE)
usable <- names(which(tapply(reg$status[reg$set %in% fut_sets] != "missing",
                             reg$var[reg$set %in% fut_sets], all)))
writeLines(usable, file.path(LOG_DIR, "vars_available_in_all_futures.txt"))
log_msg("Tüm gelecek setlerinde bulunan değişkenler: ", paste(usable, collapse = ", "))

b1 <- file.path(CLIM_DIR, "current", "bio1.tif")
if (file.exists(b1)) {
  m <- global(rast(b1), "mean", na.rm = TRUE)[1, 1]
  log_msg("Birim denetimi — yerel bio1 ortalaması: ", round(m, 2),
          if (m > 200) "  (Kelvin×10 görünüyor; ölçek/ofset uygulanmamış olabilir)" else "  (°C)")
}
log_msg("ADIM 1 tamamlandı.")
