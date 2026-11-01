# =============================================================================
#  ADIM 1 — WorldClim 2.1 iklim değişkenleri (30″ ≈ 1 km)
#  Güncel (1970–2000): geodata::worldclim_tile ile 30°'lik karo indirilir, pencere kesilir
#     bio1–bio19; WC_MONTHLY ise tavg ve prec karolarından tas_spr, pr_spr (Mart–Mayıs)
#  Gelecek (CMIP6): WC_GCMS × WC_SSPS × WC_PERIODS; dosyalar uzaktan (/vsicurl/) okunup
#     yalnız pencere kesilir; olmazsa geodata::cmip6_tile ile karo indirilir
#     bio1–bio19 ("bioc"); WC_MONTHLY ise tmin, tmax, prec → tas_spr, pr_spr
#  Pencere: çalışma alanı sınırı + CLIM_BUFFER_KM (ölçek küçültme regresyonu için)
#  Çıktı: data/climate/current/, data/climate/future_<dönem>_<ssp>_<GCM>/  (<var>.tif)
#  Kaynak: Fick & Hijmans (2017) Int J Climatol 37:4302–4315
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
if (!requireNamespace("geodata", quietly = TRUE)) install.packages("geodata", repos = "https://cloud.r-project.org")
options(timeout = 7200)

bb  <- clim_bbox(); e <- bbox_ext(bb)
lon0 <- mean(bb[c("xmin", "xmax")]); lat0 <- mean(bb[c("ymin", "ymax")])
wc_raw <- file.path(RAW_DIR, "worldclim"); dir.create(wc_raw, FALSE, TRUE)
reg <- list(); R <- function(set, var, st) { reg[[length(reg) + 1]] <<- data.frame(set = set, var = var, status = st) }
wr  <- function(x, f) writeRaster(x, f, overwrite = TRUE, wopt = WOPT)

# ---- 1. Güncel dönem --------------------------------------------------------
od <- file.path(CLIM_DIR, "current"); dir.create(od, FALSE, TRUE)
if (!all(file.exists(file.path(od, paste0("bio", 1:19, ".tif"))))) {
  b <- geodata::worldclim_tile("bio", lon = lon0, lat = lat0, path = wc_raw)
  b <- crop(b, e, snap = "out")
  for (i in 1:19) { x <- b[[i]]; names(x) <- paste0("bio", i); wr(x, file.path(od, paste0("bio", i, ".tif"))) }
}
for (i in 1:19) R("current", paste0("bio", i), "ok")
if (WC_MONTHLY && !file.exists(file.path(od, "tas_spr.tif"))) {
  ta <- crop(geodata::worldclim_tile("tavg", lon = lon0, lat = lat0, path = wc_raw), e, snap = "out")
  pr <- crop(geodata::worldclim_tile("prec", lon = lon0, lat = lat0, path = wc_raw), e, snap = "out")
  x <- mean(ta[[FRUITING_MONTHS]]); names(x) <- "tas_spr"; wr(x, file.path(od, "tas_spr.tif"))
  x <- sum(pr[[FRUITING_MONTHS]]);  names(x) <- "pr_spr";  wr(x, file.path(od, "pr_spr.tif"))
}
if (WC_MONTHLY) { R("current", "tas_spr", "ok"); R("current", "pr_spr", "ok") }
log_msg("WorldClim güncel dönem hazır: ", od)

# ---- 2. Gelecek -------------------------------------------------------------
cmip_url <- function(var, g, s, p) sprintf("%s/%s/%s/wc2.1_30s_%s_%s_%s_%s.tif", WC_CMIP6_BASE, g, s, var, g, s, p)
get_cmip <- function(var, g, s, p) {
  x <- tryCatch(crop(rast(paste0("/vsicurl/", cmip_url(var, g, s, p))), e, snap = "out"), error = function(err) NULL)
  if (!is.null(x)) return(list(r = x, how = "vsicurl"))
  x <- tryCatch(crop(geodata::cmip6_tile(lon0, lat0, model = g, ssp = sub("ssp", "", s), time = p,
                                          var = var, path = wc_raw), e, snap = "out"), error = function(err) NULL)
  if (!is.null(x)) return(list(r = x, how = "tile"))
  NULL
}
for (p in WC_PERIODS) for (s in WC_SSPS) for (g in WC_GCMS) {
  tg <- sprintf("future_%s_%s_%s", p, s, g)
  od <- file.path(CLIM_DIR, tg); dir.create(od, FALSE, TRUE)
  if (!all(file.exists(file.path(od, paste0("bio", 1:19, ".tif"))))) {
    b <- get_cmip("bioc", g, s, p)
    if (is.null(b)) { log_msg(tg, ": bioc bulunamadı"); R(tg, "bioc", "missing"); next }
    for (i in 1:19) { x <- b$r[[i]]; names(x) <- paste0("bio", i); wr(x, file.path(od, paste0("bio", i, ".tif"))) }
    log_msg(tg, ": bio1–19 (", b$how, ")")
  }
  for (i in 1:19) R(tg, paste0("bio", i), "ok")
  if (WC_MONTHLY && !file.exists(file.path(od, "tas_spr.tif"))) {
    tn <- get_cmip("tmin", g, s, p); tx <- get_cmip("tmax", g, s, p); pr <- get_cmip("prec", g, s, p)
    if (!is.null(tn) && !is.null(tx)) {
      x <- mean((tn$r[[FRUITING_MONTHS]] + tx$r[[FRUITING_MONTHS]]) / 2); names(x) <- "tas_spr"
      wr(x, file.path(od, "tas_spr.tif")); R(tg, "tas_spr", "ok")
    } else R(tg, "tas_spr", "missing")
    if (!is.null(pr)) { x <- sum(pr$r[[FRUITING_MONTHS]]); names(x) <- "pr_spr"; wr(x, file.path(od, "pr_spr.tif")); R(tg, "pr_spr", "ok") }
    else R(tg, "pr_spr", "missing")
  }
}

# Güncel tas_spr tavg'den, gelecekteki (tmin+tmax)/2'den gelir; WorldClim'de tavg zaten
# (tmin+tmax)/2 olarak üretildiğinden iki tanım tutarlıdır.
rg <- do.call(rbind, reg)
write.csv(rg, file.path(LOG_DIR, "climate_download_registry.csv"), row.names = FALSE)
fut <- unique(rg$set[grepl("^future_", rg$set)])
ok_all <- names(which(tapply(rg$status[rg$set %in% fut] == "ok", rg$var[rg$set %in% fut], all)))
writeLines(ok_all, file.path(LOG_DIR, "vars_available_in_all_futures.txt"))
log_msg("ADIM 1 (WorldClim) tamamlandı. Gelecek setleri: ", length(fut), "; eksik kayıt: ", sum(rg$status != "ok"))
