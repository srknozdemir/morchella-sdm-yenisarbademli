# =============================================================================
#  Morchella esculenta s.l. — Yenişarbademli (Isparta)
#  Yerel ölçekte (30 m) mikroiklim, topoğrafya ve toprak temelli dağılım modeli
#  Ortak ayarlar: her betik bunu source() ile okur. Yalnızca burayı düzenleyin.
# =============================================================================

# ---- 1. Klasörler ----------------------------------------------------------
if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- if (exists("SCRIPT_DIR")) dirname(SCRIPT_DIR) else dirname(getwd())
}
DATA_DIR  <- file.path(PROJECT_DIR, "data")
CLIM_DIR  <- file.path(DATA_DIR, "climate")      # CHELSA 1 km kesitleri
CLIM30_DIR<- file.path(DATA_DIR, "climate_30m")  # 30 m'ye indirgenmiş iklim
TOPO_DIR  <- file.path(DATA_DIR, "topo")         # topoğrafik ve kanopi türevleri (30 m)
SOIL_DIR  <- file.path(DATA_DIR, "soil")         # SoilGrids türevleri (30 m)
EXTRA_DIR <- file.path(DATA_DIR, "extra")        # elle eklenecek katmanlar (litoloji, meşcere)
RAW_DIR   <- file.path(DATA_DIR, "raw")          # ham veriler (DTM/DSM, DEM karoları)
OCC_DIR   <- file.path(DATA_DIR, "occurrences")
LOG_DIR   <- file.path(DATA_DIR, "logs")
for (d in c(CLIM_DIR, CLIM30_DIR, TOPO_DIR, SOIL_DIR, EXTRA_DIR, RAW_DIR, OCC_DIR, LOG_DIR))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# ---- 2. Saha kayıtları (yalnız var) ----------------------------------------
OCC_FILE     <- file.path(PROJECT_DIR, "kuzu.csv")   # Sp;Enlem;Boylam, WGS84
SPECIES_NAME <- "Morchella esculenta"
# Tekleme: her 30 m model hücresinde bir kayıt. Ek mesafe inceltmesi istenirse > 0 verin.
# Mekânsal otokorelasyon inceltmeyle değil, mekânsal blok çapraz doğrulamayla ele alınır.
THIN_DIST_M  <- 0
N_BACKGROUND <- 10000
BIAS_SIGMA_M <- 5000     # örnekleme yanlılığı çekirdeği (Gauss, m): örneklenen peyzaj ölçeği (~12 × 14 km)

# ---- 3. Alanlar ------------------------------------------------------------
# Çalışma alanı: kullanıcının sınır dosyası. Yoksa kayıtlar + LOCAL_BUFFER_KM kullanılır.
STUDY_AREA_FILE <- file.path(PROJECT_DIR, "serkan_sam_sym", "yenisar.shp")
LOCAL_BUFFER_KM <- 5
CLIM_BUFFER_KM  <- 20    # iklim indirme penceresi = alan + tampon (ölçek küçültme regresyonu için)
LOCAL_BBOX_DEFAULT <- c(xmin = 31.19, xmax = 31.47, ymin = 37.64, ymax = 37.84)
LAKE_MAX_ELEV   <- 1130  # Beyşehir Gölü yüzeyi (düz ve bu kotun altı) arka plandan çıkarılır

# ---- 4. Projeksiyon ve çözünürlük -----------------------------------------
CRS_LOCAL <- "EPSG:32636"   # WGS 84 / UTM 36N
RES_LOCAL <- 30             # model çözünürlüğü (m)

# HGM ham verileri (02a_prepare_terrain.R bunları birleştirip aşağıdaki dtm/dsm dosyalarını üretir)
HGM_SAM_PATH <- file.path(PROJECT_DIR, "yanisar_sam")                 # SAM 12 m, ArcGRID klasörü
HGM_SYM_DIRS <- c(file.path(PROJECT_DIR, "serkan_sam_sym", "sym_5m"),  # SYM 5 m paftaları
                  file.path(PROJECT_DIR, "sym_5m"))
COREG_MAX_SHIFT_M <- 15  # DTM'nin DSM'ye yatay hizalanmasında aranan en büyük kayma (m)
DSM_RES <- 10             # SYM 5 m paftaları bu çözünürlüğe ortalama ile toplanır (bellek için)

# Yüksek çözünürlüklü arazi modelleri (varsa). Herhangi bir terra okunabilir biçim.
# DTM: çıplak arazi; DSM: yüzey (kanopi dahil). İkisi birlikte kanopi yüksekliği verir.
DTM_FILE <- file.path(RAW_DIR, "dtm.tif")
DSM_FILE <- file.path(RAW_DIR, "dsm.tif")
# DTM yoksa Copernicus GLO-30 kullanılır (dikkat: GLO-30 bir yüzey modelidir).
TOPO_CALC_RES <- 10   # türevlerin hesaplandığı en ince çözünürlük; DTM daha kabaysa (ör. 12 m) DTM'ninki kullanılır
CANOPY_MIN_H  <- 3    # kanopi örtüsü eşiği (m); DTM DSM'den kabaysa 2 m yerine 3 m daha güvenlidir

# ---- 5. İklim: CHELSA v2.1 -------------------------------------------------
CHELSA_BASE <- c(
  "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies",
  "https://os.unil.cloud.switch.ch/chelsa02/chelsa/global/climatologies"
)
CHELSA_BASELINE <- "1981-2010"
CHELSA_PERIODS  <- c("2041-2070", "2071-2100")
CHELSA_SSPS     <- c("ssp126", "ssp370", "ssp585")
CHELSA_GCMS     <- c("GFDL-ESM4", "IPSL-CM6A-LR", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
CHELSA_VARS <- c(paste0("bio", 1:19), "gdd5", "gsl", "gsp", "gst", "scd", "swe", "fcf")
CHELSA_MONTHLY  <- TRUE
FRUITING_MONTHS <- 3:5

# Ölçek küçültmede yükseltiyle ilişkilendirilecek değişkenler.
# "lin": doğrusal; "log": log dönüşümlü (yağış ve kar); diğerleri bilineer enterpolasyon.
DOWNSCALE_RULES <- list(
  lin = c("bio1", "bio5", "bio6", "bio8", "bio9", "bio10", "bio11", "tas_spr", "gst",
          "gdd5", "gsl", "fcf", "bio2", "bio3", "bio4", "bio7"),
  log = c("bio12", "bio13", "bio14", "bio16", "bio17", "bio18", "bio19", "pr_spr", "gsp", "scd", "swe")
)
# Regresyon yerel model alanı (+3 km) içindeki CHELSA hücreleri üzerinden tek eşitlikle kurulur.

# ---- 5b. İklim: WorldClim 2.1 (kullanılan kaynak) ---------------------------
CLIM_SOURCE   <- "worldclim"
WC_BASELINE   <- "1970-2000"
WC_PERIODS    <- c("2041-2060", "2081-2100")
WC_SSPS       <- c("ssp126", "ssp245", "ssp370", "ssp585")
WC_GCMS       <- c("IPSL-CM6A-LR", "MIROC6", "MPI-ESM1-2-HR", "MRI-ESM2-0", "UKESM1-0-LL")
WC_CMIP6_BASE <- "https://geodata.ucdavis.edu/cmip6/30s"
WC_MONTHLY    <- TRUE      # Mart–Mayıs sıcaklık ve yağışı (tas_spr, pr_spr)

# ---- 6. Toprak: SoilGrids 2.0 ---------------------------------------------
SOIL_VARS   <- c("phh2o", "soc", "clay", "sand", "silt", "cfvo", "bdod", "nitrogen", "cec")
SOIL_DEPTHS <- c("0-5cm", "5-15cm", "15-30cm")

# ---- 7. Değişken eleme ve modelleme ---------------------------------------
SEL_DIR    <- file.path(DATA_DIR, "selection")
MODEL_DIR  <- file.path(DATA_DIR, "model_inputs")
RESULT_DIR <- file.path(PROJECT_DIR, "results")
COR_MAX    <- 0.70     # Spearman |r| eşiği (Dormann et al. 2013)
VIF_MAX    <- 5
MAX_VARS   <- 10       # M3 için üst sınır (221 kayıt / 10 ≈ 22 kayıt/değişken)
# Ekolojik öncelik: aynı kümede tek değişkenli performans farkı < 0.01 ise öncelikli olan seçilir
ECO_PRIORITY <- c("tas_spr", "pr_spr", "bio1", "bio12", "scd", "gst", "gsp", "bio2", "bio7", "bio17",
                  "insol_spr", "hli", "northness", "twi", "tpi300", "canopy_cover", "chm_mean",
                  "soil_phh2o", "soil_soc", "soil_clay", "soil_cfvo")
EXCLUDE_VARS <- c("elev")          # modele girmeyecek katmanlar
# 30 m'ye anlamlı biçimde indirgenemeyen iklim değişkenleri aday dışı: yükselti regresyonu
# R² < DOWNSCALE_MIN_R2 olanlar ve yalnız enterpole edilenler (30 m deseni 1 km'nin
# yumuşatılmış kopyasıdır; yerel kalibrasyon yapılamaz)
DOWNSCALE_MIN_R2 <- 0.30
CLAMP_TO_CALIBRATION <- TRUE       # projeksiyonda değişkenler kalibrasyon aralığına sıkıştırılır
# Değişken grupları (M1 ⊂ M2 ⊂ M3)
group_of <- function(v) {
  if (grepl("^soil_", v)) return("toprak")
  if (v %in% c("chm_mean", "chm_sd", "canopy_cover")) return("kanopi")
  if (v %in% c("slope", "slope_sd", "northness", "eastness", "tri", "hli", "twi", "insol_spr",
               "tpi100", "tpi300", "tpi1000")) return("topografya")
  if (file.exists(file.path(EXTRA_DIR, paste0(v, ".tif")))) return("ek")
  "iklim"
}
MODEL_SETS <- list(M1 = c("iklim"),
                   M2 = c("iklim", "topografya"),
                   M3 = c("iklim", "topografya", "kanopi", "toprak", "ek"))
SDM_MODELS <- c("GLM", "GAM", "GBM", "RFd", "MAXNET")
ENS_MIN_AUC <- 0.60        # ansambl için blok CV eşikleri
ENS_MIN_TSS <- 0.30
N_CORES     <- 2
PA_N   <- 2000
PA_REP <- 3
# Set bazlı ayarlar: M1 (yalnız iklim) blok testinde hiçbir algoritmada eşiği geçmediği için
# aşırı uyuma yatkın ağaç yöntemleri (GBM, RFd) çıkarılır ve eşik uygulanmaz.
MODEL_OVERRIDES <- list(M1 = list(models = c("GLM", "GAM", "MAXNET"), ens_min_auc = 0.5, ens_min_tss = 0))
FUTURE_SETS <- c("M1", "M3")   # gelecek projeksiyonu yapılacak setler (M2 yalnız güncel karşılaştırma)
PROJECT_FULL_MODELS <- TRUE    # ansambl/projeksiyon tam veri modelleriyle; CV koşuları değerlendirme için

# ---- 8. Yardımcı fonksiyonlar ---------------------------------------------
suppressPackageStartupMessages(library(terra))
terraOptions(progress = 0, memfrac = 0.4)
Sys.setenv(GDAL_HTTP_MAX_RETRY = "4", GDAL_HTTP_RETRY_DELAY = "5",
           GDAL_DISABLE_READDIR_ON_OPEN = "EMPTY_DIR",
           CPL_VSIL_CURL_ALLOWED_EXTENSIONS = ".tif,.vrt,.ovr",
           GDAL_HTTP_TIMEOUT = "120")
WOPT <- list(gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2"))
for (d in c(SEL_DIR, MODEL_DIR, RESULT_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_msg <- function(..., file = file.path(LOG_DIR, "pipeline.log")) {
  line <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", paste0(...))
  cat(line, "\n"); cat(line, "\n", file = file, append = TRUE)
}

read_occ <- function(f = OCC_FILE) {
  if (!file.exists(f)) return(NULL)
  first <- readLines(f, n = 1, warn = FALSE)
  sep <- if (grepl(";", first)) ";" else ","
  d <- read.csv(f, sep = sep, dec = ".", stringsAsFactors = FALSE, check.names = FALSE,
                fileEncoding = "UTF-8-BOM")
  nm <- tolower(names(d))
  ix <- which(nm %in% c("lon", "long", "longitude", "boylam", "x"))[1]
  iy <- which(nm %in% c("lat", "latitude", "enlem", "y"))[1]
  if (is.na(ix) || is.na(iy)) stop("Kayıt dosyasında lon/lat (Boylam/Enlem) sütunları bulunamadı: ", f)
  d$lon <- as.numeric(d[[ix]]); d$lat <- as.numeric(d[[iy]])
  d[is.finite(d$lon) & is.finite(d$lat), ]
}

study_area <- function() {
  if (file.exists(STUDY_AREA_FILE)) project(vect(STUDY_AREA_FILE), "EPSG:4326") else NULL
}
occ_bbox <- function(buffer_km) {
  sa <- study_area()
  if (!is.null(sa)) {
    e <- as.vector(ext(sa)); lat <- mean(e[c("ymin", "ymax")])
    dlat <- buffer_km / 111; dlon <- buffer_km / (111 * cos(lat * pi / 180))
    return(c(xmin = e[["xmin"]] - dlon, xmax = e[["xmax"]] + dlon,
             ymin = e[["ymin"]] - dlat, ymax = e[["ymax"]] + dlat))
  }
  occ <- read_occ()
  if (is.null(occ)) {
    b <- LOCAL_BBOX_DEFAULT; lat <- mean(b[c("ymin", "ymax")])
    return(b + c(-1, 1, 0, 0) * buffer_km / (111 * cos(lat * pi / 180)) + c(0, 0, -1, 1) * buffer_km / 111)
  }
  dlat <- buffer_km / 111
  dlon <- buffer_km / (111 * cos(mean(occ$lat) * pi / 180))
  c(xmin = min(occ$lon) - dlon, xmax = max(occ$lon) + dlon,
    ymin = min(occ$lat) - dlat, ymax = max(occ$lat) + dlat)
}
local_bbox <- function() occ_bbox(if (file.exists(STUDY_AREA_FILE)) 0 else LOCAL_BUFFER_KM)
clim_bbox  <- function() occ_bbox(CLIM_BUFFER_KM + if (file.exists(STUDY_AREA_FILE)) 0 else LOCAL_BUFFER_KM)
# Sınır poligonu dışındaki hücreleri NA yapan maske (UTM şablonunda)
study_mask <- function(tmpl) {
  sa <- study_area(); if (is.null(sa)) { m <- tmpl; values(m) <- 1; return(m) }
  rasterize(project(sa, crs(tmpl)), tmpl, field = 1)
}
bbox_ext   <- function(b) ext(b[["xmin"]], b[["xmax"]], b[["ymin"]], b[["ymax"]])

# Model şablon ızgarası (UTM 36N, RES_LOCAL m, tam sayıya hizalı)
local_template <- function(res = RES_LOCAL) {
  pv <- project(vect(bbox_ext(local_bbox()), crs = "EPSG:4326"), CRS_LOCAL)
  e  <- ext(pv)
  e  <- ext(floor(e$xmin / res) * res, ceiling(e$xmax / res) * res,
            floor(e$ymin / res) * res, ceiling(e$ymax / res) * res)
  rast(e, resolution = res, crs = CRS_LOCAL)
}

.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
