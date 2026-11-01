# =============================================================================
#  ADIM 11 — Makale şekilleri (results/figures/*.png, 300 dpi)
#  Fig1 çalışma alanı ve kayıtlar; Fig2 güncel uygunluk (M3); Fig3 2071–2100 değişim
#  sınıfları (M3, üç SSP); Fig4 korunan uygun alan M1–M3; Fig5 değişken önemi ve
#  tepki eğrileri (M3)
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages(library(ggplot2))
FIG <- file.path(RESULT_DIR, "figures"); dir.create(FIG, FALSE, TRUE)
SUM_DIR <- file.path(RESULT_DIR, "summary")
th <- theme_minimal(base_size = 9) + theme(panel.grid = element_blank(), axis.title = element_blank(),
                                           legend.key.height = unit(0.4, "cm"))
to_df <- function(r, fact = 3) {
  if (fact > 1) r <- aggregate(r, fact, mean, na.rm = TRUE)
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; d
}
find1 <- function(dir, pat) { f <- list.files(dir, pattern = pat, full.names = TRUE); if (length(f)) f[1] else NA }
km_axis <- list(scale_x_continuous(labels = function(x) x / 1000), scale_y_continuous(labels = function(x) x / 1000))

# ---- Fig 1 ---------------------------------------------------------------------
z  <- rast(file.path(TOPO_DIR, "elev.tif"))
hs <- shade(terrain(z, "slope", unit = "radians"), terrain(z, "aspect", unit = "radians"), 40, 315)
occ <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv"))
lake <- ifel(z <= LAKE_MAX_ELEV & rast(file.path(TOPO_DIR, "slope.tif")) < 0.5, 1, NA)
p1 <- ggplot() +
  geom_raster(data = to_df(z, 4), aes(x, y, fill = v)) +
  scale_fill_gradientn("Elevation (m)", colours = terrain.colors(10)) +
  geom_raster(data = to_df(hs, 4), aes(x, y, alpha = v), fill = "black", show.legend = FALSE) +
  scale_alpha(range = c(0.35, 0)) +
  geom_raster(data = to_df(lake, 4), aes(x, y), fill = "#9ecae1") +
  geom_point(data = occ, aes(x, y), shape = 21, size = 0.9, fill = "red", colour = "black", stroke = 0.2) +
  coord_equal() + km_axis + th + labs(caption = "UTM 36N (km)")
ggsave(file.path(FIG, "Fig1_study_area.png"), p1, width = 16, height = 13, units = "cm", dpi = 300, bg = "white")

# ---- Fig 2: güncel uygunluk -------------------------------------------------------
rd3 <- file.path(RESULT_DIR, "M3", "proj")
fc <- file.path(rd3, "current_suitability.tif"); if (!file.exists(fc)) fc <- NA
if (!is.na(fc)) {
  s0 <- rast(fc)
  p2 <- ggplot() + geom_raster(data = to_df(s0), aes(x, y, fill = v)) +
    scale_fill_viridis_c("Suitability", limits = c(0, 1)) +
    geom_point(data = occ, aes(x, y), size = 0.3, colour = "white") +
    coord_equal() + km_axis + th
  ggsave(file.path(FIG, "Fig2_current_suitability_M3.png"), p2, width = 16, height = 13, units = "cm", dpi = 300, bg = "white")
}

# ---- Fig 3: değişim sınıfları 2071–2100 ------------------------------------------
fb <- file.path(rd3, "current_binary.tif"); if (!file.exists(fb)) fb <- NA
pa <- list.files(rd3, pattern = "^gcm_2071-2100_.*_presence_agreement\\.tif$", full.names = TRUE)
if (!is.na(fb) && length(pa)) {
  b0 <- rast(fb)
  dd <- do.call(rbind, lapply(pa, function(f) {
    ssp <- sub("^gcm_2071-2100_(ssp[0-9]+)_.*$", "\\1", basename(f))
    fu <- rast(f) >= 0.8
    cl <- ifel(b0 == 1 & fu, 2, ifel(b0 == 1, 1, ifel(fu, 3, NA)))
    d <- as.data.frame(aggregate(cl, 3, "modal", na.rm = TRUE), xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"
    d$ssp <- toupper(sub("ssp(\\d)(\\d)(\\d)", "SSP\\1-\\2.\\3", ssp)); d
  }))
  dd$v <- factor(dd$v, levels = 1:3, labels = c("Lost", "Retained (refugia)", "Gained"))
  landdf <- as.data.frame(aggregate(!is.na(b0), 3, "max"), xy = TRUE); landdf <- landdf[landdf[[3]] == 1, 1:2]
  p3 <- ggplot(dd, aes(x, y, fill = v)) + geom_raster(data = landdf, aes(x, y), fill = "grey92", inherit.aes = FALSE) +
    geom_raster() + facet_wrap(~ssp, nrow = 1) +
    scale_fill_manual(NULL, values = c("#d73027", "#1a9850", "#4575b4")) +
    coord_equal() + km_axis + th + theme(legend.position = "bottom")
  ggsave(file.path(FIG, "Fig3_change_2071_2100_M3.png"), p3, width = 18, height = 8.5, units = "cm", dpi = 300, bg = "white")
}

# ---- Fig 4: korunan alan M1 vs M3 -------------------------------------------------
fa <- file.path(SUM_DIR, "suitable_area_change.csv")
if (file.exists(fa)) {
  a <- read.csv(fa); a <- a[a$ssp != "current", ]
  a$ssp <- sub("ssp(\\d)(\\d)(\\d)", "SSP\\1-\\2.\\3", a$ssp)
  a$model_set <- factor(a$model_set, levels = c("M1", "M3"), labels = c("Climate only (M1)", "Full model (M3)"))
  p4 <- ggplot(a, aes(ssp, retained_pct, fill = model_set)) + geom_col(position = position_dodge(0.75), width = 0.7) +
    facet_wrap(~period) + scale_fill_manual(NULL, values = c("grey60", "#1a9850")) +
    labs(x = NULL, y = "Currently suitable habitat retained (%)") +
    theme_minimal(base_size = 9) + theme(legend.position = "bottom")
  ggsave(file.path(FIG, "Fig4_retained_habitat.png"), p4, width = 14, height = 9, units = "cm", dpi = 300, bg = "white")
}

# ---- Fig 5: değişken önemi ve tepki eğrileri (M3) ----------------------------------
vi <- file.path(RESULT_DIR, "M3", "var_importance_ensemble_avg.csv")
rc <- file.path(RESULT_DIR, "M3", "response_curves", "ensemble_response_curves_data.csv")
if (file.exists(vi)) {
  v <- read.csv(vi); nm <- names(v)
  vc <- nm[grepl("var|expl", nm, ignore.case = TRUE)][1]; ic <- nm[sapply(v, is.numeric)][1]
  v <- v[order(v[[ic]]), ]; v[[vc]] <- factor(v[[vc]], levels = v[[vc]])
  p5 <- ggplot(v, aes(.data[[ic]], .data[[vc]])) + geom_col(fill = "#1a9850") +
    labs(x = "Variable importance", y = NULL) + theme_minimal(base_size = 9)
  ggsave(file.path(FIG, "Fig5_variable_importance_M3.png"), p5, width = 9, height = 8, units = "cm", dpi = 300, bg = "white")
}
log_msg("ADIM 11 tamamlandı: ", FIG)
