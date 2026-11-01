# =============================================================================
#  ADIM 12 — Makale şekilleri (yayın biçimi) ve tablolar
#  Fig1 çalışma alanı + konum; Fig2 güncel uygunluk M1/M3; Fig3 önem + tepki eğrileri;
#  Fig4 gelecek uygunluk (GCM ortalaması); Fig5 değişim sınıfları; Fig6 korunan alan,
#  mikrosığınak arazisi, belirsizlik. Ek: FigS1 değişken taraması, FigS2 GCM SD, FigS3 MESS.
#  Çıktı: results/figures_ms (PNG 300 dpi + PDF), results/tables_ms (CSV)
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
invisible(Sys.setlocale("LC_CTYPE", "C.UTF-8"))
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(biomod2); library(mgcv) })
FIG <- file.path(RESULT_DIR, "figures_ms"); dir.create(FIG, FALSE, TRUE)
TAB <- file.path(RESULT_DIR, "tables_ms"); dir.create(TAB, FALSE, TRUE)
SUM_DIR <- file.path(RESULT_DIR, "summary")
rd3 <- file.path(RESULT_DIR, "M3", "proj"); rd1 <- file.path(RESULT_DIR, "M1", "proj")

save_fig <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white", device = ragg::agg_png)
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, units = "mm", bg = "white", device = cairo_pdf)
}
base <- theme_bw(base_size = 8) + theme(panel.grid = element_blank(), strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 8), plot.tag = element_text(face = "bold", size = 9),
  legend.key.size = unit(3, "mm"), legend.title = element_text(size = 7), legend.text = element_text(size = 7))
map_th <- base + theme(axis.title = element_blank(), axis.text = element_text(size = 6))
km <- list(scale_x_continuous(labels = function(x) x / 1000, expand = c(0, 0)),
           scale_y_continuous(labels = function(x) x / 1000, expand = c(0, 0)))
to_df <- function(r, fact = 3, fun = "mean") {
  if (fact > 1) r <- aggregate(r, fact, fun, na.rm = TRUE)
  d <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(d)[3] <- "v"; d
}
ssp_lab <- function(s) sub("ssp(\\d)(\\d)(\\d)", "SSP\\1-\\2.\\3", s)
lab_var <- c(gsp = "GSP (mm)", insol_spr = "Radiation (kWh m⁻²)",
             tpi100 = "TPI 100 m (m)", canopy_cover = "Canopy cover", soil_phh2o = "pH (H₂O)",
             soil_clay = "Clay (%)", soil_cfvo = "Coarse fragments (%)", soil_bdod = "Bulk density (g cm⁻³)",
             soil_cec = "CEC (cmol(+) kg⁻¹)")
lab_long <- c(gsp = "Growing-season precipitation", insol_spr = "Spring direct radiation", tpi100 = "Topographic position (100 m)",
              canopy_cover = "Canopy cover", soil_phh2o = "Soil pH", soil_clay = "Clay", soil_cfvo = "Coarse fragments",
              soil_bdod = "Bulk density", soil_cec = "Cation exchange capacity")
grp_col <- c(Climate = "#3b6fb6", Terrain = "#8c6d46", Canopy = "#1a9850", Soil = "#d95f02")
grp_of <- function(v) ifelse(v == "gsp", "Climate", ifelse(v %in% c("insol_spr", "tpi100"), "Terrain",
                       ifelse(v == "canopy_cover", "Canopy", "Soil")))

z <- rast(file.path(TOPO_DIR, "elev.tif")); sl <- rast(file.path(TOPO_DIR, "slope.tif"))
sa <- project(study_area(), crs(z))
lake <- ifel(z <= LAKE_MAX_ELEV & sl < 0.5, 1, NA)
land <- mask(ifel(is.na(lake), 1, NA), sa)
occ <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv"))
sa_df <- geom(sa, df = TRUE)
outline <- geom_polygon(data = sa_df, aes(x, y, group = interaction(part, hole)), fill = NA, colour = "grey25", linewidth = 0.25)
lake_df <- to_df(lake, 4)
lake_l <- geom_raster(data = lake_df, aes(x, y), fill = "#c6dbef", inherit.aes = FALSE)
land_l <- geom_raster(data = to_df(land, 4), aes(x, y), fill = "grey90", inherit.aes = FALSE)
e <- ext(z)
sbar <- function(len = 10000) list(
  annotate("segment", x = e$xmax - 4000 - len, xend = e$xmax - 4000, y = e$ymin + 4000, yend = e$ymin + 4000, linewidth = 0.6),
  annotate("text", x = e$xmax - 4000 - len / 2, y = e$ymin + 6500, label = paste(len / 1000, "km"), size = 2.2))

# ---- Fig 1 ------------------------------------------------------------------------
hs <- shade(terrain(z, "slope", unit = "radians"), terrain(z, "aspect", unit = "radians"), 40, 315)
zd <- to_df(mask(z, sa), 4); hd <- to_df(mask(hs, sa), 4)
p1 <- ggplot() +
  geom_raster(data = zd, aes(x, y, fill = v)) +
  scale_fill_gradientn("Elevation (m)", colours = c("#2c7a3f", "#8fbf6a", "#e9e3a6", "#c9a36b", "#8c6d46", "#f2f2f2"),
                       breaks = c(1000, 1500, 2000, 2500)) +
  geom_raster(data = hd, aes(x, y, alpha = v), fill = "black", show.legend = FALSE) + scale_alpha(range = c(0.45, 0)) +
  geom_raster(data = lake_df, aes(x, y), fill = "#9ecae1") + outline +
  geom_point(data = occ, aes(x, y), shape = 21, size = 1, fill = "#e41a1c", colour = "black", stroke = 0.2) +
  annotate("text", x = 372000, y = 4192000, label = "Lake\nBeyşehir", size = 2.3, fontface = "italic", colour = "#08519c") +
  annotate("segment", x = e$xmin + 5000, xend = e$xmin + 5000, y = e$ymax - 12000, yend = e$ymax - 5000,
           arrow = arrow(length = unit(1.5, "mm")), linewidth = 0.4) +
  annotate("text", x = e$xmin + 5000, y = e$ymax - 3500, label = "N", size = 2.6, fontface = "bold") +
  sbar() + coord_equal() + km + map_th + labs(caption = "UTM zone 36N (km)") +
  theme(legend.position = c(0.99, 0.99), legend.justification = c(1, 1), legend.background = element_rect(fill = alpha("white", 0.8), colour = NA))
tr <- map_data("world", region = "Turkey")
sa_ll <- geom(project(sa, "EPSG:4326"), df = TRUE)
inset <- ggplot() + geom_polygon(data = tr, aes(long, lat, group = group), fill = "grey85", colour = "grey40", linewidth = 0.2) +
  geom_polygon(data = sa_ll, aes(x, y, group = part), fill = "#e41a1c", colour = "#e41a1c", linewidth = 0.3) +
  coord_quickmap(xlim = c(25.6, 44.9), ylim = c(35.8, 42.2), expand = FALSE) +
  theme_void() + theme(panel.border = element_rect(fill = NA, colour = "grey30", linewidth = 0.3),
                       panel.background = element_rect(fill = "white", colour = NA))
f1 <- p1 + inset_element(inset, left = 0.0, bottom = 0.0, right = 0.36, top = 0.2, align_to = "panel")
save_fig(f1, "Fig1_study_area", 174, 160)

# ---- Fig 2: güncel uygunluk M1 ve M3 ------------------------------------------------
cur <- function(rd, lab) { d <- to_df(rast(file.path(rd, "current_suitability.tif"))); d$m <- lab; d }
a1 <- read.csv(file.path(SUM_DIR, "suitable_area_change.csv"))
cur_km <- function(M) round(a1$suitable_km2[a1$model_set == M & a1$ssp == "current"])
d2 <- rbind(cur(rd1, sprintf("Climate only (M1): %d km² suitable", cur_km("M1"))),
            cur(rd3, sprintf("Full model (M3): %d km² suitable", cur_km("M3"))))
p2 <- ggplot(d2, aes(x, y, fill = v)) + lake_l + geom_raster() + outline + facet_wrap(~m) +
  scale_fill_viridis_c("Habitat\nsuitability", limits = c(0, 1), option = "D") +
  coord_equal() + km + map_th + theme(legend.position = "right")
save_fig(p2, "Fig2_current_suitability", 174, 80)

# ---- Fig 3: önem ve tepki eğrileri (M3) ----------------------------------------------
sp <- gsub(" ", ".", SPECIES_NAME)
mdir <- file.path(RESULT_DIR, "M3", "biomod_wd", sp, "models", "binary")
algos <- read.csv(file.path(RESULT_DIR, "M3", "ensemble_candidates_after_thresholds.csv"))$model
w <- read.csv(file.path(rd3, "ensemble_weights.csv")); w_alg <- setNames(w$weight, w$algo)
files <- unlist(lapply(algos, function(a) list.files(mdir, pattern = paste0("^", sp, "_PA[0-9]+_allRun_", a, "$"), full.names = TRUE)))
mods <- lapply(files, function(f) get(load(f))); m_alg <- sub(".*_", "", basename(files))
vars <- sub("\\.tif$", "", list.files(file.path(MODEL_DIR, "M3", "current"), "\\.tif$"))
env0 <- rast(file.path(MODEL_DIR, "M3", "current", paste0(vars, ".tif"))); names(env0) <- vars
bgx <- read.csv(file.path(OCC_DIR, "bg_local_bias.csv"))
pres <- as.data.frame(extract(env0, as.matrix(occ[, c("x", "y")]))); bgv <- as.data.frame(extract(env0, as.matrix(bgx[, c("x", "y")])))
calv <- na.omit(rbind(pres, bgv)); med <- sapply(calv, median)
thr <- read.csv(file.path(rd3, "threshold_current.csv"))$threshold_maxTSS
rc <- do.call(rbind, lapply(vars, function(v) {
  g <- seq(quantile(calv[[v]], 0.01), quantile(calv[[v]], 0.99), length.out = 120)
  nd <- as.data.frame(t(replicate(length(g), med))); nd[[v]] <- g
  P <- sapply(mods, function(m) as.numeric(predict(m, newdata = nd, on_0_1000 = FALSE)))
  A <- sapply(names(w_alg), function(a) rowMeans(P[, m_alg == a, drop = FALSE]))
  data.frame(var = v, x = g, ens = as.numeric(A %*% w_alg), lo = apply(A, 1, min), hi = apply(A, 1, max))
}))
write.csv(rc, file.path(TAB, "response_curves_M3.csv"), row.names = FALSE)
vi <- read.csv(file.path(RESULT_DIR, "M3", "var_importance_ensemble_avg.csv"))
vi$grp <- grp_of(vi$variable); vi$lab <- lab_long[vi$variable]
vi$lab <- factor(vi$lab, levels = vi$lab[order(vi$importance_percent)])
pa <- ggplot(vi, aes(importance_percent, lab, fill = grp)) + geom_col(width = 0.7) +
  scale_fill_manual(NULL, values = grp_col) + scale_x_continuous(expand = expansion(c(0, 0.05))) +
  labs(x = "Relative importance (%)", y = NULL) + base + theme(legend.position = c(0.72, 0.22), legend.background = element_blank())
ord <- as.character(vi$variable[order(-vi$importance_percent)])
rc$var <- factor(rc$var, levels = ord, labels = lab_var[ord])
rug <- do.call(rbind, lapply(vars, function(v) data.frame(var = v, x = pres[[v]])))
rug <- rug[!is.na(rug$x), ]; rug$var <- factor(rug$var, levels = ord, labels = lab_var[ord])
rug <- do.call(rbind, lapply(split(rug, rug$var), function(d) { r <- range(rc$x[rc$var == d$var[1]]); d[d$x >= r[1] & d$x <= r[2], ] }))
pb <- ggplot(rc, aes(x, ens)) + geom_ribbon(aes(ymin = lo, ymax = hi), fill = "grey80") +
  geom_hline(yintercept = thr, linetype = 2, linewidth = 0.3, colour = "grey40") +
  geom_line(linewidth = 0.5, colour = "#1a9850") + geom_rug(data = rug, aes(x = x), inherit.aes = FALSE, sides = "b", linewidth = 0.15, alpha = 0.5) +
  facet_wrap(~var, scales = "free_x", ncol = 3) + scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Predicted suitability") + base + theme(strip.text = element_text(size = 6.5, face = "bold"), axis.text = element_text(size = 6), panel.spacing.x = unit(4, "mm"))
f3 <- pa + pb + plot_layout(widths = c(0.9, 2.8)) + plot_annotation(tag_levels = "a")
save_fig(f3, "Fig3_importance_response", 174, 118)

# ---- Fig 4: gelecek uygunluk (GCM ortalaması, M3) -------------------------------------
fm <- list.files(rd3, "^gcm_.*_mean\\.tif$", full.names = TRUE)
d4 <- do.call(rbind, lapply(fm, function(f) {
  k <- regmatches(basename(f), regexec("gcm_(.+?)_(ssp[0-9]+)_mean", basename(f)))[[1]]
  d <- to_df(rast(f)); d$period <- k[2]; d$ssp <- ssp_lab(k[3]); d }))
p4 <- ggplot(d4, aes(x, y, fill = v)) + lake_l + geom_raster() + outline + facet_grid(period ~ ssp) +
  scale_fill_viridis_c("Habitat\nsuitability\n(GCM mean)", limits = c(0, 1)) + coord_equal() + km + map_th
save_fig(p4, "Fig4_future_suitability", 174, 112)

# ---- Fig 5: değişim sınıfları (her iki dönem) --------------------------------------------
b0 <- rast(file.path(rd3, "current_binary.tif"))
pa_f <- list.files(rd3, "^gcm_.*_presence_agreement\\.tif$", full.names = TRUE)
d5 <- do.call(rbind, lapply(pa_f, function(f) {
  k <- regmatches(basename(f), regexec("gcm_(.+?)_(ssp[0-9]+)_presence", basename(f)))[[1]]
  fu <- rast(f) >= 0.8
  cl <- ifel(b0 == 1 & fu, 2, ifel(b0 == 1, 1, ifel(fu, 3, NA)))
  d <- to_df(cl, 3, "modal"); d$period <- k[2]; d$ssp <- ssp_lab(k[3]); d }))
d5$v <- factor(round(d5$v), levels = 1:3, labels = c("Lost", "Retained", "Gained"))
p5 <- ggplot(d5, aes(x, y, fill = v)) + land_l + lake_l + geom_raster() + outline + facet_grid(period ~ ssp) +
  scale_fill_manual(NULL, values = c(Lost = "#D55E00", Retained = "#0072B2", Gained = "#E69F00")) +
  coord_equal() + km + map_th + theme(legend.position = "bottom")
save_fig(p5, "Fig5_change_classes", 174, 120)

# ---- Fig 6: korunan alan, mikrosığınak arazisi, belirsizlik --------------------------
a <- a1[a1$ssp != "current", ]; a$ssp <- sub("SSP", "", ssp_lab(a$ssp))
a$model_set <- factor(a$model_set, c("M1", "M3"), c("Climate only (M1)", "Full model (M3)"))
f6a <- ggplot(a, aes(ssp, retained_pct, fill = model_set)) + geom_col(position = position_dodge(0.75), width = 0.7) +
  facet_wrap(~period) + scale_fill_manual(NULL, values = c("grey65", "#0072B2")) +
  scale_y_continuous(expand = expansion(c(0, 0.05))) + labs(x = "SSP", y = "Current habitat retained (%)") +
  base + theme(legend.position = "bottom", axis.text.x = element_text(size = 6))
ag <- rast(file.path(rd3, "gcm_2071-2100_ssp370_presence_agreement.tif"))
tv <- c(northness = "Northness", slope = "Slope (°)", twi = "TWI", elev = "Elevation (m)")
st <- rast(file.path(TOPO_DIR, paste0(names(tv), ".tif"))); names(st) <- names(tv)
cls <- ifel(b0 == 1 & ag >= 0.8, 1, ifel(b0 == 1, 2, NA))
set.seed(4); cs <- spatSample(c(cls, st), 6000, "random", na.rm = TRUE); names(cs)[1] <- "cl"
cs$cl <- factor(cs$cl, 1:2, c("Refugia", "Lost"))
cl <- do.call(rbind, lapply(names(tv), function(v) data.frame(var = tv[[v]], cl = cs$cl, val = cs[[v]])))
cl$var <- factor(cl$var, levels = tv)
f6b <- ggplot(cl, aes(cl, val, fill = cl)) + geom_boxplot(outlier.shape = NA, width = 0.6, linewidth = 0.3) +
  facet_wrap(~var, scales = "free_y", nrow = 1) + scale_fill_manual(NULL, values = c(Refugia = "#0072B2", Lost = "#D55E00")) +
  labs(x = NULL, y = NULL) + base + theme(legend.position = "none")
up <- read.csv(file.path(SUM_DIR, "uncertainty_partitioning_2071_2100.csv"))
up$component <- factor(up$component, levels = rev(up$component),
  labels = rev(c("Algorithm", "GCM", "SSP", "Algorithm × GCM", "Algorithm × SSP", "GCM × SSP", "Residual")))
f6c <- ggplot(up, aes(mean_share * 100, component)) + geom_col(fill = "grey45", width = 0.7) +
  geom_errorbarh(aes(xmin = pmax(mean_share - sd_share, 0) * 100, xmax = pmin(mean_share + sd_share, 1) * 100), height = 0.25, linewidth = 0.3) +
  labs(x = "Share of variance (%)", y = NULL) + scale_x_continuous(expand = expansion(c(0, 0.05))) + base
f6 <- (f6a | f6c) / f6b + plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "a")
save_fig(f6, "Fig6_refugia_uncertainty", 174, 130)

# ---- Ek şekiller -------------------------------------------------------------------------
vs <- read.csv(file.path(DATA_DIR, "selection", "variable_selection.csv"))
grp_en <- c(iklim = "Climate", topografya = "Terrain", kanopi = "Canopy", toprak = "Soil", ek = "Terrain")
vs$grp <- grp_en[vs$group]
vs$status <- ifelse(vs$selected, "Selected", "Not selected")
vs <- vs[!is.na(vs$auc_block), ]; vs$var <- factor(vs$var, levels = vs$var[order(vs$auc_block)])
s1 <- ggplot(vs, aes(auc_block, var, colour = grp, shape = status)) + geom_vline(xintercept = 0.5, linetype = 3, linewidth = 0.3) +
  geom_point(size = 1.6) + scale_shape_manual(NULL, values = c(Selected = 16, `Not selected` = 1)) +
  scale_colour_manual(NULL, values = grp_col) + labs(x = "Univariate spatial-block AUC", y = NULL) +
  base + theme(axis.text.y = element_text(size = 5.5))
save_fig(s1, "FigS1_variable_screening", 120, 150)

fsd <- list.files(rd3, "^gcm_2071-2100_.*_sd\\.tif$", full.names = TRUE)
ds <- do.call(rbind, lapply(fsd, function(f) { d <- to_df(rast(f)); d$ssp <- ssp_lab(sub(".*_(ssp[0-9]+)_sd.*", "\\1", basename(f))); d }))
s2 <- ggplot(ds, aes(x, y, fill = v)) + lake_l + geom_raster() + outline + facet_wrap(~ssp) +
  scale_fill_viridis_c("SD among\nGCMs", option = "magma", direction = -1) + coord_equal() + km + map_th
save_fig(s2, "FigS2_gcm_sd_2071_2100", 174, 62)

mf <- list.files(file.path(SUM_DIR, "mess"), "2071-2100.*\\.tif$", full.names = TRUE)
if (length(mf)) {
  dm <- do.call(rbind, lapply(split(mf, sub(".*_(ssp[0-9]+)_.*", "\\1", basename(mf))), function(ff) {
    r <- app(mask(rast(ff), sa), mean); d <- to_df(r); d$ssp <- ssp_lab(sub(".*_(ssp[0-9]+)_.*", "\\1", basename(ff[1]))); d }))
  dm$v <- pmax(dm$v, -50)
  s3 <- ggplot(dm, aes(x, y, fill = v)) + geom_raster() + outline + facet_wrap(~ssp) +
    scale_fill_gradient2("MESS", low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0) +
    coord_equal() + km + map_th
  save_fig(s3, "FigS3_mess_2071_2100", 174, 62)
}

# ---- Tablolar -------------------------------------------------------------------------
sel <- vs[vs$selected, c("var", "grp", "auc_block", "vif")]
sel$importance <- vi$importance_percent[match(sel$var, vi$variable)]
write.csv(sel[order(-sel$importance), ], file.path(TAB, "Table1_selected_predictors.csv"), row.names = FALSE)
write.csv(vs[order(vs$grp, -vs$auc_block), c("var", "grp", "auc_block", "cluster", "selected", "vif")],
          file.path(TAB, "TableS1_candidate_predictors.csv"), row.names = FALSE)
log_msg("ADIM 12 tamamlandı: ", FIG)
