# =============================================================================
#  ADIM 10 — Model sonrası analizler
#  1. Performans özeti: M1–M3 için algoritma bazında blok CV AUC, TSS, Boyce
#  2. Ekstrapolasyon: MESS (Elith et al. 2010) her gelecek seti için; MESS < 0
#     yeni (kalibrasyon dışı) iklim koşullarını gösterir
#  3. Uygun alan ve değişim: güncel ve gelecek (GCM uyumu ≥ %80) uygun alan, km²
#  4. Topoğrafik tamponlama: güncel uygun alanın ne kadarının korunduğu,
#     M1 (yalnız iklim) ile M3 (iklim + mikrotopoğrafya + kanopi + toprak) arasında
#  5. Mikrosığınakların topoğrafik özellikleri (M3): korunan ve kaybedilen hücreler
#  6. Belirsizlik ayrıştırma (M3, 2071–2100): algoritma, GCM ve SSP'nin
#     uygunluk varyansındaki payı (Buisson et al. 2010)
#  Çıktılar: results/summary/*.csv, results/summary/*.tif
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
set.seed(11)
SUM_DIR <- file.path(RESULT_DIR, "summary"); dir.create(SUM_DIR, FALSE, TRUE)
AGREE <- 0.8
cell_km2 <- (RES_LOCAL / 1000)^2

# ---- 1. Performans ----------------------------------------------------------
perf <- do.call(rbind, lapply(names(MODEL_SETS), function(M) {
  f <- file.path(RESULT_DIR, M, "eval_model.csv"); if (!file.exists(f)) return(NULL)
  d <- read.csv(f); d <- d[d$value.run != "allRun" & !is.na(d$value.validation), ]
  d <- d[!duplicated(d[, c("value.full.name", "value.metric.eval")]), ]
  a <- aggregate(value.validation ~ value.algo + value.metric.eval, data = d,
                 FUN = function(x) c(mean = mean(x), sd = sd(x)))
  data.frame(model_set = M, algo = a$value.algo, metric = a$value.metric.eval,
             mean = round(a$value.validation[, "mean"], 3), sd = round(a$value.validation[, "sd"], 3))
}))
if (!is.null(perf)) write.csv(perf, file.path(SUM_DIR, "performance_blockCV.csv"), row.names = FALSE)

# Model setleri arası eşleştirilmiş karşılaştırma: aynı algoritma, PA seti ve blok
pair_tab <- NULL
ev_all <- lapply(names(MODEL_SETS), function(M) {
  f <- file.path(RESULT_DIR, M, "eval_model.csv"); if (!file.exists(f)) return(NULL)
  d <- read.csv(f); d <- d[d$value.run != "allRun" & !is.na(d$value.validation), ]
  d <- d[!duplicated(d[, c("value.full.name", "value.metric.eval")]), ]
  data.frame(set = M, algo = d$value.algo, pa = d$value.PA, run = d$value.run, metric = d$value.metric.eval, v = d$value.validation)
})
ev_all <- do.call(rbind, ev_all)
if (!is.null(ev_all)) for (cmp in list(c("M1", "M2"), c("M2", "M3"), c("M1", "M3"))) for (mt in c("AUCroc", "TSS", "BOYCE")) {
  a <- ev_all[ev_all$set == cmp[1] & ev_all$metric == mt, ]; b <- ev_all[ev_all$set == cmp[2] & ev_all$metric == mt, ]
  m <- merge(a, b, by = c("algo", "pa", "run"), suffixes = c("_a", "_b"))
  m <- m[m$algo %in% c("GLM", "GAM", "MAXNET"), ]
  if (nrow(m) < 5) next
  wt <- suppressWarnings(wilcox.test(m$v_b, m$v_a, paired = TRUE))
  pair_tab <- rbind(pair_tab, data.frame(comparison = paste(cmp, collapse = " vs "), metric = mt, n_pairs = nrow(m),
    mean_a = round(mean(m$v_a), 3), mean_b = round(mean(m$v_b), 3), mean_diff = round(mean(m$v_b - m$v_a), 3),
    p_wilcoxon = signif(wt$p.value, 3)))
}
if (!is.null(pair_tab)) write.csv(pair_tab, file.path(SUM_DIR, "performance_paired_comparison.csv"), row.names = FALSE)

# Değişken grubu katkısı (M3, ansambl değişken önemi)
vi <- file.path(RESULT_DIR, "M3", "var_importance_ensemble_avg.csv")
if (file.exists(vi)) {
  v <- read.csv(vi); v$group <- vapply(v$variable, group_of, "")
  g <- aggregate(importance_percent ~ group, v, sum); g$importance_percent <- round(g$importance_percent, 1)
  write.csv(g[order(-g$importance_percent), ], file.path(SUM_DIR, "importance_by_group_M3.csv"), row.names = FALSE)
}

# ---- yardımcılar -------------------------------------------------------------
cur_bin  <- function(M) { f <- file.path(RESULT_DIR, M, "proj", "current_binary.tif"); if (file.exists(f)) f else NA }
cur_cont <- function(M) { f <- file.path(RESULT_DIR, M, "proj", "current_suitability.tif"); if (file.exists(f)) f else NA }
tmpl <- local_template(); smask <- study_mask(tmpl)

# ---- 2. MESS ------------------------------------------------------------------
mess_fun <- function(X, ref) {
  out <- rep(Inf, nrow(X))
  for (j in seq_len(ncol(X))) {
    r <- ref[, j]; mn <- min(r); mx <- max(r); Fn <- ecdf(r)
    p <- X[, j]; f <- Fn(p) * 100
    s <- ifelse(f == 0, (p - mn) / (mx - mn) * 100,
         ifelse(f <= 50, 2 * f, ifelse(f < 100, 2 * (100 - f), (mx - p) / (mx - mn) * 100)))
    out <- pmin(out, s)
  }
  out
}
mess_tab <- NULL
M <- "M3"
in_dir <- file.path(MODEL_DIR, M)
if (dir.exists(in_dir)) {
  vars <- sub("\\.tif$", "", list.files(file.path(in_dir, "current"), "\\.tif$"))
  bg <- read.csv(file.path(OCC_DIR, "bg_local_bias.csv")); pr <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv"))
  cur <- rast(file.path(in_dir, "current", paste0(vars, ".tif")))
  ref <- na.omit(as.matrix(extract(cur, as.matrix(rbind(bg[, c("x", "y")], pr[, c("x", "y")])))))
  fsets <- grep("^future_", list.dirs(in_dir, recursive = FALSE, full.names = FALSE), value = TRUE)
  mess_dir <- file.path(SUM_DIR, "mess"); dir.create(mess_dir, FALSE)
  for (fs in fsets) {
    fr <- rast(file.path(in_dir, fs, paste0(vars, ".tif")))
    out <- file.path(mess_dir, paste0("mess_", fs, ".tif"))
    if (!file.exists(out)) {
      v <- values(fr); r <- rep(NA_real_, nrow(v)); ok <- complete.cases(v)   # vektörel (hücre başına çağrı yok)
      r[ok] <- mess_fun(v[ok, , drop = FALSE], ref); rm(v)
      m <- rast(fr[[1]]); values(m) <- r; names(m) <- "mess"; rm(r)
      m <- writeRaster(m, out, overwrite = TRUE, gdal = "COMPRESS=DEFLATE"); invisible(gc())
    } else m <- rast(out)
    m <- mask(m, smask)
    info <- regmatches(fs, regexec("^future_(.+?)_(ssp[0-9]{3})_(.+)$", fs))[[1]]
    novel <- global(m < 0, "mean", na.rm = TRUE)[1, 1]
    mess_tab <- rbind(mess_tab, data.frame(period = info[2], ssp = info[3], gcm = info[4], novel_share = round(novel, 4)))
  }
  if (!is.null(mess_tab)) {
    write.csv(mess_tab, file.path(SUM_DIR, "mess_by_gcm.csv"), row.names = FALSE)
    agg <- aggregate(novel_share ~ period + ssp, mess_tab, function(x) c(mean = mean(x), min = min(x), max = max(x)))
    write.csv(do.call(data.frame, agg), file.path(SUM_DIR, "mess_summary.csv"), row.names = FALSE)
  }
}

# ---- 3–4. Uygun alan, değişim ve tamponlama ----------------------------------
area_tab <- NULL; keep_rasters <- list()
for (M in intersect(FUTURE_SETS, names(MODEL_SETS))) {
  rd <- file.path(RESULT_DIR, M, "proj"); fb <- cur_bin(M); if (is.na(fb)) next
  b0 <- mask(rast(fb), smask); a0 <- global(b0 == 1, "sum", na.rm = TRUE)[1, 1] * cell_km2
  area_tab <- rbind(area_tab, data.frame(model_set = M, period = "1981-2010", ssp = "current",
                                         suitable_km2 = round(a0, 2), retained_km2 = NA, retained_pct = NA,
                                         gain_km2 = NA, net_change_pct = NA))
  pa <- list.files(rd, pattern = "^gcm_.*_presence_agreement\\.tif$", full.names = TRUE)
  for (f in pa) {
    key <- sub("^gcm_(.*)_presence_agreement\\.tif$", "\\1", basename(f))
    p <- strsplit(key, "_")[[1]]
    ag <- mask(rast(f), smask)
    fut <- ag >= AGREE
    ret <- global((b0 == 1) & fut, "sum", na.rm = TRUE)[1, 1] * cell_km2
    gain <- global((b0 == 0) & fut, "sum", na.rm = TRUE)[1, 1] * cell_km2
    area_tab <- rbind(area_tab, data.frame(model_set = M, period = p[1], ssp = p[2], suitable_km2 = round(ret + gain, 2),
      retained_km2 = round(ret, 2), retained_pct = round(100 * ret / a0, 1), gain_km2 = round(gain, 2),
      net_change_pct = round(100 * (ret + gain - a0) / a0, 1)))
    if (M == "M3") keep_rasters[[key]] <- ((b0 == 1) & fut)
  }
}
if (!is.null(area_tab)) write.csv(area_tab, file.path(SUM_DIR, "suitable_area_change.csv"), row.names = FALSE)
if (!is.null(area_tab) && all(c("M1", "M3") %in% area_tab$model_set)) {
  w <- merge(area_tab[area_tab$model_set == "M1" & area_tab$ssp != "current", c("period", "ssp", "retained_pct")],
             area_tab[area_tab$model_set == "M3" & area_tab$ssp != "current", c("period", "ssp", "retained_pct")],
             by = c("period", "ssp"), suffixes = c("_M1", "_M3"))
  w$buffering_pp <- w$retained_pct_M3 - w$retained_pct_M1
  write.csv(w, file.path(SUM_DIR, "topographic_buffering.csv"), row.names = FALSE)
}

# ---- 5. Mikrosığınakların topoğrafik özellikleri (M3) -------------------------
if (length(keep_rasters)) {
  fb <- cur_bin("M3"); b0 <- mask(rast(fb), smask)
  key <- grep("2071-2100_ssp370|2071-2100_ssp585", names(keep_rasters), value = TRUE)[1]
  if (is.na(key)) key <- names(keep_rasters)[1]
  ref_r <- keep_rasters[[key]]
  writeRaster(ref_r, file.path(SUM_DIR, paste0("microrefugia_M3_", key, ".tif")), overwrite = TRUE)
  cls <- ifel(b0 == 1 & ref_r, 2, ifel(b0 == 1, 1, NA))       # 2 sığınak, 1 kaybedilen
  tv <- c("northness", "eastness", "insol_spr", "hli", "twi", "tpi300", "slope", "canopy_cover", "elev")
  tv <- tv[file.exists(file.path(TOPO_DIR, paste0(tv, ".tif")))]
  st <- c(cls, rast(file.path(TOPO_DIR, paste0(tv, ".tif")))); names(st) <- c("cls", tv)
  sm <- spatSample(st, 60000, method = "random", na.rm = TRUE)
  if (nrow(sm) > 50) {
    tab <- do.call(rbind, lapply(tv, function(v) {
      x2 <- sm[[v]][sm$cls == 2]; x1 <- sm[[v]][sm$cls == 1]
      w <- tryCatch(wilcox.test(x2, x1), error = function(e) NULL)
      data.frame(variable = v, refugia_median = round(median(x2), 3), lost_median = round(median(x1), 3),
                 n_refugia = length(x2), n_lost = length(x1), p_wilcox = if (is.null(w)) NA else signif(w$p.value, 3))
    }))
    tab$scenario <- key
    write.csv(tab, file.path(SUM_DIR, "microrefugia_topography.csv"), row.names = FALSE)
  }
}

# ---- 6. Belirsizlik ayrıştırma (M3, 2071–2100) ------------------------------
fp <- file.path(RESULT_DIR, "M3", "proj", "probe_algorithm_predictions.csv")
if (file.exists(fp)) {
  pr <- read.csv(fp); pr <- pr[grepl("^future_2071-2100_", pr$scenario), ]
  algs <- setdiff(names(pr), c("scenario", "cell", "ensemble"))
  long <- do.call(rbind, lapply(algs, function(a) data.frame(cell = pr$cell, scenario = pr$scenario, algo = a, y = pr[[a]])))
  long$ssp <- sub("^future_[^_]+_(ssp[0-9]+)_.*$", "\\1", long$scenario)
  long$gcm <- sub("^future_[^_]+_ssp[0-9]+_", "", long$scenario)
  long <- long[is.finite(long$y), ]
  ss <- function(d) {
    s <- summary(aov(y ~ algo + gcm + ssp + algo:gcm + algo:ssp + gcm:ssp, data = d))[[1]]
    v <- s[, "Sum Sq"]; names(v) <- trimws(rownames(s)); v / sum(v)
  }
  cl <- sample(unique(long$cell), min(2000, length(unique(long$cell))))
  sp_long <- split(long[long$cell %in% cl, ], long$cell[long$cell %in% cl])
  props <- do.call(rbind, lapply(sp_long, function(d) if (sd(d$y) > 0) ss(d) else NULL))
  vp <- data.frame(component = colnames(props), mean_share = round(colMeans(props, na.rm = TRUE), 4),
                   sd_share = round(apply(props, 2, sd, na.rm = TRUE), 4))
  write.csv(vp, file.path(SUM_DIR, "uncertainty_partitioning_2071_2100.csv"), row.names = FALSE)
}
log_msg("ADIM 10 tamamlandı: ", SUM_DIR)
