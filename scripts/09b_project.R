# =============================================================================
#  ADIM 9b — Projeksiyon ve GCM ansamblı (runner dışında, bellek ve hız için)
#  Runner (skip_projection = TRUE) modelleri kalibre edip değerlendirdikten sonra:
#   1. Ansambla giren algoritmaların tam veri ("allRun") modelleri yüklenir
#   2. Ağırlık: algoritmanın mekânsal blok CV TSS ortalaması (bağımsız doğrulama)
#   3. GAM eklemeli olduğundan her düzgünleştirici terim, tüm senaryoların değer
#      aralığını kapsayan 2000 noktalık ızgarada bir kez hesaplanır ve doğrusal
#      enterpolasyonla uygulanır (mgcv TPRS tahmininin 30 m ızgarada saatler süren
#      yükünü kaldırır; ızgara hatası < 1e-4)
#   4. Güncel ve her GCM × SSP × dönem için ansambl uygunluğu; maxTSS eşiği güncel
#      tahminden (kayıtlar ve yanlılık ağırlıklı arka plan) bir kez hesaplanır
#   5. Dönem × SSP için GCM ortalaması, SD ve varlık uyumu (GCM'lerin payı)
#   6. Belirsizlik ayrıştırma için 20 000 sabit hücrede algoritma düzeyinde tahminler
#  Çıktılar: results/<M>/proj/*.tif, results/<M>/proj/*.csv
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
suppressPackageStartupMessages({ library(biomod2); library(mgcv) })
set.seed(99)
sets <- commandArgs(trailingOnly = TRUE); if (!length(sets)) sets <- FUTURE_SETS
sp <- gsub(" ", ".", SPECIES_NAME)
tmpl <- local_template(); smask <- study_mask(tmpl)
z <- rast(file.path(TOPO_DIR, "elev.tif")); sl <- rast(file.path(TOPO_DIR, "slope.tif"))
land <- mask(ifel(z <= LAKE_MAX_ELEV & sl < 0.5, NA, 1), smask)
cells <- which(!is.na(values(land, mat = FALSE)))
probe <- sort(sample(cells, 20000))

for (M in sets) {
  rd <- file.path(RESULT_DIR, M); wd <- file.path(rd, "biomod_wd"); od <- file.path(rd, "proj")
  dir.create(od, FALSE, TRUE)
  mdir <- file.path(wd, sp, "models", "binary")
  cfg <- jsonlite::read_json(file.path(rd, "config.json"))
  algos <- tryCatch(read.csv(file.path(rd, "ensemble_candidates_after_thresholds.csv"))$model, error = function(e) character())
  if (!length(algos)) algos <- unlist(cfg$models)
  ev <- read.csv(file.path(rd, "eval_model.csv"))
  ev <- ev[ev$value.run != "allRun" & ev$value.metric.eval == "TSS" & !is.na(ev$value.validation), ]
  tss <- tapply(ev$value.validation, ev$value.algo, mean)[algos]
  w_alg <- pmax(tss, 0.01); w_alg <- w_alg / sum(w_alg)
  files <- unlist(lapply(algos, function(a) list.files(mdir, pattern = paste0("^", sp, "_PA[0-9]+_allRun_", a, "$"), full.names = TRUE)))
  mods <- lapply(files, function(f) get(load(f)))
  m_alg <- sub(".*_", "", basename(files))
  log_msg(M, ": ", length(files), " tam veri modeli; ağırlıklar ", paste(names(w_alg), round(w_alg, 3), sep = "=", collapse = ", "))
  write.csv(data.frame(algo = names(w_alg), block_TSS = round(tss, 3), weight = round(w_alg, 3)),
            file.path(od, "ensemble_weights.csv"), row.names = FALSE)

  in_dir <- file.path(MODEL_DIR, M)
  vars <- sub("\\.tif$", "", list.files(file.path(in_dir, "current"), "\\.tif$"))
  scen <- c("current", grep("^future_", list.dirs(in_dir, recursive = FALSE, full.names = FALSE), value = TRUE))
  read_env <- function(s) { r <- rast(file.path(in_dir, s, paste0(vars, ".tif"))); names(r) <- vars; r }

  # GAM terim ızgaraları (tüm senaryoların aralığı)
  # Kalibrasyon aralığı (kayıtlar + arka plan adayları); projeksiyonda değerler bu aralığa
  # sıkıştırılır (clamping), aralık dışı koşullar ayrıca MESS ile raporlanır (ADIM 10)
  prx <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv")); bgx <- read.csv(file.path(OCC_DIR, "bg_local_bias.csv"))
  calv <- extract(read_env("current"), as.matrix(rbind(prx[, c("x", "y")], bgx[, c("x", "y")])))
  lo <- sapply(vars, function(v) min(calv[[v]], na.rm = TRUE)); hi <- sapply(vars, function(v) max(calv[[v]], na.rm = TRUE))
  if (!CLAMP_TO_CALIBRATION) {
    rr <- do.call(rbind, lapply(scen, function(s) { g <- global(read_env(s), c("min", "max"), na.rm = TRUE); data.frame(var = vars, min = g$min, max = g$max) }))
    lo <- tapply(rr$min, rr$var, min)[vars]; hi <- tapply(rr$max, rr$var, max)[vars]
  }
  write.csv(data.frame(var = vars, cal_min = lo, cal_max = hi), file.path(od, "calibration_range.csv"), row.names = FALSE)
  gam_fun <- lapply(seq_along(mods), function(k) {
    if (m_alg[k] != "GAM") return(NULL)
    g <- get_formal_model(mods[[k]])
    med <- as.data.frame(lapply(vars, function(v) rep((lo[[v]] + hi[[v]]) / 2, 2000))); names(med) <- vars
    terms <- list()
    for (v in vars) {
      gd <- med; gd[[v]] <- seq(lo[[v]], hi[[v]], length.out = 2000)
      tt <- predict(g, newdata = gd, type = "terms")
      col <- grep(paste0("\\(", v, "\\)$"), colnames(tt))
      terms[[v]] <- list(x = gd[[v]], y = if (length(col)) tt[, col] else rep(0, 2000))
    }
    list(const = attr(predict(g, newdata = med[1:2, , drop = FALSE], type = "terms"), "constant"), terms = terms, link = g$family$linkinv)
  })
  pred_one <- function(k, df) {
    if (m_alg[k] == "GAM") {
      gf <- gam_fun[[k]]; eta <- rep(gf$const, nrow(df))
      for (v in vars) eta <- eta + approx(gf$terms[[v]]$x, gf$terms[[v]]$y, df[[v]], rule = 2)$y
      return(gf$link(eta))
    }
    as.numeric(predict(mods[[k]], newdata = df, on_0_1000 = FALSE))
  }

  probe_long <- list(); thr <- NA
  for (s in scen) {
    out <- file.path(od, paste0(s, "_suitability.tif"))
    env <- read_env(s)
    X <- as.data.frame(env[cells])                         # yalnız kara hücreleri
    ok <- complete.cases(X)
    if (CLAMP_TO_CALIBRATION) for (v in vars) X[[v]] <- pmin(pmax(X[[v]], lo[[v]]), hi[[v]])
    idx <- which(ok)
    run_k <- function(k) {
      out <- rep(NA_real_, nrow(X))
      for (ch in split(idx, ceiling(seq_along(idx) / 5e5))) out[ch] <- pred_one(k, X[ch, , drop = FALSE])
      out
    }
    cols <- parallel::mclapply(seq_along(mods), run_k, mc.cores = N_CORES)
    bad <- which(!vapply(cols, function(x) is.numeric(x) && length(x) == nrow(X), logical(1)))
    for (k in bad) cols[[k]] <- run_k(k)                 # başarısız çocuk süreçler sırayla yeniden
    P <- do.call(cbind, cols); rm(cols); invisible(gc())
    alg_mean <- sapply(names(w_alg), function(a) rowMeans(P[, m_alg == a, drop = FALSE]))
    ens <- as.numeric(alg_mean %*% w_alg)
    r <- rast(tmpl); r[cells] <- ens; names(r) <- "suitability"
    writeRaster(r, out, overwrite = TRUE, wopt = WOPT)
    pi <- match(probe, cells)
    probe_long[[s]] <- data.frame(scenario = s, cell = probe, alg_mean[pi, , drop = FALSE], ensemble = ens[pi])
    if (s == "current") {
      pr <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv")); bg <- read.csv(file.path(OCC_DIR, "bg_local_bias.csv"))
      sp_ <- extract(r, as.matrix(pr[, c("x", "y")]))[, 1]; sb <- extract(r, as.matrix(bg[, c("x", "y")]))[, 1]
      sp_ <- sp_[is.finite(sp_)]; sb <- sb[is.finite(sb)]
      th <- seq(0.001, 0.999, by = 0.001)
      tssv <- sapply(th, function(t) mean(sp_ >= t) + mean(sb < t) - 1)
      thr <- th[which.max(tssv)]
      auc_tr <- (sum(rank(c(sp_, sb))[seq_along(sp_)]) - length(sp_) * (length(sp_) + 1) / 2) / (length(sp_) * length(sb))
      write.csv(data.frame(threshold_maxTSS = thr, TSS_at_threshold = round(max(tssv), 3),
                           sensitivity = round(mean(sp_ >= thr), 3), specificity_bg = round(mean(sb < thr), 3),
                           AUC_presence_vs_background = round(auc_tr, 3)),
                file.path(od, "threshold_current.csv"), row.names = FALSE)
      writeRaster(r >= thr, file.path(od, "current_binary.tif"), overwrite = TRUE, datatype = "INT1U", wopt = WOPT)
    }
    log_msg(M, " ", s, ": tamam")
    rm(X, P, alg_mean, ens); invisible(gc())
  }
  write.csv(do.call(rbind, probe_long), file.path(od, "probe_algorithm_predictions.csv"), row.names = FALSE)

  # GCM ansamblı
  fs <- grep("^future_", scen, value = TRUE)
  info <- do.call(rbind, lapply(fs, function(s) { m <- regmatches(s, regexec("^future_(.+?)_(ssp[0-9]{3})_(.+)$", s))[[1]]
                                                 data.frame(set = s, period = m[2], ssp = m[3], gcm = m[4]) }))
  for (key in unique(paste(info$period, info$ssp, sep = "_"))) {
    ss <- info$set[paste(info$period, info$ssp, sep = "_") == key]
    st <- rast(file.path(od, paste0(ss, "_suitability.tif")))
    writeRaster(mean(st), file.path(od, paste0("gcm_", key, "_mean.tif")), overwrite = TRUE, wopt = WOPT)
    writeRaster(stdev(st), file.path(od, paste0("gcm_", key, "_sd.tif")), overwrite = TRUE, wopt = WOPT)
    writeRaster(mean(st >= thr), file.path(od, paste0("gcm_", key, "_presence_agreement.tif")), overwrite = TRUE, wopt = WOPT)
  }
  write.csv(info, file.path(od, "scenario_registry.csv"), row.names = FALSE)
  log_msg(M, ": projeksiyon ve GCM ansamblı tamamlandı (eşik ", thr, ")")
}
