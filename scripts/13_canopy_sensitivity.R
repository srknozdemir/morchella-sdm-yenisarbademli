# =============================================================================
#  ADIM 13 — Kanopi duyarlılık analizi (M3)
#  Gelecek projeksiyonlarında kanopi örtüsü bugünkü değerinde sabit tutulur. Bu adım,
#  kanopi örtüsü %20 ve %40 azaldığında ve tamamen kalktığında (alt sınır) uygun alanın
#  ve mikrosığınakların nasıl değiştiğini sınar. Hız için sistematik 1/9 hücre örneği
#  (her üç satır ve sütundan biri) kullanılır; alanlar örnek oranıyla ölçeklenir.
#  Çıktı: results/summary/canopy_sensitivity.csv
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
sets <- "M3"
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


  thr <- read.csv(file.path(od, "threshold_current.csv"))$threshold_maxTSS
  rc <- rowColFromCell(tmpl, cells); sub <- cells[rc[, 1] %% 3 == 1 & rc[, 2] %% 3 == 1]
  cell_km2 <- prod(res(tmpl)) / 1e6 * length(cells) / length(sub)
  b0 <- extract(rast(file.path(od, "current_binary.tif")), sub)[, 1] == 1
  ens_pred <- function(X) {
    ok <- complete.cases(X)
    for (v in vars) X[[v]] <- pmin(pmax(X[[v]], lo[[v]]), hi[[v]])
    cols <- parallel::mclapply(seq_along(mods), function(k) { o <- rep(NA_real_, nrow(X)); o[ok] <- pred_one(k, X[ok, , drop = FALSE]); o }, mc.cores = N_CORES)
    bad <- which(!vapply(cols, function(x) is.numeric(x) && length(x) == nrow(X), logical(1)))
    for (k in bad) { o <- rep(NA_real_, nrow(X)); o[ok] <- pred_one(k, X[ok, , drop = FALSE]); cols[[k]] <- o }
    P <- do.call(cbind, cols)
    A <- sapply(names(w_alg), function(a) rowMeans(P[, m_alg == a, drop = FALSE]))
    as.numeric(A %*% w_alg)
  }
  levels_f <- c(1, 0.8, 0.6, 0)
  fsets <- grep("^future_2071-2100_", scen, value = TRUE)
  ssp_of <- sub("^future_2071-2100_(ssp[0-9]{3})_.*$", "\\1", fsets)
  res_tab <- NULL
  # bugünkü iklim, azaltılmış kanopi
  Xc <- as.data.frame(read_env("current")[sub])
  for (f in levels_f) {
    X <- Xc; X$canopy_cover <- X$canopy_cover * f
    s <- ens_pred(X) >= thr
    res_tab <- rbind(res_tab, data.frame(period = "1981-2010", ssp = "current", canopy_factor = f,
      suitable_km2 = round(sum(s, na.rm = TRUE) * cell_km2, 1), retained_pct = round(100 * mean(s[b0], na.rm = TRUE), 1),
      refugia_retained_pct = NA))
  }
  for (sp in unique(ssp_of)) {
    ss <- fsets[ssp_of == sp]
    Xs <- lapply(ss, function(s) as.data.frame(read_env(s)[sub]))
    agree <- list()
    for (f in levels_f) {
      B <- sapply(Xs, function(X0) { X <- X0; X$canopy_cover <- X$canopy_cover * f; ens_pred(X) >= thr })
      agree[[as.character(f)]] <- rowMeans(B, na.rm = TRUE) >= 0.8
    }
    ref1 <- b0 & agree[["1"]]
    for (f in levels_f) {
      a <- agree[[as.character(f)]]
      res_tab <- rbind(res_tab, data.frame(period = "2071-2100", ssp = sp, canopy_factor = f,
        suitable_km2 = round(sum(a, na.rm = TRUE) * cell_km2, 1), retained_pct = round(100 * mean(a[b0], na.rm = TRUE), 1),
        refugia_retained_pct = round(100 * mean(a[ref1], na.rm = TRUE), 1)))
    }
    log_msg("Kanopi duyarlılığı ", sp, ": tamam")
  }
  write.csv(res_tab, file.path(RESULT_DIR, "summary", "canopy_sensitivity.csv"), row.names = FALSE)
  print(res_tab)
}
