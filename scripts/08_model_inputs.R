# =============================================================================
#  ADIM 8 — Model girdi klasörleri ve biomod2-sdm-runner yapılandırmaları
#  Her model seti (M1, M2, M3) için:
#    data/model_inputs/<M>/current/<var>.tif
#    data/model_inputs/<M>/future_<dönem>_<ssp>_<GCM>/<var>.tif
#      (iklim değişkenleri GCM'e göre değişir; topoğrafya, kanopi, toprak sabit)
#    results/<M>/config.json    (runner girdisi)
#  Model bir kez güncel iklimle kalibre edilir, her GCM × SSP × dönem için ayrı
#  projekte edilir; GCM'ler projeksiyon düzeyinde birleştirilir (Araújo & New 2007).
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))

sets <- readRDS(file.path(SEL_DIR, "model_sets.rds"))
src_of <- function(v, set_dir = file.path(CLIM30_DIR, "current")) {
  g <- group_of(v)
  if (g == "iklim") return(file.path(set_dir, paste0(v, ".tif")))
  d <- switch(g, topografya = TOPO_DIR, kanopi = TOPO_DIR, toprak = SOIL_DIR, ek = EXTRA_DIR)
  file.path(d, paste0(v, ".tif"))
}
fut_dirs <- list.dirs(CLIM30_DIR, recursive = FALSE)
fut_dirs <- fut_dirs[grepl("^future_", basename(fut_dirs))]
parse_fut <- function(b) {
  m <- regmatches(b, regexec("^future_(.+?)_(ssp[0-9]{3})_(.+)$", b))[[1]]
  list(period = m[2], ssp = m[3], gcm = m[4])
}

for (M in names(sets)) {
  vars <- sets[[M]]
  if (!length(vars)) { log_msg(M, ": değişken yok, atlandı."); next }
  root <- file.path(MODEL_DIR, M)
  unlink(root, recursive = TRUE); dir.create(file.path(root, "current"), recursive = TRUE)
  for (v in vars) file.copy(src_of(v), file.path(root, "current", paste0(v, ".tif")), overwrite = TRUE)

  has_clim <- any(vapply(vars, group_of, "") == "iklim") && M %in% FUTURE_SETS
  fut_blocks <- list()
  if (has_clim && length(fut_dirs)) {
    info <- lapply(basename(fut_dirs), parse_fut)
    periods <- unique(vapply(info, `[[`, "", "period"))
    for (p in periods) {
      scs <- list()
      for (k in which(vapply(info, `[[`, "", "period") == p)) {
        fd <- file.path(root, basename(fut_dirs[k])); dir.create(fd, showWarnings = FALSE)
        for (v in vars) file.copy(src_of(v, fut_dirs[k]), file.path(fd, paste0(v, ".tif")), overwrite = TRUE)
        scs[[length(scs) + 1]] <- list(label = paste0(info[[k]]$ssp, "_", info[[k]]$gcm),
                                       asc_dir = normalizePath(fd), ssp = info[[k]]$ssp, gcm = info[[k]]$gcm)
      }
      fut_blocks[[length(fut_blocks) + 1]] <- list(period_label = p, scenarios = scs)
    }
  }

  out <- file.path(RESULT_DIR, M); dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cfg <- list(
    run_id = paste0("morchella_", M), asc_dir = normalizePath(file.path(root, "current")),
    occ_file = normalizePath(file.path(OCC_DIR, "occ_local_thinned.csv")),
    bg_file  = normalizePath(file.path(OCC_DIR, "bg_local_bias.csv")),
    sp_name = SPECIES_NAME, out_dir = normalizePath(out),
    resp_kind = "binary", data_type = "binary", filter_raster = TRUE,
    pa_strategy = "random", pa_rep = PA_REP, pa_n = PA_N, pa_dist_min = 0L, pa_dist_max = 0L,
    cv_strategy = "block", cv_k = 4L, cv = list(strategy = "block", k = 4, rep = 1, perc = 70),
    cpu_n = 1L, models = as.list(SDM_MODELS), opt_strategy = "adaptive",
    ens_min_auc = ENS_MIN_AUC, ens_min_tss = ENS_MIN_TSS, ens_near_one = 0.999,
    out_scale = "0_1", rc_em_pick = "EMwmean", change_maps = TRUE,
    change_em = "EMwmean", change_agree = 0.8, change_tol_frac = 0.05, binary_consensus = 0.5,
    project_full_models = PROJECT_FULL_MODELS, terra_memfrac = 0.4, clamping_mask = FALSE, skip_projection = TRUE,
    future_n = length(fut_blocks), future = fut_blocks
  )
  if (!is.null(MODEL_OVERRIDES[[M]])) for (k in names(MODEL_OVERRIDES[[M]]))
    cfg[[k]] <- if (k == "models") as.list(MODEL_OVERRIDES[[M]][[k]]) else MODEL_OVERRIDES[[M]][[k]]
  jsonlite::write_json(cfg, file.path(out, "config.json"), auto_unbox = TRUE, pretty = TRUE, null = "null")
  log_msg(M, ": ", length(vars), " değişken, ", sum(lengths(lapply(fut_blocks, `[[`, "scenarios"))),
          " gelecek seti → ", file.path(out, "config.json"))
}
log_msg("ADIM 8 tamamlandı.")
