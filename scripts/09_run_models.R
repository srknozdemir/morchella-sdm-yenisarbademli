# =============================================================================
#  ADIM 9 — Modellerin çalıştırılması (biomod2-sdm-runner v1.2, başsız mod)
#  Her model seti için Rscript run_biomod.R <config.json> çağrılır. Runner kendi
#  klasöründe çalışmalıdır (yardımcı dosyaları oradan okur).
#  Kullanım: Rscript 09_run_models.R [M1 M2 M3]
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
RUNNER_DIR <- file.path(SCRIPT_DIR, "runner")
if (!file.exists(file.path(RUNNER_DIR, "run_biomod.R"))) stop("runner/run_biomod.R yok: ", RUNNER_DIR)

which_sets <- commandArgs(trailingOnly = TRUE)
if (!length(which_sets)) which_sets <- names(MODEL_SETS)
rscript <- file.path(R.home("bin"), "Rscript")
for (M in which_sets) {
  cfgf <- file.path(RESULT_DIR, M, "config.json")
  if (!file.exists(cfgf)) { log_msg(M, ": config yok, atlandı."); next }
  log_msg(M, ": model çalıştırılıyor …")
  t0 <- Sys.time()
  # Her set kendi çalışma klasöründe: biomod2 model dosyalarını tür adıyla çalışma
  # klasörüne yazdığı için setler paralel çalıştırıldığında çakışmazlar.
  wd <- file.path(RESULT_DIR, M, "biomod_wd"); dir.create(wd, FALSE, TRUE)
  file.copy(list.files(RUNNER_DIR, pattern = "\\.R$", full.names = TRUE), wd, overwrite = TRUE)
  old <- setwd(wd)
  st <- system2(rscript, c("run_biomod.R", shQuote(normalizePath(cfgf))),
                stdout = file.path(RESULT_DIR, M, "run.log"), stderr = file.path(RESULT_DIR, M, "run.log"))
  setwd(old)
  log_msg(M, ": bitti (çıkış kodu ", st, ", ", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " dk)")
}
