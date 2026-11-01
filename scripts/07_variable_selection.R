# =============================================================================
#  ADIM 7 — Değişken eleme
#  1. Aday katmanlar: climate_30m/current (tüm gelecek setlerinde bulunanlar), topo,
#     soil, extra; EXCLUDE_VARS dışarıda
#  2. Değerler: 221 kayıt + yanlılık ağırlıklı arka plan noktaları
#  3. Neredeyse sabit değişkenler çıkarılır
#  4. Tek değişkenli önem: GLM (x + x²), 4 mekânsal blokta (Muscarella et al. 2014
#     düzeninde) test AUC'si
#  5. Spearman |r| ≥ COR_MAX kümeleri (Dormann et al. 2013); her kümeden en yüksek
#     blok AUC'li değişken (fark < 0.01 ise ECO_PRIORITY önceliklidir)
#  6. VIF < VIF_MAX (adımsal)
#  7. İç içe model setleri M1 ⊂ M2 ⊂ M3 (MODEL_SETS) ve gelecekte korelasyon kayması
#     (Braunisch et al. 2013)
#  Çıktılar: data/selection/*.csv, *.png
# =============================================================================
.script_dir <- function() {
  for (f in rev(sys.frames())) if (!is.null(f$ofile)) return(dirname(normalizePath(f$ofile)))
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) return(dirname(normalizePath(sub("^--file=", "", a[1]))))
  getwd()
}
SCRIPT_DIR <- .script_dir()
source(file.path(SCRIPT_DIR, "00_config.R"))
set.seed(2027)

# ---- 1. Aday katmanlar ------------------------------------------------------
cur <- file.path(CLIM30_DIR, "current")
clim_files <- list.files(cur, pattern = "\\.tif$", full.names = TRUE)
fut_sets <- list.dirs(CLIM30_DIR, recursive = FALSE)
fut_sets <- fut_sets[grepl("^future_", basename(fut_sets))]
if (length(fut_sets)) {
  in_all <- Reduce(intersect, lapply(fut_sets, function(d) sub("\\.tif$", "", list.files(d, "\\.tif$"))))
  clim_files <- clim_files[sub("\\.tif$", "", basename(clim_files)) %in% in_all]
}
files <- c(clim_files,
           list.files(TOPO_DIR,  pattern = "\\.tif$", full.names = TRUE),
           list.files(SOIL_DIR,  pattern = "\\.tif$", full.names = TRUE),
           list.files(EXTRA_DIR, pattern = "\\.tif$", full.names = TRUE))
nm <- sub("\\.tif$", "", basename(files))
dc <- file.path(LOG_DIR, "downscale_coefficients.csv")
weak <- character()
if (file.exists(dc)) {
  co <- read.csv(dc); weak <- co$var[co$rule == "interp" | (!is.na(co$r2) & co$r2 < DOWNSCALE_MIN_R2)]
  log_msg("İndirgeme R² < ", DOWNSCALE_MIN_R2, " ya da yalnız enterpolasyon, aday dışı: ", paste(weak, collapse = ", "))
}
keep <- !nm %in% c(EXCLUDE_VARS, weak) & !duplicated(nm)
files <- files[keep]; nm <- nm[keep]
env <- rast(files); names(env) <- nm
grp <- vapply(nm, group_of, "")
log_msg("Aday değişken: ", length(nm), " (", paste(names(table(grp)), table(grp), sep = "=", collapse = ", "), ")")

# ---- 2. Değerler -------------------------------------------------------------
pres <- read.csv(file.path(OCC_DIR, "occ_local_thinned.csv"))
bg   <- read.csv(file.path(OCC_DIR, "bg_local_bias.csv"))
P <- extract(env, as.matrix(pres[, c("x", "y")]))
B <- extract(env, as.matrix(bg[, c("x", "y")]))
okP <- complete.cases(P); okB <- complete.cases(B)
P <- P[okP, , drop = FALSE]; B <- B[okB, , drop = FALSE]
xyP <- as.matrix(pres[okP, c("x", "y")]); xyB <- as.matrix(bg[okB, c("x", "y")])
log_msg("Kayıt: ", nrow(P), "; arka plan: ", nrow(B))

# ---- 3. Neredeyse sabit değişkenler -----------------------------------------
allv <- rbind(P, B)
nzv  <- vapply(allv, function(x) sd(x) > 0 && length(unique(round(x, 6))) >= 10 &&
                 (sd(x) / max(abs(mean(x)), 1e-9)) > 1e-3, logical(1))
dropped_nzv <- names(nzv)[!nzv]
if (length(dropped_nzv)) log_msg("Sabite yakın, çıkarıldı: ", paste(dropped_nzv, collapse = ", "))
P <- P[, nzv, drop = FALSE]; B <- B[, nzv, drop = FALSE]

# ---- 4. Tek değişkenli blok AUC ---------------------------------------------
block4 <- function(xyP, xyB) {
  lon <- median(xyP[, 1])
  gP <- integer(nrow(xyP)); gB <- integer(nrow(xyB))
  for (side in 0:1) {
    iP <- if (side == 0) xyP[, 1] <= lon else xyP[, 1] > lon
    iB <- if (side == 0) xyB[, 1] <= lon else xyB[, 1] > lon
    lat <- median(xyP[iP, 2])
    gP[iP] <- side * 2 + ifelse(xyP[iP, 2] <= lat, 1, 2)
    gB[iB] <- side * 2 + ifelse(xyB[iB, 2] <= lat, 1, 2)
  }
  list(P = gP, B = gB)
}
auc <- function(sp, sb) {
  r <- rank(c(sp, sb)); np <- length(sp); nb <- length(sb)
  (sum(r[seq_len(np)]) - np * (np + 1) / 2) / (np * nb)
}
blk <- block4(xyP, xyB)
uni <- data.frame(var = names(P), group = vapply(names(P), group_of, ""), auc_block = NA_real_,
                  auc_train = NA_real_, stringsAsFactors = FALSE)
for (v in names(P)) {
  d <- data.frame(y = c(rep(1, nrow(P)), rep(0, nrow(B))), x = scale(c(P[[v]], B[[v]]))[, 1],
                  g = c(blk$P, blk$B))
  w <- ifelse(d$y == 1, 1, nrow(P) / nrow(B))
  a <- c()
  for (k in 1:4) {
    tr <- d$g != k
    if (sum(d$y[!tr]) < 5) next
    f <- suppressWarnings(glm(y ~ x + I(x^2), family = binomial, data = d[tr, ], weights = w[tr]))
    p <- predict(f, d[!tr, ], type = "response")
    a <- c(a, auc(p[d$y[!tr] == 1], p[d$y[!tr] == 0]))
  }
  f <- suppressWarnings(glm(y ~ x + I(x^2), family = binomial, data = d, weights = w))
  p <- fitted(f)
  uni$auc_block[uni$var == v] <- mean(a)
  uni$auc_train[uni$var == v] <- auc(p[d$y == 1], p[d$y == 0])
}
uni <- uni[order(-uni$auc_block), ]

# ---- 5. Korelasyon kümeleri -------------------------------------------------
R  <- cor(B, method = "spearman")
write.csv(round(R, 3), file.path(SEL_DIR, "spearman_all_candidates.csv"))
hc <- hclust(as.dist(1 - abs(R)), method = "complete")
cl <- cutree(hc, h = 1 - COR_MAX)
uni$cluster <- cl[uni$var]
pri <- function(v) { i <- match(v, ECO_PRIORITY); ifelse(is.na(i), 999, i) }
chosen <- character()
for (k in unique(uni$cluster)) {
  s <- uni[uni$cluster == k, ]
  best <- s$auc_block[1]
  cand <- s[s$auc_block >= best - 0.01, ]
  chosen <- c(chosen, cand$var[order(pri(cand$var), -cand$auc_block)][1])
}
png(file.path(SEL_DIR, "dendrogram.png"), 1800, 900, res = 150)
plot(hc, main = sprintf("Spearman |r| kümeleri (kesim |r| = %.2f)", COR_MAX), xlab = "", sub = "", ylab = "1 − |r|")
abline(h = 1 - COR_MAX, col = "red", lty = 2)
dev.off()

# ---- 6. VIF ------------------------------------------------------------------
vif_of <- function(X) vapply(names(X), function(v) {
  r2 <- summary(lm(X[[v]] ~ ., data = X[, setdiff(names(X), v), drop = FALSE]))$r.squared
  1 / (1 - r2) }, numeric(1))
# Mekânsal blok testinde rastgeleden kötü olan (AUC < 0.5) değişkenler aday dışı bırakılır
chosen <- chosen[uni$auc_block[match(chosen, uni$var)] >= 0.5]
sel <- chosen[order(-uni$auc_block[match(chosen, uni$var)])]
if (length(sel) > MAX_VARS) sel <- sel[seq_len(MAX_VARS)]
repeat {
  if (length(sel) < 3) break
  vf <- vif_of(B[, sel, drop = FALSE])
  if (max(vf) < VIF_MAX) break
  hi_vif <- names(vf)[vf >= VIF_MAX]                 # VIF eşiğini aşanlar arasından
  drop <- hi_vif[which.min(uni$auc_block[match(hi_vif, uni$var)])]   # en zayıf blok AUC'li çıkar
  sel <- setdiff(sel, drop)
}
vf <- if (length(sel) >= 2) vif_of(B[, sel, drop = FALSE]) else setNames(1, sel)
uni$selected <- uni$var %in% sel
uni$vif <- round(vf[uni$var], 2)

# ---- 7. Model setleri ve korelasyon kayması ----------------------------------
sets <- lapply(MODEL_SETS, function(g) sel[vapply(sel, group_of, "") %in% g])
for (m in names(sets)) uni[[m]] <- uni$var %in% sets[[m]]
write.csv(uni, file.path(SEL_DIR, "variable_selection.csv"), row.names = FALSE)
writeLines(unlist(lapply(names(sets), function(m) paste0(m, ": ", paste(sets[[m]], collapse = ", ")))),
           file.path(SEL_DIR, "model_sets.txt"))
saveRDS(sets, file.path(SEL_DIR, "model_sets.rds"))

clim_sel <- sets$M1
if (length(clim_sel) >= 2 && length(fut_sets)) {
  idx <- sample(nrow(xyB), min(5000, nrow(xyB)))
  r0 <- cor(B[idx, clim_sel], method = "spearman")
  shift <- do.call(rbind, lapply(fut_sets, function(d) {
    Fv <- extract(rast(file.path(d, paste0(clim_sel, ".tif"))), xyB[idx, , drop = FALSE])
    names(Fv) <- clim_sel
    r1 <- cor(Fv, method = "spearman", use = "complete.obs")
    data.frame(set = basename(d), max_abs_delta_r = max(abs(r1 - r0)[upper.tri(r0)], na.rm = TRUE))
  }))
  write.csv(shift, file.path(SEL_DIR, "correlation_shift_future.csv"), row.names = FALSE)
  log_msg("Gelecekte korelasyon kayması: en büyük |Δr| = ", round(max(shift$max_abs_delta_r), 3))
}

print(uni[, c("var", "group", "auc_block", "cluster", "selected", "vif")], row.names = FALSE, digits = 3)
for (m in names(sets)) log_msg(m, ": ", paste(sets[[m]], collapse = ", "))
log_msg("ADIM 7 tamamlandı.")
