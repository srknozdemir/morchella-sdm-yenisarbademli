# Morchella esculenta s.l. — Yenişarbademli: yerel ölçek (30 m) SDM iş akışı

Yalnız saha kayıtları (kuzu.csv) kullanılır. Model alanı: çalışma alanı sınırı
(serkan_sam_sym/yenisar.shp, 4560 km²), 30 m, UTM 36N. İklim kesitleri sınır + 20 km.
Tüm ayarlar 00_config.R içindedir.

## Çalıştırma sırası

    setwd(".../kuzugobegi/sdm_pipeline")
    source("01_download_climate.R")   # CHELSA v2.1, 1981–2010 + 2 dönem × 3 SSP × 5 GCM (1 km)
    source("04_soil.R")               # SoilGrids 2.0, 0–30 cm
    source("02a_prepare_terrain.R")   # HGM SAM + SYM paftaları: birleştirme, hizalama
    system("python3 -I 02b_twi.py ..")# TWI (pysheds, D-infinity); whitebox varsa gerekmez
    source("02_dem_topography.R")     # arazi, ışınım, kanopi (30 m)
    source("03_occurrences.R")        # kayıt temizliği, 30 m hücre tekleme, yanlılık ağırlıklı arka plan
    source("05_downscale_climate.R")  # iklimi 30 m'ye indirger (regresyon + delta)
    source("07_variable_selection.R") # blok AUC + Spearman kümeleri + VIF → M1, M2, M3
    source("08_model_inputs.R")       # model klasörleri ve runner config.json dosyaları
    source("09_run_models.R")         # biomod2-sdm-runner (runner/, v1.2) ile kalibrasyon
    source("09b_project.R")           # tam veri modelleriyle projeksiyon + GCM ansamblı
    source("10_post_analyses.R")      # performans karşılaştırması, MESS, alan değişimi,
                                      # tamponlama, mikrosığınak, belirsizlik
    source("11_figures.R")            # makale şekilleri (results/figures)

Terminal: Rscript 09_run_models.R M3 ; Rscript 09b_project.R M3 M1

RUN_DOWNLOADS.R yalnız indirme adımlarını (01 + 04) çalıştırır. WorldClim 2.1 alternatifi
01_download_climate_worldclim.R içindedir; makalede CHELSA kullanıldı.

## Modelleme kurgusu
- M1 iklim, M2 + topoğrafya, M3 + kanopi ve toprak (iç içe setler).
- GLM, GAM, GBM, RFd, MAXNET; 3 × 2000 sözde-yok (yanlılık ağırlıklı arka plandan, σ = 5 km);
  4 mekânsal blok CV (Muscarella ve ark. 2014; runner'da user.defined tablo).
- Topluluk: blok AUC ≥ 0.60 ve TSS ≥ 0.30 olan algoritmaların tam veri modelleri,
  blok TSS ile ağırlıklı. M1'de eşik yok (GLM + GAM).
- Eşik: güncel tahminde kayıtlar ve arka plan arasında maxTSS; tüm senaryolarda sabit.
- Projeksiyon: değişkenler kalibrasyon aralığına sıkıştırılır; GCM'ler tahmin düzeyinde
  birleştirilir (ortalama, SD, GCM uyumu; uyum ≥ 0.8 = uygun).
- runner/run_biomod.R GitHub sürümüne göre eklemeler içerir (v1.2, yerel): cv_strategy
  "block", bg_file, project_full_models, skip_projection, terra_memfrac, clamping_mask.

## Çıktılar
  data/climate_30m/      30 m'ye indirgenmiş iklim
  data/topo/, data/soil/ 30 m arazi, kanopi ve toprak katmanları
  data/selection/        değişken seçimi tabloları, dendrogram, model_sets.txt
  results/<M>/           runner çıktıları (değerlendirme, önem, tepki eğrileri)
  results/<M>/proj/      current_suitability/binary, <senaryo>_suitability,
                         gcm_<dönem>_<ssp>_{mean,sd,presence_agreement}.tif, eşik ve ağırlıklar
  results/summary/       performance_*.csv, suitable_area_change.csv, mess_*.csv,
                         topographic_buffering.csv, microrefugia_*, uncertainty_partitioning_*.csv
  results/figures/       Fig1–Fig5 (300 dpi)

## Notlar
- logs/downscale_coefficients.csv: yükseltiyle azalma oranları ve R² (ek tablo).
- DTM'nin DSM'ye yatay hizalanması kodda yapılır (logs/terrain_prep.csv).
