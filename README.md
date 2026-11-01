# Fine-scale habitat modelling of *Morchella esculenta* s.l. in Yenişarbademli (Isparta, Türkiye)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23861480.svg)](https://doi.org/10.5281/zenodo.23861480)

Code and outputs for the manuscript:

> Kaya, C., Özdemir, S., & Aday Kaya, A. G. Soil and canopy define the habitat of the yellow morel (*Morchella esculenta* s.l.) but do not buffer its projected decline under climate change. *In preparation.*

The workflow builds a 30-m species distribution model from field presence records (2022–2023), CHELSA v2.1 climate downscaled to 30 m, terrain and canopy variables from the national terrain (SAM) and surface (SYM) models of the General Directorate of Mapping (HGM), and SoilGrids 2.0 soil properties. Models are fitted with biomod2 through [biomod2-sdm-runner](https://github.com/srknozdemir/biomod2-sdm-runner) using spatial block cross-validation and a sampling-bias-weighted background, and projected to five CMIP6 GCMs under SSP1-2.6, SSP3-7.0 and SSP5-8.5 for 2041–2070 and 2071–2100.

## Repository layout

| Path | Content |
|---|---|
| `scripts/00_config.R` | Paths, study area, resolution and all modelling settings |
| `scripts/01_download_climate.R` | CHELSA v2.1 baseline and CMIP6 climatologies (cropped download) |
| `scripts/01_download_climate_worldclim.R` | Alternative: WorldClim 2.1 (not used in the paper) |
| `scripts/02a_prepare_terrain.R` | HGM terrain (SAM) and surface (SYM) models: mosaicking, alignment |
| `scripts/02b_twi.py` | Topographic wetness index with D-infinity routing (pysheds) |
| `scripts/02_dem_topography.R` | Terrain, radiation and canopy variables at 30 m |
| `scripts/03_occurrences.R` | Record cleaning, one-per-cell filtering, bias-weighted background |
| `scripts/04_soil.R` | SoilGrids 2.0, thickness-weighted 0–30 cm means |
| `scripts/05_downscale_climate.R` | Elevation-regression downscaling (1 km → 30 m) and delta method |
| `scripts/07_variable_selection.R` | Block-AUC screening, Spearman clustering, VIF; nested model sets |
| `scripts/08_model_inputs.R` | Model inputs and runner configuration files |
| `scripts/09_run_models.R` | Calibration with biomod2-sdm-runner |
| `scripts/09b_project.R` | Ensemble projection, clamping and GCM aggregation |
| `scripts/10_post_analyses.R` | Performance tests, MESS, area change, buffering, microrefugia, uncertainty |
| `scripts/11_figures.R`, `scripts/12_manuscript_figures.R` | Figures and tables |
| `scripts/13_canopy_sensitivity.R` | Sensitivity of the projections to a decline in canopy cover |
| `scripts/runner/` | biomod2-sdm-runner v1.2 used for the analyses |
| `patches/` | Changes of runner v1.2 relative to the public v1.1 (spatial block CV, background file, full-data models) |
| `data/occurrences_1km_generalised.csv` | Presence records generalised to 1-km cells (see Data note) |
| `results/` | Summary tables, manuscript tables and figures (PNG, PDF) |

## Running the workflow

Requirements: R ≥ 4.3 with terra, biomod2 (4.3-4-8), mgcv, gbm, randomForest, maxnet, ggplot2, patchwork, ragg; Python 3 with pysheds for the wetness index.

1. Place the record file (`kuzu.csv`, columns `Sp;Enlem;Boylam`), the study-area polygon and the HGM terrain and surface models as described in `scripts/00_config.R`.
2. Run the scripts in numerical order from the `scripts/` directory, e.g. `Rscript 01_download_climate.R`. Downloads resume if interrupted.
3. Calibration and projection: `Rscript 09_run_models.R M1 M2 M3`, then `Rscript 09b_project.R M3 M1`.
4. Post-hoc analyses and figures: `Rscript 10_post_analyses.R`, `Rscript 12_manuscript_figures.R`.

A Turkish description of the workflow is given in `scripts/BENIOKU_tr.md`.

## Data note

The exact locations of harvested morel stands are commercially sensitive for local communities. The repository therefore provides the records generalised to 1-km cells (UTM 36N cell centres with the number of 30-m records per cell). Exact coordinates are available from the corresponding author (serkanozdemir@isparta.edu.tr) on reasonable request. The HGM terrain and surface models are licensed data and are not redistributed; CHELSA and SoilGrids layers are downloaded by the scripts.

## Citation

Please cite the article and this repository:

> Kaya, C., Özdemir, S., Aday Kaya, A.G., 2026. morchella-sdm-yenisarbademli: fine-scale habitat modelling of *Morchella esculenta* s.l. in south-western Anatolia (v1.0.0) [Software]. Zenodo. https://doi.org/10.5281/zenodo.23861480

## Licence

Code: MIT licence. Outputs in `results/`: CC BY 4.0.
