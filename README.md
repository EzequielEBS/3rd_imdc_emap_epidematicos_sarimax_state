# 3rd Infodengue-Mosqlimate Dengue Challenge (IMDC) 2026 — State-level SARIMAX

Submission repository for the **2026 3rd Infodengue-Mosqlimate Dengue Challenge**: state-level (UF) weekly probable-case forecasts for **dengue** and **chikungunya** in Brazil, produced with SARIMAX models fit per state.

## 1. Team and Contributors

**Team:** Epidemáticos

**Affiliation:** School of Applied Mathematics, Fundação Getulio Vargas (FGV/EMAp), Rio de Janeiro, Brazil

**Contributors:**

- Eduardo Adame, M.Sc. — FGV/EMAp
- Ezequiel Braga, M.Sc. — FGV/EMAp
- Iara Cristina, M.Sc. — FGV/EMAp
- Isaque Pim, Ph.D. — FGV/EMAp

## 2. Repository Structure

```
.
├── README.md                   # this file
├── LICENSE                     # GNU GPLv3
├── pyproject.toml              # Python dependencies (Poetry), used only by sub_pred.r via reticulate
├── .gitattributes              # marks *.csv / *.gz as Git LFS objects
├── .gitignore
│
├── data_prep/                  # data preparation pipeline (R), run in this order
│   ├── prep_climate.r              # 1. fills climate gaps for island municipalities
│   ├── merge_dengue.r               # 2a. joins raw dengue cases with climate/env/ocean/population covariates
│   ├── merge_chikungunya.r          # 2b. same as above, for chikungunya
│   ├── handle_na.r                  # 3. interpolates missing ENSO/IOD/PDO values
│   ├── agg_data_uf.r                # 4. aggregates city-level data to state (UF) level + lag/rolling features
│   ├── sel_cities.r                 # optional: extracts a few focal cities for diagnostics
│   ├── summarise_sel_cities.r       # placeholder, currently empty / unused
│   └── README.md                    # details on each script above
│
├── raw_data/                   # original input datasets (Git LFS)
│   ├── dengue.csv.gz, chikungunya.csv.gz        # case counts by geocode/epiweek + train/target split flags
│   ├── climate.csv.gz                            # weather variables by geocode/epiweek
│   ├── environ_vars.csv.gz                       # Köppen climate classification, biome (static, by geocode)
│   ├── ocean_climate_oscillations.csv.gz         # ENSO, IOD, PDO indices by epiweek (national scale)
│   ├── datasus_population_2001_2025.csv.gz       # DATASUS population by geocode/year
│   ├── map_regional_health.csv                   # geocode → regional/macroregional health mapping
│   ├── shape_muni.gpkg, shape_regional_health.gpkg, shape_macroregional_health.gpkg  # spatial boundaries (currently unused as model inputs)
│   └── forecasting_climate.csv.gz                # forecasted climate variables (currently unused by the pipeline)
│
├── processed_data/              # output of data_prep/ — the modeling-ready tables
│   ├── climate/climate.csv.gz
│   ├── dengue/
│   │   ├── dengue_merged.csv.gz             # city-level, all covariates merged
│   │   ├── dengue_<UF>_agg.csv.gz           # one file per state — read directly by the SARIMAX scripts
│   │   └── sel_cities/                      # focal-city extracts (diagnostics only)
│   └── chikungunya/                         # same structure as dengue/
│
└── sarimax/                     # modeling code and results
    ├── README.md                            # modeling workflow and how-to-run instructions
    ├── src/
    │   ├── utils.r                          # shared helpers: fitting, covariate selection, metrics, PCA
    │   ├── model_sel.r                      # cross-validated covariate + SARIMAX-order selection (training)
    │   ├── fit.r                            # refits the best model per state and generates forecasts
    │   └── sub_pred.r                       # uploads forecasts to the Mosqlimate platform
    └── results/
        ├── concluded_states_dengue.csv, concluded_states_chikungunya.csv  # resumability checkpoints
        ├── metrics/                          # per-state/disease CV metrics, incl. best_wis_*_all_states.csv
        └── preds/                            # final forecast CSVs, one per state x target window
```

## 3. Libraries and Dependencies

**R** (used throughout `data_prep/` and `sarimax/src/`):

| Package | Used for |
|---|---|
| `tidyverse` (`dplyr`, `tidyr`, `purrr`, `readr`) | data wrangling |
| `lubridate`, `aweek` | date and epiweek handling |
| `zoo` | linear interpolation of missing ocean-climate index values |
| `runner` | rolling-window (3/6/9/12-month) covariate means |
| `sf`, `geobr` | municipality geometries/centroids, used to gap-fill island climate data |
| `forecast` | `Arima()` / `forecast()` — the SARIMAX engine |
| `scoringutils` | Weighted Interval Score (WIS) and interval-coverage scoring |
| `parallel`, `pbapply` | parallelized, progress-tracked grid search |
| `reticulate`, `data.table` | Python interop and fast I/O for the submission step |
| `dotenv` | loads the Mosqlimate API key from a local `.env` file (never hardcoded) |

**Python** (declared in `pyproject.toml`, only used via `reticulate` in `sarimax/src/sub_pred.r` to submit forecasts):

- `python = ">=3.10,<3.13"`
- `mosqlient==1.5.2` — Mosqlimate Predictions Registry client
- `epiweeks`, `python-dotenv`
- `jupyter`, `pmdarima` are declared but not currently invoked anywhere in this repository's code path — they are leftovers from the challenge template and can be removed if unused in your own workflow.

Install the R packages with, e.g., `install.packages(c("tidyverse","lubridate","aweek","zoo","runner","sf","geobr","forecast","scoringutils","pbapply","reticulate","data.table","dotenv"))`, and the Python side with `poetry install` (or `pip install mosqlient==1.5.2 epiweeks python-dotenv`).

## 4. Data and Variables

### Datasets

All disease and most covariate data are distributed by the Mosqlimate platform (Infodengue feed, probable cases — the `casprov`/`target_*` columns already reflect this); population comes from DATASUS.

- **Case data** — `raw_data/dengue.csv.gz`, `raw_data/chikungunya.csv.gz`: weekly probable case counts by municipality (`geocode`) and `epiweek`, with the official `train_1..4` / `target_1..4` split indicators used for both retrospective validation and the real submission (see Section 6).
- **Climate** — `raw_data/climate.csv.gz`: minimum/mean/maximum temperature, precipitation, atmospheric pressure, relative humidity, thermal range, and rainy-day count, by municipality and epiweek.
- **Environmental attributes** — `raw_data/environ_vars.csv.gz`: Köppen climate classification and biome, static per municipality.
- **Ocean-climate indices** — `raw_data/ocean_climate_oscillations.csv.gz`: ENSO, IOD, and PDO, indexed by epiweek (national scale, not municipality-specific).
- **Population** — `raw_data/datasus_population_2001_2025.csv.gz`: municipality population by year (2026 obtained by carrying the 2025 value forward — see `data_prep/merge_dengue.r`).
- **Regional-health mapping and shapefiles** — `raw_data/map_regional_health.csv` and the `shape_*.gpkg` files are read in by the merge scripts / available for spatial work but are **not** used as direct model covariates in the current pipeline; `sf`/`geobr` are used only in `data_prep/prep_climate.r` to fill climate gaps for island municipalities by borrowing their nearest mainland neighbor's series.
- `raw_data/forecasting_climate.csv.gz` is included but currently unused by any script.

### Pre-processing

Implemented in `data_prep/` (see `data_prep/README.md` for full detail), in order:

1. `prep_climate.r` — fills missing climate series for island municipalities.
2. `merge_dengue.r` / `merge_chikungunya.r` — left-joins case data with climate, environmental, ocean-climate, and population covariates.
3. `handle_na.r` — linearly interpolates the remaining gaps in ENSO/IOD/PDO.
4. `agg_data_uf.r` — sums cases and averages weather covariates from city to state (UF) level, and engineers 4/8/12/16-week lagged and 3/6/9/12-month rolling-mean versions of every weather covariate. Writes the final `processed_data/<disease>/<disease>_<UF>_agg.csv.gz` files used for modeling. Espírito Santo (ES) is excluded.

Inside model fitting itself, covariates are additionally standardized (z-scored) using the **training split's mean/SD only**, to avoid leaking information from the forecast window into the scaling — see `fit_sarimax()` in `sarimax/src/utils.r` (around the "Build & standardize regressor matrices" step) and `fit_sarimax_epiweek()` in `sarimax/src/fit.r`.

### Variable selection

Variable selection is fully automated per state/disease inside `run_model_selection()` (`sarimax/src/model_sel.r`), using only training-split rows at every filtering step to avoid leakage:

1. **`get_candidates()`** (`sarimax/src/utils.r`) enumerates the contemporaneous, lagged, and rolling-mean versions of the 8 base covariates (`temp_med_mean`, `precip_med_mean`, `rel_humid_med_mean`, `thermal_range_mean`, `rainy_days_mean`, `enso`, `iod`, `pdo`).
2. **`filter_low_variance()`** drops near-constant covariates (SD below a threshold).
3. **`filter_by_correlation()`** drops covariates whose absolute correlation with `log1p(cases)` is below a minimum threshold.
4. **Dimensionality reduction** — by default (`pca = TRUE`): **`pca_all()`** reduces the surviving local-weather covariates to the smallest number of principal components that explain ≥ 90% of variance in the most restrictive of the four CV folds. ENSO/IOD/PDO are deliberately excluded from the PCA and tested separately with **`filter_redundant_indices()`**, which keeps an index only if it still correlates with the response after the local-weather PCs' linear effect is removed — so ocean-climate signal isn't silently absorbed (or discarded) by the PCA. Candidate formulas use the leading 1..k PCs, with and without the surviving indices.
   An alternative, non-PCA path (`pca = FALSE`) is also implemented: **`select_best_per_variable()`** keeps only the best lag/rolling version of each base variable, and **`build_covariate_combinations()`** enumerates de-correlated covariate subsets (pairwise |correlation| below a threshold) up to a maximum size, from which a random sample is grid-searched.
5. **`run_grid_search()`** (`sarimax/src/utils.r`) cross-validates every candidate formula against a grid of SARIMAX `(p,d,q)(P,D,Q)` orders across the four train/target splits, and **`run_model_selection()`** picks the formula/order combination with the lowest mean cross-validated score (WIS by default).

## 5. Model Training

Each state and disease is modeled independently with 