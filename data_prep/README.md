# Data Preparation Pipeline

Scripts in this folder turn the raw Mosqlimate/DATASUS extracts in `raw_data/` into the state-level, modeling-ready tables in `processed_data/` that `sarimax/src/model_sel.r` and `sarimax/src/fit.r` read directly. **Run them in order** from the repository root (each step writes the input the next step expects):

| Order | Script | Reads | Writes | Purpose |
|---|---|---|---|---|
| 1 | `prep_climate.r` | `raw_data/climate.csv.gz` | `processed_data/climate/climate.csv.gz` | Gap-fills climate series for island municipalities (no weather coverage of their own) by copying their nearest mainland neighbor's series under the island's geocode (`sf`/`geobr` centroid distance). |
| 2a | `merge_dengue.r` | `raw_data/dengue.csv.gz`, `processed_data/climate/...`, `raw_data/environ_vars.csv.gz`, `raw_data/ocean_climate_oscillations.csv.gz`, `raw_data/datasus_population_2001_2025.csv.gz`, `raw_data/map_regional_health.csv` | `processed_data/dengue/dengue_merged.csv.gz` | Left-joins dengue cases with climate (by geocode+epiweek), environmental attributes (by geocode+uf_code), ocean-climate indices (by epiweek), and population (by geocode+year, extrapolated to 2026). |
| 2b | `merge_chikungunya.r` | same as 2a, chikungunya | `processed_data/chikungunya/chikungunya_merged.csv.gz` | Identical procedure to `merge_dengue.r`, for chikungunya. |
| 3 | `handle_na.r` | the two `*_merged.csv.gz` files | overwrites them in place | Linearly interpolates remaining gaps in ENSO/IOD/PDO (`zoo::na.approx`); these are the only systematic missing values left after the merge step. |
| 4 | `agg_data_uf.r` | the two `*_merged.csv.gz` files | `processed_data/<disease>/<disease>_<UF>_agg.csv.gz` (one per state, both diseases) | Aggregates city-level rows up to state (UF) level (sums cases, averages weather, takes the modal Köppen/biome), then engineers 4/8/12/16-week lagged and 3/6/9/12-month rolling-mean versions of every weather covariate. **This is the final modeling input.** Espírito Santo (`ES`) is excluded. |
| optional | `sel_cities.r` | the two `*_merged.csv.gz` files | `processed_data/<disease>/sel_cities/<disease>_<geocode>.csv.gz` | Extracts a fixed list of focal cities for city-level diagnostics; not part of the state-level submission data path. |
| n/a | `summarise_sel_cities.r` | — | — | Placeholder, currently empty / not implemented. Nothing else in the repository depends on it. |

## Key columns produced

- `cases`: weekly probable case count (response variable), summed to state level.
- `*_mean`: state-level average of each climate variable (temperature, precipitation, pressure, relative humidity, thermal range, rainy days).
- `*_mean_lag{4,8,12,16}`: the same variables lagged by 4–16 weeks.
- `*_mean_mean_{3,6,9,12}mo`: rolling means of the same variables over 3–12 months, lagged by one period so no current-week information leaks in.
- `enso`, `iod`, `pdo`: ocean-climate oscillation indices (national scale).
- `koppen`, `biome`: most common climate classification / biome among the state's cities (categorical; not currently used as model regressors, only carried through for reference).
- `train_1..4`, `target_1..4`: boolean indicators defining the four train/forecast splits used throughout `sarimax/`. See the root `README.md`, Section 6 ("Data Usage Restriction"), for exactly how these map to the EW25→EW41/EW40 challenge rule.

## Re-running the pipeline

```r
source("data_prep/prep_climate.r")
source("data_prep/merge_dengue.r")
source("data_prep/merge_chikungunya.r")
source("data_prep/handle_na.r")
source("data_prep/agg_data_uf.r")
source("data_prep/sel_cities.r")  # optional
```

Required R packages: `tidyverse`, `lubridate`, `aweek`, `zoo`, `runner`, `sf`, `geobr` (see the root README's "Libraries and Dependencies" section).
