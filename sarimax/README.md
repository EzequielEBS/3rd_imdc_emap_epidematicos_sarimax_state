# SARIMAX Models — Dengue & Chikungunya (State Level)

Code in this folder fits, selects, and submits the state-level SARIMAX models used for this challenge entry. It reads the per-state files produced by `data_prep/agg_data_uf.r` (`processed_data/dengue/dengue_<UF>_agg.csv.gz`, `processed_data/chikungunya/chikungunya_<UF>_agg.csv.gz`) and writes everything under `sarimax/results/`.

## Files

| File | Role |
|---|---|
| `src/utils.r` | Shared library: `fit_sarimax()` (core SARIMAX fit + forecast), covariate-selection helpers (`get_candidates()`, `filter_low_variance()`, `filter_by_correlation()`, `filter_redundant_indices()`, `select_best_per_variable()`, `build_covariate_combinations()`, `pca_all()`), the grid-search engine (`run_grid_search()`), and the scoring helper (`compute_metrics()`, computing WIS/MAE/RMSE/MAPE/coverage). Sourced by every other script. |
| `src/model_sel.r` | **Training / model selection.** Defines `run_model_selection()` and loops it over all states (and a few focal cities) for both diseases, writing per-state CV metrics and a `best_wis_<disease>_all_states.csv` summary. |
| `src/fit.r` | **Final fit & forecast generation.** Refits each state's selected best model on all four train/target splits and writes one forecast CSV per state x target window to `results/preds/`. Includes the retry loop for chikungunya states whose initial fit produced convergence/NaN-SE warnings, and `fit_sarimax_epiweek()`, the variant used for the real submission window (split 4 — see the root README's "Data Usage Restriction" section). |
| `src/sub_pred.r` | **Submission.** Uploads the forecast CSVs in `results/preds/` to the Mosqlimate Predictions Registry via the `mosqlient` Python package (through `reticulate`). |
| `results/concluded_states_*.csv` | Checkpoints so `model_sel.r` can be safely interrupted and resumed without redoing finished states. |
| `results/metrics/` | One `metrics_all_formulas_<disease>_<state>.csv` per state/disease (every formula x order combination tried, ranked by mean CV metric), plus `best_wis_<disease>_all_states.csv` (the winning combination per state) and `mean_wis_best_model_<disease>_<state>.csv` (final per-split metrics of the chosen model). |
| `results/preds/` | `pred_<disease>_<UF>_target_<n>.csv` — the actual forecast files (date, point prediction, and 50/80/90/95% interval bounds), ready to submit. |

## How to run

From the repository root, with R and the packages listed in the root README's "Libraries and Dependencies" section installed:

```r
# Training / model selection for every state (and focal city), both diseases.
# Long-running (grid search x 26 states x 2 diseases). Safe to interrupt and
# re-run: already-finished states are skipped via results/concluded_states_*.csv.
source("sarimax/src/model_sel.r")

# Refit the winning model per state and generate the 4 forecast windows.
source("sarimax/src/fit.r")

# Submit the forecasts to Mosqlimate (see "Credentials" below).
source("sarimax/src/sub_pred.r")
```

To experiment with a single state without running the full loop, call `run_model_selection()` directly after sourcing `src/utils.r` and `src/model_sel.r`'s function definition, e.g.:

```r
source("sarimax/src/utils.r")
source("sarimax/src/model_sel.r")  # defines run_model_selection(); the script's
                                    # own loop will also start running unless you
                                    # source only the function definition block
concluded_states <- data.frame(state = character(0))
result <- run_model_selection(state = "SP", concluded_states = concluded_states)
```

## Credentials for `sub_pred.r`

`sub_pred.r` never hardcodes the Mosqlimate API key. Create a `.env` file in the **repository root** (already excluded from Git via `.gitignore`) with:

```
MOSQLIMATE_API_KEY=your-key-here
```

Get your key from the "Auth" section of your Mosqlimate profile at [mosqlimate.org](https://mosqlimate.org/). Also update the `repository` and `commit` constants at the top of `sub_pred.r` to match the GitHub repository and commit hash you are actually submitting.

## Modeling summary

- Response: `log1p(cases)`, back-transformed with `expm1()` and clipped at zero.
- Covariates: PCA components of de-correlated weather variables (contemporaneous, lagged, and rolling-mean), plus ENSO/IOD/PDO when they add signal beyond the retained weather components — see the root README's "Data and Variables" section for the full selection procedure.
- SARIMAX order: grid-searched per state, bounded by `max_order` (default `p≤2, d≤1, q≤2, P≤1, D=0, Q≤1`), fit with `method = "CSS-ML"`, `optim.method = "BFGS"`.
- Selection metric: mean Weighted Interval Score (WIS) across 4 cross-validation splits (configurable via the `metric` argument of `run_model_selection()`).
- Prediction intervals: simulated (`bootstrap = TRUE`, `npaths = 1000`) at 50/80/90/95% coverage.

See the root `README.md` for the full Team, Data, Model Training, Data Usage Restriction, Predictive Uncertainty, and References sections required by the IMDC submission format.
