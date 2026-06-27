# SARIMAX Dengue Models

This folder contains code to fit SARIMAX models for the state-level dengue
aggregate files in `processed_data/dengue`.

Run all states and validation windows with the full shared-model grid:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' sarimax/model_sel.R
```

Fast smoke test for one state and one validation split:

```powershell
& 'C:\Program Files\R\R-4.4.1\bin\Rscript.exe' sarimax/fit_sarimax.R --states SP --splits 3 --quick
```

Outputs are written to `sarimax/results`:

- `sarimax_metrics.csv`: one row per state and validation split, including MAE,
  RMSE, MAPE, bias, globally chosen SARIMAX orders, and globally selected
  covariates.
- `sarimax_predictions.csv`: one row per validation epiweek with observed,
  predicted, error, and absolute error.
- `sarimax_model_grid.csv`: one row per shared model specification evaluated
  during global selection.

The script first selects a single global, de-correlated covariate ranking using
only training rows across the requested states and splits. It then evaluates
shared SARIMAX specifications using a general validation metric, `mean_mae` by
default, and writes final metrics/predictions for the best shared specification.

Covariates are standardized separately within each state/split using training
rows only. The response is modeled as `log1p(cases)`, and forecasts are converted
back to case counts with `expm1`.

Useful options:

```powershell
--metric mean_mae
--metric median_mae
--metric total_absolute_error
--covariate-counts 0 4 8 12
--max-covariates 12
--corr-threshold 0.75
```
