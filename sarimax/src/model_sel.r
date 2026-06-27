#' model_sel.r — Cross-validated covariate and SARIMAX-order selection
#'
#' Orchestrates the full model-selection pipeline (covariate screening →
#' dimensionality reduction → grid search → best-model evaluation) for every
#' state (and, in the driver code at the bottom of this file, every
#' selected city) and disease. Outputs are written to `sarimax/results/`:
#' one `metrics_all_formulas_<disease>_<state>.csv` per state/disease
#' (every formula x order combination tried, ranked by mean CV metric) and
#' a `best_wis_<disease>_all_states.csv` summary used downstream by
#' `fit.r` to refit the winning model and generate the submission forecasts.
source("sarimax/src/utils.r")

#' Select the best covariate set and SARIMAX order for one state or city
#'
#' For the given disease's aggregated table, this: (1) lists candidate
#' covariates with `get_candidates()`; (2) drops near-constant
#' (`filter_low_variance()`) and weakly-correlated
#' (`filter_by_correlation()`) ones, using only `train_1` rows to avoid
#' leakage; (3) either reduces the remainder to a handful of PCA components
#' (`pca_all()`, re-testing ocean-climate indices with
#' `filter_redundant_indices()`) when `pca = TRUE`, or falls back to an
#' explicit, de-correlated combinatorial search
#' (`select_best_per_variable()` + `build_covariate_combinations()`) when
#' `pca = FALSE`; (4) runs `run_grid_search()` over the resulting formulas
#' and a grid of SARIMAX orders, scored by cross-validated `metric` (WIS by
#' default) across the four train/target splits; and (5) refits the overall
#' best formula/order combination to report final per-split metrics. Skips
#' states/cities already present in `concluded_states`/`concluded_cities` so
#' a long run can be safely resumed.
#'
#' @param state, city                Exactly one should be supplied: a UF
#'                                    code (state-level) or an IBGE geocode
#'                                    (city-level) to process.
#' @param concluded_states, concluded_cities  Data frames of already-
#'                                    processed state/city codes, used to
#'                                    skip work that's already done.
#' @param disease                    "dengue" or "chikungunya"; selects
#'                                    which `processed_data/<disease>/...`
#'                                    file to read.
#' @param metric                     Column name in the metrics tibble used
#'                                    to rank formula x order combinations
#'                                    (default "wis", giving `mean_wis`).
#' @param sample_size                When `pca = FALSE`, number of covariate
#'                                    combinations to randomly sample for the
#'                                    grid search (full combinatorial search
#'                                    is otherwise too large).
#' @param threshold_low_variance, threshold_cor, min_cor  Thresholds passed
#'                                    to `filter_low_variance()`,
#'                                    `build_covariate_combinations()`, and
#'                                    `filter_by_correlation()` respectively.
#' @param max_size_covariates        Max covariates per combination when
#'                                    `pca = FALSE`.
#' @param levels                     Prediction-interval coverage levels
#'                                    (percent) used throughout.
#' @param max_order                  Upper bound on each SARIMAX order
#'                                    component for the grid search.
#' @param n_cores                    Parallel workers for `run_grid_search()`.
#' @param method, lambda             Passed to `forecast::Arima()`.
#' @param pca, pca_var_threshold, k  Whether to use PCA for dimensionality
#'                                    reduction, the variance threshold for
#'                                    choosing the number of components, and
#'                                    the max number of leading PCs to try
#'                                    (formulas are built for 1..k PCs).
#' @param index_cor_threshold        Passed to `filter_redundant_indices()`.
#' @param fixed_stat_par             If TRUE, lock d/D to the
#'                                    `determine_d()` recommendation instead
#'                                    of searching them in the order grid.
#' @param bootstrap, npaths          Passed to `fit_sarimax()` /
#'                                    `run_grid_search()` for simulation-based
#'                                    prediction intervals.
#' @param concluded_states_path, concluded_cities_path  Checkpoint files to
#'                                    append to after a successful run
#'                                    (overridable so test runs don't touch
#'                                    the real checkpoint).
#'
#' @return A list with `best_order`, `best_formula`, `best_pred` (predictions
#'   for the winning formula/order), `final_metrics` (per-split metrics for
#'   the winning model), and `all_metrics` (every formula x order combination
#'   tried, ranked by mean CV metric). Returns NULL (and skips all work) if
#'   the state/city was already concluded.
run_model_selection <- function(
  state = NULL,
  concluded_states = NULL,
  city = NULL,
  concluded_cities = NULL,
  disease = "dengue",
  metric = "wis",
  sample_size = 10,
  threshold_low_variance = 0.01,
  threshold_cor = 0.6,
  min_cor = 0.1,
  max_size_covariates = 3,
  levels = c(50, 80, 90, 95),
  max_order = c(p = 2, d = 1, q = 2, P = 1, D = 0, Q = 1),
  n_cores = parallel::detectCores() - 1,
  method = "CSS-ML",
  lambda = NULL,
  pca = T,
  pca_var_threshold = 0.90,
  k = 5,
  index_cor_threshold = 0.3,  # min |partial cor| for enso/iod/pdo to be kept (see filter_redundant_indices)
  fixed_stat_par = F,         # search d/D via CV WIS instead of locking to a KPSS/OCSB pick
  bootstrap = TRUE,           # simulate forecast paths for better-calibrated intervals
  npaths = 1000,
  concluded_states_path = "sarimax/results/concluded_states_dengue.csv",  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
  concluded_cities_path = "sarimax/results/concluded_cities_dengue.csv"  # overridable so test
                                                                   # runs don't touch the
                                                                   # real checkpoint file
) {
  if (!is.null(state) | !is.null(city)) {
    if (!is.null(state)) {
      if (state %in% concluded_states$state) {
        message("State ", state, " already processed. Skipping.")
        return(NULL)
      }
      message("Processing state: ", state)
      file_name <- paste0("processed_data/", disease,"/", disease, "_", state, "_agg.csv.gz")
    } else if (!is.null(city)) {
      if (city %in% concluded_cities$city) {
        message("City ", city, " already processed. Skipping.")
        return(NULL)
      }
      message("Processing city: ", city)
      file_name <- paste0("processed_data/", disease,"/sel_cities/", disease, "_", city, ".csv.gz")
    }
  }
  
  dengue_state <- read_csv(file_name, show_col_types = FALSE)

  train_ids  <- paste0("train_", 1:4)
  
  candidates <- get_candidates(dengue_state)
  candidates <- filter_low_variance(
    dengue_state[
      dengue_state[[train_ids[1]]] == 1, 
    ],
    candidates, 
    threshold = threshold_low_variance
  )
  candidates <- filter_by_correlation(
    dengue_state[
      dengue_state[[train_ids[1]]] == 1, 
    ],
    candidates, 
    min_cor = min_cor
  )

  if (pca) {
    pca_result <- pca_all(
      data = dengue_state,
      candidates = candidates[!grepl("enso|iod|pdo", candidates)],
      var_threshold = pca_var_threshold
    )
    dengue_state <- pca_result$data
    pcs <- pca_result$variables
    # Only go up to the components selected by var_threshold
    max_k <- min(k, length(pcs))

    # ENSO/IOD/PDO were deliberately excluded from the PCA matrix above so they
    # aren't blended away by an unsupervised projection. Test whether they carry
    # signal beyond what the retained local-weather PCs already explain, and add
    # back any that do, instead of dropping them outright.
    climate_indices <- intersect(c("enso", "iod", "pdo"), candidates)
    kept_indices <- character(0)
    if (length(climate_indices) > 0) {
      retained <- filter_redundant_indices(
        data       = dengue_state,
        covariates = c(pcs, climate_indices),
        indices    = climate_indices,
        threshold  = index_cor_threshold
      )
      kept_indices <- intersect(climate_indices, retained)
    }

    formulas <- lapply(seq_len(max_k), function(i) {
      reformulate(pcs[1:i], response = "cases")
    })
    if (length(kept_indices) > 0) {
      formulas <- c(formulas, lapply(seq_len(max_k), function(i) {
        reformulate(c(pcs[1:i], kept_indices), response = "cases")
      }))
    }
    candidates <- c(pcs, kept_indices)
  } else {
    candidates <- select_best_per_variable(
      dengue_state[
        dengue_state[[train_ids[1]]] == 1, 
      ],
      candidates
    )  
    combos <- build_covariate_combinations(
      data       = dengue_state,
      covariates = candidates,
      threshold  = threshold_cor,
      max_size   = max_size_covariates
    )
    formulas   <- lapply(combos, \(vars) reformulate(vars, response = "cases"))
  }
  
  if (pca) {
    sample_formulas <- formulas
  } else {
    sample_formulas <- formulas[sample(length(formulas), min(sample_size, length(formulas)))]
  }
  
  order_screen <- run_grid_search(
    data      = dengue_state,
    formulas  = sample_formulas,
    levels    = levels,
    max_order = max_order,
    n_cores   = n_cores,
    method    = method,
    lambda    = lambda,
    fixed_stat_par = fixed_stat_par,
    bootstrap = bootstrap,
    npaths    = npaths
  )

  metric_name <- paste0("mean_", metric)
  best_par <- order_screen$metrics |>
    group_by(formula_id, order) |>
    filter(all(!failed)) |>
    summarise(mean_metric = mean(get(metric)), .groups = "drop") |>
    arrange(mean_metric) |>
    slice(1)
  best_order <- best_par$order
  ord <- parse_order(best_order)

  if (pca) {
    best_formula <- best_par$formula_id
    metrics <- order_screen$metrics |>
      group_by(formula_id, order) |>
      filter(all(!failed)) |>
      summarise(mean_metric = mean(get(metric)), .groups = "drop") |>
      arrange(mean_metric)
    best_pred <- order_screen$predictions |>
      filter(formula_id == best_formula, order == best_order)
  } else {
    formula_results <- run_grid_search(
      data      = dengue_state,
      formulas  = formulas,
      levels    = levels,
      fixed_order = c(p = ord$order[1], d = ord$order[2], q = ord$order[3],
                      P = ord$seasonal_order[1], D = ord$seasonal_order[2], Q = ord$seasonal_order[3]),
      n_cores   = n_cores,
      method    = method,
      lambda    = lambda,
      bootstrap = bootstrap,
      npaths    = npaths
    )
    metrics <- formula_results$metrics |>
      group_by(formula_id, order) |>
      filter(all(!failed)) |>        # keep only groups where ALL splits succeeded
      summarise(mean_metric = mean(get(metric)), .groups = "drop") |>
      arrange(mean_metric)

    best_formula <- metrics$formula_id[1]
    best_pred <- formula_results$predictions |>
      filter(formula_id == best_formula, order == best_order)
  }

  train_ids <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)

  final_metrics <- lapply(seq_along(train_ids), function(i) {
    actual    <- dengue_state$cases[dengue_state[[target