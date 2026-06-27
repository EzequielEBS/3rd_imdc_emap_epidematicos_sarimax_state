source("sarimax/src/utils.r")

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
    actual    <- dengue_state$cases[dengue_state[[target_ids[i]]] == 1]
    pred     <- best_pred |> filter(train_id == train_ids[i])
    pred_df  <- pred |> select(date, pred, starts_with("lower_"), starts_with("upper_"))
    compute_metrics(pred, actual)
  }) |> bind_rows() |> mutate(target_id = target_ids)
  # save files
  write_csv(final_metrics, file.path("sarimax/results/metrics/", paste0(metric_name, "_best_model_", disease, "_", state, ".csv")))
  write_csv(metrics, file.path("sarimax/results/metrics/", paste0("metrics_all_formulas_", disease, "_", state, ".csv")))

  # update concluded states
  if (!is.null(state)) {
    concluded_states <- rbind(concluded_states, data.frame(state = state, stringsAsFactors = FALSE))
    write_csv(concluded_states, concluded_states_path)
  } else if (!is.null(city)) {
    concluded_cities <- rbind(concluded_cities, data.frame(city = city, stringsAsFactors = FALSE))
    write_csv(concluded_cities, concluded_cities_path)
  }
  return(list(
    best_order = best_order,
    best_formula = best_formula,
    best_pred = best_pred,
    final_metrics = final_metrics,
    all_metrics = metrics
  ))
}


states <- c("AC", "AL", "AM", "AP", "BA", "CE", "DF", "GO", "MA",
            "MG", "MS", "MT", "PA", "PB", "PE", "PI", "PR", "RJ",
            "RN", "RO", "RR", "RS", "SC", "SE", "SP", "TO")


#-----------------------------------------------------------------------------------
# Run model selection for each state and save results
#-----------------------------------------------------------------------------------

results_state_dengue <- lapply(states, function(st) {
  concluded_states <- tryCatch({
    read_csv("sarimax/results/concluded_states_dengue.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(state = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    state = st,
    concluded_states = concluded_states
  )
}
)

# get summary
best_wis_df_dengue <- lapply(states, function(st) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_dengue_", st, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(state = st, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_dengue, "sarimax/results/metrics/best_wis_dengue_all_states.csv")

results_state_chikungunya <- lapply(states, function(st) {
  concluded_states <- tryCatch({
    read_csv("sarimax/results/concluded_states_chikungunya.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(state = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    state = st,
    concluded_states = concluded_states,
    disease = "chikungunya",
    concluded_states_path = "sarimax/results/concluded_states_chikungunya.csv"
  )
}
)

best_wis_df_chikungunya <- lapply(states, function(st) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_chikungunya", st, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(state = st, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_chikungunya, "sarimax/results/metrics/best_wis_chikungunya_all_states.csv")

#-----------------------------------------------------------------------------------
# Run model selection for each city and save results
#-----------------------------------------------------------------------------------

cities_dengue <- c(
  2931350,
  2933307,
  2302503,
  3119401,
  3549805,
  3541406,
  1200401,
  1200203,
  1716109,
  4113700,
  4103701,
  4104808,
  5201405,
  5102637,
  5215231
)

cities_chikungunya <- c(
  2211001,
  2931350,
  3143302,
  3119401,
  1721000,
  1716109,
  4104808,
  4219507,
  5103403,
  5102637
)

results_city_dengue <- lapply(cities_dengue, function(city) {
  concluded_cities <- tryCatch({
    read_csv("sarimax/results/concluded_cities_dengue.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(city = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    city = city,
    concluded_cities = concluded_cities
  )
}
)

best_wis_df_dengue_cities <- lapply(cities_dengue, function(city) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_dengue", city, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(city = city, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_dengue_cities, "sarimax/results/metrics/best_wis_dengue_all_cities.csv")

results_city_chikungunya <- lapply(cities_chikungunya, function(city) {
  concluded_cities <- tryCatch({
    read_csv("sarimax/results/concluded_cities_chikungunya.csv", show_col_types = FALSE)
  }, error = function(e) {
    data.frame(city = character(0), stringsAsFactors = FALSE)
  })
  run_model_selection(
    city = city,
    concluded_cities = concluded_cities,
    disease = "chikungunya"
  )
}
)

best_wis_df_chikungunya_cities <- lapply(cities_chikungunya, function(city) {
  file_name <- paste0("sarimax/results/metrics/metrics_all_formulas_chikungunya", city, ".csv")
  data <- read_csv(file_name, show_col_types = FALSE)
  data[1, ] |> mutate(city = city, .before = 1)
}) |> bind_rows()
write_csv(best_wis_df_chikungunya_cities, "sarimax/results/metrics/best_wis_chikungunya_all_cities.csv")