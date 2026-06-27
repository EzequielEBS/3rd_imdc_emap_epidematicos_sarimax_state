source("sarimax/src/utils.r")

# ── Helper: does a given year have 53 epiweeks? ───────────────────────────
has_53_weeks <- function(year) {
  as.integer(format(as.Date(paste0(year, "-12-28")), "%V")) == 53
}

# ── Helper: enumerate all epiweeks between two YYYYWW integers ────────────
enumerate_epiweeks <- function(start, end) {
  start_year <- start %/% 100L
  start_week <- start  %% 100L
  end_year   <- end   %/% 100L
  end_week   <- end    %% 100L

  epiweeks <- integer(0)
  yr <- start_year
  wk <- start_week

  repeat {
    epiweeks <- c(epiweeks, yr * 100L + wk)
    if (yr == end_year && wk == end_week) break
    n_weeks <- if (has_53_weeks(yr)) 53L else 52L
    if (wk < n_weeks) {
      wk <- wk + 1L
    } else {
      yr <- yr + 1L
      wk <- 1L
    }
  }
  epiweeks
}

# ── Helper: given YYYYWW, return the previous year's equivalent epiweek ───
# If the previous year does not have that week (week 53 in a 52-week year),
# shift forward by one week.
prev_year_epiweek <- function(yw) {
  yr <- yw %/% 100L
  wk <- yw  %% 100L
  prev_yr      <- yr - 1L
  prev_n_weeks <- if (has_53_weeks(prev_yr)) 53L else 52L

  if (wk <= prev_n_weeks) {
    list(epiweek = prev_yr * 100L + wk, adjusted = FALSE)
  } else {
    # week 53 doesn't exist in prev year → shift to week 1 of current year
    list(epiweek = yr * 100L + 1L, adjusted = TRUE)
  }
}

# ── Helper: get the Monday date of a given YYYYWW epiweek (ISO 8601) ──────
epiweek_to_date <- function(yw) {
  yr <- yw %/% 100L
  wk <- yw  %% 100L
  # ISO week date: Monday of week WW of year YYYY
  # Jan 4th is always in week 1
  jan4     <- as.Date(paste0(yr, "-01-04"))
  monday_w1 <- jan4 - as.integer(format(jan4, "%u")) + 1L
  monday_w1 + (wk - 1L) * 7L
}

fit_sarimax_epiweek <- function(
                        data,
                        formula,
                        train_start,
                        train_end,
                        forecast_start,
                        forecast_end,
                        method        = "CSS-ML",
                        lambda        = NULL,
                        optim.control = list(maxit = 500),
                        optim.method  = "BFGS",
                        levels        = c(50, 80, 90, 95),
                        order         = c(1, 1, 1),
                        seasonal      = list(order = c(1, 0, 1), period = 52),
                        bootstrap     = TRUE,
                        npaths        = 1000) {
 
  stopifnot(is.data.frame(data))
  stopifnot("epiweek" %in% names(data))
  stopifnot("cases"   %in% names(data))
 
  # ── 1. Training rows ───────────────────────────────────────────────────────
  train_rows <- data[data$epiweek >= train_start & data$epiweek <= train_end, ]
  if (nrow(train_rows) == 0) {
    stop("No training rows found for epiweeks ", train_start, "–", train_end)
  }
 
  # ── 2. Enumerate forecast epiweeks & build forecast dates/prev-year info ──
  fc_epiweeks <- enumerate_epiweeks(forecast_start, forecast_end)
  n_fc        <- length(fc_epiweeks)
 
  prev_epiweeks  <- integer(n_fc)
  forecast_dates <- as.Date(rep(NA, n_fc))
  adjusted_idx   <- logical(n_fc)
 
  for (i in seq_len(n_fc)) {
    res              <- prev_year_epiweek(fc_epiweeks[i])
    prev_epiweeks[i] <- res$epiweek
    adjusted_idx[i]  <- res$adjusted
    # Forecast date = Monday of the forecast epiweek
    # (Sunday = Monday + 6 if your dates are Sundays)
    forecast_dates[i] <- epiweek_to_date(fc_epiweeks[i]) - 1L  # Sunday
  }
 
  if (any(adjusted_idx)) {
    adj_info <- paste0(
      fc_epiweeks[adjusted_idx], " (prev: ", prev_epiweeks[adjusted_idx], ")",
      collapse = "; "
    )
    warning(
      "53-week year mismatch: the following forecast epiweeks have no direct ",
      "previous-year equivalent and were shifted forward by one week:\n  ",
      adj_info,
      call. = FALSE
    )
  }
 
  # ── 3. Log-transform response ──────────────────────────────────────────────
  y <- if (is.null(lambda)) log1p(train_rows$cases) else train_rows$cases
 
  # ── 4. Build & standardize regressor matrices ──────────────────────────────
  rhs_terms <- {
    if (is.null(formula)) {
      character(0)
    } else if (is.character(formula)) {
      if (length(formula) == 0) character(0)
      else {
        fo <- stats::as.formula(paste("~", paste(formula, collapse = "+")))
        attr(stats::terms(fo), "term.labels")
      }
    } else {
      attr(stats::terms(formula), "term.labels")
    }
  }
 
  if (length(rhs_terms) == 0) {
    xreg_train  <- NULL
    xreg_future <- NULL
  } else {
    raw_train <- as.matrix(train_rows[, rhs_terms, drop = FALSE])
 
    col_means <- colMeans(raw_train, na.rm = TRUE)
    col_sds   <- apply(raw_train, 2, sd, na.rm = TRUE)
    col_sds[col_sds == 0] <- 1
 
    xreg_train <- scale(raw_train, center = col_means, scale = col_sds)
 
    # Look up previous-year rows by epiweek
    prev_rows  <- data[match(prev_epiweeks, data$epiweek), rhs_terms, drop = FALSE]
 
    if (any(is.na(prev_rows))) {
      missing_ew <- prev_epiweeks[apply(is.na(prev_rows), 1, any)]
      stop(
        "Could not find previous-year covariate data for epiweeks: ",
        paste(missing_ew, collapse = ", ")
      )
    }
 
    raw_future  <- as.matrix(prev_rows)
    xreg_future <- scale(raw_future, center = col_means, scale = col_sds)
  }
 
  # ── 5. Fit SARIMAX ─────────────────────────────────────────────────────────
  y_ts         <- ts(y, frequency = seasonal$period)
  fit_warnings <- character(0)
 
  fit <- withCallingHandlers(
    forecast::Arima(
      y_ts,
      order         = order,
      seasonal      = seasonal,
      xreg          = xreg_train,
      method        = method,
      lambda        = lambda,
      optim.control = optim.control,
      optim.method  = optim.method
    ),
    error   = function(e) stop("Arima() failed: ", conditionMessage(e)),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
 
  # ── 6. Detect NaN standard errors ─────────────────────────────────────────
  nan_se_params <- names(which(is.nan(sqrt(diag(fit$var.coef)))))
  if (length(nan_se_params) > 0) {
    fit_warnings <- c(
      fit_warnings,
      paste0(
        "NaN standard errors for: ",
        paste(nan_se_params, collapse = ", "),
        ". Prediction intervals may be unreliable. ",
        "Consider reducing model order or using auto.arima()."
      )
    )
  }
 
  # ── 7. Forecast & back-transform ───────────────────────────────────────────
  levels <- sort(unique(levels))
  fc_warnings <- character(0)
 
  fc <- withCallingHandlers(
    forecast::forecast(
      fit,
      h         = n_fc,
      xreg      = xreg_future,
      level     = levels,
      bootstrap = bootstrap,
      npaths    = npaths
    ),
    warning = function(w) {
      fc_warnings <<- c(fc_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
 
  bt <- if (is.null(lambda)) expm1 else identity
 
  # ── 8. Assemble output tibble ──────────────────────────────────────────────
  out <- tibble::tibble(
    date    = forecast_dates,
    pred    = bt(as.numeric(fc$mean))
  )
 
  for (lv in levels) {
    lv_char <- paste0(lv, "%")
    out[[paste0("lower_", lv)]] <- bt(as.numeric(fc$lower[, lv_char]))
    out[[paste0("upper_", lv)]] <- bt(as.numeric(fc$upper[, lv_char]))
  }
 
  out <- out |>
    dplyr::mutate(dplyr::across(
      c(pred, dplyr::starts_with("lower_"), dplyr::starts_with("upper_")),
      \(x) pmax(x, 0)
    ))
 
  # ── 9. Attach metadata ─────────────────────────────────────────────────────
  all_warnings <- c(fit_warnings, fc_warnings)
  attr(out, "warnings") <- if (length(all_warnings) > 0) all_warnings else NULL
  attr(out, "fit")      <- fit
 
  out
}

best_wis_df_dengue_state <- read_csv("sarimax/results/metrics/best_wis_dengue_all_states.csv", show_col_types = FALSE)

preds_dengue_state <- lapply(best_wis_df_dengue_state$state, function(st) {
  file_name <- paste0("processed_data/dengue/dengue_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)
  train_rows <- d$train_1 == 1

  candidates <- get_candidates(d)
  candidates <- filter_low_variance(d, candidates, threshold = 0.01)
  candidates <- filter_by_correlation(
    d[train_rows, ],
    candidates,
    min_cor = 0.1
  )
    

  pca_result <- pca_all(
    data = d,
    candidates = candidates[!grepl("enso|iod|pdo", candidates)],
    var_threshold = 0.9
  )
  d <- pca_result$data

  best_order <- best_wis_df_dengue_state |> filter(state == st) |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_wis_df_dengue_state |> filter(state == st) |> pull(formula_id)
  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)
  pred_target <- lapply(4, function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    if (train_id == "train_4") {
      fit <- fit_sarimax_epiweek(
        data = d,
        formula = best_formula,
        train_start = 201001,
        train_end = 202525,
        forecast_start = 202541,
        forecast_end = 202640,
        order = c(ord$order[1], ord$order[2], ord$order[3]),
        seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52)
      )
    } else {
      fit <- fit_sarimax(
        data = d,
        formula = best_formula,
        order = c(ord$order[1], ord$order[2], ord$order[3]),
        seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
        train_id = train_id
      )
    }
    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_dengue_", st, "_", target_id, ".csv")))
  })
  names(pred_target) <- target_ids
  bind_rows(pred_target, .id = "target_id") |>
    mutate(state = st)
})

best_wis_df_chikungunya_state <- read_csv("sarimax/results/metrics/best_wis_chikungunya_all_states.csv", show_col_types = FALSE)

preds_chikungunya_state <- lapply(best_wis_df_chikungunya_state$state, function(st) {
  file_name <- paste0("processed_data/chikungunya/chikungunya_", st, "_agg.csv.gz")
  d <- read_csv(file_name, show_col_types = FALSE)
  train_rows <- d$train_1 == 1 | d$train_2 == 1 | d$train_3 == 1 | d$train_4 == 1

  candidates <- get_candidates(d)
  candidates <- filter_low_variance(d, candidates, threshold = 0.01)
  candidates <- filter_by_correlation(
    d[train_rows, ],
    candidates,
    min_cor = 0.1
  )

  pca_result <- pca_all(
    data = d,
    candidates = candidates[!grepl("enso|iod|pdo", candidates)],
    var_threshold = 0.9
  )
  d <- pca_result$data

  best_order <- best_wis_df_chikungunya_state |> filter(state == st) |> pull(order)
  ord <- parse_order(best_order)
  best_formula <- best_wis_df_chikungunya_state |> filter(state == st) |> pull(formula_id)
  train_ids  <- paste0("train_", 1:4)
  target_ids <- paste0("target_", 1:4)
  pred_target <- lapply(seq_along(train_ids), function(i) {
    train_id <- train_ids[i]
    target_id <- target_ids[i]
    fit <- fit_sarimax(
      data = d,
      formula = best_formula,
      order = c(ord$order[1], ord$order[2], ord$order[3]),
      seasonal = list(order = c(ord$seasonal_order[1], ord$seasonal_order[2], ord$seasonal_order[3]), period = 52),
      train_id = train_id
    )
    write_csv(fit, file.path("sarimax/results/preds/", paste0("pred_chikungunya_", st, "_", target_id, ".csv")))
  })
  names(pred_target) <- target_ids
  bind_rows(pred_target, .id = "target_id") |>
    mutate(state = st)
})