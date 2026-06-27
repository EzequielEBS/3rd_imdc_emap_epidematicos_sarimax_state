library(forecast)
library(dplyr)
library(tidyverse)
library(tidyr)
library(purrr)
library(lubridate)
library(scoringutils)
library(parallel)
library(pbapply)

#' Fit a SARIMAX model and return forecasts with prediction intervals
#'
#' @param data        A data frame containing all columns referenced in formula,
#'                    plus `epiweek`, `cases`, and train/target indicator columns.
#' @param formula     A formula with `cases` as LHS and covariates as RHS,
#'                    e.g. ~ temp_med_mean + precip_med_mean + enso
#' @param train_id    Character string: one of "train_1","train_2","train_3","train_4"
#' @param quantiles   Numeric vector of coverage levels, e.g. c(0.50, 0.80, 0.95).
#'                    Each value q produces lower_{q*100} and upper_{q*100} columns.
#' @param order       ARIMA (p,d,q) order. Default c(1,1,1).
#' @param seasonal    Seasonal ARIMA list, e.g. list(order=c(1,0,1), period=52).
#'                    Default list(order=c(1,0,1), period=52) for weekly epiweek data.
#'
#' @return A tibble with columns: date, pred, lower_*, upper_*
fit_sarimax <- function(data,
                        formula,
                        train_id,
                        method = "CSS-ML",
                        lambda = NULL,
                        optim.control = list(maxit = 500),   # more iterations
                        optim.method  = "BFGS",
                        levels = c(50, 80, 90, 95),
                        order      = c(1, 1, 1),
                        seasonal   = list(order = c(1, 0, 1), period = 52),
                        bootstrap  = TRUE,   # simulate forecast paths instead of
                        npaths     = 1000) { # assuming Gaussian normal-theory intervals

  stopifnot(is.data.frame(data))
  stopifnot(train_id %in% paste0("train_", 1:4))

  target_id <- sub("train_", "target_", train_id)
  stopifnot(target_id %in% names(data))

  # ── 1. Split train / target rows ─────────────────────────────────────────
  train_rows  <- data[data[[train_id]]  == 1, ]
  target_rows <- data[data[[target_id]] == 1, ]

  if (nrow(train_rows)  == 0) stop("No training rows found for ", train_id)
  if (nrow(target_rows) == 0) stop("No target rows found for ",  target_id)

  # ── 2. Log-transform response (log1p to handle zero counts) ──────────────
  if (is.null(lambda)) {
    y <- log1p(train_rows$cases)
  } else {
    y <- train_rows$cases
  }
  

  # ── 3. Build & standardize regressor matrices ─────────────────────────────
  if (is.null(formula)) {
    rhs_terms <- character(0)
  } else if (is.character(formula)) {
    if (length(formula) == 0) {
      rhs_terms <- character(0)
    } else {
      rhs <- paste(formula, collapse = "+")
      formula_obj <- stats::as.formula(paste("~", rhs))
      rhs_terms <- attr(stats::terms(formula_obj), "term.labels")
    }
  } else {
    rhs_terms <- attr(stats::terms(formula), "term.labels")
  }

  if (length(rhs_terms) == 0) {
    xreg_train  <- NULL
    xreg_future <- NULL
  } else {
    raw_train  <- as.matrix(train_rows[,  rhs_terms, drop = FALSE])
    raw_future <- as.matrix(target_rows[, rhs_terms, drop = FALSE])

    # Compute mean and sd from training data only (no data leakage)
    col_means <- colMeans(raw_train, na.rm = TRUE)
    col_sds   <- apply(raw_train, 2, sd, na.rm = TRUE)

    # Avoid division by zero for constant columns
    col_sds[col_sds == 0] <- 1

    xreg_train  <- scale(raw_train,  center = col_means, scale = col_sds)
    xreg_future <- scale(raw_future, center = col_means, scale = col_sds)
  }

  # ── 4. Fit SARIMAX — capture warnings ─────────────────────────────────────
  y_ts <- ts(y, frequency = seasonal$period)
  fit_warnings <- character(0)

  fit <- withCallingHandlers(
    Arima(y_ts,
          order    = order,
          seasonal = seasonal,
          xreg     = xreg_train,
          method   = method,
          lambda   = lambda,
          optim.control = optim.control,
          optim.method  = optim.method
        ),
    error = function(e) stop("Arima() failed: ", conditionMessage(e)),
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # ── 5. Detect NaN standard errors ─────────────────────────────────────────
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

  # ── 6. Forecast & back-transform ──────────────────────────────────────────
  h      <- nrow(target_rows)
  levels <- sort(unique(levels))

  fc_warnings <- character(0)

  fc <- withCallingHandlers(
    forecast(fit, h = h, xreg = xreg_future, level = levels,
             bootstrap = bootstrap, npaths = npaths),
    warning = function(w) {
      fc_warnings <<- c(fc_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  # Back-transform from log1p scale: expm1(x) = exp(x) - 1
  bt <- expm1  # shorthand

  # ── 7. Access date ─────────────────────────────────────────────
  dates <- data$date[data[[target_id]] == 1]

  # ── 8. Assemble output tibble ─────────────────────────────────────────────
  if (is.null(lambda)) {
    out <- tibble(
      date = dates,
      pred = bt(as.numeric(fc$mean))
    )
  } else {
    out <- tibble(
      date = dates,
      pred = as.numeric(fc$mean)
    )
  }
  

  for (lv in levels) {
    lv_char <- paste0(lv, "%")
    if (is.null(lambda)) {
      out[[paste0("lower_", lv)]] <- bt(as.numeric(fc$lower[, lv_char]))
      out[[paste0("upper_", lv)]] <- bt(as.numeric(fc$upper[, lv_char]))
    } else {
      out[[paste0("lower_", lv)]] <- as.numeric(fc$lower[, lv_char])
      out[[paste0("upper_", lv)]] <- as.numeric(fc$upper[, lv_char])
    }
  }

  # Clip to zero (back-transform should already ensure this, but be safe)
  out <- out |>
    mutate(across(c(pred, starts_with("lower_"), starts_with("upper_")),
                  \(x) pmax(x, 0)))
  
  # ── 9. Attach warnings as attribute ───────────────────────────────────────
  all_warnings <- c(fit_warnings, fc_warnings)

  attr(out, "warnings") <- if (length(all_warnings) > 0) all_warnings else NULL
  attr(out, "fit")      <- fit

  out
}


#' Build all covariate combinations where no pair exceeds a correlation threshold
#'
#' @param data        A data frame containing the candidate covariate columns.
#' @param covariates  Character vector of candidate covariate names.
#' @param threshold   Absolute Pearson correlation threshold (default 0.7).
#'                    Any pair with |r| >= threshold is considered collinear.
#' @param min_size    Minimum number of covariates per combination (default 1).
#' @param max_size    Maximum number of covariates per combination (default Inf).
#'
#' @return A list of character vectors, each a valid (uncorrelated) covariate set.
build_covariate_combinations <- function(data,
                                         covariates,
                                         threshold = 0.7,
                                         min_size  = 1,
                                         max_size  = Inf) {

  stopifnot(is.data.frame(data))
  stopifnot(all(covariates %in% names(data)))
  stopifnot(threshold > 0 & threshold <= 1)

  # ── 1. Correlation matrix ─────────────────────────────────────────────────
  cor_mat  <- cor(data[, covariates, drop = FALSE], use = "pairwise.complete.obs")
  collinear <- abs(cor_mat) >= threshold
  diag(collinear) <- FALSE

  # ── 2. Pre-compute total combinations for progress tracking ───────────────
  max_size  <- min(max_size, length(covariates))
  sizes     <- seq(min_size, max_size)
  combos_per_size <- sapply(sizes, \(k) choose(length(covariates), k))
  total     <- sum(combos_per_size)

  # ── 3. Progress bar setup ─────────────────────────────────────────────────
  start_time  <- proc.time()["elapsed"]
  checked     <- 0L
  valid       <- vector("list", 0)

  fmt_time <- function(secs) {
    if (secs < 60)  return(sprintf("%.0fs",       secs))
    if (secs < 3600) return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
    sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
  }

  draw_progress <- function(checked, total, n_valid, start_time) {
    elapsed  <- proc.time()["elapsed"] - start_time
    pct      <- checked / total
    rate     <- if (elapsed > 0) checked / elapsed else 0
    eta      <- if (rate > 0) (total - checked) / rate else NA

    bar_width <- 30L
    filled    <- round(pct * bar_width)
    bar       <- paste0(
      "[", strrep("=", filled),
      if (filled < bar_width) ">" else "",
      strrep(" ", max(0L, bar_width - filled - 1L)),
      "]"
    )

    cat(sprintf(
      "\r%s %5.1f%%  checked: %d/%d  valid: %d  elapsed: %s  ETA: %s     ",
      bar, pct * 100, checked, total, n_valid,
      fmt_time(elapsed),
      if (is.na(eta)) "---" else fmt_time(eta)
    ))
    flush.console()
  }

  # ── 4. Iterate combinations with live progress ────────────────────────────
  for (k in sizes) {
    combos_k <- combn(covariates, k, simplify = FALSE)

    for (combo in combos_k) {
      sub_collinear <- collinear[combo, combo, drop = FALSE]
      if (!any(sub_collinear)) valid <- c(valid, list(combo))

      checked <- checked + 1L
      if (checked %% 500L == 0L || checked == total)
        draw_progress(checked, total, length(valid), start_time)
    }
  }

  # Final completed bar
  draw_progress(total, total, length(valid), start_time)
  cat("\n")

  elapsed_total <- proc.time()["elapsed"] - start_time
  message(sprintf(
    "Done. %d valid combination(s) from %d candidates | threshold = %.2f | time: %s",
    length(valid), length(covariates), threshold, fmt_time(elapsed_total)
  ))

  valid
}

get_candidates <- function(data,
                           groups       = c("contemporaneous", "lagged", "rolling"),
                           vars         = c("temp_med_mean", "precip_med_mean",
                                            "rel_humid_med_mean", "thermal_range_mean",
                                            "rainy_days_mean", "enso", "iod", "pdo"),
                           lag_select   = NULL,
                           roll_select  = NULL
                          ) # single rolling window to keep
{
  all_cols <- names(data)
  selected <- character(0)

  if ("contemporaneous" %in% groups)
    selected <- c(selected, intersect(vars, all_cols))

  if ("lagged" %in% groups) {
    if (!is.null(lag_select)) {
      lag_pat  <- sprintf("_lag%d$", lag_select)
      selected <- c(selected, grep(lag_pat, all_cols, value = TRUE))
    } else {
      selected <- c(selected, grep("_lag\\d+$",    all_cols, value = TRUE))
    }
  }

  if ("rolling" %in% groups) {
    if (!is.null(roll_select)) {
      roll_pat <- sprintf("_mean_%s$", roll_select)
      selected <- c(selected, grep(roll_pat, all_cols, value = TRUE))
    } else {
      selected <- c(selected, grep("_mean_\\d+mo$", all_cols, value = TRUE))
    }
  }

  unique(selected)
}

filter_low_variance <- function(data, covariates, threshold = 0.01) {
  vars_sd <- sapply(covariates, \(v) sd(data[[v]], na.rm = TRUE))
  keep    <- names(vars_sd[vars_sd > threshold])
  dropped <- setdiff(covariates, keep)
  if (length(dropped) > 0)
    message("Dropped low-variance columns: ", paste(dropped, collapse = ", "))
  keep
}

# Scoringutils WIS helper ────────────────────────────────────────────
# scoringutils expects a long data frame with quantile forecasts
compute_metrics <- function(
  pred_df,
  actual,
  levels = c(50, 80, 90, 95)
) {
  levels     <- sort(levels)
  lower_cols <- paste0("lower_", levels)
  upper_cols <- paste0("upper_", levels)
  alphas     <- 1 - levels / 100

  n <- nrow(pred_df)

  make_quantile_rows <- function(predicted, quantile_level) {
    tibble(
      time           = seq_len(n),
      observed       = actual,
      predicted      = predicted,
      quantile_level = quantile_level,
      model          = "sarimax"
    )
  }

  su <- bind_rows(
    # Lower bounds: alpha/2  (e.g. level=95 -> 0.025, level=50 -> 0.25)
    lapply(seq_along(levels), function(i)
      make_quantile_rows(pred_df[[lower_cols[i]]], alphas[i] / 2)),
    # Upper bounds: 1 - alpha/2
    lapply(seq_along(levels), function(i)
      make_quantile_rows(pred_df[[upper_cols[i]]], 1 - alphas[i] / 2)),
    # Median — always added; alpha/2 never equals 0.5 so no duplication risk
    make_quantile_rows(pred_df$pred, 0.5)
  ) |>
    arrange(time, quantile_level) |>
    as_forecast_quantile(
      observed       = "observed",
      predicted      = "predicted",
      quantile_level = "quantile_level",
      model          = "model"
    )

  sc <- withCallingHandlers(
    score(su) |> summarise_scores(by = "model"),
    warning = function(w) {
      if (grepl("interval_coverage", conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )

  coverage <- setNames(
    sapply(seq_along(levels), function(i) {
      mean(actual >= pred_df[[lower_cols[i]]] &
           actual <= pred_df[[upper_cols[i]]])
    }),
    paste0("coverage_", levels)
  )

  y    <- actual
  yhat <- pred_df$pred

  c(
    list(
      wis  = sc$wis,
      mae  = mean(abs(y - yhat)),
      mse  = mean((y - yhat)^2),
      rmse = sqrt(mean((y - yhat)^2)),
      mape = mean(abs((y - yhat) / pmax(y, 1))) * 100
    ),
    as.list(coverage)
  )
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Grid search over ARIMA orders and covariate formulas with CV evaluation
#'
#' @param data       Data frame with all covariates and train/target indicators.
#' @param formulas   List of formulas (from build_covariate_combinations output).
#' @param quantiles  Numeric vector of coverage levels, e.g. c(0.80, 0.95).
#' @param seasonal   Seasonal period list passed to fit_sarimax.
#' @param max_order  Named integer vector with upper bounds for each order
#'                   component: c(p=2, d=2, q=2, P=1, D=1, Q=1).
#' @param n_cores    Number of parallel workers for the order loop.
#'                   Defaults to all available cores minus 1.
#'
#' @return A list with two tibbles:
#'   $predictions : one row per formula x order x train_id x date
#'   $metrics     : one row per formula x order x train_id, with all scores
run_grid_search <- function(data,
                            formulas,
                            levels = c(50, 80, 90, 95),
                            method = "CSS",
                            lambda = NULL,
                            optim.control = list(maxit = 500),   # more iterations
                            optim.method  = "BFGS",
                            seasonal  = list(order = c(1, 0, 1), period = 52),
                            max_order = c(p = 2, d = 2, q = 2,
                                          P = 1, D = 1, Q = 1),
                            fixed_order = F,
                            # FALSE searches d (and D) over the grid via CV WIS instead of
                            # locking it to a single KPSS/OCSB-recommended value: on short,
                            # irregular epidemic series the auto test can pick a worse d than
                            # what actually scores best out-of-sample.
                            fixed_stat_par = F,
                            bootstrap = TRUE,
                            npaths    = 1000,
                            n_cores   = max(1L, detectCores() - 1L)) {

  # ── 1. Build order grid ───────────────────────────────────────────────────
  if (fixed_order) {
    # Single row grid from fixed order — no filtering needed
    order_grid <- data.frame(
      p = fixed_order["p"], d = fixed_order["d"], q = fixed_order["q"],
      P = fixed_order["P"], D = fixed_order["D"], Q = fixed_order["Q"],
      row.names = NULL
    )
  } else {
    required_names <- c("p", "d", "q", "P", "D", "Q")
    stopifnot(all(required_names %in% names(max_order)))
    stopifnot(all(max_order >= 0))

    if (fixed_stat_par) {
      diff_par <- determine_d(log1p(data$cases))
    }

    order_grid <- expand.grid(
      p = seq(0, max_order["p"]),
      d = if (fixed_stat_par) diff_par$d else seq(0, max_order["d"]),
      q = seq(0, max_order["q"]),
      P = seq(0, max_order["P"]),
      D = if (fixed_stat_par) diff_par$D else seq(0, max_order["D"]),
      Q = seq(0, max_order["Q"])
    ) |> filter(!(p == 0 & q == 0), !(P == 0 & Q == 0))
  }

  train_ids  <- paste0("train_", 1:4)
  n_orders   <- nrow(order_grid)
  n_formulas <- length(formulas)

  # ── 3. Build flat grid: every (formula_idx, order_idx) pair ───────────────
  job_grid <- expand.grid(
    fi = seq_len(n_formulas),
    oi = seq_len(n_orders)
  )
  total_jobs <- nrow(job_grid)

  message(sprintf(
    "Grid: %d formula(s) x %d order(s) x %d splits = %d total fits.",
    n_formulas, n_orders, length(train_ids),
    total_jobs * length(train_ids)
  ))

  # ── 4. Precompute formula metadata (outside workers) ──────────────────────
  formula_ids <- vapply(formulas, FUN.VALUE = character(1), function(f) {
    if (is.character(f)) {
      if (length(f) > 1) {
        paste(trimws(f), collapse = "+")
      } else { # length == 1
        single <- f[[1]]
        if (grepl("~", single)) {
          rhs <- sub(".*~", "", single)
          terms <- strsplit(rhs, "\\+")[[1]]
          paste(trimws(terms), collapse = "+")
        } else {
          paste(trimws(single), collapse = "+")
        }
      }
    } else if (inherits(f, "formula")) {
      paste(attr(stats::terms(f), "term.labels"), collapse = "+")
    } else {
      stop("Each element of 'formulas' must be a formula or character string/vector")
    }
  })
  # ── 5. Parallel backend ───────────────────────────────────────────────────
  cl <- makeCluster(n_cores)
  on.exit(stopCluster(cl), add = TRUE)          # no progress_file to unlink

  clusterExport(cl, varlist = c(
    "fit_sarimax", "compute_metrics",
    "data", "formulas", "formula_ids",
    "order_grid", "job_grid",
    "seasonal", "levels", "train_ids",
    "levels", "method", "lambda", "optim.control", "optim.method",
    "bootstrap", "npaths"
  ), envir = environment())

  clusterEvalQ(cl, {
    library(forecast)
    library(dplyr)
    library(tidyr)
    library(scoringutils)
    library(pbapply)
  })

  message(sprintf("Running %d jobs on %d core(s)...", total_jobs, n_cores))
  start_time <- proc.time()["elapsed"]

  # ── 6. Single flat parallel loop over all (formula, order) pairs ──────────
  op <- pbapply::pboptions(type="timer") 
  pbapply::pboptions(op)
  all_results <- pbapply::pblapply(
    seq_len(total_jobs),
    function(job) {

      fi <- job_grid$fi[job]
      oi <- job_grid$oi[job]

      formula    <- formulas[[fi]]
      formula_id <- formula_ids[[fi]]
      og         <- order_grid[oi, ]
      ord        <- c(og$p, og$d, og$q)
      sea        <- list(order = c(og$P, og$D, og$Q), period = seasonal$period)
      order_str  <- sprintf("(%d,%d,%d)(%d,%d,%d)",
                            og$p, og$d, og$q, og$P, og$D, og$Q)

      split_results <- lapply(train_ids, function(train_id) {
        target_id <- sub("train_", "target_", train_id)
        actual    <- data$cases[data[[target_id]] == 1]

        preds <- tryCatch(
          fit_sarimax(data,
                      formula   = formula,
                      train_id  = train_id,
                      levels    = levels,
                      method    = method,
                      order     = ord,
                      seasonal  = sea,
                      lambda    = lambda,
                      optim.control = optim.control,
                      optim.method  = optim.method,
                      bootstrap = bootstrap,
                      npaths    = npaths
                    ),
          error = function(e) NULL
        )

        failed  <- is.null(preds) || nrow(preds) != length(actual)
        metrics <- if (!failed) {
          tryCatch(
            as_tibble(compute_metrics(preds, actual)),
            error = function(e) NULL
          )
        } else NULL

        list(
          predictions = if (!failed) {
            preds |> mutate(formula_id = formula_id,
                            order      = order_str,
                            train_id   = train_id,
                            failed     = FALSE)
          } else {
            tibble(formula_id = formula_id, order = order_str,
                  train_id = train_id, failed = TRUE)
          },
          metrics = if (!is.null(metrics)) {
            metrics |> mutate(formula_id = formula_id,
                              order      = order_str,
                              train_id   = train_id,
                              failed     = FALSE)
          } else {
            tibble(formula_id = formula_id, order = order_str,
                  train_id = train_id, failed = TRUE)
          }
        )
      })

      split_results
    },
    cl = cl   # PSOCK cluster — pblapply handles distribution automatically
  )

  # ── 7. Progress summary ───────────────────────────────────────────────────
  elapsed_total <- proc.time()["elapsed"] - start_time
  fmt_time <- function(secs) {
    if (secs < 60)    return(sprintf("%.0fs", secs))
    if (secs < 3600)  return(sprintf("%.0fm %.0fs", secs %/% 60, secs %% 60))
    sprintf("%.0fh %.0fm", secs %/% 3600, (secs %% 3600) %/% 60)
  }

  # ── 8. Flatten and bind ───────────────────────────────────────────────────
  flat <- unlist(all_results, recursive = FALSE)

  predictions <- bind_rows(lapply(flat, `[[`, "predictions"))
  metrics     <- bind_rows(lapply(flat, `[[`, "metrics"))

  n_failed <- sum(metrics$failed, na.rm = TRUE)
  message(sprintf("Done. Total time: %s. Failed fits: %d / %d.",
                  fmt_time(elapsed_total), n_failed,
                  total_jobs * length(train_ids)))

  list(predictions = predictions, metrics = metrics)
}

select_best_per_variable <- function(data, covariates, response = "cases") {
  # Extract base variable name by stripping lag/rolling suffixes
  base_var <- gsub("_lag\\d+$|_mean_\\d+mo$", "", covariates)

  y <- log1p(data[[response]])

  # For each base variable, keep the candidate most correlated with response
  result <- tapply(covariates, base_var, function(candidates) {
    cors <- sapply(candidates, function(v)
      abs(cor(data[[v]], y, use = "pairwise.complete.obs")))
    candidates[which.max(cors)]
  })

  unname(unlist(result))
}

filter_by_correlation <- function(data, covariates,
                                  response  = "cases",
                                  min_cor   = 0.1) {
  y    <- log1p(data[[response]])
  cors <- sapply(covariates, function(v)
    abs(cor(data[[v]], y, use = "pairwise.complete.obs")))

  keep    <- names(cors[cors >= min_cor])
  dropped <- setdiff(covariates, keep)

  if (length(dropped) > 0)
    message(sprintf("Dropped %d weak predictor(s): %s",
                    length(dropped), paste(dropped, collapse = ", ")))
  keep
}

filter_redundant_indices <- function(data, covariates,
                                     indices   = c("enso", "iod", "pdo"),
                                     threshold = 0.3) {
  # Keep an index only if its partial correlation with response
  # (after removing linear effect of local weather) exceeds threshold
  local_weather <- setdiff(covariates, indices)
  present_indices <- intersect(indices, covariates)

  if (length(local_weather) == 0 || length(present_indices) == 0)
    return(covariates)

  y       <- log1p(data$cases)
  X_local <- as.matrix(data[, local_weather, drop = FALSE])

  # Residualise response on local weather
  y_resid <- residuals(lm(y ~ X_local))

  keep_indices <- Filter(function(idx) {
    abs(cor(data[[idx]], y_resid, use = "pairwise.complete.obs")) >= threshold
  }, present_indices)

  dropped <- setdiff(present_indices, keep_indices)
  if (length(dropped) > 0)
    message(sprintf("Dropped redundant index/indices: %s",
                    paste(dropped, collapse = ", ")))

  c(local_weather, keep_indices)
}

parse_order <- function(order_str) {
  nums <- as.integer(regmatches(order_str, gregexpr("[0-9]", order_str))[[1]])
  list(order = nums[1:3], seasonal_order = nums[4:6])
}

determine_d <- function(y, max_d = 2) {
  d  <- ndiffs(y,  test = "kpss", max.d = max_d)
  D  <- nsdiffs(y, test = "ocsb", m = 52)
  message(sprintf("Recommended: d = %d, D = %d", d, D))
  list(d = d, D = D)
}

pca_all <- function(data,
                    candidates,
                    train_cols    = paste0("train_", 1:4),
                    var_threshold = 0.90) {

  mat    <- as.matrix(data[, candidates, drop = FALSE])
  n_rows <- nrow(data)

  # ── 1. Pre-compute n_comp consistently across all folds ───────────────────
  n_comps <- sapply(seq_along(train_cols), function(i) {
    train_rows <- data[[train_cols[i]]] == 1
    pca_fit    <- prcomp(mat[train_rows, , drop = FALSE],
                         center = TRUE, scale. = TRUE)
    var_exp    <- cumsum(pca_fit$sdev^2) / sum(pca_fit$sdev^2)
    max(1L, which(var_exp >= var_threshold)[1])
  })

  # Use minimum across folds so all folds produce the same number of components
  n_comp     <- min(n_comps)
  comp_names <- paste0("PC", seq_len(n_comp))
  scores_out <- matrix(NA_real_, nrow = n_rows, ncol = n_comp,
                       dimnames = list(NULL, comp_names))

  message(sprintf("Using %d components (min across folds: %s).",
                  n_comp, paste(n_comps, collapse = ", ")))

  # ── 2. Fit per fold and write scores ──────────────────────────────────────
  pca_fit_last <- NULL

  for (i in seq_along(train_cols)) {
    train_col  <- train_cols[i]
    target_col <- sub("train_", "target_", train_col)
    train_rows <- data[[train_col]] == 1
    fold_rows  <- data[[train_col]] == 1 | data[[target_col]] == 1

    message(sprintf("Fold %d: fitting PCA on %d training rows.", i, sum(train_rows)))

    pca_fit <- prcomp(mat[train_rows, , drop = FALSE],
                      center = TRUE, scale. = TRUE)

    mat_scaled <- scale(mat[fold_rows, , drop = FALSE],
                        center = pca_fit$center,
                        scale  = pca_fit$scale)

    # Always take exactly n_comp components
    scores <- mat_scaled %*% pca_fit$rotation[, seq_len(n_comp), drop = FALSE]

    scores_out[fold_rows, ] <- scores
    pca_fit_last <- pca_fit
  }

  # ── 3. Variance table (based on consistent n_comp) ────────────────────────
  var_table <- tibble(
    component = comp_names,
    var_exp   = pca_fit_last$sdev[seq_len(n_comp)]^2 / sum(pca_fit_last$sdev^2),
    cum_var   = cumsum(pca_fit_last$sdev[seq_len(n_comp)]^2 /
                       sum(pca_fit_last$sdev^2))
  )
  print(var_table)

  list(
    data      = bind_cols(data, as_tibble(scores_out)),
    variables = comp_names,
    var_table = var_table
  )
}