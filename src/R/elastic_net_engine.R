# Configuration-driven Elastic Net fitting with nested cross-validation.

# Fit one inner cross-validation model while capturing warnings for the audit trail.
capture_cv_glmnet <- function(x, y, fold_id, alpha, penalty_factor, maxit) {
  warning_messages <- character()
  fit <- withCallingHandlers(
    glmnet::cv.glmnet(
      x = x,
      y = y,
      family = "binomial",
      alpha = alpha,
      foldid = fold_id,
      type.measure = "auc",
      standardize = TRUE,
      intercept = TRUE,
      penalty.factor = penalty_factor,
      grouped = TRUE,
      keep = FALSE,
      parallel = FALSE,
      maxit = maxit
    ),
    warning = function(condition) {
      warning_messages <<- c(warning_messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warning_messages))
}

# Tune alpha and lambda entirely within the outer-training partition and predict the held-out fold.
fit_elastic_net_model <- function(
  model_name,
  predictors,
  context,
  parameters
) {
  data <- context$analysis_data
  training_data <- data[context$training_indices, , drop = FALSE]
  test_data <- data[context$test_indices, , drop = FALSE]

  x_training <- encode_predictor_matrix(training_data, predictors, context$binary_mapping)
  x_test <- encode_predictor_matrix(test_data, predictors, context$binary_mapping)
  y_training <- numeric_vector(training_data$hdremit.all, "hdremit.all")
  y_test <- numeric_vector(test_data$hdremit.all, "hdremit.all")

  penalty_factor <- ifelse(predictors %in% context$forced_covariates, 0, 1)
  alpha_grid <- as_numeric_list(parameters$alpha_grid)
  maxit <- as_integer_parameter(parameters, "elastic_net_maxit")
  master_seed <- as_integer_parameter(parameters, "master_seed")
  tolerance <- as_numeric_parameter(parameters, "tie_tolerance")
  epsilon <- as_numeric_parameter(parameters, "probability_epsilon")

  tuning_rows <- list()
  warning_rows <- list()

  for (alpha_index in seq_along(alpha_grid)) {
    alpha_value <- alpha_grid[[alpha_index]]
    set.seed(master_seed + alpha_index)
    cv_result <- capture_cv_glmnet(
      x_training,
      y_training,
      context$inner_fold_id,
      alpha_value,
      penalty_factor,
      maxit
    )
    cv_fit <- cv_result$fit
    finite_indices <- which(is.finite(cv_fit$cvm))
    if (length(finite_indices) == 0L) {
      stop(paste("No finite inner AUC values for", model_name, "at alpha", alpha_value))
    }
    best_auc <- max(cv_fit$cvm[finite_indices])
    tied <- finite_indices[abs(cv_fit$cvm[finite_indices] - best_auc) <= tolerance]
    best_index <- tied[which.max(cv_fit$lambda[tied])]
    tuning_rows[[length(tuning_rows) + 1L]] <- data.frame(
      model = model_name,
      alpha = alpha_value,
      selected_lambda = cv_fit$lambda[[best_index]],
      mean_inner_auc = best_auc,
      inner_auc_standard_error = cv_fit$cvsd[[best_index]],
      lambda_index = best_index,
      warning_count = length(cv_result$warnings),
      stringsAsFactors = FALSE
    )
    if (length(cv_result$warnings) > 0L) {
      for (warning_text in cv_result$warnings) {
        warning_rows[[length(warning_rows) + 1L]] <- data.frame(
          stage = "inner_tuning",
          model = model_name,
          alpha = alpha_value,
          warning = warning_text,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  tuning_table <- do.call(rbind, tuning_rows)
  maximum_auc <- max(tuning_table$mean_inner_auc)
  candidates <- which(abs(tuning_table$mean_inner_auc - maximum_auc) <= tolerance)
  largest_lambda <- max(tuning_table$selected_lambda[candidates])
  candidates <- candidates[
    abs(tuning_table$selected_lambda[candidates] - largest_lambda) <= tolerance
  ]
  selected_row <- candidates[which.min(tuning_table$alpha[candidates])]
  selected_alpha <- tuning_table$alpha[[selected_row]]
  selected_lambda <- tuning_table$selected_lambda[[selected_row]]
  tuning_table$selected <- seq_len(nrow(tuning_table)) == selected_row

  set.seed(master_seed + 5000L)
  final_warnings <- character()
  final_fit <- withCallingHandlers(
    glmnet::glmnet(
      x = x_training,
      y = y_training,
      family = "binomial",
      alpha = selected_alpha,
      lambda = selected_lambda,
      standardize = TRUE,
      intercept = TRUE,
      penalty.factor = penalty_factor,
      maxit = maxit
    ),
    warning = function(condition) {
      final_warnings <<- c(final_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  if (length(final_warnings) > 0L) {
    for (warning_text in unique(final_warnings)) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        stage = "outer_final_fit",
        model = model_name,
        alpha = selected_alpha,
        warning = warning_text,
        stringsAsFactors = FALSE
      )
    }
  }

  predictions <- as.numeric(predict(final_fit, newx = x_test, s = selected_lambda, type = "response"))
  if (length(predictions) != nrow(test_data) || anyNA(predictions) || any(!is.finite(predictions))) {
    stop(paste("Invalid held-out predictions for:", model_name))
  }
  if (any(predictions < 0 | predictions > 1)) stop("Predictions fall outside [0,1].")

  coefficient_matrix <- as.matrix(stats::coef(final_fit, s = selected_lambda))
  coefficient_table <- data.frame(
    model = model_name,
    term = rownames(coefficient_matrix),
    coefficient = as.numeric(coefficient_matrix[, 1]),
    stringsAsFactors = FALSE
  )
  coefficient_table$forced_covariate <- coefficient_table$term %in% context$forced_covariates
  coefficient_table$PRS_term <- coefficient_table$term %in% context$prs_columns
  coefficient_table$penalised <- !(coefficient_table$term %in% c("(Intercept)", context$forced_covariates))
  coefficient_table$nonzero <- abs(coefficient_table$coefficient) > tolerance

  model_metrics <- metric_row(model_name, y_test, predictions, epsilon)
  model_metrics$predictors <- length(predictors)
  model_metrics$forced_unpenalised_covariates <- sum(penalty_factor == 0)
  model_metrics$penalised_predictors <- sum(penalty_factor == 1)
  model_metrics$selected_alpha <- selected_alpha
  model_metrics$selected_lambda <- selected_lambda
  model_metrics$selected_inner_auc <- tuning_table$mean_inner_auc[[selected_row]]
  model_metrics$selected_inner_auc_standard_error <- tuning_table$inner_auc_standard_error[[selected_row]]
  model_metrics$nonzero_coefficients_including_intercept <- sum(coefficient_table$nonzero)
  model_metrics$nonzero_penalised_predictors <- sum(coefficient_table$nonzero & coefficient_table$penalised)
  model_metrics$nonzero_PRS_predictors <- sum(coefficient_table$nonzero & coefficient_table$PRS_term)

  prediction_table <- data.frame(
    participant_id = as.character(test_data$Row.names),
    outcome = as.integer(y_test),
    drug = numeric_vector(test_data$drug, "drug"),
    prediction = predictions,
    model = model_name,
    stringsAsFactors = FALSE
  )

  warning_table <- if (length(warning_rows) > 0L) {
    do.call(rbind, warning_rows)
  } else {
    data.frame(stage = character(), model = character(), alpha = numeric(), warning = character())
  }

  list(
    metrics = model_metrics,
    predictions = prediction_table,
    tuning = tuning_table,
    coefficients = coefficient_table,
    warnings = warning_table
  )
}

# Run paired comparator and combined Elastic Net fits for one outer split.
run_elastic_net_split <- function(context, parameters, comparator_label, combined_label) {
  comparator <- fit_elastic_net_model(
    comparator_label,
    context$comparator_predictors,
    context,
    parameters
  )
  combined <- fit_elastic_net_model(
    combined_label,
    context$combined_predictors,
    context,
    parameters
  )
  list(comparator = comparator, combined = combined)
}
