# Configuration-driven Random Forest fitting with nested cross-validation.

# Fit one deterministic Random Forest candidate using the supplied tuning values.
fit_random_forest <- function(training_x, training_y, mtry, nodesize, ntree, seed) {
  warnings <- character()
  set.seed(seed)
  fit <- withCallingHandlers(
    randomForest::randomForest(
      x = training_x,
      y = factor(training_y, levels = c(0, 1)),
      ntree = ntree,
      mtry = mtry,
      nodesize = nodesize,
      replace = TRUE,
      importance = FALSE,
      proximity = FALSE,
      keep.forest = TRUE,
      na.action = stats::na.fail
    ),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

# Extract the probability of remission from a fitted Random Forest.
predict_event_probability <- function(fit, new_data) {
  probabilities <- stats::predict(fit, newdata = new_data, type = "prob")
  if (!"1" %in% colnames(probabilities)) stop("Event-probability column is missing.")
  as.numeric(probabilities[, "1"])
}

# Tune Random Forest hyperparameters within the outer-training data only.
fit_random_forest_model <- function(
  model_name,
  model_position,
  predictors,
  context,
  parameters
) {
  data <- context$analysis_data
  full_matrix <- encode_predictor_matrix(data, predictors, context$binary_mapping)
  full_frame <- as.data.frame(full_matrix, check.names = FALSE)
  outcome <- numeric_vector(data$hdremit.all, "hdremit.all")
  predictor_count <- length(predictors)
  ntree <- as_integer_parameter(parameters, "random_forest_ntree")
  master_seed <- as_integer_parameter(parameters, "master_seed")
  nodesize_values <- as.integer(as_numeric_list(parameters$random_forest_nodesize))
  epsilon <- as_numeric_parameter(parameters, "probability_epsilon")

  mtry_values <- unique(pmin(
    predictor_count,
    pmax(1L, as.integer(c(
      floor(sqrt(predictor_count) / 2),
      floor(sqrt(predictor_count)),
      floor(2 * sqrt(predictor_count)),
      floor(predictor_count / 3)
    )))
  ))
  tuning_grid <- expand.grid(
    mtry = mtry_values,
    nodesize = nodesize_values,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  tuning_grid$config_id <- seq_len(nrow(tuning_grid))

  training_indices <- context$training_indices
  test_indices <- context$test_indices
  inner_rows <- list()
  warning_rows <- list()

  for (config_position in seq_len(nrow(tuning_grid))) {
    config_id <- tuning_grid$config_id[[config_position]]
    mtry <- tuning_grid$mtry[[config_position]]
    nodesize <- tuning_grid$nodesize[[config_position]]
    for (inner_fold in sort(unique(context$inner_fold_id))) {
      validation_local <- which(context$inner_fold_id == inner_fold)
      inner_training_local <- which(context$inner_fold_id != inner_fold)
      inner_training_indices <- training_indices[inner_training_local]
      validation_indices <- training_indices[validation_local]
      seed <- master_seed + context$outer_repeat * 100000L + context$outer_fold * 10000L +
        model_position * 1000L + config_id * 10L + inner_fold
      fit_result <- fit_random_forest(
        full_frame[inner_training_indices, , drop = FALSE],
        outcome[inner_training_indices],
        mtry,
        nodesize,
        ntree,
        seed
      )
      if (length(fit_result$warnings) > 0L) {
        for (warning_text in fit_result$warnings) {
          warning_rows[[length(warning_rows) + 1L]] <- data.frame(
            stage = "inner_tuning",
            model = model_name,
            config_id = config_id,
            inner_fold = inner_fold,
            warning = warning_text,
            stringsAsFactors = FALSE
          )
        }
      }
      prediction <- predict_event_probability(
        fit_result$fit,
        full_frame[validation_indices, , drop = FALSE]
      )
      y_validation <- outcome[validation_indices]
      inner_rows[[length(inner_rows) + 1L]] <- data.frame(
        model = model_name,
        predictor_count = predictor_count,
        config_id = config_id,
        mtry = mtry,
        nodesize = nodesize,
        ntree = ntree,
        inner_fold = inner_fold,
        training_participants = length(inner_training_indices),
        validation_participants = length(validation_indices),
        validation_events = sum(y_validation == 1L),
        validation_non_events = sum(y_validation == 0L),
        validation_AUC = calculate_auc(y_validation, prediction),
        validation_Brier = calculate_brier(y_validation, prediction),
        validation_log_loss = calculate_log_loss(y_validation, prediction, epsilon),
        fit_seed = seed,
        stringsAsFactors = FALSE
      )
    }
  }

  inner_results <- do.call(rbind, inner_rows)
  summary_rows <- lapply(seq_len(nrow(tuning_grid)), function(position) {
    config_id <- tuning_grid$config_id[[position]]
    rows <- inner_results[inner_results$config_id == config_id, , drop = FALSE]
    data.frame(
      model = model_name,
      predictor_count = predictor_count,
      config_id = config_id,
      mtry = tuning_grid$mtry[[position]],
      nodesize = tuning_grid$nodesize[[position]],
      ntree = ntree,
      completed_inner_folds = nrow(rows),
      mean_inner_AUC = mean(rows$validation_AUC),
      mean_inner_Brier = mean(rows$validation_Brier),
      mean_inner_log_loss = mean(rows$validation_log_loss),
      stringsAsFactors = FALSE
    )
  })
  tuning_summary <- do.call(rbind, summary_rows)
  selection_order <- order(
    -tuning_summary$mean_inner_AUC,
    tuning_summary$mean_inner_Brier,
    tuning_summary$mean_inner_log_loss,
    -tuning_summary$nodesize,
    tuning_summary$mtry,
    tuning_summary$config_id
  )
  selected_position <- selection_order[[1]]
  tuning_summary$selected <- seq_len(nrow(tuning_summary)) == selected_position
  selected <- tuning_summary[selected_position, , drop = FALSE]

  final_seed <- master_seed + context$outer_repeat * 1000000L + context$outer_fold * 100000L +
    model_position * 10000L + 999L
  final_fit <- fit_random_forest(
    full_frame[training_indices, , drop = FALSE],
    outcome[training_indices],
    as.integer(selected$mtry[[1]]),
    as.integer(selected$nodesize[[1]]),
    ntree,
    final_seed
  )
  if (length(final_fit$warnings) > 0L) {
    for (warning_text in final_fit$warnings) {
      warning_rows[[length(warning_rows) + 1L]] <- data.frame(
        stage = "outer_final_fit",
        model = model_name,
        config_id = selected$config_id[[1]],
        inner_fold = NA_integer_,
        warning = warning_text,
        stringsAsFactors = FALSE
      )
    }
  }

  prediction <- predict_event_probability(final_fit$fit, full_frame[test_indices, , drop = FALSE])
  y_test <- outcome[test_indices]
  metrics <- metric_row(model_name, y_test, prediction, epsilon)
  metrics$predictors <- predictor_count
  metrics$selected_config_id <- selected$config_id[[1]]
  metrics$selected_mtry <- selected$mtry[[1]]
  metrics$selected_nodesize <- selected$nodesize[[1]]
  metrics$ntree <- ntree
  metrics$final_fit_seed <- final_seed

  predictions <- data.frame(
    participant_id = as.character(data$Row.names[test_indices]),
    outcome = as.integer(y_test),
    drug = numeric_vector(data$drug[test_indices], "drug"),
    prediction = prediction,
    model = model_name,
    stringsAsFactors = FALSE
  )

  warnings <- if (length(warning_rows) > 0L) {
    do.call(rbind, warning_rows)
  } else {
    data.frame(stage = character(), model = character(), config_id = integer(), inner_fold = integer(), warning = character())
  }
  coefficients <- data.frame(
    model = character(), term = character(), coefficient = numeric(),
    forced_covariate = logical(), PRS_term = logical(), penalised = logical(), nonzero = logical()
  )

  list(
    metrics = metrics,
    predictions = predictions,
    tuning = tuning_summary,
    inner_results = inner_results,
    coefficients = coefficients,
    warnings = warnings
  )
}

# Run paired comparator and combined Random Forest fits for one outer split.
run_random_forest_split <- function(context, parameters, comparator_label, combined_label) {
  comparator <- fit_random_forest_model(
    comparator_label, 1L, context$comparator_predictors, context, parameters
  )
  combined <- fit_random_forest_model(
    combined_label, 2L, context$combined_predictors, context, parameters
  )
  list(comparator = comparator, combined = combined)
}
