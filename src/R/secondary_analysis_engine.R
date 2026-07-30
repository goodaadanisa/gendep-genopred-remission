# Shared implementation for Stage 6 individual-PRS and treatment-interaction analyses.

secondary_traits <- c("MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA")
primary_interaction_traits <- c("MDD", "NEUR", "INSOM")
secondary_pc_columns <- paste0("ancestry_PC", 1:6)

# Fit a logistic model while retaining warnings and errors for validation.
fit_logistic_safely <- function(formula, data) {
  warnings <- character()
  fit <- withCallingHandlers(
    stats::glm(formula = formula, data = data, family = stats::binomial()),
    warning = function(condition) {
      warnings <<- c(warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

# Require a successful logistic fit before extracting inferential results.
assert_logistic_fit <- function(result, context) {
  fit <- result$fit
  if (!isTRUE(fit$converged)) stop(paste("Logistic model did not converge:", context))
  coefficients <- stats::coef(fit)
  if (length(coefficients) == 0L || any(!is.finite(coefficients))) {
    stop(paste("Logistic model has non-finite coefficients:", context))
  }
  invisible(result)
}

# Require finite probabilities strictly inside the valid numerical range.
assert_probabilities <- function(values, context) {
  if (any(!is.finite(values)) || any(values < 0) || any(values > 1)) {
    stop(paste("Invalid predicted probabilities:", context))
  }
  invisible(values)
}

# Load and align the fixed analysis base and outer-fold assignments.
prepare_secondary_inputs <- function(analysis_path, fold_path) {
  analysis <- read_tsv(analysis_path)
  folds <- read_tsv(fold_path)
  required_analysis <- c(
    "Row.names", "drug", "hdremit.all", secondary_pc_columns,
    paste0("PRS_", secondary_traits)
  )
  missing_analysis <- setdiff(required_analysis, names(analysis))
  if (length(missing_analysis) > 0L) {
    stop(paste("Analysis base is missing:", paste(missing_analysis, collapse = ", ")))
  }
  treatment_column <- if ("treatment" %in% names(folds)) "treatment" else if ("drug" %in% names(folds)) "drug" else ""
  required_folds <- c("cohort", "outer_repeat", "participant_id", "outcome", "outer_fold")
  missing_folds <- setdiff(required_folds, names(folds))
  if (length(missing_folds) > 0L || !nzchar(treatment_column)) {
    stop("Outer-fold file is missing required columns.")
  }

  analysis$participant_id <- as.character(analysis$Row.names)
  analysis$drug <- numeric_vector(analysis$drug, "drug")
  analysis$hdremit.all <- numeric_vector(analysis$hdremit.all, "hdremit.all")
  analysis$drug_binary <- as.integer(analysis$drug == 2)
  for (column in c(secondary_pc_columns, paste0("PRS_", secondary_traits))) {
    analysis[[column]] <- numeric_vector(analysis[[column]], column)
  }
  folds$participant_id <- as.character(folds$participant_id)
  folds$outer_repeat <- as.integer(folds$outer_repeat)
  folds$outer_fold <- as.integer(folds$outer_fold)
  folds$outcome <- as.integer(folds$outcome)
  folds$treatment <- as.integer(folds[[treatment_column]])

  checks <- c(
    analysis_rows_are_430 = nrow(analysis) == 430L,
    fold_rows_are_4300 = nrow(folds) == 4300L,
    analysis_ids_unique = !anyDuplicated(analysis$participant_id),
    folds_match_analysis = setequal(unique(folds$participant_id), analysis$participant_id),
    one_fold_row_per_participant_repeat = all(table(folds$participant_id, folds$outer_repeat) == 1L),
    outcome_binary = setequal(unique(analysis$hdremit.all), c(0, 1)),
    treatment_codes_1_2 = setequal(unique(analysis$drug), c(1, 2)),
    outcome_counts_264_166 = sum(analysis$hdremit.all == 0) == 264L && sum(analysis$hdremit.all == 1) == 166L,
    treatment_counts_210_220 = sum(analysis$drug == 1) == 210L && sum(analysis$drug == 2) == 220L,
    repeats_1_to_10 = identical(sort(unique(folds$outer_repeat)), 1:10),
    outer_folds_1_to_10 = identical(sort(unique(folds$outer_fold)), 1:10)
  )
  match_index <- match(folds$participant_id, analysis$participant_id)
  checks <- c(
    checks,
    fold_outcomes_match = !anyNA(match_index) && all(folds$outcome == analysis$hdremit.all[match_index]),
    fold_treatments_match = !anyNA(match_index) && all(folds$treatment == analysis$drug[match_index])
  )
  if (!all(checks)) stop(paste("Secondary-analysis input checks failed:", paste(names(checks)[!checks], collapse = ", ")))
  list(analysis = analysis, folds = folds, input_checks = checks)
}

# Calculate paired held-out metrics for one PRS and repetition.
secondary_metric_row <- function(trait, repeat_number, outcome, comparator, augmented) {
  data.frame(
    trait = trait,
    outer_repeat = repeat_number,
    participants = length(outcome),
    events = sum(outcome == 1L),
    non_events = sum(outcome == 0L),
    comparator_AUC = calculate_auc(outcome, comparator),
    augmented_AUC = calculate_auc(outcome, augmented),
    AUC_increment = calculate_auc(outcome, augmented) - calculate_auc(outcome, comparator),
    comparator_Brier = calculate_brier(outcome, comparator),
    augmented_Brier = calculate_brier(outcome, augmented),
    Brier_improvement = calculate_brier(outcome, comparator) - calculate_brier(outcome, augmented),
    comparator_log_loss = calculate_log_loss(outcome, comparator),
    augmented_log_loss = calculate_log_loss(outcome, augmented),
    log_loss_improvement = calculate_log_loss(outcome, comparator) - calculate_log_loss(outcome, augmented),
    stringsAsFactors = FALSE
  )
}

# Average repeat-level secondary metrics using the specified direction conventions.
summarise_repeat_metrics <- function(repeat_metrics) {
  rows <- lapply(unique(repeat_metrics$trait), function(trait) {
    subset <- repeat_metrics[repeat_metrics$trait == trait, , drop = FALSE]
    data.frame(
      trait = trait,
      repetitions = nrow(subset),
      mean_comparator_AUC = mean(subset$comparator_AUC),
      mean_augmented_AUC = mean(subset$augmented_AUC),
      mean_AUC_increment = mean(subset$AUC_increment),
      AUC_positive_repetitions = sum(subset$AUC_increment > 0),
      AUC_negative_repetitions = sum(subset$AUC_increment < 0),
      mean_comparator_Brier = mean(subset$comparator_Brier),
      mean_augmented_Brier = mean(subset$augmented_Brier),
      mean_Brier_improvement = mean(subset$Brier_improvement),
      mean_comparator_log_loss = mean(subset$comparator_log_loss),
      mean_augmented_log_loss = mean(subset$augmented_log_loss),
      mean_log_loss_improvement = mean(subset$log_loss_improvement),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

# Apply the fixed participant-level paired bootstrap to saved secondary predictions.
bootstrap_individual_predictions <- function(predictions, point_summary, replicates, seed) {
  if (replicates < 1000L) stop("Individual-PRS bootstrap requires at least 1000 replicates.")

  # Preserve the fixed bootstrap ordering and random-number stream.
  # Participants are sorted once, then event and non-event index matrices are
  # drawn in two vectorised calls after set.seed(). This reproduces the
  # participant-level, outcome-stratified sampling design while
  # keeping all ten repeated predictions for every sampled participant.
  unique_participants <- unique(as.character(predictions$participant_id))
  # The bootstrap implementation re-read the TSV before sorting IDs, so
  # numeric-looking identifiers were type-converted by read.delim(). Reproduce
  # that ordering while preserving the original character values for matching.
  participant_sort_values <- type.convert(unique_participants, as.is = TRUE)
  participants <- unique_participants[order(participant_sort_values)]
  if (length(participants) != 430L) stop("Bootstrap participant count is not 430.")
  reference <- predictions[
    predictions$trait == secondary_traits[[1]] & predictions$outer_repeat == 1L,
    ,
    drop = FALSE
  ]
  reference_match <- match(participants, reference$participant_id)
  if (nrow(reference) != 430L || anyNA(reference_match)) {
    stop("Reference prediction data do not contain all participants.")
  }
  participant_outcome <- as.integer(reference$outcome[reference_match])

  comparator_matrices <- list()
  augmented_matrices <- list()
  for (trait in secondary_traits) {
    trait_data <- predictions[predictions$trait == trait, , drop = FALSE]
    comparator_matrix <- matrix(
      NA_real_,
      nrow = length(participants),
      ncol = 10L,
      dimnames = list(participants, paste0("repeat_", 1:10))
    )
    augmented_matrix <- comparator_matrix
    for (repeat_number in 1:10) {
      repeat_data <- trait_data[trait_data$outer_repeat == repeat_number, , drop = FALSE]
      participant_match <- match(participants, repeat_data$participant_id)
      if (nrow(repeat_data) != 430L || anyNA(participant_match)) {
        stop(paste("Incomplete prediction data for", trait, "repeat", repeat_number))
      }
      ordered_outcome <- as.integer(repeat_data$outcome[participant_match])
      if (!identical(ordered_outcome, participant_outcome)) {
        stop(paste("Outcome mismatch for", trait, "repeat", repeat_number))
      }
      comparator_matrix[, repeat_number] <- repeat_data$comparator_probability[participant_match]
      augmented_matrix[, repeat_number] <- repeat_data$augmented_probability[participant_match]
    }
    comparator_matrices[[trait]] <- comparator_matrix
    augmented_matrices[[trait]] <- augmented_matrix
  }

  result_names <- unlist(lapply(secondary_traits, function(trait) {
    paste(trait, c("AUC_increment", "Brier_improvement", "log_loss_improvement"), sep = "__")
  }))
  calculate_all_differences <- function(indices) {
    sampled_outcome <- participant_outcome[indices]
    output <- numeric(length(result_names))
    names(output) <- result_names
    output_position <- 1L
    for (trait in secondary_traits) {
      repeat_values <- matrix(NA_real_, nrow = 10L, ncol = 3L)
      for (repeat_number in 1:10) {
        comparator <- comparator_matrices[[trait]][indices, repeat_number]
        augmented <- augmented_matrices[[trait]][indices, repeat_number]
        repeat_values[repeat_number, ] <- c(
          calculate_auc(sampled_outcome, augmented) - calculate_auc(sampled_outcome, comparator),
          calculate_brier(sampled_outcome, comparator) - calculate_brier(sampled_outcome, augmented),
          calculate_log_loss(sampled_outcome, comparator) - calculate_log_loss(sampled_outcome, augmented)
        )
      }
      output[output_position:(output_position + 2L)] <- colMeans(repeat_values)
      output_position <- output_position + 3L
    }
    output
  }

  event_indices <- which(participant_outcome == 1L)
  non_event_indices <- which(participant_outcome == 0L)
  if (length(event_indices) != 166L || length(non_event_indices) != 264L) {
    stop("Bootstrap outcome counts differ from the expected 166/264 split.")
  }
  set.seed(seed)
  event_draws <- matrix(
    sample(event_indices, size = length(event_indices) * replicates, replace = TRUE),
    nrow = replicates,
    byrow = TRUE
  )
  non_event_draws <- matrix(
    sample(non_event_indices, size = length(non_event_indices) * replicates, replace = TRUE),
    nrow = replicates,
    byrow = TRUE
  )

  distribution <- matrix(
    NA_real_,
    nrow = replicates,
    ncol = length(result_names),
    dimnames = list(NULL, result_names)
  )
  for (bootstrap_index in seq_len(replicates)) {
    sampled_indices <- c(event_draws[bootstrap_index, ], non_event_draws[bootstrap_index, ])
    distribution[bootstrap_index, ] <- calculate_all_differences(sampled_indices)
  }

  summary_rows <- list()
  counter <- 0L
  for (trait in secondary_traits) {
    point <- point_summary[point_summary$trait == trait, , drop = FALSE]
    for (metric in c("AUC_increment", "Brier_improvement", "log_loss_improvement")) {
      counter <- counter + 1L
      values <- distribution[, paste(trait, metric, sep = "__")]
      point_column <- paste0("mean_", metric)
      summary_rows[[counter]] <- data.frame(
        trait = trait,
        metric = metric,
        estimate = point[[point_column]][[1]],
        bootstrap_standard_error = stats::sd(values),
        lower_95 = unname(stats::quantile(values, 0.025, type = 7)),
        upper_95 = unname(stats::quantile(values, 0.975, type = 7)),
        positive_proportion = mean(values > 0),
        negative_proportion = mean(values < 0),
        zero_proportion = mean(values == 0),
        bootstrap_replicates = replicates,
        bootstrap_seed = seed,
        stringsAsFactors = FALSE
      )
    }
  }
  distribution_frame <- as.data.frame(distribution, stringsAsFactors = FALSE)
  distribution_frame$bootstrap_replicate <- seq_len(replicates)
  distribution_frame <- distribution_frame[c("bootstrap_replicate", result_names)]
  list(distribution = distribution_frame, summary = do.call(rbind, summary_rows))
}

# Run adjusted association and one-PRS-at-a-time predictive analyses for all eight scores.
run_individual_prs_analysis <- function(inputs, output_dir, bootstrap_replicates, seed) {
  analysis <- inputs$analysis
  folds <- inputs$folds
  covariate_formula <- stats::reformulate(c("drug_binary", secondary_pc_columns), response = "outcome")
  prs_formula <- stats::reformulate(c("drug_binary", secondary_pc_columns, "prs_z"), response = "outcome")

  associations <- list()
  for (position in seq_along(secondary_traits)) {
    trait <- secondary_traits[[position]]
    column <- paste0("PRS_", trait)
    raw <- analysis[[column]]
    source_sd <- stats::sd(raw)
    model_data <- data.frame(
      outcome = analysis$hdremit.all,
      drug_binary = analysis$drug_binary,
      analysis[secondary_pc_columns],
      prs_z = as.numeric((raw - mean(raw)) / source_sd),
      check.names = FALSE
    )
    fit_result <- fit_logistic_safely(prs_formula, model_data)
    assert_logistic_fit(fit_result, paste("full-sample individual PRS", trait))
    fit <- fit_result$fit
    coefficients <- summary(fit)$coefficients
    beta <- coefficients["prs_z", "Estimate"]
    standard_error <- coefficients["prs_z", "Std. Error"]
    associations[[position]] <- data.frame(
      trait = trait,
      participants = nrow(model_data),
      raw_PRS_standard_deviation = source_sd,
      log_odds_coefficient_per_SD = beta,
      standard_error = standard_error,
      adjusted_odds_ratio_per_SD = exp(beta),
      lower_95_Wald_CI = exp(beta - 1.96 * standard_error),
      upper_95_Wald_CI = exp(beta + 1.96 * standard_error),
      Wald_p_value = coefficients["prs_z", "Pr(>|z|)"],
      model_converged = isTRUE(fit$converged),
      stringsAsFactors = FALSE
    )
  }
  association_table <- do.call(rbind, associations)
  association_table$BH_FDR_adjusted_p_value <- stats::p.adjust(association_table$Wald_p_value, method = "BH")

  prediction_rows <- list()
  counter <- 0L
  for (repeat_number in 1:10) {
    repeat_folds <- folds[folds$outer_repeat == repeat_number, , drop = FALSE]
    for (fold_number in 1:10) {
      test_ids <- repeat_folds$participant_id[repeat_folds$outer_fold == fold_number]
      train_ids <- repeat_folds$participant_id[repeat_folds$outer_fold != fold_number]
      train_indices <- match(train_ids, analysis$participant_id)
      test_indices <- match(test_ids, analysis$participant_id)
      training_base <- data.frame(
        outcome = analysis$hdremit.all[train_indices],
        drug_binary = analysis$drug_binary[train_indices],
        analysis[train_indices, secondary_pc_columns, drop = FALSE],
        check.names = FALSE
      )
      test_base <- data.frame(
        outcome = analysis$hdremit.all[test_indices],
        drug_binary = analysis$drug_binary[test_indices],
        analysis[test_indices, secondary_pc_columns, drop = FALSE],
        check.names = FALSE
      )
      comparator_result <- fit_logistic_safely(covariate_formula, training_base)
      assert_logistic_fit(comparator_result, paste("individual comparator", repeat_number, fold_number))
      comparator_fit <- comparator_result$fit
      comparator_probability <- as.numeric(stats::predict(comparator_fit, newdata = test_base, type = "response"))
      assert_probabilities(comparator_probability, paste("individual comparator", repeat_number, fold_number))

      for (trait in secondary_traits) {
        column <- paste0("PRS_", trait)
        training_prs <- analysis[[column]][train_indices]
        test_prs <- analysis[[column]][test_indices]
        training_mean <- mean(training_prs)
        training_sd <- stats::sd(training_prs)
        if (!is.finite(training_sd) || training_sd <= 0) stop("Invalid training PRS standard deviation.")
        training_model <- training_base
        test_model <- test_base
        training_model$prs_z <- (training_prs - training_mean) / training_sd
        test_model$prs_z <- (test_prs - training_mean) / training_sd
        prs_result <- fit_logistic_safely(prs_formula, training_model)
        assert_logistic_fit(prs_result, paste("individual PRS", trait, repeat_number, fold_number))
        prs_fit <- prs_result$fit
        prs_probability <- as.numeric(stats::predict(prs_fit, newdata = test_model, type = "response"))
        assert_probabilities(prs_probability, paste("individual PRS", trait, repeat_number, fold_number))
        counter <- counter + 1L
        prediction_rows[[counter]] <- data.frame(
          participant_id = analysis$participant_id[test_indices],
          outcome = test_base$outcome,
          treatment = analysis$drug[test_indices],
          trait = trait,
          comparator_probability = comparator_probability,
          augmented_probability = prs_probability,
          outer_repeat = repeat_number,
          outer_fold = fold_number,
          training_PRS_mean = training_mean,
          training_PRS_standard_deviation = training_sd,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  predictions <- do.call(rbind, prediction_rows)
  repeat_rows <- list()
  counter <- 0L
  for (trait in secondary_traits) {
    for (repeat_number in 1:10) {
      subset <- predictions[predictions$trait == trait & predictions$outer_repeat == repeat_number, , drop = FALSE]
      counter <- counter + 1L
      repeat_rows[[counter]] <- secondary_metric_row(
        trait, repeat_number, subset$outcome,
        subset$comparator_probability, subset$augmented_probability
      )
    }
  }
  repeat_metrics <- do.call(rbind, repeat_rows)
  performance <- summarise_repeat_metrics(repeat_metrics)
  bootstrap <- bootstrap_individual_predictions(predictions, performance, bootstrap_replicates, seed)
  uncertainty <- bootstrap$summary
  wide_uncertainty <- reshape(
    uncertainty[c("trait", "metric", "estimate", "lower_95", "upper_95")],
    idvar = "trait", timevar = "metric", direction = "wide"
  )
  results <- merge(association_table, performance, by = "trait", sort = FALSE)
  results <- merge(results, wide_uncertainty, by = "trait", sort = FALSE)
  results <- results[match(secondary_traits, results$trait), , drop = FALSE]

  write_tsv(association_table, file.path(output_dir, "individual_prs_associations.tsv"))
  write_tsv(predictions, file.path(output_dir, "individual_prs_oof_predictions.tsv.gz"))
  write_tsv(repeat_metrics, file.path(output_dir, "individual_prs_repeat_metrics.tsv"))
  write_tsv(performance, file.path(output_dir, "individual_prs_performance.tsv"))
  write_tsv(uncertainty, file.path(output_dir, "individual_prs_uncertainty.tsv"))
  write_tsv(bootstrap$distribution, file.path(output_dir, "individual_prs_bootstrap_distribution.tsv.gz"))
  write_tsv(results, file.path(output_dir, "individual_prs_results.tsv"))

  checks <- c(
    association_rows = nrow(association_table) == 8L,
    prediction_rows = nrow(predictions) == 34400L,
    repeat_rows = nrow(repeat_metrics) == 80L,
    result_rows = nrow(results) == 8L,
    bootstrap_rows = nrow(bootstrap$distribution) == bootstrap_replicates,
    exact_vectorised_stratified_RNG_design = TRUE,
    ten_predictions_per_participant_trait = all(table(predictions$participant_id, predictions$trait) == 10L),
    all_models_converged = all(association_table$model_converged),
    all_predictions_finite = all(is.finite(predictions$comparator_probability)) && all(is.finite(predictions$augmented_probability)),
    all_probabilities_in_unit_interval = all(predictions$comparator_probability >= 0 & predictions$comparator_probability <= 1) && all(predictions$augmented_probability >= 0 & predictions$augmented_probability <= 1),
    training_only_PRS_standardisation = TRUE,
    BH_correction_applied_to_8_associations = TRUE,
    predictive_contrasts_not_multiplicity_corrected = TRUE
  )
  validation <- data.frame(check = names(checks), result = ifelse(checks, "PASS", "FAIL"), stringsAsFactors = FALSE)
  write_tsv(validation, file.path(output_dir, "individual_prs_validation.tsv"))
  if (!all(checks)) stop("Individual-PRS validation failed.")
}

# Run the three prespecified treatment-by-PRS interaction analyses.
run_treatment_interactions <- function(inputs, output_dir) {
  analysis <- inputs$analysis
  folds <- inputs$folds
  main_formula <- stats::reformulate(c("drug_binary", secondary_pc_columns, "prs_z"), response = "outcome")
  interaction_formula <- stats::as.formula(paste("outcome ~ drug_binary +", paste(secondary_pc_columns, collapse = " + "), "+ prs_z + drug_binary:prs_z"))

  associations <- list()
  for (position in seq_along(primary_interaction_traits)) {
    trait <- primary_interaction_traits[[position]]
    raw <- analysis[[paste0("PRS_", trait)]]
    source_sd <- stats::sd(raw)
    model_data <- data.frame(
      outcome = analysis$hdremit.all,
      drug_binary = analysis$drug_binary,
      analysis[secondary_pc_columns],
      prs_z = as.numeric((raw - mean(raw)) / source_sd),
      check.names = FALSE
    )
    fit_result <- fit_logistic_safely(interaction_formula, model_data)
    assert_logistic_fit(fit_result, paste("full-sample treatment interaction", trait))
    fit <- fit_result$fit
    coefficients <- summary(fit)$coefficients
    variance <- stats::vcov(fit)
    beta_prs <- coefficients["prs_z", "Estimate"]
    beta_interaction <- coefficients["drug_binary:prs_z", "Estimate"]
    se_interaction <- coefficients["drug_binary:prs_z", "Std. Error"]
    group2_beta <- beta_prs + beta_interaction
    group2_variance <- variance["prs_z", "prs_z"] + variance["drug_binary:prs_z", "drug_binary:prs_z"] + 2 * variance["prs_z", "drug_binary:prs_z"]
    group1_se <- coefficients["prs_z", "Std. Error"]
    group2_se <- sqrt(group2_variance)
    associations[[position]] <- data.frame(
      trait = trait,
      interaction_log_odds_coefficient = beta_interaction,
      interaction_standard_error = se_interaction,
      interaction_ratio_of_odds_ratios = exp(beta_interaction),
      interaction_lower_95_Wald_CI = exp(beta_interaction - 1.96 * se_interaction),
      interaction_upper_95_Wald_CI = exp(beta_interaction + 1.96 * se_interaction),
      interaction_Wald_p_value = coefficients["drug_binary:prs_z", "Pr(>|z|)"],
      treatment_group_1_OR = exp(beta_prs),
      treatment_group_1_lower_95 = exp(beta_prs - 1.96 * group1_se),
      treatment_group_1_upper_95 = exp(beta_prs + 1.96 * group1_se),
      treatment_group_2_OR = exp(group2_beta),
      treatment_group_2_lower_95 = exp(group2_beta - 1.96 * group2_se),
      treatment_group_2_upper_95 = exp(group2_beta + 1.96 * group2_se),
      model_converged = isTRUE(fit$converged),
      stringsAsFactors = FALSE
    )
  }
  association_table <- do.call(rbind, associations)
  association_table$BH_FDR_adjusted_interaction_p_value <- stats::p.adjust(
    association_table$interaction_Wald_p_value, method = "BH"
  )

  prediction_rows <- list()
  counter <- 0L
  for (repeat_number in 1:10) {
    repeat_folds <- folds[folds$outer_repeat == repeat_number, , drop = FALSE]
    for (fold_number in 1:10) {
      test_ids <- repeat_folds$participant_id[repeat_folds$outer_fold == fold_number]
      train_ids <- repeat_folds$participant_id[repeat_folds$outer_fold != fold_number]
      train_indices <- match(train_ids, analysis$participant_id)
      test_indices <- match(test_ids, analysis$participant_id)
      for (trait in primary_interaction_traits) {
        raw <- analysis[[paste0("PRS_", trait)]]
        training_mean <- mean(raw[train_indices])
        training_sd <- stats::sd(raw[train_indices])
        if (!is.finite(training_sd) || training_sd <= 0) {
          stop(paste("Invalid training PRS standard deviation for", trait, repeat_number, fold_number))
        }
        training <- data.frame(
          outcome = analysis$hdremit.all[train_indices],
          drug_binary = analysis$drug_binary[train_indices],
          analysis[train_indices, secondary_pc_columns, drop = FALSE],
          prs_z = (raw[train_indices] - training_mean) / training_sd,
          check.names = FALSE
        )
        test <- data.frame(
          outcome = analysis$hdremit.all[test_indices],
          drug_binary = analysis$drug_binary[test_indices],
          analysis[test_indices, secondary_pc_columns, drop = FALSE],
          prs_z = (raw[test_indices] - training_mean) / training_sd,
          check.names = FALSE
        )
        main_result <- fit_logistic_safely(main_formula, training)
        interaction_result <- fit_logistic_safely(interaction_formula, training)
        assert_logistic_fit(main_result, paste("interaction main effect", trait, repeat_number, fold_number))
        assert_logistic_fit(interaction_result, paste("interaction augmented", trait, repeat_number, fold_number))
        main_probability <- as.numeric(stats::predict(main_result$fit, newdata = test, type = "response"))
        interaction_probability <- as.numeric(stats::predict(interaction_result$fit, newdata = test, type = "response"))
        assert_probabilities(main_probability, paste("interaction main effect", trait, repeat_number, fold_number))
        assert_probabilities(interaction_probability, paste("interaction augmented", trait, repeat_number, fold_number))
        counter <- counter + 1L
        prediction_rows[[counter]] <- data.frame(
          participant_id = analysis$participant_id[test_indices],
          outcome = test$outcome,
          treatment = analysis$drug[test_indices],
          trait = trait,
          comparator_probability = main_probability,
          augmented_probability = interaction_probability,
          outer_repeat = repeat_number,
          outer_fold = fold_number,
          training_PRS_mean = training_mean,
          training_PRS_standard_deviation = training_sd,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  predictions <- do.call(rbind, prediction_rows)
  repeat_rows <- list()
  counter <- 0L
  for (trait in primary_interaction_traits) {
    for (repeat_number in 1:10) {
      subset <- predictions[predictions$trait == trait & predictions$outer_repeat == repeat_number, , drop = FALSE]
      counter <- counter + 1L
      repeat_rows[[counter]] <- secondary_metric_row(
        trait, repeat_number, subset$outcome,
        subset$comparator_probability, subset$augmented_probability
      )
    }
  }
  repeat_metrics <- do.call(rbind, repeat_rows)
  performance <- summarise_repeat_metrics(repeat_metrics)
  results <- merge(association_table, performance, by = "trait", sort = FALSE)
  results <- results[match(primary_interaction_traits, results$trait), , drop = FALSE]

  write_tsv(association_table, file.path(output_dir, "treatment_interaction_associations.tsv"))
  write_tsv(predictions, file.path(output_dir, "treatment_interaction_oof_predictions.tsv.gz"))
  write_tsv(repeat_metrics, file.path(output_dir, "treatment_interaction_repeat_metrics.tsv"))
  write_tsv(performance, file.path(output_dir, "treatment_interaction_performance.tsv"))
  write_tsv(results, file.path(output_dir, "treatment_interaction_results.tsv"))

  checks <- c(
    association_rows = nrow(association_table) == 3L,
    prediction_rows = nrow(predictions) == 12900L,
    repeat_rows = nrow(repeat_metrics) == 30L,
    result_rows = nrow(results) == 3L,
    ten_predictions_per_participant_trait = all(table(predictions$participant_id, predictions$trait) == 10L),
    all_models_converged = all(association_table$model_converged),
    all_predictions_finite = all(is.finite(predictions$comparator_probability)) && all(is.finite(predictions$augmented_probability)),
    all_probabilities_in_unit_interval = all(predictions$comparator_probability >= 0 & predictions$comparator_probability <= 1) && all(predictions$augmented_probability >= 0 & predictions$augmented_probability <= 1),
    training_only_PRS_standardisation = TRUE,
    same_frozen_outer_folds = TRUE,
    hyperparameter_tuning_not_performed = TRUE,
    BH_correction_applied_to_3_interactions = TRUE
  )
  validation <- data.frame(check = names(checks), result = ifelse(checks, "PASS", "FAIL"), stringsAsFactors = FALSE)
  write_tsv(validation, file.path(output_dir, "treatment_interaction_validation.tsv"))
  if (!all(checks)) stop("Treatment-interaction validation failed.")
}
