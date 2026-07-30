#!/usr/bin/env Rscript
# Estimate paired participant-level uncertainty from fixed repeated-CV predictions.
#
# Purpose:
#   Resample participants within remission strata while retaining all ten held-out
#   predictions for each sampled participant. Models are never refitted or retuned.
#
# Inputs:
#   Aggregated fixed held-out prediction table for one configured paired analysis.
#
# Outputs:
#   Metric point estimates, 10,000-replicate percentile intervals and validation
#   records under the requested output directory.
#
# Example:
#   Rscript scripts/bootstrap_uncertainty.R --analysis primary \
#     --predictions work/modelling/aggregated/primary/all_held_out_predictions.tsv.gz \
#     --output-dir work/modelling/uncertainty/primary

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_arg, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else "scripts/bootstrap_uncertainty.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "src", "R", "model_common.R"))

options <- parse_cli(args)
analysis <- required_cli(options, "analysis")
prediction_path <- required_cli(options, "predictions")
output_dir <- required_cli(options, "output-dir")
analysis_config <- options[["analysis-config"]] %||% file.path(root, "config", "model_analyses.tsv")
parameter_config <- options[["parameter-config"]] %||% file.path(root, "config", "model_parameters.tsv")
analysis_spec <- read_analysis_spec(analysis_config, analysis)
parameters <- read_parameters(parameter_config)
replicates <- suppressWarnings(as.integer(options[["replicates"]] %||% parameters$bootstrap_replicates))
seed <- suppressWarnings(as.integer(options[["seed"]] %||% parameters$bootstrap_seed))
include_calibration <- isTRUE(options[["include-calibration"]])
save_distribution <- isTRUE(options[["save-distribution"]])
epsilon <- as_numeric_parameter(parameters, "probability_epsilon")
outer_repeats <- as_integer_parameter(parameters, "outer_repeats")

if (is.na(replicates) || replicates < 1L) stop("--replicates must be a positive integer.")
if (is.na(seed)) stop("Bootstrap seed must be an integer.")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

predictions <- read_tsv(prediction_path)
required <- c(
  "outer_repeat", "participant_id", "outcome",
  "comparator_probability", "combined_probability"
)
missing <- setdiff(required, names(predictions))
if (length(missing) > 0L) stop(paste("Predictions missing:", paste(missing, collapse = ", ")))
predictions$participant_id <- as.character(predictions$participant_id)
predictions$outcome <- as.integer(predictions$outcome)
predictions$outer_repeat <- as.integer(predictions$outer_repeat)
if (anyNA(predictions[, required])) stop("Predictions contain missing values.")
if (!setequal(unique(predictions$outcome), c(0L, 1L))) stop("Outcome must contain both classes.")
if (length(unique(predictions$outer_repeat)) != outer_repeats) stop("Unexpected number of outer repeats.")

participant_outcomes <- unique(predictions[, c("participant_id", "outcome"), drop = FALSE])
participant_outcomes <- participant_outcomes[order(participant_outcomes$participant_id), , drop = FALSE]
rownames(participant_outcomes) <- NULL
if (anyDuplicated(participant_outcomes$participant_id)) stop("Participant outcome changes across repeats.")
expected_n <- as.integer(analysis_spec$expected_participants)
if (nrow(participant_outcomes) != expected_n) stop("Unexpected participant count.")
repeat_sizes <- table(predictions$outer_repeat)
if (!all(repeat_sizes == expected_n)) stop("Every repeat must contain one prediction per participant.")
if (anyDuplicated(predictions[, c("outer_repeat", "participant_id")])) stop("Duplicate repeat-participant prediction.")

repeat_frames <- lapply(seq_len(outer_repeats), function(repeat_id) {
  frame <- predictions[predictions$outer_repeat == repeat_id, , drop = FALSE]
  frame[match(participant_outcomes$participant_id, frame$participant_id), , drop = FALSE]
})

repeat_metrics_for_ids <- function(sampled_ids, include_calibration_metrics = FALSE) {
  rows <- vector("list", outer_repeats)
  for (repeat_id in seq_len(outer_repeats)) {
    frame <- repeat_frames[[repeat_id]]
    sampled <- frame[match(sampled_ids, frame$participant_id), , drop = FALSE]
    y <- sampled$outcome
    comparator <- sampled$comparator_probability
    combined <- sampled$combined_probability
    row <- c(
      comparator_AUC = calculate_auc(y, comparator),
      combined_AUC = calculate_auc(y, combined),
      AUC_increment = calculate_auc(y, combined) - calculate_auc(y, comparator),
      comparator_Brier = calculate_brier(y, comparator),
      combined_Brier = calculate_brier(y, combined),
      Brier_increment = calculate_brier(y, comparator) - calculate_brier(y, combined),
      comparator_log_loss = calculate_log_loss(y, comparator, epsilon),
      combined_log_loss = calculate_log_loss(y, combined, epsilon),
      log_loss_increment = calculate_log_loss(y, comparator, epsilon) - calculate_log_loss(y, combined, epsilon)
    )
    if (include_calibration_metrics) {
      comparator_cal <- calculate_calibration(y, comparator, epsilon)
      combined_cal <- calculate_calibration(y, combined, epsilon)
      row <- c(
        row,
        comparator_calibration_intercept = comparator_cal[["calibration_intercept"]],
        combined_calibration_intercept = combined_cal[["calibration_intercept"]],
        calibration_intercept_difference = combined_cal[["calibration_intercept"]] - comparator_cal[["calibration_intercept"]],
        comparator_calibration_slope = comparator_cal[["calibration_slope"]],
        combined_calibration_slope = combined_cal[["calibration_slope"]],
        calibration_slope_difference = combined_cal[["calibration_slope"]] - comparator_cal[["calibration_slope"]]
      )
    }
    rows[[repeat_id]] <- row
  }
  colMeans(do.call(rbind, rows))
}

all_ids <- participant_outcomes$participant_id
non_event_ids <- participant_outcomes$participant_id[participant_outcomes$outcome == 0L]
event_ids <- participant_outcomes$participant_id[participant_outcomes$outcome == 1L]
point_metrics <- repeat_metrics_for_ids(all_ids, include_calibration_metrics = FALSE)

performance_names <- c(
  "comparator_AUC", "combined_AUC", "AUC_increment",
  "comparator_Brier", "combined_Brier", "Brier_increment",
  "comparator_log_loss", "combined_log_loss", "log_loss_increment"
)
case_indices <- which(participant_outcomes$outcome == 1L)
control_indices <- which(participant_outcomes$outcome == 0L)
set.seed(seed)
case_samples <- matrix(
  sample(case_indices, length(case_indices) * replicates, replace = TRUE),
  nrow = replicates,
  ncol = length(case_indices),
  byrow = TRUE
)
control_samples <- matrix(
  sample(control_indices, length(control_indices) * replicates, replace = TRUE),
  nrow = replicates,
  ncol = length(control_indices),
  byrow = TRUE
)
performance_distribution <- matrix(
  NA_real_, nrow = replicates, ncol = length(performance_names),
  dimnames = list(NULL, performance_names)
)
for (bootstrap_index in seq_len(replicates)) {
  sampled_indices <- c(case_samples[bootstrap_index, ], control_samples[bootstrap_index, ])
  sampled_ids <- participant_outcomes$participant_id[sampled_indices]
  performance_distribution[bootstrap_index, ] <- repeat_metrics_for_ids(sampled_ids, FALSE)[performance_names]
}

summarise_distribution <- function(names, distribution, point) {
  rows <- lapply(names, function(metric_name) {
    interval <- stats::quantile(distribution[, metric_name], probs = c(0.025, 0.975), names = FALSE, type = 7)
    data.frame(
      analysis = analysis,
      metric = metric_name,
      estimate = point[[metric_name]],
      lower_95 = interval[[1]],
      upper_95 = interval[[2]],
      bootstrap_replicates = nrow(distribution),
      seed = seed,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

performance_summary <- summarise_distribution(performance_names, performance_distribution, point_metrics)
write_tsv(
  data.frame(analysis = analysis, metric = names(point_metrics[performance_names]), estimate = as.numeric(point_metrics[performance_names])),
  file.path(output_dir, "performance_point_estimates.tsv")
)
write_tsv(performance_summary, file.path(output_dir, "performance_uncertainty.tsv"))
if (save_distribution) {
  write_tsv(
    data.frame(bootstrap_replicate = seq_len(replicates), performance_distribution, check.names = FALSE),
    file.path(output_dir, "performance_bootstrap_distribution.tsv.gz")
  )
}

calibration_summary <- data.frame()
if (include_calibration) {
  stop(
    paste(
      "Formal calibration uses NumPy's fixed participant bootstrap and is run by",
      "gendep model-diagnostics. Omit --include-calibration here."
    )
  )
}

checks <- data.frame(
  check = c(
    "participant_count", "outer_repeats", "rows_per_repeat",
    "one_prediction_per_participant_per_repeat", "performance_replicates",
    "performance_estimates_are_finite", "calibration_delegated_to_model_diagnostics"
  ),
  observed = c(
    nrow(participant_outcomes), length(unique(predictions$outer_repeat)), min(repeat_sizes),
    max(table(predictions$outer_repeat, predictions$participant_id)), nrow(performance_distribution),
    all(is.finite(performance_summary$estimate)) && all(is.finite(performance_summary$lower_95)) && all(is.finite(performance_summary$upper_95)),
    TRUE
  ),
  expected = c(expected_n, outer_repeats, expected_n, 1L, replicates, TRUE, TRUE),
  stringsAsFactors = FALSE
)
checks$result <- ifelse(as.character(checks$observed) == as.character(checks$expected), "PASS", "FAIL")
write_tsv(checks, file.path(output_dir, "bootstrap_validation.tsv"))
if (any(checks$result == "FAIL")) stop("Bootstrap validation failed.")

writeLines(
  c(
    "status=PASS",
    paste0("analysis=", analysis),
    paste0("participants=", expected_n),
    paste0("outer_repeats=", outer_repeats),
    paste0("bootstrap_replicates=", replicates),
    paste0("performance_seed=", seed),
    "calibration_included=FALSE",
    "calibration_command=gendep model-diagnostics",
    "models_refitted=FALSE",
    "hyperparameters_retuned=FALSE"
  ),
  file.path(output_dir, "_SUCCESS")
)
cat(sprintf("PAIRED_BOOTSTRAP=PASS analysis=%s replicates=%d output=%s\n", analysis, replicates, output_dir))
