#!/usr/bin/env Rscript
# Run one configured outer split for Elastic Net or Random Forest.
#
# Purpose:
#   Fit paired comparator and comparator-plus-PRS models using fixed nested-CV
#   assignments. The same command supports every Stage 5 analysis specification.
#
# Inputs:
#   Controlled analysis data, release-safe predictor policy, deterministic outer
#   and inner fold manifests, and the fixed model configuration tables.
#
# Outputs:
#   Held-out predictions, split metrics, tuning records, coefficients where
#   applicable, warnings, validation checks and a completion marker.
#
# Example:
#   Rscript scripts/run_models.R --analysis primary \
#     --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
#     --predictor-policy-dir work/analysis/predictor_policy \
#     --resampling-dir work/analysis/resampling \
#     --output-dir work/modelling --repeat 1 --fold 1

args <- commandArgs(trailingOnly = TRUE)
script_arg <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_arg, value = TRUE)
script_path <- if (length(file_arg) == 1L) sub("^--file=", "", file_arg) else "scripts/run_models.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

source(file.path(root, "src", "R", "model_common.R"))
source(file.path(root, "src", "R", "elastic_net_engine.R"))
source(file.path(root, "src", "R", "random_forest_engine.R"))

options <- parse_cli(args)
analysis <- required_cli(options, "analysis")
analysis_base <- required_cli(options, "analysis-base")
predictor_policy_dir <- required_cli(options, "predictor-policy-dir")
resampling_dir <- required_cli(options, "resampling-dir")
output_root <- required_cli(options, "output-dir")
analysis_config <- options[["analysis-config"]] %||% file.path(root, "config", "model_analyses.tsv")
parameter_config <- options[["parameter-config"]] %||% file.path(root, "config", "model_parameters.tsv")

analysis_spec <- read_analysis_spec(analysis_config, analysis)
parameters <- read_parameters(parameter_config)
outer_repeats <- as_integer_parameter(parameters, "outer_repeats")
outer_folds <- as_integer_parameter(parameters, "outer_folds")
split <- resolve_repeat_fold(options, outer_repeats, outer_folds)
outer_repeat <- unname(split[["outer_repeat"]])
outer_fold <- unname(split[["outer_fold"]])

split_dir <- file.path(
  output_root,
  analysis,
  sprintf("repeat_%02d", outer_repeat),
  sprintf("fold_%02d", outer_fold)
)
if (dir.exists(split_dir) && file.exists(file.path(split_dir, "_SUCCESS")) && !isTRUE(options[["overwrite"]])) {
  stop(paste("Completed split already exists. Use --overwrite to replace:", split_dir))
}
dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)

context <- prepare_split_context(
  analysis_base,
  predictor_policy_dir,
  resampling_dir,
  analysis_spec,
  parameters,
  outer_repeat,
  outer_fold
)

algorithm <- as.character(analysis_spec$algorithm)
comparator_label <- as.character(analysis_spec$comparator_label)
combined_label <- as.character(analysis_spec$combined_label)

if (algorithm == "elastic_net") {
  if (!requireNamespace("glmnet", quietly = TRUE)) stop("R package 'glmnet' is required.")
  engine_result <- run_elastic_net_split(context, parameters, comparator_label, combined_label)
  package_name <- "glmnet"
} else if (algorithm == "random_forest") {
  if (!requireNamespace("randomForest", quietly = TRUE)) stop("R package 'randomForest' is required.")
  engine_result <- run_random_forest_split(context, parameters, comparator_label, combined_label)
  package_name <- "randomForest"
} else {
  stop(paste("Unsupported algorithm:", algorithm))
}

comparator_prediction <- engine_result$comparator$predictions
combined_prediction <- engine_result$combined$predictions
if (!identical(comparator_prediction$participant_id, combined_prediction$participant_id)) {
  stop("Comparator and combined participant predictions are not aligned.")
}
if (!identical(comparator_prediction$outcome, combined_prediction$outcome)) {
  stop("Comparator and combined outcomes are not aligned.")
}

predictions <- data.frame(
  analysis = analysis,
  algorithm = algorithm,
  outer_repeat = outer_repeat,
  outer_fold = outer_fold,
  participant_id = comparator_prediction$participant_id,
  outcome = comparator_prediction$outcome,
  drug = comparator_prediction$drug,
  comparator_probability = comparator_prediction$prediction,
  combined_probability = combined_prediction$prediction,
  combined_minus_comparator_probability = combined_prediction$prediction - comparator_prediction$prediction,
  stringsAsFactors = FALSE
)

metrics <- rbind(engine_result$comparator$metrics, engine_result$combined$metrics)
metrics$analysis <- analysis
metrics$algorithm <- algorithm
metrics$outer_repeat <- outer_repeat
metrics$outer_fold <- outer_fold
metrics <- metrics[, c("analysis", "algorithm", "outer_repeat", "outer_fold", setdiff(names(metrics), c("analysis", "algorithm", "outer_repeat", "outer_fold"))), drop = FALSE]

comparator_metrics <- metrics[metrics$model == comparator_label, , drop = FALSE]
combined_metrics <- metrics[metrics$model == combined_label, , drop = FALSE]
comparison <- data.frame(
  analysis = analysis,
  algorithm = algorithm,
  outer_repeat = outer_repeat,
  outer_fold = outer_fold,
  comparator_AUC = comparator_metrics$AUC,
  combined_AUC = combined_metrics$AUC,
  combined_minus_comparator_AUC = combined_metrics$AUC - comparator_metrics$AUC,
  comparator_Brier = comparator_metrics$Brier,
  combined_Brier = combined_metrics$Brier,
  comparator_minus_combined_Brier = comparator_metrics$Brier - combined_metrics$Brier,
  comparator_log_loss = comparator_metrics$log_loss,
  combined_log_loss = combined_metrics$log_loss,
  comparator_minus_combined_log_loss = comparator_metrics$log_loss - combined_metrics$log_loss,
  comparator_calibration_intercept = comparator_metrics$calibration_intercept,
  combined_calibration_intercept = combined_metrics$calibration_intercept,
  combined_minus_comparator_calibration_intercept = combined_metrics$calibration_intercept - comparator_metrics$calibration_intercept,
  comparator_calibration_slope = comparator_metrics$calibration_slope,
  combined_calibration_slope = combined_metrics$calibration_slope,
  combined_minus_comparator_calibration_slope = combined_metrics$calibration_slope - comparator_metrics$calibration_slope,
  stringsAsFactors = FALSE
)

tuning <- rbind(engine_result$comparator$tuning, engine_result$combined$tuning)
tuning$analysis <- analysis
tuning$outer_repeat <- outer_repeat
tuning$outer_fold <- outer_fold
tuning <- tuning[, c("analysis", "outer_repeat", "outer_fold", setdiff(names(tuning), c("analysis", "outer_repeat", "outer_fold"))), drop = FALSE]

coefficients <- rbind(engine_result$comparator$coefficients, engine_result$combined$coefficients)
if (nrow(coefficients) > 0L) {
  coefficients$analysis <- analysis
  coefficients$outer_repeat <- outer_repeat
  coefficients$outer_fold <- outer_fold
  coefficients <- coefficients[, c("analysis", "outer_repeat", "outer_fold", setdiff(names(coefficients), c("analysis", "outer_repeat", "outer_fold"))), drop = FALSE]
}

warnings <- rbind(engine_result$comparator$warnings, engine_result$combined$warnings)
if (nrow(warnings) > 0L) {
  warnings$analysis <- analysis
  warnings$outer_repeat <- outer_repeat
  warnings$outer_fold <- outer_fold
  warnings <- warnings[, c("analysis", "outer_repeat", "outer_fold", setdiff(names(warnings), c("analysis", "outer_repeat", "outer_fold"))), drop = FALSE]
}

expected_package_version <- parameters[[paste0("expected_", package_name, "_version")]] %||% ""
observed_package_version <- as.character(utils::packageVersion(package_name))
checks <- data.frame(
  check = c(
    "analysis_participant_count",
    "outer_training_and_test_disjoint",
    "outer_split_covers_cohort",
    "inner_assignment_matches_outer_training",
    "comparator_is_combined_subset",
    "combined_adds_exactly_eight_prs",
    "prediction_rows_equal_test_rows",
    "predictions_are_finite_probabilities",
    "outcomes_match",
    "software_version_matches_frozen_runtime"
  ),
  observed = c(
    nrow(context$analysis_data),
    length(intersect(context$training_ids, context$test_ids)) == 0L,
    length(unique(c(context$training_ids, context$test_ids))),
    length(context$inner_fold_id),
    all(context$comparator_predictors %in% context$combined_predictors),
    length(context$prs_columns),
    nrow(predictions),
    all(is.finite(predictions$comparator_probability)) && all(is.finite(predictions$combined_probability)) &&
      all(predictions$comparator_probability >= 0 & predictions$comparator_probability <= 1) &&
      all(predictions$combined_probability >= 0 & predictions$combined_probability <= 1),
    identical(comparator_prediction$outcome, combined_prediction$outcome),
    if (nzchar(expected_package_version)) {
      normalise_package_version(observed_package_version) == normalise_package_version(expected_package_version) || isTRUE(options[["allow-version-drift"]])
    } else TRUE
  ),
  expected = c(
    as.integer(analysis_spec$expected_participants),
    TRUE,
    as.integer(analysis_spec$expected_participants),
    length(context$training_ids),
    TRUE,
    8L,
    length(context$test_ids),
    TRUE,
    TRUE,
    TRUE
  ),
  stringsAsFactors = FALSE
)
checks$result <- ifelse(as.character(checks$observed) == as.character(checks$expected), "PASS", "FAIL")

write_tsv(predictions, file.path(split_dir, "predictions.tsv.gz"))
write_tsv(metrics, file.path(split_dir, "model_metrics.tsv"))
write_tsv(comparison, file.path(split_dir, "comparison.tsv"))
write_tsv(tuning, file.path(split_dir, "tuning.tsv.gz"))
write_tsv(coefficients, file.path(split_dir, "coefficients.tsv.gz"))
write_tsv(warnings, file.path(split_dir, "warnings.tsv"))
write_tsv(checks, file.path(split_dir, "split_validation.tsv"))

runtime <- data.frame(
  field = c(
    "analysis", "algorithm", "outer_repeat", "outer_fold", "R_version",
    "package", "package_version", "expected_package_version", "prs_version"
  ),
  value = c(
    analysis, algorithm, outer_repeat, outer_fold, R.version.string,
    package_name, observed_package_version, expected_package_version,
    as.character(analysis_spec$prs_version)
  ),
  stringsAsFactors = FALSE
)
write_tsv(runtime, file.path(split_dir, "runtime.tsv"))

if (any(checks$result == "FAIL")) {
  stop(paste("Split validation failed:", paste(checks$check[checks$result == "FAIL"], collapse = ", ")))
}
writeLines(
  c(
    "status=PASS",
    paste0("analysis=", analysis),
    paste0("algorithm=", algorithm),
    paste0("outer_repeat=", outer_repeat),
    paste0("outer_fold=", outer_fold),
    paste0("package_version=", observed_package_version)
  ),
  file.path(split_dir, "_SUCCESS")
)
cat(sprintf("MODEL_SPLIT=PASS analysis=%s repeat=%d fold=%d output=%s\n", analysis, outer_repeat, outer_fold, split_dir))
