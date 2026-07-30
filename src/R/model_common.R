# Shared input, configuration, encoding and metric helpers for Stage 5 models.

`%||%` <- function(value, fallback) {
  if (is.null(value) || length(value) == 0L || is.na(value) || !nzchar(value)) fallback else value
}

# Read a tab-separated workflow table and fail clearly when the file is absent.
read_tsv <- function(path) {
  if (!file.exists(path)) stop(paste("Missing file:", path))
  connection <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  read.delim(
    connection,
    header = TRUE,
    sep = "\t",
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

# Write a tab-separated output after creating its parent directory.
write_tsv <- function(data, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  connection <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(connection), add = TRUE)
  write.table(
    data,
    file = connection,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = TRUE,
    na = ""
  )
}

# Parse --key value command-line options used by the R entry points.
parse_cli <- function(args) {
  result <- list()
  index <- 1L
  while (index <= length(args)) {
    token <- args[[index]]
    if (!startsWith(token, "--")) stop(paste("Unexpected argument:", token))
    key <- sub("^--", "", token)
    if (index == length(args) || startsWith(args[[index + 1L]], "--")) {
      result[[key]] <- TRUE
      index <- index + 1L
    } else {
      result[[key]] <- args[[index + 1L]]
      index <- index + 2L
    }
  }
  result
}

# Return a required CLI option or stop with a clear message.
required_cli <- function(options, key) {
  value <- options[[key]]
  if (is.null(value) || identical(value, TRUE) || !nzchar(value)) {
    stop(paste("Missing required option --", key, sep = ""))
  }
  value
}

# Load the fixed modelling parameters as a named list.
read_parameters <- function(path) {
  table <- read_tsv(path)
  if (!all(c("parameter", "value") %in% names(table))) {
    stop("Parameter configuration must contain parameter and value columns.")
  }
  values <- as.list(setNames(as.character(table$value), as.character(table$parameter)))
  values
}

# Load one analysis row from the configuration table.
read_analysis_spec <- function(path, analysis) {
  table <- read_tsv(path)
  required <- c(
    "analysis", "algorithm", "cohort_name", "expected_participants",
    "comparator_predictor_set", "combined_predictor_set",
    "comparator_label", "combined_label", "prs_version"
  )
  missing <- setdiff(required, names(table))
  if (length(missing) > 0L) stop(paste("Analysis configuration missing:", paste(missing, collapse = ", ")))
  row <- table[table$analysis == analysis, , drop = FALSE]
  if (nrow(row) != 1L) stop(paste("Expected one analysis configuration for:", analysis))
  as.list(row[1, , drop = FALSE])
}

# Read and validate one integer modelling parameter.
as_integer_parameter <- function(parameters, name) {
  value <- suppressWarnings(as.integer(parameters[[name]]))
  if (is.na(value)) stop(paste("Parameter is not an integer:", name))
  value
}

# Read and validate one numeric modelling parameter.
as_numeric_parameter <- function(parameters, name) {
  value <- suppressWarnings(as.numeric(parameters[[name]]))
  if (is.na(value) || !is.finite(value)) stop(paste("Parameter is not numeric:", name))
  value
}

# Parse a comma-separated numeric parameter vector.
as_numeric_list <- function(value) {
  result <- suppressWarnings(as.numeric(strsplit(as.character(value), ",", fixed = TRUE)[[1]]))
  if (anyNA(result) || any(!is.finite(result))) stop(paste("Invalid numeric list:", value))
  result
}

# Convert a predictor column to numeric and reject non-finite values.
numeric_vector <- function(values, variable) {
  result <- suppressWarnings(as.numeric(values))
  if (anyNA(result) || any(!is.finite(result))) {
    stop(paste("Non-numeric or non-finite values found for:", variable))
  }
  result
}

# Load the comparator and combined predictor sets in their fixed order.
read_predictor_set <- function(policy_dir, set_name, combined_set_name = NULL) {
  path <- file.path(policy_dir, "model_predictor_sets.tsv")
  table <- read_tsv(path)
  required <- c("model", "predictor_order", "variable")
  if (!all(required %in% names(table))) stop("model_predictor_sets.tsv has an invalid schema.")

  if (identical(set_name, "__remove_prs__")) {
    if (is.null(combined_set_name)) stop("combined_set_name is required for __remove_prs__.")
    rows <- table[table$model == combined_set_name, , drop = FALSE]
    rows <- rows[!grepl("^PRS_", rows$variable), , drop = FALSE]
  } else {
    rows <- table[table$model == set_name, , drop = FALSE]
  }

  if (nrow(rows) == 0L) stop(paste("No predictors found for set:", set_name))
  rows <- rows[order(as.integer(rows$predictor_order)), , drop = FALSE]
  predictors <- as.character(rows$variable)
  if (anyDuplicated(predictors)) stop(paste("Duplicate predictors in set:", set_name))
  predictors
}

# Load training-independent binary encoding rules.
read_binary_mapping <- function(policy_dir) {
  path <- file.path(policy_dir, "binary_source_code_mappings.tsv")
  table <- read_tsv(path)
  required <- c("variable", "source_level_mapped_to_zero", "source_level_mapped_to_one")
  if (!all(required %in% names(table))) stop("Binary mapping file has an invalid schema.")
  table
}

# Apply the fixed semantic encoding policy without full-cohort scaling.
encode_predictor_matrix <- function(data, predictors, binary_mapping) {
  missing <- setdiff(predictors, names(data))
  if (length(missing) > 0L) stop(paste("Missing predictors:", paste(missing, collapse = ", ")))

  result <- matrix(
    NA_real_,
    nrow = nrow(data),
    ncol = length(predictors),
    dimnames = list(NULL, predictors)
  )
  mapped <- as.character(binary_mapping$variable)

  for (column_index in seq_along(predictors)) {
    variable <- predictors[[column_index]]
    values <- numeric_vector(data[[variable]], variable)
    if (variable %in% mapped) {
      row <- binary_mapping[binary_mapping$variable == variable, , drop = FALSE]
      if (nrow(row) != 1L) stop(paste("Expected one binary mapping for:", variable))
      zero <- as.numeric(row$source_level_mapped_to_zero[[1]])
      one <- as.numeric(row$source_level_mapped_to_one[[1]])
      encoded <- rep(NA_real_, length(values))
      encoded[abs(values - zero) <= 1e-12] <- 0
      encoded[abs(values - one) <= 1e-12] <- 1
      if (anyNA(encoded)) stop(paste("Unexpected binary source levels for:", variable))
      values <- encoded
    }
    result[, column_index] <- values
  }

  if (anyNA(result) || any(!is.finite(result))) stop("Encoded predictor matrix contains invalid values.")
  result
}

# Calculate ROC AUC from ranks, including average ranks for ties.
calculate_auc <- function(outcome, prediction) {
  outcome <- as.integer(outcome)
  prediction <- as.numeric(prediction)
  positive_count <- sum(outcome == 1L)
  negative_count <- sum(outcome == 0L)
  if (positive_count == 0L || negative_count == 0L) stop("AUC requires both outcome classes.")
  ranks <- rank(prediction, ties.method = "average")
  (sum(ranks[outcome == 1L]) - positive_count * (positive_count + 1) / 2) /
    (positive_count * negative_count)
}

# Calculate the Brier score for held-out probabilities.
calculate_brier <- function(outcome, prediction) {
  mean((as.numeric(prediction) - as.numeric(outcome))^2)
}

# Calculate bounded binary logarithmic loss.
calculate_log_loss <- function(outcome, prediction, epsilon = 1e-15) {
  bounded <- pmin(pmax(as.numeric(prediction), epsilon), 1 - epsilon)
  outcome <- as.numeric(outcome)
  -mean(outcome * log(bounded) + (1 - outcome) * log(1 - bounded))
}

# Convert bounded probabilities to logits for calibration models.
bounded_logit <- function(prediction, epsilon = 1e-15) {
  bounded <- pmin(pmax(as.numeric(prediction), epsilon), 1 - epsilon)
  log(bounded / (1 - bounded))
}

# Estimate calibration-in-the-large and calibration slope from held-out predictions.
calculate_calibration <- function(outcome, prediction, epsilon = 1e-15) {
  outcome <- as.numeric(outcome)
  logits <- bounded_logit(prediction, epsilon)
  intercept_fit <- suppressWarnings(glm(outcome ~ 1, family = binomial(), offset = logits))
  slope_fit <- suppressWarnings(glm(outcome ~ logits, family = binomial()))
  c(
    calibration_intercept = unname(coef(intercept_fit)[[1]]),
    calibration_slope = unname(coef(slope_fit)[[2]])
  )
}

# Create one standardised metric record for a fitted model.
metric_row <- function(model, outcome, prediction, epsilon = 1e-15) {
  calibration <- calculate_calibration(outcome, prediction, epsilon)
  data.frame(
    model = model,
    participants = length(outcome),
    events = sum(as.integer(outcome) == 1L),
    non_events = sum(as.integer(outcome) == 0L),
    AUC = calculate_auc(outcome, prediction),
    Brier = calculate_brier(outcome, prediction),
    log_loss = calculate_log_loss(outcome, prediction, epsilon),
    calibration_intercept = calibration[["calibration_intercept"]],
    calibration_slope = calibration[["calibration_slope"]],
    stringsAsFactors = FALSE
  )
}

# Resolve and validate the requested outer repetition and fold.
resolve_repeat_fold <- function(options, outer_repeats, outer_folds) {
  repeat_value <- options[["repeat"]]
  fold_value <- options[["fold"]]

  if (is.null(repeat_value) && is.null(fold_value)) {
    task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
    if (!nzchar(task)) stop("Provide --repeat and --fold, or set SLURM_ARRAY_TASK_ID.")
    task_id <- suppressWarnings(as.integer(task))
    if (is.na(task_id) || task_id < 1L || task_id > outer_repeats * outer_folds) {
      stop("SLURM_ARRAY_TASK_ID is outside the configured split range.")
    }
    repeat_value <- ((task_id - 1L) %/% outer_folds) + 1L
    fold_value <- ((task_id - 1L) %% outer_folds) + 1L
  }

  outer_repeat <- suppressWarnings(as.integer(repeat_value))
  outer_fold <- suppressWarnings(as.integer(fold_value))
  if (is.na(outer_repeat) || !outer_repeat %in% seq_len(outer_repeats)) stop("Invalid --repeat.")
  if (is.na(outer_fold) || !outer_fold %in% seq_len(outer_folds)) stop("Invalid --fold.")
  c(outer_repeat = outer_repeat, outer_fold = outer_fold)
}

# Assemble one outer split while enforcing sample, predictor and fold invariants.
prepare_split_context <- function(
  analysis_path,
  predictor_policy_dir,
  resampling_dir,
  analysis_spec,
  parameters,
  outer_repeat,
  outer_fold
) {
  analysis_data <- read_tsv(analysis_path)
  if (!"Row.names" %in% names(analysis_data)) stop("Analysis data must contain Row.names.")
  analysis_data$participant_id <- as.character(analysis_data$Row.names)
  if (anyDuplicated(analysis_data$participant_id)) stop("Analysis participant IDs are duplicated.")

  expected_n <- as.integer(analysis_spec$expected_participants)
  if (nrow(analysis_data) != expected_n) {
    stop(paste("Expected", expected_n, "participants; found", nrow(analysis_data)))
  }

  comparator_predictors <- read_predictor_set(
    predictor_policy_dir,
    as.character(analysis_spec$comparator_predictor_set),
    as.character(analysis_spec$combined_predictor_set)
  )
  combined_predictors <- read_predictor_set(
    predictor_policy_dir,
    as.character(analysis_spec$combined_predictor_set)
  )
  if (!all(comparator_predictors %in% combined_predictors)) {
    stop("Comparator predictors are not a subset of combined predictors.")
  }
  prs_difference <- setdiff(combined_predictors, comparator_predictors)
  if (length(prs_difference) != 8L || !all(grepl("^PRS_", prs_difference))) {
    stop("The combined-minus-comparator difference must be exactly eight PRS.")
  }

  forbidden <- c("hdremit.all", "subjectid", "Row.names", "bloodsampleid.x", "participant_id", "model_based_EUR_keep")
  if (any(forbidden %in% combined_predictors)) stop("Forbidden variable found in predictor set.")

  required_columns <- unique(c("participant_id", "hdremit.all", "drug", combined_predictors))
  missing <- setdiff(required_columns, names(analysis_data))
  if (length(missing) > 0L) stop(paste("Analysis data missing:", paste(missing, collapse = ", ")))
  if (anyNA(analysis_data[, unique(c("hdremit.all", combined_predictors)), drop = FALSE])) {
    stop("Outcome or predictors contain missing values.")
  }

  cohort <- as.character(analysis_spec$cohort_name)
  outer_path <- file.path(resampling_dir, paste0(cohort, "_outer_folds.tsv.gz"))
  inner_path <- file.path(resampling_dir, paste0(cohort, "_inner_folds.tsv.gz"))
  outer_manifest <- read_tsv(outer_path)
  inner_manifest <- read_tsv(inner_path)
  outer_manifest$participant_id <- as.character(outer_manifest$participant_id)
  inner_manifest$participant_id <- as.character(inner_manifest$participant_id)

  selected_outer <- outer_manifest[outer_manifest$outer_repeat == outer_repeat, , drop = FALSE]
  if (nrow(selected_outer) != expected_n) stop("Outer manifest does not contain the expected cohort.")
  if (!setequal(selected_outer$participant_id, analysis_data$participant_id)) stop("Outer manifest IDs do not match analysis data.")

  test_ids <- selected_outer$participant_id[selected_outer$outer_fold == outer_fold]
  training_ids <- selected_outer$participant_id[selected_outer$outer_fold != outer_fold]
  if (length(intersect(test_ids, training_ids)) != 0L) stop("Outer train/test leakage detected.")
  if (!setequal(c(test_ids, training_ids), analysis_data$participant_id)) stop("Outer split does not cover the cohort.")

  selected_inner <- inner_manifest[
    inner_manifest$outer_repeat == outer_repeat &
      inner_manifest$held_out_outer_fold == outer_fold,
    ,
    drop = FALSE
  ]
  if (!setequal(selected_inner$participant_id, training_ids)) stop("Inner manifest does not match outer training IDs.")
  configured_inner <- as_integer_parameter(parameters, "inner_folds")
  if (!identical(sort(unique(as.integer(selected_inner$inner_fold))), seq_len(configured_inner))) {
    stop("Inner fold labels do not match the configured design.")
  }

  training_indices <- match(training_ids, analysis_data$participant_id)
  test_indices <- match(test_ids, analysis_data$participant_id)
  inner_fold_id <- selected_inner$inner_fold[match(training_ids, selected_inner$participant_id)]
  if (anyNA(c(training_indices, test_indices, inner_fold_id))) stop("Participant matching failed.")

  list(
    analysis_data = analysis_data,
    comparator_predictors = comparator_predictors,
    combined_predictors = combined_predictors,
    binary_mapping = read_binary_mapping(predictor_policy_dir),
    training_indices = training_indices,
    test_indices = test_indices,
    inner_fold_id = as.integer(inner_fold_id),
    training_ids = training_ids,
    test_ids = test_ids,
    outer_repeat = outer_repeat,
    outer_fold = outer_fold,
    forced_covariates = c("drug", paste0("ancestry_PC", 1:6)),
    prs_columns = prs_difference
  )
}

normalise_package_version <- function(value) {
  gsub("-", ".", trimws(as.character(value)), fixed = TRUE)
}
