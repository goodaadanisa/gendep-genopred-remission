#!/usr/bin/env Rscript
# Purpose:
# Merge the eight participant-level GenoPred PRS profiles in validated PSAM order.
#
# Inputs:
# GenoPred output directory, validated PSAM and config/traits.tsv.
#
# Outputs:
# Controlled participant-level 430 x 8 PRS matrix and aggregate validation tables.
# The remission outcome is neither required nor read by this script.
#
# Usage:
# Rscript scripts/collect_prs.R --genopred-output work/prs/genopred/pipeline_output \
#   --psam work/genotype/target_validation/pfiles/GENDEP.chr1.psam \
#   --traits config/traits.tsv --output-dir work/prs/final --audit-dir work/prs/validation

suppressPackageStartupMessages(library(data.table))

parse_args <- function(values) {
  result <- list()
  index <- 1L
  while (index <= length(values)) {
    key <- values[[index]]
    if (!startsWith(key, "--")) stop(paste("Unexpected argument:", key))
    if (index == length(values)) stop(paste("Missing value for", key))
    result[[substring(key, 3L)]] <- values[[index + 1L]]
    index <- index + 2L
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
required <- c("genopred-output", "psam", "traits", "output-dir", "audit-dir")
missing_args <- setdiff(required, names(args))
if (length(missing_args) > 0L) {
  stop(paste("Missing required argument(s):", paste(missing_args, collapse = ", ")))
}

`%||%` <- function(left, right) if (is.null(left)) right else left
expected_participants <- as.integer(args[["expected-participants"]] %||% "430")
genopred_output <- normalizePath(args[["genopred-output"]], mustWork = TRUE)
psam_path <- normalizePath(args[["psam"]], mustWork = TRUE)
traits_path <- normalizePath(args[["traits"]], mustWork = TRUE)
output_dir <- args[["output-dir"]]
audit_dir <- args[["audit-dir"]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)

clean_names <- function(values) sub("^#", "", values)
read_table <- function(path) {
  result <- fread(path, data.table = TRUE, check.names = FALSE)
  setnames(result, clean_names(names(result)))
  result
}
find_column <- function(columns, requested) {
  matched <- columns[toupper(columns) == toupper(requested)]
  if (length(matched) != 1L) stop(paste("Expected exactly one", requested, "column"))
  matched[[1L]]
}

traits_table <- fread(traits_path)
if (!"trait" %in% names(traits_table)) stop("traits.tsv must contain a trait column")
traits <- as.character(traits_table$trait)
if (length(traits) != 8L || anyDuplicated(traits)) stop("Exactly eight unique traits are required")

psam <- read_table(psam_path)
fid_psam <- find_column(names(psam), "FID")
iid_psam <- find_column(names(psam), "IID")
expected <- data.table(FID = as.character(psam[[fid_psam]]), IID = as.character(psam[[iid_psam]]))
expected[, participant_key := paste(FID, IID, sep = "::")]
if (nrow(expected) != expected_participants) stop("Unexpected PSAM participant count")
if (anyDuplicated(expected$participant_key)) stop("PSAM participant identifiers are duplicated")

combined <- expected[, .(FID, IID)]
validation <- list()
summary <- list()

for (trait in traits) {
  profile_path <- file.path(
    genopred_output, "GENDEP", "pgs", "TRANS", "prscs", trait,
    paste0("GENDEP-", trait, "-TRANS.profiles")
  )
  if (!file.exists(profile_path) || file.info(profile_path)$size <= 0) {
    stop(paste("Missing participant-level PRS profile:", profile_path))
  }
  profile <- read_table(profile_path)
  fid <- find_column(names(profile), "FID")
  iid <- find_column(names(profile), "IID")
  score_columns <- setdiff(names(profile), c(fid, iid))
  if (length(score_columns) != 1L) {
    stop(paste(trait, "profile must contain exactly one score column"))
  }
  score_column <- score_columns[[1L]]
  profile[, participant_key := paste(as.character(get(fid)), as.character(get(iid)), sep = "::")]
  values <- suppressWarnings(as.numeric(profile[[score_column]]))
  match_index <- match(expected$participant_key, profile$participant_key)
  aligned <- values[match_index]
  pass <- (
    nrow(profile) == expected_participants &&
    uniqueN(profile$participant_key) == expected_participants &&
    !anyDuplicated(profile$participant_key) &&
    !anyNA(match_index) &&
    all(is.finite(aligned)) &&
    is.finite(sd(aligned)) && sd(aligned) > 0
  )
  validation[[trait]] <- data.table(
    trait = trait,
    profile_file = normalizePath(profile_path),
    score_column = score_column,
    participant_rows = nrow(profile),
    unique_participants = uniqueN(profile$participant_key),
    exact_psam_order = identical(profile$participant_key, expected$participant_key),
    missing_scores = sum(is.na(aligned)),
    nonfinite_scores = sum(!is.finite(aligned)),
    standard_deviation = sd(aligned),
    validation = ifelse(pass, "PASS", "FAIL")
  )
  summary[[trait]] <- data.table(
    trait = trait,
    participants = length(aligned),
    minimum = min(aligned),
    first_quartile = quantile(aligned, 0.25, names = FALSE),
    median = median(aligned),
    mean = mean(aligned),
    third_quartile = quantile(aligned, 0.75, names = FALSE),
    maximum = max(aligned),
    standard_deviation = sd(aligned)
  )
  combined[, (paste0("PRS_", trait)) := aligned]
}

validation_table <- rbindlist(validation, use.names = TRUE)
summary_table <- rbindlist(summary, use.names = TRUE)
score_columns <- paste0("PRS_", traits)
correlations <- cor(combined[, ..score_columns], method = "pearson")
correlation_long <- as.data.table(as.table(correlations))
setnames(correlation_long, c("trait_1", "trait_2", "pearson_r"))

matrix_path <- file.path(output_dir, "gendep_prscs_auto_8trait_prs.tsv.gz")
fwrite(combined, matrix_path, sep = "\t")
fwrite(validation_table, file.path(audit_dir, "participant_prs_validation.tsv"), sep = "\t")
fwrite(summary_table, file.path(audit_dir, "participant_prs_summary.tsv"), sep = "\t")
fwrite(correlation_long, file.path(audit_dir, "participant_prs_correlations.tsv"), sep = "\t")

overall_pass <- (
  nrow(combined) == expected_participants &&
  ncol(combined) == 10L &&
  all(validation_table$validation == "PASS") &&
  all(is.finite(as.matrix(combined[, ..score_columns])))
)
record <- c(
  paste0("traits=", paste(traits, collapse = ",")),
  paste0("participant_rows=", nrow(combined)),
  paste0("prs_columns=", length(score_columns)),
  paste0("final_matrix=", normalizePath(matrix_path, mustWork = FALSE)),
  paste0("validation=", ifelse(overall_pass, "PASS", "FAIL"))
)
writeLines(record, file.path(audit_dir, "participant_prs_collection_validation.txt"))
cat(paste(record, collapse = "\n"), "\n")
if (!overall_pass) stop("Participant-level PRS collection failed")
