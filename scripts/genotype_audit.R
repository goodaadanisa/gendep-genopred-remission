#!/usr/bin/env Rscript

# Purpose:
# Audit the supplied rsID-only genotype matrix and export the variant list and
# aggregate 0/1/2 frequencies required for reference reconstruction.
#
# Inputs:
# An RData file containing one two-dimensional genotype object.
#
# Outputs:
# Source variant IDs, aggregate genotype counts/frequencies, an audit summary
# and R session information. No participant-level genotype table is exported.
#
# Usage:
# Rscript scripts/genotype_audit.R \
#   --input /authorised/path/dataGen.RData \
#   --output-dir work/genotype/source_audit \
#   --object data.gen --expected-samples 430 --expected-variants 524876

options(stringsAsFactors = FALSE)

parse_arguments <- function(arguments) {
  defaults <- list(
    object = "data.gen",
    expected_samples = NA_integer_,
    expected_variants = NA_integer_,
    chunk_size = 4000L
  )

  values <- defaults
  index <- 1L

  while (index <= length(arguments)) {
    key <- arguments[[index]]

    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }

    if (index == length(arguments)) {
      stop("Missing value for argument: ", key)
    }

    value <- arguments[[index + 1L]]
    name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    values[[name]] <- value
    index <- index + 2L
  }

  if (is.null(values$input) || is.null(values$output_dir)) {
    stop("Required arguments: --input and --output-dir")
  }

  values$expected_samples <- suppressWarnings(as.integer(values$expected_samples))
  values$expected_variants <- suppressWarnings(as.integer(values$expected_variants))
  values$chunk_size <- suppressWarnings(as.integer(values$chunk_size))

  if (is.na(values$chunk_size) || values$chunk_size < 1L) {
    stop("--chunk-size must be a positive integer")
  }

  values
}

write_metric_table <- function(path, metrics) {
  output <- data.frame(
    metric = names(metrics),
    value = unname(unlist(metrics)),
    stringsAsFactors = FALSE
  )

  write.table(
    output,
    file = path,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
input_path <- normalizePath(arguments$input, mustWork = TRUE)
output_dir <- normalizePath(arguments$output_dir, mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_environment <- new.env(parent = baseenv())
loaded_objects <- load(input_path, envir = load_environment)

if (!arguments$object %in% loaded_objects) {
  stop(
    "Expected object '", arguments$object, "' was not found. Loaded objects: ",
    paste(loaded_objects, collapse = ", ")
  )
}

genotypes <- load_environment[[arguments$object]]

if (is.null(dim(genotypes)) || length(dim(genotypes)) != 2L) {
  stop("The genotype object is not two-dimensional")
}

sample_count <- nrow(genotypes)
variant_count <- ncol(genotypes)

if (!is.na(arguments$expected_samples) && sample_count != arguments$expected_samples) {
  stop(
    "Expected ", arguments$expected_samples, " samples, but found ", sample_count
  )
}

if (!is.na(arguments$expected_variants) && variant_count != arguments$expected_variants) {
  stop(
    "Expected ", arguments$expected_variants, " variants, but found ", variant_count
  )
}

variant_ids <- colnames(genotypes)

if (is.null(variant_ids)) {
  stop("The genotype object does not have variant column names")
}

if (anyNA(variant_ids) || any(trimws(variant_ids) == "")) {
  stop("The genotype object contains missing or empty variant IDs")
}

if (anyDuplicated(variant_ids) > 0L) {
  stop("The genotype object contains duplicate variant IDs")
}

row_names <- rownames(genotypes)
if (is.null(row_names)) {
  row_names <- as.character(seq_len(sample_count))
}

n_code_0 <- integer(variant_count)
n_code_1 <- integer(variant_count)
n_code_2 <- integer(variant_count)
n_missing <- integer(variant_count)
code2_frequency <- numeric(variant_count)

chunk_starts <- seq.int(1L, variant_count, by = arguments$chunk_size)

for (chunk_number in seq_along(chunk_starts)) {
  start_index <- chunk_starts[[chunk_number]]
  end_index <- min(start_index + arguments$chunk_size - 1L, variant_count)
  indexes <- start_index:end_index

  chunk <- genotypes[, indexes, drop = FALSE]

  if (is.data.frame(chunk)) {
    numeric_columns <- vapply(
      chunk,
      function(column) is.numeric(column) || is.integer(column),
      logical(1)
    )

    if (!all(numeric_columns)) {
      stop("A non-numeric genotype column was found in chunk ", chunk_number)
    }

    chunk <- as.matrix(chunk)
  }

  storage.mode(chunk) <- "double"

  if (!all(is.na(chunk) | chunk %in% c(0, 1, 2))) {
    stop("A genotype value outside 0, 1 and 2 was found in chunk ", chunk_number)
  }

  nonmissing <- colSums(!is.na(chunk))

  if (any(nonmissing == 0L)) {
    stop("A variant with no observed genotypes was found in chunk ", chunk_number)
  }

  n_code_0[indexes] <- colSums(chunk == 0, na.rm = TRUE)
  n_code_1[indexes] <- colSums(chunk == 1, na.rm = TRUE)
  n_code_2[indexes] <- colSums(chunk == 2, na.rm = TRUE)
  n_missing[indexes] <- colSums(is.na(chunk))
  code2_frequency[indexes] <- colSums(chunk, na.rm = TRUE) / (2 * nonmissing)

  message(
    "Processed genotype chunk ", chunk_number, " of ", length(chunk_starts)
  )

  rm(chunk)
  invisible(gc())
}

frequency_table <- data.frame(
  SNP = variant_ids,
  N_CODE_0 = n_code_0,
  N_CODE_1 = n_code_1,
  N_CODE_2 = n_code_2,
  N_MISSING = n_missing,
  GENDEP_CODE2_FREQ = code2_frequency,
  stringsAsFactors = FALSE
)

if (any(frequency_table$N_CODE_0 + frequency_table$N_CODE_1 +
        frequency_table$N_CODE_2 + frequency_table$N_MISSING != sample_count)) {
  stop("Aggregate genotype counts do not close to the sample count")
}

writeLines(variant_ids, file.path(output_dir, "source_variant_ids.txt"))

frequency_connection <- gzfile(
  file.path(output_dir, "source_genotype_frequencies.tsv.gz"),
  open = "wt"
)
write.table(
  frequency_table,
  file = frequency_connection,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
close(frequency_connection)

write.table(
  data.frame(
    row_order = seq_len(sample_count),
    source_row_name = row_names,
    stringsAsFactors = FALSE
  ),
  file = file.path(output_dir, "source_row_order.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

metrics <- c(
  genotype_objects_loaded = length(loaded_objects),
  samples = sample_count,
  variants = variant_count,
  unique_variant_ids = length(unique(variant_ids)),
  variants_with_missing_calls = sum(n_missing > 0L),
  total_missing_calls = sum(n_missing),
  observed_minimum_code2_frequency = min(code2_frequency),
  observed_maximum_code2_frequency = max(code2_frequency),
  source_row_names_are_sequence = identical(row_names, as.character(seq_len(sample_count)))
)

write_metric_table(file.path(output_dir, "source_genotype_audit.tsv"), metrics)
capture.output(sessionInfo(), file = file.path(output_dir, "R_session_info.txt"))

message("Genotype source audit completed successfully: ", output_dir)
