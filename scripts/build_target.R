#!/usr/bin/env Rscript

# Purpose:
# Convert the source 0/1/2 genotype matrix to chromosome-specific VCF files
# using the fixed target-orientation table.
#
# Inputs:
# The authorised genotype RData file, final target markers and a validated
# sample-ID list in exact source-row order.
#
# Outputs:
# One VCF and variant manifest per autosome, plus build summaries. PLINK2 import
# and exact round-trip validation are performed by validate_target.py.
#
# Usage:
# Rscript scripts/build_target.R \
#   --input /authorised/path/dataGen.RData \
#   --markers work/genotype/reconstruction/final_target_markers.tsv.gz \
#   --samples /authorised/path/validated_sample_ids.txt \
#   --output-dir work/genotype/target

options(stringsAsFactors = FALSE)

parse_arguments <- function(arguments) {
  values <- list(
    object = "data.gen",
    expected_samples = NA_integer_,
    expected_variants = NA_integer_
  )

  index <- 1L
  while (index <= length(arguments)) {
    key <- arguments[[index]]
    if (!startsWith(key, "--")) {
      stop("Unexpected positional argument: ", key)
    }
    if (index == length(arguments)) {
      stop("Missing value for argument: ", key)
    }
    name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
    values[[name]] <- arguments[[index + 1L]]
    index <- index + 2L
  }

  required <- c("input", "markers", "samples", "output_dir")
  missing <- required[!vapply(required, function(name) !is.null(values[[name]]), logical(1))]
  if (length(missing) > 0L) {
    stop("Missing required arguments: ", paste(paste0("--", gsub("_", "-", missing)), collapse = ", "))
  }

  values$expected_samples <- suppressWarnings(as.integer(values$expected_samples))
  values$expected_variants <- suppressWarnings(as.integer(values$expected_variants))
  values
}

read_marker_table <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(connection))
  read.delim(
    connection,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE
  )
}

parse_logical <- function(value) {
  toupper(trimws(as.character(value))) %in% c("TRUE", "T", "1")
}

arguments <- parse_arguments(commandArgs(trailingOnly = TRUE))
input_path <- normalizePath(arguments$input, mustWork = TRUE)
marker_path <- normalizePath(arguments$markers, mustWork = TRUE)
sample_path <- normalizePath(arguments$samples, mustWork = TRUE)
output_dir <- normalizePath(arguments$output_dir, mustWork = FALSE)

vcf_dir <- file.path(output_dir, "vcf")
manifest_dir <- file.path(output_dir, "manifests")
summary_dir <- file.path(output_dir, "summaries")
for (directory in c(output_dir, vcf_dir, manifest_dir, summary_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

load_environment <- new.env(parent = baseenv())
loaded_objects <- load(input_path, envir = load_environment)
if (!arguments$object %in% loaded_objects) {
  stop("Expected genotype object '", arguments$object, "' was not found")
}

genotypes <- load_environment[[arguments$object]]
if (is.null(dim(genotypes)) || length(dim(genotypes)) != 2L) {
  stop("The genotype object is not two-dimensional")
}

sample_count <- nrow(genotypes)
source_variant_count <- ncol(genotypes)
if (!is.na(arguments$expected_samples) && sample_count != arguments$expected_samples) {
  stop("Expected ", arguments$expected_samples, " samples, but found ", sample_count)
}

source_ids <- colnames(genotypes)
if (is.null(source_ids) || anyDuplicated(source_ids) > 0L) {
  stop("Source genotype variant IDs are missing or duplicated")
}

samples <- trimws(readLines(sample_path, warn = FALSE))
samples <- samples[nzchar(samples)]
if (length(samples) != sample_count) {
  stop("Sample list contains ", length(samples), " IDs, but the genotype matrix has ", sample_count, " rows")
}
if (anyDuplicated(samples) > 0L) {
  stop("Sample IDs are not unique")
}
if (any(!grepl("^[A-Za-z0-9_.-]+$", samples))) {
  stop("A sample ID contains a VCF-unsafe character")
}

markers <- read_marker_table(marker_path)
required_columns <- c(
  "SNP", "GENOPRED_CHR", "GENOPRED_BP", "GENOPRED_REF", "GENOPRED_ALT",
  "GENDEP_CODE2_FREQ", "EUR_ALT_FREQ", "DIFFERENCE_TO_EUR_MAJOR_FREQ",
  "CODE2_SIDE", "CODE0_SIDE", "CODE2_ALLELE", "CODE0_ALLELE",
  "ORIENTATION_DECISION"
)
missing_columns <- setdiff(required_columns, names(markers))
if (length(missing_columns) > 0L) {
  stop("Marker table is missing columns: ", paste(missing_columns, collapse = ", "))
}

if ("KEEP_FOR_TARGET" %in% names(markers)) {
  markers <- markers[parse_logical(markers$KEEP_FOR_TARGET), , drop = FALSE]
} else {
  if (is.na(arguments$expected_variants)) {
    stop(
      "Marker table does not contain KEEP_FOR_TARGET; --expected-variants is ",
      "required when supplying an already-filtered final marker table"
    )
  }

  if (nrow(markers) != arguments$expected_variants) {
    stop(
      "Marker table does not contain KEEP_FOR_TARGET and has ", nrow(markers),
      " rows, but --expected-variants is ", arguments$expected_variants
    )
  }

  markers$KEEP_FOR_TARGET <- TRUE
}

markers <- markers[, c(required_columns, "KEEP_FOR_TARGET"), drop = FALSE]
markers$GENOPRED_CHR <- as.integer(markers$GENOPRED_CHR)
markers$GENOPRED_BP <- as.integer(markers$GENOPRED_BP)
markers$GENOPRED_REF <- toupper(markers$GENOPRED_REF)
markers$GENOPRED_ALT <- toupper(markers$GENOPRED_ALT)
markers$CODE2_SIDE <- toupper(markers$CODE2_SIDE)
markers$CODE0_SIDE <- toupper(markers$CODE0_SIDE)
markers$CODE2_ALLELE <- toupper(markers$CODE2_ALLELE)
markers$CODE0_ALLELE <- toupper(markers$CODE0_ALLELE)

if (!is.na(arguments$expected_variants) && nrow(markers) != arguments$expected_variants) {
  stop("Expected ", arguments$expected_variants, " target markers, but found ", nrow(markers))
}
if (anyNA(markers$GENOPRED_CHR) || anyNA(markers$GENOPRED_BP)) {
  stop("Target markers contain missing coordinates")
}
if (any(!markers$GENOPRED_CHR %in% 1:22) || any(markers$GENOPRED_BP <= 0L)) {
  stop("Target markers contain invalid autosomal coordinates")
}
valid_bases <- c("A", "C", "G", "T")
if (any(!markers$GENOPRED_REF %in% valid_bases) || any(!markers$GENOPRED_ALT %in% valid_bases)) {
  stop("Target markers contain non-SNV alleles")
}
if (any(markers$GENOPRED_REF == markers$GENOPRED_ALT)) {
  stop("Target markers contain identical REF and ALT alleles")
}
if (any(!markers$CODE2_SIDE %in% c("REF", "ALT")) || any(!markers$CODE0_SIDE %in% c("REF", "ALT"))) {
  stop("Target markers contain invalid code-to-allele sides")
}
if (any(markers$CODE2_SIDE == markers$CODE0_SIDE)) {
  stop("A target marker assigns code 0 and code 2 to the same allele")
}
if (anyDuplicated(markers$SNP) > 0L) {
  stop("Target marker rsIDs are duplicated")
}
variant_key <- paste(markers$GENOPRED_CHR, markers$GENOPRED_BP, markers$GENOPRED_REF, markers$GENOPRED_ALT, sep = ":")
if (anyDuplicated(variant_key) > 0L) {
  stop("Target marker chromosome-position-REF-ALT keys are duplicated")
}

expected_code2 <- ifelse(markers$CODE2_SIDE == "REF", markers$GENOPRED_REF, markers$GENOPRED_ALT)
expected_code0 <- ifelse(markers$CODE0_SIDE == "REF", markers$GENOPRED_REF, markers$GENOPRED_ALT)
if (any(markers$CODE2_ALLELE != expected_code2) || any(markers$CODE0_ALLELE != expected_code0)) {
  stop("Recorded code alleles disagree with the REF/ALT side")
}

source_indexes <- match(markers$SNP, source_ids)
if (anyNA(source_indexes)) {
  stop("Target markers are absent from the source matrix: ", paste(head(markers$SNP[is.na(source_indexes)], 10L), collapse = ", "))
}

writeLines(samples, file.path(output_dir, "validated_sample_ids.txt"))
summary_rows <- vector("list", 22L)

for (chromosome in 1:22) {
  chromosome_markers <- markers[markers$GENOPRED_CHR == chromosome, , drop = FALSE]
  chromosome_markers <- chromosome_markers[order(chromosome_markers$GENOPRED_BP, chromosome_markers$SNP), , drop = FALSE]
  chromosome_indexes <- match(chromosome_markers$SNP, source_ids)

  vcf_path <- file.path(vcf_dir, paste0("GENDEP.chr", chromosome, ".vcf"))
  temporary_path <- paste0(vcf_path, ".tmp")
  manifest_path <- file.path(manifest_dir, paste0("GENDEP.chr", chromosome, ".variant_manifest.tsv"))

  expected_alt_count <- integer(nrow(chromosome_markers))
  observed_alt_frequency <- numeric(nrow(chromosome_markers))
  expected_alt_frequency <- ifelse(
    chromosome_markers$CODE2_SIDE == "ALT",
    chromosome_markers$GENDEP_CODE2_FREQ,
    1 - chromosome_markers$GENDEP_CODE2_FREQ
  )

  connection <- file(temporary_path, open = "wt")
  writeLines(
    c(
      "##fileformat=VCFv4.2",
      "##source=gendep-genopred-remission",
      "##reference=GRCh37",
      paste0("##contig=<ID=", chromosome, ">"),
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Unphased genotype\">"
    ),
    connection
  )
  writeLines(
    paste(c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", samples), collapse = "\t"),
    connection
  )

  for (variant_index in seq_len(nrow(chromosome_markers))) {
    source_values <- as.numeric(genotypes[, chromosome_indexes[[variant_index]], drop = TRUE])
    if (anyNA(source_values) || any(!source_values %in% c(0, 1, 2))) {
      close(connection)
      unlink(temporary_path)
      stop("Invalid or missing genotype calls for ", chromosome_markers$SNP[[variant_index]])
    }

    alt_dosage <- if (chromosome_markers$CODE2_SIDE[[variant_index]] == "ALT") source_values else 2 - source_values
    expected_alt_count[[variant_index]] <- sum(alt_dosage)
    observed_alt_frequency[[variant_index]] <- expected_alt_count[[variant_index]] / (2 * sample_count)
    genotype_strings <- c("0/0", "0/1", "1/1")[alt_dosage + 1L]

    fields <- c(
      chromosome,
      chromosome_markers$GENOPRED_BP[[variant_index]],
      chromosome_markers$SNP[[variant_index]],
      chromosome_markers$GENOPRED_REF[[variant_index]],
      chromosome_markers$GENOPRED_ALT[[variant_index]],
      ".", "PASS", ".", "GT", genotype_strings
    )
    writeLines(paste(fields, collapse = "\t"), connection)
  }

  close(connection)
  if (file.exists(vcf_path)) unlink(vcf_path)
  if (!file.rename(temporary_path, vcf_path)) {
    unlink(temporary_path)
    stop("Could not move completed VCF into place: ", vcf_path)
  }

  difference <- abs(expected_alt_frequency - observed_alt_frequency)
  if (length(difference) > 0L && any(difference > 1e-12)) {
    stop("Chromosome ", chromosome, " VCF frequencies do not reproduce the source coding")
  }

  manifest <- data.frame(
    CHROM = chromosome,
    POS = chromosome_markers$GENOPRED_BP,
    ID = chromosome_markers$SNP,
    REF = chromosome_markers$GENOPRED_REF,
    ALT = chromosome_markers$GENOPRED_ALT,
    CODE2_SIDE = chromosome_markers$CODE2_SIDE,
    CODE2_ALLELE = chromosome_markers$CODE2_ALLELE,
    CODE0_SIDE = chromosome_markers$CODE0_SIDE,
    CODE0_ALLELE = chromosome_markers$CODE0_ALLELE,
    EXPECTED_ALT_COUNT = expected_alt_count,
    EXPECTED_OBS_CT = 2 * sample_count,
    EXPECTED_ALT_FREQ = expected_alt_frequency,
    RECODED_ALT_FREQ = observed_alt_frequency,
    ABSOLUTE_RECODING_DIFFERENCE = difference,
    EUR_ALT_FREQ = chromosome_markers$EUR_ALT_FREQ,
    DIFFERENCE_TO_EUR_MAJOR_FREQ = chromosome_markers$DIFFERENCE_TO_EUR_MAJOR_FREQ,
    ORIENTATION_DECISION = chromosome_markers$ORIENTATION_DECISION,
    stringsAsFactors = FALSE
  )
  write.table(manifest, manifest_path, sep = "\t", quote = FALSE, row.names = FALSE)

  summary_rows[[chromosome]] <- data.frame(
    chromosome = chromosome,
    participants = sample_count,
    variants = nrow(chromosome_markers),
    code2_is_ref = sum(chromosome_markers$CODE2_SIDE == "REF"),
    code2_is_alt = sum(chromosome_markers$CODE2_SIDE == "ALT"),
    missing_genotypes = 0,
    maximum_internal_alt_frequency_difference = if (length(difference) == 0L) 0 else max(difference),
    stringsAsFactors = FALSE
  )
  write.table(
    summary_rows[[chromosome]],
    file.path(summary_dir, paste0("GENDEP.chr", chromosome, ".build_summary.tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
  message("Built chromosome ", chromosome, ": ", nrow(chromosome_markers), " variants")
}

summary <- do.call(rbind, summary_rows)
if (sum(summary$variants) != nrow(markers)) {
  stop("Chromosome summaries do not close to the target marker count")
}
write.table(summary, file.path(output_dir, "target_build_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(
  data.frame(
    metric = c("source_samples", "source_variants", "target_variants", "chromosomes", "missing_genotypes"),
    value = c(sample_count, source_variant_count, nrow(markers), 22L, 0L),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "target_build_totals.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
capture.output(sessionInfo(), file = file.path(output_dir, "R_session_info.txt"))
message("Target VCF construction completed successfully: ", output_dir)
