#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

repository <- normalizePath(".", mustWork = TRUE)
script <- normalizePath(
  file.path(repository, "scripts", "build_target.R"),
  mustWork = TRUE
)
rscript <- file.path(R.home("bin"), "Rscript")

if (!file.exists(rscript)) {
  stop("Rscript was not found at: ", rscript)
}

work <- tempfile("build_target_regression_")
dir.create(work, recursive = TRUE)
on.exit(unlink(work, recursive = TRUE), add = TRUE)

variant_count <- 22L
variant_ids <- paste0("rs", seq_len(variant_count))

data.gen <- matrix(
  rep(c(2, 0), variant_count),
  nrow = 2,
  ncol = variant_count
)
colnames(data.gen) <- variant_ids

input_path <- file.path(work, "dataGen.RData")
save(data.gen, file = input_path)

sample_path <- file.path(work, "samples.txt")
writeLines(c("S1", "S2"), sample_path)

markers <- data.frame(
  SNP = variant_ids,
  GENOPRED_CHR = seq_len(variant_count),
  GENOPRED_BP = seq_len(variant_count) * 100L,
  GENOPRED_REF = rep("A", variant_count),
  GENOPRED_ALT = rep("G", variant_count),
  GENDEP_CODE2_FREQ = rep(0.5, variant_count),
  EUR_ALT_FREQ = rep(0.5, variant_count),
  DIFFERENCE_TO_EUR_MAJOR_FREQ = rep(0, variant_count),
  CODE2_SIDE = rep("ALT", variant_count),
  CODE0_SIDE = rep("REF", variant_count),
  CODE2_ALLELE = rep("G", variant_count),
  CODE0_ALLELE = rep("A", variant_count),
  ORIENTATION_DECISION = rep("DIRECT", variant_count)
)

run_build <- function(marker_path, output_dir, expected_variants = NULL) {
  arguments <- c(
    script,
    "--input", input_path,
    "--markers", marker_path,
    "--samples", sample_path,
    "--output-dir", output_dir,
    "--object", "data.gen",
    "--expected-samples", "2"
  )

  if (!is.null(expected_variants)) {
    arguments <- c(
      arguments,
      "--expected-variants",
      as.character(expected_variants)
    )
  }

  output <- suppressWarnings(
    system2(
      rscript,
      arguments,
      stdout = TRUE,
      stderr = TRUE
    )
  )

  status <- attr(output, "status")
  if (is.null(status)) {
    status <- 0L
  }

  list(status = status, output = output)
}

final_marker_path <- file.path(work, "final_markers.tsv")
write.table(
  markers,
  final_marker_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

final_output <- file.path(work, "final_output")
final_result <- run_build(
  final_marker_path,
  final_output,
  expected_variants = variant_count
)

if (final_result$status != 0L) {
  stop(
    "Already-filtered marker-table test failed:\n",
    paste(final_result$output, collapse = "\n")
  )
}

totals <- read.delim(
  file.path(final_output, "target_build_totals.tsv"),
  check.names = FALSE
)
observed <- setNames(as.character(totals$value), totals$metric)

if (observed[["source_samples"]] != "2") {
  stop("Unexpected source sample count")
}
if (
  observed[["target_variants"]] != as.character(variant_count)
) {
  stop("Unexpected target variant count")
}
if (observed[["chromosomes"]] != "22") {
  stop("Unexpected chromosome count")
}

vcf_files <- list.files(
  file.path(final_output, "vcf"),
  pattern = "^GENDEP[.]chr[0-9]+[.]vcf$",
  full.names = TRUE
)
if (length(vcf_files) != 22L) {
  stop("Expected 22 chromosome VCF files")
}

missing_expected_output <- file.path(work, "missing_expected_output")
missing_expected_result <- run_build(
  final_marker_path,
  missing_expected_output
)

if (missing_expected_result$status == 0L) {
  stop(
    "An already-filtered table without --expected-variants ",
    "should have been rejected"
  )
}

if (!any(grepl(
  "--expected-variants is required",
  missing_expected_result$output,
  fixed = TRUE
))) {
  stop("The expected missing --expected-variants error was not produced")
}

filterable_markers <- markers
filterable_markers$KEEP_FOR_TARGET <- TRUE

unused_marker <- filterable_markers[1, , drop = FALSE]
unused_marker$SNP <- "rs_unused"
unused_marker$GENOPRED_BP <- 150L
unused_marker$KEEP_FOR_TARGET <- FALSE

filterable_markers <- rbind(
  filterable_markers,
  unused_marker
)

filterable_marker_path <- file.path(work, "filterable_markers.tsv")
write.table(
  filterable_markers,
  filterable_marker_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

filterable_output <- file.path(work, "filterable_output")
filterable_result <- run_build(
  filterable_marker_path,
  filterable_output,
  expected_variants = variant_count
)

if (filterable_result$status != 0L) {
  stop(
    "KEEP_FOR_TARGET filtering test failed:\n",
    paste(filterable_result$output, collapse = "\n")
  )
}

cat("BUILD_TARGET_REGRESSION=PASS\n")
