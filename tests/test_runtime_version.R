#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", script_args, value = TRUE)

script_path <- if (length(file_arg) == 1L) {
  sub("^--file=", "", file_arg)
} else {
  "tests/test_runtime_version.R"
}

root <- normalizePath(
  file.path(dirname(script_path), ".."),
  mustWork = TRUE
)

source(file.path(root, "src", "R", "model_common.R"))

stopifnot(
  identical(normalise_package_version("4.1-8"), "4.1.8"),
  identical(normalise_package_version("4.1.8"), "4.1.8"),
  identical(
    normalise_package_version(" 4.7-1.1 "),
    "4.7.1.1"
  ),
  normalise_package_version("4.1-8") ==
    normalise_package_version("4.1.8")
)

cat("RUNTIME_VERSION_NORMALISATION_REGRESSION=PASS\n")
