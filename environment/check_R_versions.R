#!/usr/bin/env Rscript
# Check the R packages used directly by the modelling workflow.

required <- c(
  glmnet = "4.1.8",
  randomForest = "4.7.1.1"
)

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("data.table is required but not installed", call. = FALSE)
}

failures <- character()
for (package in names(required)) {
  if (!requireNamespace(package, quietly = TRUE)) {
    failures <- c(failures, paste(package, "is not installed"))
    next
  }
  observed <- as.character(utils::packageVersion(package))
  if (observed != required[[package]]) {
    failures <- c(
      failures,
      sprintf("%s version mismatch: expected %s, observed %s", package, required[[package]], observed)
    )
  }
}

if (length(failures)) {
  stop(paste(failures, collapse = "\n"), call. = FALSE)
}
cat("R_ENVIRONMENT_CHECK=PASS\n")
