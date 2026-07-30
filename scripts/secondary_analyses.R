#!/usr/bin/env Rscript
# Purpose:
# Run either the eight individual-PRS analyses or the three prespecified
# PRS-by-treatment-group interaction analyses through one documented command.
#
# Inputs:
# Fixed 430-participant analysis base and primary outer-fold manifest.
#
# Outputs:
# Full-sample association estimates, repeated held-out predictions, metric
# summaries and validation records. Individual PRS analyses additionally use
# 10,000 paired participant-level bootstrap replicates by default.
#
# Usage:
# Rscript scripts/secondary_analyses.R --analysis individual_prs \
#   --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
#   --outer-fold-file work/analysis/resampling/primary_430_outer_folds.tsv.gz \
#   --output-dir work/secondary/individual_prs
#
# Rscript scripts/secondary_analyses.R --analysis treatment_interactions \
#   --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
#   --outer-fold-file work/analysis/resampling/primary_430_outer_folds.tsv.gz \
#   --output-dir work/secondary/treatment_interactions

script_path <- normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
repository_root <- dirname(dirname(script_path))
source(file.path(repository_root, "src", "R", "model_common.R"))
source(file.path(repository_root, "src", "R", "secondary_analysis_engine.R"))

options <- parse_cli(commandArgs(trailingOnly = TRUE))
analysis_name <- required_cli(options, "analysis")
analysis_base <- required_cli(options, "analysis-base")
outer_fold_file <- required_cli(options, "outer-fold-file")
output_dir <- required_cli(options, "output-dir")
bootstrap_replicates <- suppressWarnings(as.integer(options[["bootstrap-replicates"]] %||% "10000"))
bootstrap_seed <- suppressWarnings(as.integer(options[["bootstrap-seed"]] %||% "20260715"))
if (is.na(bootstrap_replicates) || is.na(bootstrap_seed)) stop("Invalid bootstrap settings.")
if (!analysis_name %in% c("individual_prs", "treatment_interactions")) {
  stop("--analysis must be individual_prs or treatment_interactions")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
inputs <- prepare_secondary_inputs(analysis_base, outer_fold_file)
if (analysis_name == "individual_prs") {
  run_individual_prs_analysis(inputs, output_dir, bootstrap_replicates, bootstrap_seed)
} else {
  run_treatment_interactions(inputs, output_dir)
}
writeLines(
  c(
    paste0("analysis=", analysis_name),
    "participants=430",
    "outer_repeats=10",
    "outer_folds=10",
    "training_only_PRS_standardisation=TRUE",
    paste0("bootstrap_replicates=", if (analysis_name == "individual_prs") bootstrap_replicates else 0L),
    paste0("bootstrap_seed=", if (analysis_name == "individual_prs") bootstrap_seed else "NA"),
    "analysis_status=PASS"
  ),
  file.path(output_dir, "_SUCCESS")
)
cat("SECONDARY_ANALYSIS=PASS\n")
cat(paste0("analysis=", analysis_name, "\n"))
cat(paste0("output_dir=", normalizePath(output_dir), "\n"))
