# Reproducibility

## Environments

Python dependencies are pinned in `pyproject.toml`. Recorded software versions are listed under `environment/`; `environment/check_R_versions.R` verifies the R packages used directly by the workflow. PLINK2, Snakemake, GenoPred, PRS-CS and reference resources must be supplied in the authorised computing environment.

## Configuration

Copy `config/project.example.yml` to the ignored `config/project.yml` and supply machine-specific paths. Runtime configuration is deliberately separated by responsibility:

- `project.yml`: authorised paths and executable locations;
- `genopred.yml`, `gwas.yml` and `traits.tsv`: fixed PRS inputs and settings;
- `analyses.yml`: deterministic repeated nested-CV design;
- `model_analyses.tsv` and `model_parameters.tsv`: model specifications and parameters;
- `expected_results.yml`: aggregate validation expectations.

Predictor lists and clinical schema decisions are generated from the fixed implementation and validated against the controlled source workbook rather than duplicated in decorative configuration files.

## Randomness and leakage control

The analysis uses fixed seeds and stored resampling assignments. Data-dependent preprocessing and model tuning are confined to training partitions. Paired uncertainty operates on held-out predictions and does not refit models.

## Public reporting

Release-safe figure inputs are committed under `results/figure_data/`, rather than hidden in a nested archive. The participant-level ancestry probability input for Supplementary Figure S4 remains controlled. Exact table values are stored in `results/table_data.yml`; `results/tables.md` is a rounded human-readable rendering. Reporting commands never access controlled participant-level inputs.

## Automated checks

GitHub Actions runs Python compilation, the synthetic test suite, public release validation, both R regression tests, R syntax parsing and shell syntax checking on every push and pull request.

## Validation coverage

Controlled-environment validation records cover 63,369,100 target-genotype comparisons, all 3,440 participant-by-trait PRS assignments, both integrated analysis matrices, every fixed predictor list, all nested-CV fold manifests and one complete primary Elastic Net outer split. Public tests additionally cover synthetic genotype processing, PRS preparation, clinical integration, resampling, model metrics, aggregation, sensitivity analyses and reporting structure.

The public repository excludes controlled participant-level inputs and therefore requires authorised data and external resources for full workflow execution.
