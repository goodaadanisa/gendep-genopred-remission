# Novel Application of GenoPred Polygenic Risk Scores to Predict Antidepressant Treatment Response Using AI

Reproducible code and release-safe outputs for an MSc Applied Bioinformatics project evaluating whether eight Bayesian polygenic risk scores improve machine-learning prediction of binary antidepressant remission beyond clinical variables and ancestry principal components.

> **Repository scope:** This is a controlled-data, source-available assessment repository. The public release supports methodological inspection, automated testing and regeneration of release-safe outputs. Full end-to-end execution requires authorised access to GENDEP participant-level data and the specified external resources. The repository is not distributed as open-source software.

## Repository at a glance

Eight executable files expose the workflow; reusable implementation code lives under `src/`.

```text
scripts/
├── gendep.py                 # Compatibility wrapper for the installed Python CLI
├── genotype_audit.R          # Audit the supplied rsID-only genotype matrix
├── build_target.R            # Build chromosome-specific target VCFs
├── collect_prs.R             # Collect and validate eight participant PRS
├── run_models.R              # Run one configured nested-CV outer split
├── bootstrap_uncertainty.R   # Paired participant-level uncertainty analysis
├── secondary_analyses.R      # Individual PRS and treatment interactions
└── run_genopred.sh           # Configure and launch GenoPred/PRS-CS-auto
```

The installed `gendep` command provides named Python subcommands; implementation modules are not intended to be executed directly.

## Scientific workflow

```text
Clinical and genotype audit
        ↓
Variant metadata annotation and allele-orientation inference
        ↓
Validated PLINK2 target and Bayesian PRS generation
        ↓
Clinical, PRS and ancestry-PC integration
        ↓
Repeated nested Elastic Net and Random Forest analyses
        ↓
Paired uncertainty, secondary analyses and sensitivity checks
        ↓
Tables and figures from release-safe aggregate outputs
```

The final target contained 147,370 variants, of which 141,038 were compatible with the selected PRS-CS LD resource. Mean held-out AUC was 0.694 for the clinical comparator and 0.688 after adding eight PRS.

## Installation and public checks

```bash
python -m pip install -e '.[test]'
gendep --help
python -m pytest
python -m compileall -q scripts src tests
gendep validate
```

The compatibility wrapper `python scripts/gendep.py ...` remains available when the package is not installed.

The public reporting and validation commands are cross-platform. Controlled-data GenoPred and orientation-sensitivity stages target Linux/HPC. On Windows, the symlink-dependent synthetic integration test is skipped when the operating system does not grant symlink privileges; the full test remains active in Linux CI.

## Results and reporting inputs

- [`results/tables.md`](results/tables.md) contains five main and 15 supplementary tables.
- [`results/figures/`](results/figures/) contains five main and eight supplementary figures as PNG files.
- [`results/figure_data/`](results/figure_data/) contains 27 shareable figure-source files: 26 regenerate 12 figures, while one aggregate summary accompanies Supplementary Figure S4, whose participant-level plotting input remains controlled.
- [`results/table_data.yml`](results/table_data.yml) preserves the unrounded machine-readable source for all 20 tables.

```bash
gendep figures
gendep tables
```

Reporting does not access participant identifiers, generate PRS, fit models or perform bootstrap analyses.

## Controlled data

Participant-level clinical data, genotype calls, PRS, ancestry probabilities, resampling assignments and held-out predictions are not included. See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md).

## Documentation

- [Controlled execution workflow](docs/workflow.md)
- [Methods-to-code map](docs/code_map.md)
- [Reproducibility and validation boundary](docs/reproducibility.md)
- [Interpretation and limitations](docs/interpretation.md)

## Validation boundary

The repository includes automated unit, integration, syntax and privacy checks, together with controlled-environment validation records for the genotype target, PRS collection, analysis-data integration, predictor policy, resampling and representative model outputs. See [`results/validation_summary.md`](results/validation_summary.md).
