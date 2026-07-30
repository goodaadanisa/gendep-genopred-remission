# Repository validation

- Thirteen final dissertation figures are included as PNG files.
- Five main and 15 supplementary tables are included in `tables.md`.
- Twenty-six shareable inputs regenerate 12 figures. An additional aggregate summary accompanies Supplementary Figure S4, whose participant-level probability input is excluded.
- Exact table sources are stored in `table_data.yml`.
- Synthetic tests cover genotype processing, PRS preparation, clinical integration, resampling, model metrics, aggregation, orientation sensitivity and reporting structure.
- GitHub Actions checks Python compilation, Python tests, public-release privacy, R parsing, both R regression tests and shell syntax.

## Controlled workflow checks

Controlled-environment validation records include:

- 430 participants, 147,370 target variants and 63,369,100 genotype comparisons with zero genotype or variant-metadata mismatches;
- 3,440 participant-by-trait PRS values checked after identifier alignment;
- primary and strict-EUR analysis matrices covering 130,592 cells;
- fixed predictor lists and deterministic nested-CV fold manifests; and
- a complete primary Elastic Net outer split covering predictions, tuning values, coefficients and performance metrics.

## Public release boundary

Participant-level clinical data, genotype calls, PRS, ancestry probabilities, resampling assignments, held-out predictions and model objects are excluded. Full execution requires authorised data, external genetic resources and the recorded software environments.
