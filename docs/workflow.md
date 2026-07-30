# Controlled execution workflow

This page gives the execution order for an authorised environment. Copy `config/project.example.yml` to the ignored `config/project.yml`, replace every placeholder path, install the Python package, and load the fixed R/PLINK/GenoPred runtime before starting.

```bash
cp config/project.example.yml config/project.yml
python -m pip install -e '.[test]'
```

The commands below write participant-level intermediates beneath `work/`; that directory must remain outside version control.

## 1. Audit the supplied sources

```bash
Rscript scripts/genotype_audit.R \
  --input /authorised/path/dataGen.RData \
  --output-dir work/genotype/source_audit \
  --object data.gen \
  --expected-samples 430 \
  --expected-variants 524876

gendep audit-clinical \
  --clinical /authorised/path/data-rem.xlsx \
  --output-dir work/analysis/clinical_audit \
  --expected-rows 430
```

The genotype audit exports only variant-level counts and frequencies. The clinical audit exports aggregate checks and predictor metadata, not participant rows.

## 2. Prepare, build and validate the genotype target

```bash
gendep reconstruct-target \
  --config config/project.yml \
  --expected config/expected_results.yml

Rscript scripts/build_target.R \
  --input /authorised/path/dataGen.RData \
  --markers work/genotype/reconstruction/final_target_markers.tsv.gz \
  --samples /authorised/path/validated_sample_ids.txt \
  --output-dir work/genotype/target \
  --expected-samples 430 \
  --expected-variants 147370

gendep validate-target \
  --target-dir work/genotype/target \
  --output-dir work/genotype/target_validation \
  --plink2 /path/to/plink2 \
  --threads 4 \
  --expected config/expected_results.yml
```

This stage is the authoritative allele-aware round-trip check. Do not substitute a file-existence check for `gendep validate-target`.

## 3. Standardise GWAS inputs and run GenoPred/PRS-CS-auto

```bash
gendep prepare-gwas \
  --config config/project.yml \
  --gwas-config config/gwas.yml

bash scripts/run_genopred.sh \
  --config config/project.yml \
  --cores 4
```

The run script verifies the fixed GenoPred and GenoUtils commits unless `--skip-commit-check` is explicitly supplied. That override is for diagnosis only and is not part of the specified workflow.

## 4. Collect the eight PRS and build analysis matrices

```bash
Rscript scripts/collect_prs.R \
  --genopred-output work/prs/genopred/pipeline_output \
  --psam work/genotype/target_validation/pfiles/GENDEP.chr1.psam \
  --traits config/traits.tsv \
  --output-dir work/prs/final \
  --audit-dir work/prs/validation \
  --expected-participants 430

gendep define-predictors \
  --clinical /authorised/path/data-rem.xlsx \
  --output-dir work/analysis/predictor_policy

gendep build-analysis \
  --clinical /authorised/path/data-rem.xlsx \
  --prs work/prs/final/gendep_prscs_auto_8trait_prs.tsv.gz \
  --pcs /authorised/path/GENDEP-TRANS.profiles \
  --ancestry /authorised/path/GENDEP.Ancestry.model_pred \
  --eur-keep /authorised/path/EUR.keep \
  --output-dir work/analysis/datasets \
  --audit-dir work/analysis/integration_audit
```

All participant joins are by validated identifier. The clinical workbook remains the master row order.

## 5. Generate the fixed repeated nested-CV assignments

```bash
gendep resampling \
  --cohort primary_430=work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz:0:430 \
  --cohort EUR_418=work/analysis/datasets/gendep_clinical_prs_pc_analysis_base_eur418.tsv.gz:50000000:418 \
  --config config/analyses.yml \
  --output-dir work/analysis/resampling
```

The manifests are deterministic and jointly stratified by remission outcome and treatment group.

## 6. Fit configured outer splits

Each call fits one outer split. The complete primary analysis consists of 10 repeats × 10 folds.

```bash
for repeat in $(seq 1 10); do
  for fold in $(seq 1 10); do
    Rscript scripts/run_models.R \
      --analysis primary \
      --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
      --predictor-policy-dir work/analysis/predictor_policy \
      --resampling-dir work/analysis/resampling \
      --output-dir work/modelling \
      --repeat "$repeat" \
      --fold "$fold"
  done
done
```

Other predefined analyses use the names in `config/model_analyses.tsv` and the corresponding primary or strict-EUR analysis base. Do not tune on held-out outer-fold outcomes.

## 7. Aggregate models and calculate uncertainty

```bash
gendep aggregate-models \
  --analysis primary \
  --split-root work/modelling \
  --output-dir work/modelling/aggregated/primary

gendep model-diagnostics \
  --predictions work/modelling/aggregated/primary/all_held_out_predictions.tsv.gz \
  --coefficients work/modelling/aggregated/primary/all_coefficients.tsv.gz \
  --tuning work/modelling/aggregated/primary/all_tuning_results.tsv.gz \
  --output-dir work/modelling/diagnostics/primary

Rscript scripts/bootstrap_uncertainty.R \
  --analysis primary \
  --predictions work/modelling/aggregated/primary/all_held_out_predictions.tsv.gz \
  --output-dir work/modelling/uncertainty/primary
```

Aggregation requires all configured splits. Paired uncertainty resamples participants from held-out predictions and does not refit models.

## 8. Run predefined secondary and orientation-sensitivity analyses

```bash
Rscript scripts/secondary_analyses.R \
  --analysis individual_prs \
  --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
  --outer-fold-file work/analysis/resampling/primary_430_outer_folds.tsv.gz \
  --output-dir work/secondary/individual_prs

Rscript scripts/secondary_analyses.R \
  --analysis treatment_interactions \
  --analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz \
  --outer-fold-file work/analysis/resampling/primary_430_outer_folds.tsv.gz \
  --output-dir work/secondary/treatment_interactions

gendep orientation-sensitivity --help
```

The orientation-sensitivity command has several controlled runtime paths; use its `--help` output together with the primary score root and GenoPred configuration.

## 9. Regenerate release-safe reporting outputs

```bash
gendep tables
gendep figures
gendep validate
```

These commands use release-safe aggregate inputs. They do not access controlled participant-level inputs, generate PRS, fit predictive models or perform bootstrap analyses.
