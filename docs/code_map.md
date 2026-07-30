# Methods-to-code map

| Dissertation component | Public command | Main implementation |
|---|---|---|
| Supplied genotype audit | `Rscript scripts/genotype_audit.R` | R script |
| Variant metadata annotation and orientation | `gendep reconstruct-target` | `gendep.genotype` |
| Target VCF construction | `Rscript scripts/build_target.R` | R script |
| PLINK2 round-trip validation | `gendep validate-target` | `gendep.genotype.validation` |
| GWAS standardisation | `gendep prepare-gwas` | `gendep.prs.gwas` |
| GenoPred and PRS-CS-auto | `bash scripts/run_genopred.sh` | fixed configuration files |
| Participant PRS collection | `Rscript scripts/collect_prs.R` | R script |
| Clinical audit and predictor policy | `gendep audit-clinical`, `gendep define-predictors` | `gendep.clinical` |
| Clinical-genetic integration | `gendep build-analysis` | `gendep.clinical.integration` |
| Nested resampling | `gendep resampling` | `gendep.resampling` |
| Elastic Net and Random Forest | `Rscript scripts/run_models.R` | `src/R/*_engine.R` |
| Model aggregation and calibration | `gendep aggregate-models`, `gendep model-diagnostics` | `gendep.modeling` |
| Paired uncertainty | `Rscript scripts/bootstrap_uncertainty.R` | R script |
| Individual PRS and interactions | `Rscript scripts/secondary_analyses.R` | `secondary_analysis_engine.R` |
| Conservative orientation sensitivity | `gendep orientation-sensitivity` | `gendep.sensitivity.orientation` |
| Tables and figures | `gendep tables`, `gendep figures` | `gendep.figures` and release-safe aggregate data |
| Public release checks | `gendep validate` | structural and privacy validation |
