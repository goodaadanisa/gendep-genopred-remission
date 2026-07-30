# Software environments

The workflow used component-specific environments. The exact recorded components are listed in `component_versions.tsv`.

## R packages

Principal R dependencies included `data.table`, `dplyr`, `glmnet`, `pROC`, `randomForest`, `readxl` and `tibble`.

## Python release environment

The public workflow and reporting dependencies are pinned in `pyproject.toml`. Install them with:

```bash
python -m pip install -e '.[test]'
```

## External executables and resources

Full execution requires compatible installations of R, PLINK/PLINK2, Snakemake, GenoPred and PRS-CS plus the fixed external genetic reference files and eight GWAS inputs. These are not redistributed.
