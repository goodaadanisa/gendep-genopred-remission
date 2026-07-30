"""Independent checks for controlled primary and strict-EUR analysis bases."""

from __future__ import annotations

from collections import Counter

import numpy as np
import pandas as pd

from .audit import canonical_vector
from .io import canonical_identifier
from .schema import (
    CLINICAL_COLUMNS,
    EUR_INDICATOR,
    EXPECTED_MODEL_COUNTS,
    IDENTIFIER_COLUMNS,
    MODEL_PREDICTOR_SETS,
    OUTCOME_COLUMN,
    PC_COLUMNS,
    PRS_COLUMNS,
    TREATMENT_COLUMN,
)


def _column_values_equal(left: pd.Series, right: pd.Series) -> bool:
    """Compare a TSV round trip without mistaking equivalent float text for change."""
    left_numeric = pd.to_numeric(left, errors="coerce")
    right_numeric = pd.to_numeric(right, errors="coerce")
    if left_numeric.notna().all() and right_numeric.notna().all():
        return np.allclose(
            left_numeric.to_numpy(dtype=float),
            right_numeric.to_numpy(dtype=float),
            rtol=0.0,
            atol=1e-12,
            equal_nan=False,
        )
    return canonical_vector(left) == canonical_vector(right)


def validate_analysis_bases(
    primary: pd.DataFrame,
    strict_eur: pd.DataFrame,
    *,
    expected_primary_rows: int = 430,
    expected_eur_rows: int = 418,
    source_clinical: pd.DataFrame | None = None,
) -> list[dict[str, object]]:
    """Validate primary and strict-EUR analysis bases against linkage, schema and completeness invariants."""
    checks: list[dict[str, object]] = []

    def check(metric: str, observed: object, expected: object, stage: str = "analysis_base") -> None:
        """Append one analysis-integration validation result."""
        checks.append({
            "stage": stage,
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "result": "PASS" if observed == expected else "FAIL",
        })

    expected_header = (*CLINICAL_COLUMNS, *PRS_COLUMNS, *PC_COLUMNS, EUR_INDICATOR)
    check("primary_rows", len(primary), expected_primary_rows)
    check("strict_eur_rows", len(strict_eur), expected_eur_rows)
    check("primary_columns", len(primary.columns), 154)
    check("header_matches_fixed_schema", tuple(primary.columns) == expected_header, True)
    check("strict_eur_header_matches_primary", tuple(strict_eur.columns) == tuple(primary.columns), True)
    check("primary_missing_cells", int(primary.isna().sum().sum()), 0)
    check("strict_eur_missing_cells", int(strict_eur.isna().sum().sum()), 0)

    ids = primary["Row.names"].map(canonical_identifier)
    check("participant_ids_unique", int(ids.nunique()), len(primary))
    check("Row.names_equals_bloodsampleid.x", canonical_vector(primary["Row.names"]) == canonical_vector(primary["bloodsampleid.x"]), True)

    outcome_counts = Counter(pd.to_numeric(primary[OUTCOME_COLUMN]).astype(int))
    treatment_counts = Counter(pd.to_numeric(primary[TREATMENT_COLUMN]).astype(int))
    if expected_primary_rows == 430:
        check("non_remitters", outcome_counts[0], 264)
        check("remitters", outcome_counts[1], 166)
        check("treatment_group_1", treatment_counts[1], 210)
        check("treatment_group_2", treatment_counts[2], 220)

    for block_name, columns in {"PRS": PRS_COLUMNS, "PC": PC_COLUMNS}.items():
        values = primary.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
        check(f"{block_name}_missing_or_nonfinite", int((~np.isfinite(values)).sum()), 0)
    indicator = pd.to_numeric(primary[EUR_INDICATOR], errors="coerce")
    check("EUR_indicator_binary", set(indicator.astype(int).unique()) <= {0, 1}, True)
    check("EUR_indicator_sum", int(indicator.sum()), expected_eur_rows)
    subset = primary.loc[indicator == 1].reset_index(drop=True)
    check("strict_eur_is_exact_indicator_subset", subset.equals(strict_eur.reset_index(drop=True)), True)

    if source_clinical is not None:
        check("source_clinical_header_match", tuple(source_clinical.columns) == CLINICAL_COLUMNS, True)
        check(
            "clinical_cells_preserved_exactly",
            all(_column_values_equal(primary[column], source_clinical[column]) for column in CLINICAL_COLUMNS),
            True,
        )

    all_columns = set(primary.columns)
    forbidden = set(IDENTIFIER_COLUMNS) | {OUTCOME_COLUMN, EUR_INDICATOR}
    for model, predictors in MODEL_PREDICTOR_SETS.items():
        check(f"{model}_count", len(predictors), EXPECTED_MODEL_COUNTS[model], "predictor_policy")
        check(f"{model}_columns_present", set(predictors) <= all_columns, True, "predictor_policy")
        check(f"{model}_forbidden_excluded", not (forbidden & set(predictors)), True, "predictor_policy")
        check(f"{model}_unique", len(predictors), len(set(predictors)), "predictor_policy")

    return checks
