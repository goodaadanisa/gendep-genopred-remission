"""Generate and validate the outcome-blind predictor and encoding policy."""

from __future__ import annotations

from collections import Counter

import pandas as pd

from .audit import canonical_vector
from .schema import (
    BINARY_CATEGORICAL,
    CLINICAL_COLUMNS,
    COMPREHENSIVE_CLINICAL,
    EUR_INDICATOR,
    EXPECTED_MODEL_COUNTS,
    EXPECTED_SEMANTIC_COUNTS,
    EXACT_DUPLICATE_GROUPS,
    IDENTIFIER_COLUMNS,
    MODEL_PREDICTOR_SETS,
    OUTCOME_COLUMN,
    PC_COLUMNS,
    REMOVED_EXACT_DUPLICATES,
    SUMMARY_CLINICAL,
    expected_storage_class,
    final_disposition,
    model_role,
    semantic_type,
    source_section,
)


def ordered_levels(series: pd.Series) -> list[str]:
    """Return non-missing categorical levels in deterministic order."""
    values = [value for value in canonical_vector(series) if value != "<NA>"]
    unique = sorted(set(values), key=lambda value: (float(value), value))
    return unique


def build_predictor_policy(frame: pd.DataFrame) -> dict[str, object]:
    """Build the fixed outcome-blind predictor, duplicate-removal and encoding policy."""
    checks: list[dict[str, object]] = []

    def check(metric: str, observed: object, expected: object) -> None:
        """Append one predictor-policy validation result."""
        checks.append({
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "result": "PASS" if observed == expected else "FAIL",
        })

    check("comprehensive_clinical_predictors", len(COMPREHENSIVE_CLINICAL), 128)
    check("summary_clinical_predictors", len(SUMMARY_CLINICAL), 46)
    check("removed_exact_duplicates", len(REMOVED_EXACT_DUPLICATES), 7)
    for name, predictors in MODEL_PREDICTOR_SETS.items():
        check(f"{name}_predictors", len(predictors), EXPECTED_MODEL_COUNTS[name])
        check(f"{name}_unique", len(predictors), len(set(predictors)))
        check(f"{name}_outcome_excluded", OUTCOME_COLUMN not in predictors, True)
        check(f"{name}_identifiers_excluded", not (set(IDENTIFIER_COLUMNS) & set(predictors)), True)
        check(f"{name}_eur_indicator_excluded", EUR_INDICATOR not in predictors, True)

    for group in EXACT_DUPLICATE_GROUPS:
        reference = canonical_vector(frame[group[0]])
        check(
            f"duplicate_group_{group[0]}",
            all(canonical_vector(frame[column]) == reference for column in group[1:]),
            True,
        )

    inventory: list[dict[str, object]] = []
    for position, variable in enumerate(CLINICAL_COLUMNS, start=1):
        models = [name for name, predictors in MODEL_PREDICTOR_SETS.items() if variable in predictors]
        inventory.append({
            "position": position,
            "variable": variable,
            "source_section": source_section(variable),
            "expected_storage_class": expected_storage_class(variable),
            "final_disposition": final_disposition(variable),
            "retained_comprehensive": variable in COMPREHENSIVE_CLINICAL,
            "retained_summary": variable in SUMMARY_CLINICAL,
            "model_membership": ";".join(models),
        })

    model_rows: list[dict[str, object]] = []
    for model, predictors in MODEL_PREDICTOR_SETS.items():
        for order, variable in enumerate(predictors, start=1):
            model_rows.append({
                "model": model,
                "predictor_order": order,
                "variable": variable,
                "semantic_type": semantic_type(variable),
                "model_role": model_role(variable, model),
                "penalty_factor": 0 if variable in (*PC_COLUMNS, "drug") else 1,
            })

    combined = MODEL_PREDICTOR_SETS["primary_comprehensive_combined"]
    semantic_counts = Counter(semantic_type(variable) for variable in combined)
    semantic_rows = [
        {
            "semantic_type": kind,
            "observed_count": semantic_counts.get(kind, 0),
            "expected_count": expected,
            "result": "PASS" if semantic_counts.get(kind, 0) == expected else "FAIL",
        }
        for kind, expected in EXPECTED_SEMANTIC_COUNTS.items()
    ]
    check("all_predictors_classified", semantic_counts.get("unclassified", 0), 0)

    encoding_rows: list[dict[str, object]] = []
    binary_rows: list[dict[str, object]] = []
    for variable in combined:
        kind = semantic_type(variable)
        values = frame[variable] if variable in frame.columns else None
        levels = ordered_levels(values) if values is not None else []
        minimum = ""
        maximum = ""
        if values is not None:
            numeric = pd.to_numeric(values, errors="coerce")
            if numeric.notna().all():
                minimum = float(numeric.min())
                maximum = float(numeric.max())
        encoding_rows.append({
            "variable": variable,
            "semantic_type": kind,
            "source_levels_or_range": ";".join(levels) if len(levels) <= 20 else f"{len(levels)} unique values",
            "minimum": minimum,
            "maximum": maximum,
            "elastic_net_encoding": (
                "fixed lower-source-code=0, higher-source-code=1"
                if kind == "binary_categorical"
                else "retain numeric value; standardise within training fold"
            ),
            "random_forest_encoding": (
                "fixed lower-source-code=0, higher-source-code=1"
                if kind == "binary_categorical"
                else "retain numeric value without scaling"
            ),
        })
        if variable in BINARY_CATEGORICAL:
            check(f"{variable}_binary_levels", len(levels), 2)
            if len(levels) == 2:
                counts = Counter(canonical_vector(values))
                binary_rows.append({
                    "variable": variable,
                    "source_level_mapped_to_zero": levels[0],
                    "zero_level_count": counts[levels[0]],
                    "source_level_mapped_to_one": levels[1],
                    "one_level_count": counts[levels[1]],
                    "mapping_rule": "lower ordered source code to 0; higher ordered source code to 1",
                })

    return {
        "checks": checks,
        "inventory": inventory,
        "model_rows": model_rows,
        "semantic_rows": semantic_rows,
        "encoding_rows": encoding_rows,
        "binary_rows": binary_rows,
        "passed": all(row["result"] == "PASS" for row in checks)
            and all(row["result"] == "PASS" for row in semantic_rows),
    }
