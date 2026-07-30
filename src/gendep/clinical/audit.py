"""Release-safe audit of the supplied clinical workbook."""

from __future__ import annotations

import hashlib
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .schema import (
    BASELINE_CLINICAL_CANDIDATES,
    CLINICAL_COLUMNS,
    EXACT_DUPLICATE_GROUPS,
    IDENTIFIER_COLUMNS,
    OUTCOME_COLUMN,
    TREATMENT_COLUMN,
    expected_storage_class,
    final_disposition,
    source_section,
)


def sha256_file(path: str | Path) -> str:
    """Return a file's SHA-256 digest without loading it into memory."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def canonical_vector(series: pd.Series) -> tuple[str, ...]:
    """Convert vector to a stable canonical representation."""
    values: list[str] = []
    for value in series:
        if pd.isna(value):
            values.append("<NA>")
        elif isinstance(value, (int, np.integer)):
            values.append(str(int(value)))
        elif isinstance(value, (float, np.floating)) and float(value).is_integer():
            values.append(str(int(value)))
        else:
            values.append(str(value).strip())
    return tuple(values)


def observed_storage_class(series: pd.Series) -> str:
    """Classify a clinical series by its observed storage pattern."""
    nonmissing = series.dropna()
    numeric = pd.to_numeric(nonmissing, errors="coerce")
    if len(nonmissing) and numeric.notna().all():
        values = numeric.astype(float)
        integer = np.all(np.isclose(values, np.round(values), atol=1e-12))
        if integer and values.nunique(dropna=True) == 2:
            return "binary_integer"
        if integer:
            return "integer"
        return "continuous_numeric"
    return "text_or_categorical"


def exact_duplicate_groups(frame: pd.DataFrame) -> list[tuple[str, ...]]:
    """Return source-order groups of columns containing identical values."""
    groups: dict[tuple[str, ...], list[str]] = {}
    for column in BASELINE_CLINICAL_CANDIDATES:
        groups.setdefault(canonical_vector(frame[column]), []).append(column)
    duplicates = [tuple(columns) for columns in groups.values() if len(columns) > 1]
    return sorted(duplicates, key=lambda group: min(CLINICAL_COLUMNS.index(column) for column in group))


def audit_clinical(frame: pd.DataFrame, *, expected_rows: int = 430) -> dict[str, object]:
    """Audit the supplied clinical table without using outcome-informed selection rules."""
    checks: list[dict[str, object]] = []

    def check(metric: str, observed: object, expected: object) -> None:
        """Append one validation result to the clinical audit record."""
        checks.append({
            "metric": metric,
            "observed": observed,
            "expected": expected,
            "result": "PASS" if observed == expected else "FAIL",
        })

    check("participant_rows", len(frame), expected_rows)
    check("clinical_columns", len(frame.columns), 139)
    check("column_order_matches_fixed_schema", tuple(frame.columns) == CLINICAL_COLUMNS, True)
    check("column_names_unique", frame.columns.is_unique, True)
    check("missing_cells", int(frame.isna().sum().sum()), 0)

    for column in IDENTIFIER_COLUMNS:
        check(f"{column}_unique", int(frame[column].nunique(dropna=False)), expected_rows)
        check(f"{column}_missing", int(frame[column].isna().sum()), 0)
    check(
        "Row.names_equals_bloodsampleid.x",
        canonical_vector(frame["Row.names"]) == canonical_vector(frame["bloodsampleid.x"]),
        True,
    )

    outcome_counts = Counter(pd.to_numeric(frame[OUTCOME_COLUMN], errors="raise").astype(int))
    treatment_counts = Counter(pd.to_numeric(frame[TREATMENT_COLUMN], errors="raise").astype(int))
    check("outcome_0", outcome_counts[0], 264 if expected_rows == 430 else outcome_counts[0])
    check("outcome_1", outcome_counts[1], 166 if expected_rows == 430 else outcome_counts[1])
    check("treatment_1", treatment_counts[1], 210 if expected_rows == 430 else treatment_counts[1])
    check("treatment_2", treatment_counts[2], 220 if expected_rows == 430 else treatment_counts[2])

    observed_groups = exact_duplicate_groups(frame)
    check("exact_duplicate_group_count", len(observed_groups), len(EXACT_DUPLICATE_GROUPS))
    check("exact_duplicate_groups_match", set(observed_groups) == set(EXACT_DUPLICATE_GROUPS), True)

    inventory: list[dict[str, object]] = []
    for position, variable in enumerate(CLINICAL_COLUMNS, start=1):
        series = frame[variable]
        inventory.append({
            "position": position,
            "variable": variable,
            "source_section": source_section(variable),
            "observed_storage_class": observed_storage_class(series),
            "expected_storage_class": expected_storage_class(variable),
            "unique_values": int(series.nunique(dropna=True)),
            "missing_values": int(series.isna().sum()),
            "final_disposition": final_disposition(variable),
        })

    duplicate_rows: list[dict[str, object]] = []
    for group_index, group in enumerate(observed_groups, start=1):
        retained = group[0]
        for variable in group:
            duplicate_rows.append({
                "group": group_index,
                "group_size": len(group),
                "variable": variable,
                "source_position": CLINICAL_COLUMNS.index(variable) + 1,
                "retained_variable": retained,
                "decision": "retain" if variable == retained else "remove_later_exact_duplicate",
                "group_members": ";".join(group),
            })

    cohort = [
        {"metric": "participants", "value": len(frame)},
        {"metric": "clinical_columns", "value": len(frame.columns)},
        {"metric": "missing_cells", "value": int(frame.isna().sum().sum())},
        {"metric": "non_remitters", "value": outcome_counts[0]},
        {"metric": "remitters", "value": outcome_counts[1]},
        {"metric": "treatment_group_1", "value": treatment_counts[1]},
        {"metric": "treatment_group_2", "value": treatment_counts[2]},
    ]
    return {
        "checks": checks,
        "inventory": inventory,
        "duplicate_rows": duplicate_rows,
        "cohort": cohort,
        "passed": all(row["result"] == "PASS" for row in checks),
    }
