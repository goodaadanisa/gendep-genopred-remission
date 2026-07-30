from __future__ import annotations

# ---- test_clinical_stage.py ----

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from gendep.clinical.audit import audit_clinical, exact_duplicate_groups
from gendep.clinical.integration import integrate_analysis_data
from gendep.clinical.predictors import build_predictor_policy
from gendep.clinical.schema import (
    BINARY_CATEGORICAL,
    CLINICAL_COLUMNS,
    CONTINUOUS_CLINICAL,
    COUNT_NUMERIC,
    EXACT_DUPLICATE_GROUPS,
    MODEL_PREDICTOR_SETS,
    OUTCOME_COLUMN,
    PC_COLUMNS,
    PRS_COLUMNS,
    TREATMENT_COLUMN,
)
from gendep.clinical.validation import validate_analysis_bases


def synthetic_clinical(n: int = 430) -> pd.DataFrame:
    data: dict[str, object] = {}
    participant_ids = [f"P{index:04d}" for index in range(1, n + 1)]
    data["subjectid"] = [f"S{index:04d}" for index in range(1, n + 1)]
    data["Row.names"] = participant_ids
    data["bloodsampleid.x"] = participant_ids
    for column_index, column in enumerate(CLINICAL_COLUMNS[3:], start=3):
        rng = np.random.default_rng(20260715 + column_index * 37)
        if column == TREATMENT_COLUMN:
            values = np.array([1] * 210 + [2] * 220)
        elif column == OUTCOME_COLUMN:
            values = np.array([0] * 264 + [1] * 166)
        elif column in BINARY_CATEGORICAL:
            values = rng.integers(0, 2, size=n)
            values[0], values[1] = 0, 1
        elif column == "hamd16wk0":
            values = rng.integers(0, 4, size=n).astype(float)
            values[:10] += 0.5
        elif column in CONTINUOUS_CLINICAL:
            values = rng.normal(loc=column_index / 5, scale=2.0, size=n)
        elif column in COUNT_NUMERIC:
            values = rng.integers(0, 10, size=n)
        elif column.startswith("madrs"):
            values = rng.integers(0, 7, size=n)
        elif column.startswith("hamd"):
            values = rng.integers(0, 5, size=n)
        elif column.startswith("bdi"):
            values = rng.integers(0, 4, size=n)
        elif column.startswith("k"):
            values = rng.integers(0, 4, size=n)
        else:
            values = rng.integers(0, 8, size=n)
        data[column] = values
    frame = pd.DataFrame(data, columns=CLINICAL_COLUMNS)
    for group in EXACT_DUPLICATE_GROUPS:
        for duplicate in group[1:]:
            frame[duplicate] = frame[group[0]].copy()
    return frame


def synthetic_genetic_inputs(clinical: pd.DataFrame):
    ids = clinical["Row.names"].astype(str).tolist()
    n = len(ids)
    order = np.random.default_rng(11).permutation(n)
    prs = pd.DataFrame({"FID": ids, "IID": ids})
    for index, column in enumerate(PRS_COLUMNS, start=1):
        prs[column] = np.linspace(-2, 2, n) + index / 10
    prs = prs.iloc[order].reset_index(drop=True)

    pcs = pd.DataFrame({"FID": ids, "IID": ids})
    for index in range(1, 7):
        pcs[f"PC{index}"] = np.linspace(-1, 1, n) * index
    pcs = pcs.iloc[order[::-1]].reset_index(drop=True)

    ancestry = pd.DataFrame({"FID": ids, "IID": ids})
    ancestry["AFR"] = 0.002
    ancestry["AMR"] = 0.002
    ancestry["CSA"] = 0.002
    ancestry["EAS"] = 0.002
    ancestry["EUR"] = np.linspace(0.98, 0.995, n)
    ancestry["MID"] = 1.0 - ancestry[["AFR", "AMR", "CSA", "EAS", "EUR"]].sum(axis=1)
    ancestry = ancestry.iloc[order].reset_index(drop=True)
    keep = set(ids[:418])
    return prs, pcs, ancestry, keep


def test_static_clinical_schema_and_predictor_counts() -> None:
    assert len(CLINICAL_COLUMNS) == 139
    assert len(MODEL_PREDICTOR_SETS["primary_comprehensive_clinical"]) == 134
    assert len(MODEL_PREDICTOR_SETS["primary_comprehensive_combined"]) == 142
    assert len(MODEL_PREDICTOR_SETS["secondary_prs_focused"]) == 15
    assert len(MODEL_PREDICTOR_SETS["sensitivity_summary_clinical"]) == 52
    assert len(MODEL_PREDICTOR_SETS["sensitivity_summary_combined"]) == 60


def test_audit_policy_integration_and_independent_validation() -> None:
    clinical = synthetic_clinical()
    assert set(exact_duplicate_groups(clinical)) == set(EXACT_DUPLICATE_GROUPS)
    audit = audit_clinical(clinical)
    assert audit["passed"]
    policy = build_predictor_policy(clinical)
    assert policy["passed"]

    prs, pcs, ancestry, keep = synthetic_genetic_inputs(clinical)
    integrated = integrate_analysis_data(clinical, prs, pcs, keep, ancestry)
    primary = integrated["primary"]
    strict = integrated["strict_eur"]
    assert primary.shape == (430, 154)
    assert strict.shape == (418, 154)
    assert primary["Row.names"].tolist() == clinical["Row.names"].tolist()
    checks = validate_analysis_bases(primary, strict, source_clinical=clinical)
    assert all(row["result"] == "PASS" for row in checks)


# ---- test_resampling.py ----

from collections import Counter

from gendep.resampling import FoldDesign, Participant, generate_cohort_folds


def make_primary() -> list[Participant]:
    counts = {
        ("0", "1"): 138,
        ("0", "2"): 126,
        ("1", "1"): 72,
        ("1", "2"): 94,
    }
    participants: list[Participant] = []
    index = 1
    for (outcome, treatment), count in counts.items():
        for _ in range(count):
            participants.append(Participant(f"P{index:04d}", outcome, treatment))
            index += 1
    return participants


def test_primary_manifest_dimensions_and_balance() -> None:
    participants = make_primary()
    design = FoldDesign()
    outer, inner, balance = generate_cohort_folds(
        participants,
        cohort_name="primary_430",
        seed_offset=0,
        design=design,
    )
    assert len(outer) == 4300
    assert len(inner) == 38700
    assert list(outer[0]) == [
        "cohort",
        "outer_repeat",
        "participant_id",
        "outcome",
        "drug",
        "outer_fold",
        "outer_seed",
        "joint_stratum",
    ]
    assert list(inner[0]) == [
        "cohort",
        "outer_repeat",
        "held_out_outer_fold",
        "participant_id",
        "outcome",
        "drug",
        "inner_fold",
        "inner_seed",
        "joint_stratum",
    ]
    assert list(balance[0]) == [
        "cohort",
        "split_level",
        "outer_repeat",
        "outer_fold",
        "inner_fold",
        "role",
        "participants",
        "outcome_0",
        "outcome_1",
        "drug_1",
        "drug_2",
        "outcome_0_drug_1",
        "outcome_0_drug_2",
        "outcome_1_drug_1",
        "outcome_1_drug_2",
        "all_four_joint_strata_present",
    ]
    assert set(Counter(row["participant_id"] for row in outer).values()) == {10}
    assert set(Counter(row["participant_id"] for row in inner).values()) == {90}
    outer_valid = [row for row in balance if row["split_level"] == "outer" and row["role"] == "validation"]
    assert {row["participants"] for row in outer_valid} == {43}
    assert all(
        row["all_four_joint_strata_present"] == "TRUE"
        for row in outer_valid
    )


def test_resampling_is_deterministic() -> None:
    participants = make_primary()
    design = FoldDesign()
    first = generate_cohort_folds(participants, cohort_name="primary_430", seed_offset=0, design=design)[0]
    second = generate_cohort_folds(participants, cohort_name="primary_430", seed_offset=0, design=design)[0]
    assert first == second
