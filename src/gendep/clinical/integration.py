"""Merge clinical, PRS and ancestry outputs by validated participant ID."""

from __future__ import annotations

from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd

from .io import canonical_identifier, identifier_columns, require_unique_ids
from .schema import (
    ANCESTRY_POPULATIONS,
    CLINICAL_COLUMNS,
    EUR_INDICATOR,
    OUTCOME_COLUMN,
    PC_COLUMNS,
    PRS_COLUMNS,
    TREATMENT_COLUMN,
)


def read_keep_file(path: str | Path) -> set[str]:
    """Read keep file from disk and validate its basic structure."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep=r"\s+", header=None, dtype=str, comment="#")
    if frame.shape[1] < 2:
        raise ValueError("EUR keep file must contain FID and IID")
    frame = frame.iloc[:, :2]
    frame.columns = ["FID", "IID"]
    if len(frame) and str(frame.iloc[0, 0]).lstrip("#").upper() == "FID" and str(frame.iloc[0, 1]).upper() == "IID":
        frame = frame.iloc[1:].reset_index(drop=True)
    _, iid = require_unique_ids(frame, "EUR keep file")
    fid = frame["FID"].map(canonical_identifier)
    if not (fid == iid).all():
        raise ValueError("EUR keep file contains non-identical FID and IID")
    return set(iid)


def _aligned_by_iid(frame: pd.DataFrame, clinical_ids: pd.Series, label: str) -> pd.DataFrame:
    fid, iid = require_unique_ids(frame, label)
    if not (fid == iid).all():
        raise ValueError(f"{label} contains non-identical FID and IID")
    if set(iid) != set(clinical_ids):
        missing = sorted(set(clinical_ids) - set(iid))[:5]
        extra = sorted(set(iid) - set(clinical_ids))[:5]
        raise ValueError(f"{label} participant set mismatch; missing={missing}, extra={extra}")
    working = frame.copy()
    working["__IID_CANONICAL__"] = iid
    return working.set_index("__IID_CANONICAL__").loc[clinical_ids.tolist()].reset_index(drop=True)


def _numeric_complete(frame: pd.DataFrame, columns: list[str] | tuple[str, ...], label: str) -> pd.DataFrame:
    numeric = frame.loc[:, list(columns)].apply(pd.to_numeric, errors="coerce")
    if numeric.isna().any().any() or not np.isfinite(numeric.to_numpy(dtype=float)).all():
        raise ValueError(f"{label} contains missing or non-finite values")
    return numeric


def integrate_analysis_data(
    clinical: pd.DataFrame,
    prs: pd.DataFrame,
    pcs: pd.DataFrame,
    eur_keep_ids: set[str],
    ancestry: pd.DataFrame | None = None,
) -> dict[str, object]:
    """Merge clinical data, PRS, ancestry PCs and strict-EUR membership by validated participant identifiers."""
    if tuple(clinical.columns) != CLINICAL_COLUMNS:
        raise ValueError("Clinical columns or order do not match the fixed 139-column schema")
    clinical_ids = clinical["Row.names"].map(canonical_identifier)
    blood_ids = clinical["bloodsampleid.x"].map(canonical_identifier)
    if (clinical_ids == "").any() or clinical_ids.duplicated().any():
        raise ValueError("Clinical Row.names are missing or duplicated")
    if not (clinical_ids == blood_ids).all():
        raise ValueError("Clinical Row.names and bloodsampleid.x differ")

    prs_aligned = _aligned_by_iid(prs, clinical_ids, "PRS matrix")
    missing_prs = [column for column in PRS_COLUMNS if column not in prs_aligned.columns]
    if missing_prs:
        raise ValueError(f"PRS matrix is missing columns: {missing_prs}")
    prs_values = _numeric_complete(prs_aligned, PRS_COLUMNS, "PRS matrix")

    pcs_aligned = _aligned_by_iid(pcs, clinical_ids, "Projected-PC file")
    fid_column, iid_column = identifier_columns(pcs_aligned)
    source_pc_columns = [
        column for column in pcs_aligned.columns
        if column not in {fid_column, iid_column, "__IID_CANONICAL__"}
    ]
    if len(source_pc_columns) != 6:
        raise ValueError(f"Expected six projected PCs, found {source_pc_columns}")
    pc_values = _numeric_complete(pcs_aligned, source_pc_columns, "Projected-PC file")
    pc_values.columns = list(PC_COLUMNS)

    ancestry_summary: dict[str, object] = {}
    if ancestry is not None:
        ancestry_aligned = _aligned_by_iid(ancestry, clinical_ids, "Ancestry probabilities")
        lookup = {str(column).upper(): str(column) for column in ancestry_aligned.columns}
        missing = [population for population in ANCESTRY_POPULATIONS if population not in lookup]
        if missing:
            raise ValueError(f"Ancestry file is missing probability columns: {missing}")
        probability_columns = [lookup[population] for population in ANCESTRY_POPULATIONS]
        probabilities = _numeric_complete(ancestry_aligned, probability_columns, "Ancestry probabilities")
        sums = probabilities.sum(axis=1)
        maximum = probabilities.idxmax(axis=1).map({lookup[p]: p for p in ANCESTRY_POPULATIONS})
        ancestry_summary = {
            "probability_sum_failures": int((np.abs(sums - 1.0) > 0.01).sum()),
            "maximum_population_eur": int((maximum == "EUR").sum()),
            "mean_maximum_probability": float(probabilities.max(axis=1).mean()),
            "minimum_maximum_probability": float(probabilities.max(axis=1).min()),
            "maximum_maximum_probability": float(probabilities.max(axis=1).max()),
        }
        if ancestry_summary["probability_sum_failures"]:
            raise ValueError("Ancestry probabilities do not sum to one within 0.01")

    if not eur_keep_ids <= set(clinical_ids):
        raise ValueError("EUR keep file contains participants absent from clinical data")
    eur_indicator = clinical_ids.isin(eur_keep_ids).astype(int)

    primary = pd.concat(
        [clinical.reset_index(drop=True), prs_values.reset_index(drop=True), pc_values.reset_index(drop=True)],
        axis=1,
    )
    primary[EUR_INDICATOR] = eur_indicator.to_numpy()
    strict_eur = primary.loc[primary[EUR_INDICATOR] == 1].reset_index(drop=True)

    outcome_counts = Counter(pd.to_numeric(primary[OUTCOME_COLUMN]).astype(int))
    treatment_counts = Counter(pd.to_numeric(primary[TREATMENT_COLUMN]).astype(int))
    summary = {
        "primary_rows": len(primary),
        "strict_eur_rows": len(strict_eur),
        "analysis_columns": len(primary.columns),
        "clinical_columns": len(CLINICAL_COLUMNS),
        "prs_columns": len(PRS_COLUMNS),
        "pc_columns": len(PC_COLUMNS),
        "eur_indicator_columns": 1,
        "non_remitters": outcome_counts[0],
        "remitters": outcome_counts[1],
        "treatment_group_1": treatment_counts[1],
        "treatment_group_2": treatment_counts[2],
        **ancestry_summary,
    }
    return {"primary": primary, "strict_eur": strict_eur, "summary": summary}
