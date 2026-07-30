"""Input/output helpers that preserve clinical source order and identifiers."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd


def canonical_identifier(value: object) -> str:
    """Convert identifier to a stable canonical representation."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    text = str(value).strip()
    if re.fullmatch(r"-?\d+\.0+", text):
        return text.split(".", 1)[0]
    return text


def canonical_scalar(value: object) -> object:
    """Convert scalar to a stable canonical representation."""
    if value is None or pd.isna(value):
        return np.nan
    if isinstance(value, str):
        return value.strip()
    return value


def read_clinical_workbook(path: str | Path) -> tuple[str, pd.DataFrame]:
    """Read clinical workbook from disk and validate its basic structure."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    workbook = pd.ExcelFile(path, engine="openpyxl")
    if not workbook.sheet_names:
        raise ValueError("The clinical workbook contains no worksheets")
    sheet = workbook.sheet_names[0]
    frame = pd.read_excel(path, sheet_name=sheet, dtype=object, engine="openpyxl")
    frame.columns = [str(column).strip() for column in frame.columns]
    frame = frame.map(canonical_scalar)
    return sheet, frame


def read_table(path: str | Path) -> pd.DataFrame:
    """Read table from disk and validate its basic structure."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    compression = "gzip" if path.suffix == ".gz" else "infer"
    # sep=None handles tab/whitespace/comma while the Python engine keeps
    # participant IDs as strings.
    frame = pd.read_csv(path, sep=None, engine="python", dtype=object, compression=compression)
    frame.columns = [str(column).lstrip("#").strip() for column in frame.columns]
    return frame.map(canonical_scalar)


def identifier_columns(frame: pd.DataFrame) -> tuple[str, str]:
    """Resolve the available participant-identifier columns in source order."""
    lookup = {str(column).upper(): str(column) for column in frame.columns}
    if "IID" not in lookup:
        raise ValueError("Input table is missing IID")
    fid = lookup.get("FID", lookup["IID"])
    return fid, lookup["IID"]


def participant_keys(frame: pd.DataFrame) -> pd.Series:
    """Return stable participant keys from the configured identifier columns."""
    fid, iid = identifier_columns(frame)
    return frame[fid].map(canonical_identifier) + "::" + frame[iid].map(canonical_identifier)


def require_unique_ids(frame: pd.DataFrame, label: str) -> tuple[pd.Series, pd.Series]:
    """Require unique ids and raise a clear error when absent."""
    fid, iid = identifier_columns(frame)
    fid_values = frame[fid].map(canonical_identifier)
    iid_values = frame[iid].map(canonical_identifier)
    if (fid_values == "").any() or (iid_values == "").any():
        raise ValueError(f"{label} contains missing participant identifiers")
    keys = fid_values + "::" + iid_values
    if keys.duplicated().any():
        duplicates = keys[keys.duplicated()].unique().tolist()[:5]
        raise ValueError(f"{label} contains duplicate identifiers: {duplicates}")
    return fid_values, iid_values


def write_tsv(frame: pd.DataFrame, path: str | Path) -> None:
    """Write a deterministic TSV workflow output."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(
        path, sep="\t", index=False, compression=compression,
        lineterminator="\n", float_format="%.17g"
    )


def write_records(records: list[dict[str, object]], path: str | Path) -> None:
    """Write records to a deterministic workflow output."""
    write_tsv(pd.DataFrame.from_records(records), path)
