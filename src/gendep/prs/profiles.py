"""Independent validation helpers for GenoPred profiles and aggregate outputs."""

from __future__ import annotations

import gzip
import math
import re
from pathlib import Path
from typing import TextIO

import numpy as np
import pandas as pd


POPULATIONS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_text(path: str | Path) -> TextIO:
    """Open text while handling plain-text and compressed files."""
    file_path = Path(path)
    if _is_gzip(file_path):
        return gzip.open(file_path, "rt", encoding="utf-8", errors="replace")
    return file_path.open("rt", encoding="utf-8", errors="replace")


def read_table(path: str | Path) -> pd.DataFrame:
    """Read tab- or whitespace-delimited output while preserving IDs as strings."""
    file_path = Path(path)
    with open_text(file_path) as handle:
        first = handle.readline()
    separator = "\t" if "\t" in first else r"\s+"
    frame = pd.read_csv(file_path, sep=separator, compression="infer", dtype=str, engine="python")
    frame.columns = [re.sub(r"^#", "", str(value)) for value in frame.columns]
    return frame


def participant_keys(frame: pd.DataFrame) -> pd.Series:
    """Return FID/IID participant keys in source order."""
    columns = {name.upper(): name for name in frame.columns}
    if "FID" not in columns or "IID" not in columns:
        raise ValueError("Participant table must contain FID and IID")
    return frame[columns["FID"]].astype(str) + "::" + frame[columns["IID"]].astype(str)


def detect_score_column(frame: pd.DataFrame) -> str:
    """Detect score column and fail on ambiguous input."""
    id_columns = {name for name in frame.columns if name.upper() in {"FID", "IID"}}
    score_columns = [name for name in frame.columns if name not in id_columns]
    if len(score_columns) != 1:
        raise ValueError(f"Expected exactly one profile score column; found {score_columns}")
    return score_columns[0]


def count_data_rows(path: str | Path) -> int:
    """Count data rows without loading unnecessary participant-level data."""
    count = 0
    header_seen = False
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            if not header_seen:
                header_seen = True
                continue
            if line.startswith("#"):
                continue
            count += 1
    if not header_seen:
        raise ValueError(f"No header found: {path}")
    return count


def _split(line: str) -> list[str]:
    return line.rstrip("\r\n").split("\t") if "\t" in line else line.split()


def count_score_weights(path: str | Path) -> tuple[int, int, str]:
    """Count all and non-zero score rows without loading the full file."""
    header: list[str] | None = None
    weight_index: int | None = None
    rows = nonzero = 0
    metadata = re.compile(r"^(#?CHR|CHROM|BP|POS|POSITION|SNP|RSID|ID|A1|A2|REF|ALT|N)$", re.I)
    preferred = re.compile(r"BETA|WEIGHT|EFFECT|SCORE|PHI", re.I)
    with open_text(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("##"):
                continue
            if header is None:
                header = _split(line)
                candidates = [
                    index for index, name in enumerate(header)
                    if preferred.search(name) and not metadata.search(name)
                ]
                if not candidates:
                    candidates = [
                        index for index, name in enumerate(header)
                        if not metadata.search(name)
                    ]
                if not candidates:
                    raise ValueError(f"No candidate weight column found: {path}")
                weight_index = candidates[-1]
                continue
            fields = _split(line)
            if len(fields) != len(header):
                raise ValueError(f"Malformed score row in {path}")
            assert weight_index is not None
            try:
                value = float(fields[weight_index])
            except ValueError as error:
                raise ValueError(f"Non-numeric score weight in {path}") from error
            if not math.isfinite(value):
                raise ValueError(f"Non-finite score weight in {path}")
            rows += 1
            if value != 0:
                nonzero += 1
    if header is None or weight_index is None:
        raise ValueError(f"No score header found: {path}")
    return rows, nonzero, header[weight_index]


def validate_profile(
    profile_path: Path,
    expected_keys: pd.Series,
    merged_values: pd.Series | None = None,
) -> tuple[dict[str, object], np.ndarray]:
    """Validate one participant PRS profile and its FID/IID alignment against the target sample order."""
    profile = read_table(profile_path)
    keys = participant_keys(profile)
    score_column = detect_score_column(profile)
    values = pd.to_numeric(profile[score_column], errors="coerce").to_numpy(dtype=float)
    key_to_index = {key: index for index, key in enumerate(keys)}
    aligned_index = [key_to_index.get(key) for key in expected_keys]
    missing = sum(index is None for index in aligned_index)
    aligned = np.array(
        [values[index] if index is not None else np.nan for index in aligned_index],
        dtype=float,
    )
    matrix_match = True
    if merged_values is not None:
        merged = pd.to_numeric(merged_values, errors="coerce").to_numpy(dtype=float)
        matrix_match = bool(np.allclose(aligned, merged, rtol=0, atol=1e-12, equal_nan=False))
    record = {
        "profile_file": str(profile_path.resolve()),
        "score_column": score_column,
        "participant_rows": len(profile),
        "unique_participants": int(keys.nunique()),
        "duplicate_participants": int(keys.duplicated().sum()),
        "missing_expected_participants": int(missing),
        "unexpected_participants": int((~keys.isin(expected_keys)).sum()),
        "exact_psam_order": bool(keys.tolist() == expected_keys.tolist()),
        "missing_scores": int(np.isnan(aligned).sum()),
        "nonfinite_scores": int((~np.isfinite(aligned)).sum()),
        "standard_deviation": float(np.nanstd(aligned, ddof=1)),
        "merged_matrix_match": matrix_match,
    }
    record["validation"] = "PASS" if (
        record["participant_rows"] == len(expected_keys)
        and record["unique_participants"] == len(expected_keys)
        and record["duplicate_participants"] == 0
        and record["missing_expected_participants"] == 0
        and record["unexpected_participants"] == 0
        and record["missing_scores"] == 0
        and record["nonfinite_scores"] == 0
        and record["standard_deviation"] > 0
        and matrix_match
    ) else "FAIL"
    return record, aligned


def validate_ancestry(genopred_output: Path, expected_keys: pd.Series) -> tuple[dict[str, object], pd.DataFrame]:
    """Validate ancestry probabilities and strict population assignments for the target cohort."""
    path = genopred_output / "GENDEP" / "ancestry" / "GENDEP.Ancestry.model_pred"
    frame = read_table(path)
    keys = participant_keys(frame)
    missing_columns = sorted(set(POPULATIONS).difference(frame.columns))
    if missing_columns:
        raise ValueError(f"Ancestry output is missing: {', '.join(missing_columns)}")
    probabilities = frame.loc[:, POPULATIONS].apply(pd.to_numeric, errors="coerce")
    max_population = probabilities.idxmax(axis=1)
    max_probability = probabilities.max(axis=1)
    probability_sum = probabilities.sum(axis=1)

    keep_dir = genopred_output / "GENDEP" / "ancestry" / "keep_files" / "model_based"
    keep_counts: dict[str, int] = {}
    key_set = set(keys)
    for population in POPULATIONS:
        keep_path = keep_dir / f"{population}.keep"
        if keep_path.is_file() and keep_path.stat().st_size > 0:
            keep = pd.read_csv(keep_path, sep=r"\s+", header=None, dtype=str)
            if keep.shape[1] < 2:
                raise ValueError(f"Keep file lacks FID/IID: {keep_path}")
            keep_keys = set(keep.iloc[:, 0] + "::" + keep.iloc[:, 1])
            if not keep_keys.issubset(key_set):
                raise ValueError(f"Unexpected participant in keep file: {keep_path}")
            keep_counts[population] = len(keep_keys)
        else:
            keep_counts[population] = 0

    record: dict[str, object] = {
        "ancestry_file": str(path.resolve()),
        "participant_rows": len(frame),
        "unique_participants": int(keys.nunique()),
        "duplicate_participants": int(keys.duplicated().sum()),
        "missing_expected_participants": int((~expected_keys.isin(keys)).sum()),
        "unexpected_participants": int((~keys.isin(expected_keys)).sum()),
        "exact_psam_order": bool(keys.tolist() == expected_keys.tolist()),
        "probability_missing_values": int(probabilities.isna().sum().sum()),
        "probability_nonfinite_values": int((~np.isfinite(probabilities.to_numpy(dtype=float))).sum()),
        "probability_sum_outside_0_01": int((np.abs(probability_sum - 1) > 0.01).sum()),
        "maximum_population_eur": int((max_population == "EUR").sum()),
        "minimum_maximum_probability": float(max_probability.min()),
        "mean_maximum_probability": float(max_probability.mean()),
        "maximum_maximum_probability": float(max_probability.max()),
        "strict_eur_keep_participants": int(keep_counts["EUR"]),
    }
    # Compute the exact union separately to keep the record independent of row order.
    keep_union: set[str] = set()
    for population in POPULATIONS:
        keep_path = keep_dir / f"{population}.keep"
        if keep_path.is_file() and keep_path.stat().st_size > 0:
            keep = pd.read_csv(keep_path, sep=r"\s+", header=None, dtype=str)
            keep_union.update(set(keep.iloc[:, 0] + "::" + keep.iloc[:, 1]))
    record["participants_outside_all_model_keeps"] = len(key_set.difference(keep_union))
    summary = pd.DataFrame(
        {
            "predicted_population": max_population,
            "maximum_probability": max_probability,
        }
    ).groupby("predicted_population", as_index=False).agg(
        participants=("maximum_probability", "size"),
        mean_maximum_probability=("maximum_probability", "mean"),
        minimum_maximum_probability=("maximum_probability", "min"),
        maximum_maximum_probability=("maximum_probability", "max"),
    )
    return record, summary


def validate_projected_pcs(genopred_output: Path, expected_keys: pd.Series) -> dict[str, object]:
    """Validate the six projected ancestry principal components and participant order."""
    directory = genopred_output / "GENDEP" / "pcs" / "projected" / "TRANS"
    preferred = directory / "GENDEP-TRANS.profiles"
    candidates = [preferred] if preferred.is_file() else sorted(
        path for path in directory.rglob("*")
        if path.is_file() and re.search(r"(profiles|sscore|eigenvec)$", path.name, re.I)
    )
    if not candidates:
        raise FileNotFoundError(f"Projected-PC output was not found beneath {directory}")

    valid_candidates: list[tuple[Path, pd.DataFrame, pd.Series, list[str]]] = []
    for path in candidates:
        frame = read_table(path)
        keys = participant_keys(frame)
        id_columns = {name for name in frame.columns if name.upper() in {"FID", "IID"}}
        pc_columns = [name for name in frame.columns if name not in id_columns]
        if len(pc_columns) == 6:
            valid_candidates.append((path, frame, keys, pc_columns))
    if len(valid_candidates) != 1:
        names = [str(item[0]) for item in valid_candidates]
        raise ValueError(f"Expected exactly one six-PC output; found {names}")

    path, frame, keys, pc_columns = valid_candidates[0]
    numeric = frame[pc_columns].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    return {
        "pc_file": str(path.resolve()),
        "participant_rows": len(frame),
        "unique_participants": int(keys.nunique()),
        "duplicate_participants": int(keys.duplicated().sum()),
        "missing_expected_participants": int((~expected_keys.isin(keys)).sum()),
        "unexpected_participants": int((~keys.isin(expected_keys)).sum()),
        "exact_psam_order": bool(keys.tolist() == expected_keys.tolist()),
        "pc_columns": len(pc_columns),
        "missing_pc_values": int(np.isnan(numeric).sum()),
        "nonfinite_pc_values": int((~np.isfinite(numeric)).sum()),
    }
