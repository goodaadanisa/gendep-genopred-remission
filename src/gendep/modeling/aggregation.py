"""Aggregate fixed outer-split model outputs into repeat-level results."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from .metrics import metric_bundle


@dataclass(frozen=True)
class AnalysisSpec:
    """Represent AnalysisSpec as an immutable workflow record."""
    analysis: str
    algorithm: str
    cohort_name: str
    expected_participants: int
    comparator_label: str
    combined_label: str


def read_analysis_spec(path: Path, analysis: str) -> AnalysisSpec:
    """Read analysis spec from disk and validate its basic structure."""
    frame = pd.read_csv(path, sep="\t", dtype=str)
    selected = frame.loc[frame["analysis"] == analysis]
    if len(selected) != 1:
        raise ValueError(f"Expected one configuration row for {analysis}")
    row = selected.iloc[0]
    return AnalysisSpec(
        analysis=analysis,
        algorithm=str(row["algorithm"]),
        cohort_name=str(row["cohort_name"]),
        expected_participants=int(row["expected_participants"]),
        comparator_label=str(row["comparator_label"]),
        combined_label=str(row["combined_label"]),
    )


def read_parameters(path: Path) -> dict[str, str]:
    """Read parameters from disk and validate its basic structure."""
    frame = pd.read_csv(path, sep="\t", dtype=str)
    return dict(zip(frame["parameter"], frame["value"], strict=True))


def split_directory(root: Path, analysis: str, repeat: int, fold: int) -> Path:
    """Return the canonical output directory for one outer split."""
    return root / analysis / f"repeat_{repeat:02d}" / f"fold_{fold:02d}"


def load_split_outputs(
    split_root: Path,
    spec: AnalysisSpec,
    outer_repeats: int,
    outer_folds: int,
    *,
    allow_incomplete: bool = False,
) -> dict[str, pd.DataFrame]:
    """Load split outputs as typed workflow records."""
    collections: dict[str, list[pd.DataFrame]] = {
        "predictions": [],
        "metrics": [],
        "comparison": [],
        "tuning": [],
        "coefficients": [],
        "warnings": [],
    }
    missing: list[str] = []
    for repeat in range(1, outer_repeats + 1):
        for fold in range(1, outer_folds + 1):
            directory = split_directory(split_root, spec.analysis, repeat, fold)
            success = directory / "_SUCCESS"
            if not success.exists():
                missing.append(str(directory))
                continue
            file_map = {
                "predictions": directory / "predictions.tsv.gz",
                "metrics": directory / "model_metrics.tsv",
                "comparison": directory / "comparison.tsv",
                "tuning": directory / "tuning.tsv.gz",
                "coefficients": directory / "coefficients.tsv.gz",
                "warnings": directory / "warnings.tsv",
            }
            for key, path in file_map.items():
                if not path.exists():
                    raise FileNotFoundError(path)
                frame = pd.read_csv(path, sep="\t")
                if not frame.empty:
                    collections[key].append(frame)
    if missing and not allow_incomplete:
        raise RuntimeError(f"Missing {len(missing)} completed splits; first: {missing[0]}")
    result: dict[str, pd.DataFrame] = {}
    for key, frames in collections.items():
        result[key] = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    result["missing_splits"] = pd.DataFrame({"split_directory": missing})
    return result


def _repeat_metric_rows(predictions: pd.DataFrame, spec: AnalysisSpec, epsilon: float) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for repeat, frame in predictions.groupby("outer_repeat", sort=True):
        if len(frame) != spec.expected_participants:
            raise ValueError(f"Repeat {repeat} has {len(frame)} predictions")
        if frame["participant_id"].duplicated().any():
            raise ValueError(f"Repeat {repeat} contains duplicate participants")
        outcome = frame["outcome"].to_numpy(dtype=int)
        comparator = metric_bundle(outcome, frame["comparator_probability"].to_numpy(float), epsilon)
        combined = metric_bundle(outcome, frame["combined_probability"].to_numpy(float), epsilon)
        row: dict[str, object] = {
            "analysis": spec.analysis,
            "algorithm": spec.algorithm,
            "outer_repeat": int(repeat),
            "participants": len(frame),
            "events": int(outcome.sum()),
            "non_events": int((outcome == 0).sum()),
        }
        for name, value in comparator.items():
            row[f"comparator_{name}"] = value
        for name, value in combined.items():
            row[f"combined_{name}"] = value
        row.update(
            {
                "combined_minus_comparator_AUC": combined["AUC"] - comparator["AUC"],
                "comparator_minus_combined_Brier": comparator["Brier"] - combined["Brier"],
                "comparator_minus_combined_log_loss": comparator["log_loss"] - combined["log_loss"],
                "combined_minus_comparator_calibration_intercept": (
                    combined["calibration_intercept"] - comparator["calibration_intercept"]
                ),
                "combined_minus_comparator_calibration_slope": (
                    combined["calibration_slope"] - comparator["calibration_slope"]
                ),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def aggregate_outputs(
    raw: dict[str, pd.DataFrame],
    spec: AnalysisSpec,
    outer_repeats: int,
    outer_folds: int,
    epsilon: float,
) -> dict[str, pd.DataFrame]:
    """Aggregate complete outer-split outputs into paired held-out predictions and performance summaries."""
    predictions = raw["predictions"].copy()
    required = {
        "analysis", "algorithm", "outer_repeat", "outer_fold", "participant_id",
        "outcome", "comparator_probability", "combined_probability",
    }
    missing_columns = sorted(required - set(predictions.columns))
    if missing_columns:
        raise ValueError(f"Prediction columns missing: {missing_columns}")
    if not predictions.empty:
        predictions["participant_id"] = predictions["participant_id"].astype(str)
        predictions = predictions.sort_values(
            ["outer_repeat", "outer_fold", "participant_id"], kind="stable"
        ).reset_index(drop=True)
    expected_rows = spec.expected_participants * outer_repeats
    if len(predictions) != expected_rows:
        raise ValueError(f"Expected {expected_rows} held-out rows; found {len(predictions)}")
    if predictions[["outer_repeat", "participant_id"]].duplicated().any():
        raise ValueError("A participant has more than one held-out prediction in a repeat")
    if predictions["outer_repeat"].nunique() != outer_repeats:
        raise ValueError("Outer repeats are incomplete")
    if predictions[["outer_repeat", "outer_fold"]].drop_duplicates().shape[0] != outer_repeats * outer_folds:
        raise ValueError("Outer split grid is incomplete")

    outcome_consistency = predictions.groupby("participant_id")["outcome"].nunique()
    if (outcome_consistency != 1).any():
        raise ValueError("Participant outcomes change across repeats")

    repeat_metrics = _repeat_metric_rows(predictions, spec, epsilon)
    numeric_columns = [
        column for column in repeat_metrics.columns
        if column not in {"analysis", "algorithm", "outer_repeat", "participants", "events", "non_events"}
    ]
    performance_summary = pd.DataFrame(
        [
            {
                "analysis": spec.analysis,
                "algorithm": spec.algorithm,
                "participants": spec.expected_participants,
                "outer_repeats": outer_repeats,
                **{column: float(repeat_metrics[column].mean()) for column in numeric_columns},
            }
        ]
    )

    participant_mean = (
        predictions.groupby("participant_id", as_index=False, sort=False)
        .agg(
            outcome=("outcome", "first"),
            drug=("drug", "first"),
            comparator_probability=("comparator_probability", "mean"),
            combined_probability=("combined_probability", "mean"),
            repeated_predictions=("outer_repeat", "nunique"),
        )
    )
    participant_mean.insert(0, "analysis", spec.analysis)

    tuning = raw["tuning"].copy()
    tuning_counts = pd.DataFrame()
    if not tuning.empty and "selected" in tuning.columns:
        selected_mask = tuning["selected"].astype(str).str.lower().isin({"true", "1", "t", "yes"})
        selected = tuning.loc[selected_mask].copy()
        if spec.algorithm == "elastic_net" and "alpha" in selected:
            tuning_counts = (
                selected.groupby(["model", "alpha"], as_index=False)
                .size()
                .rename(columns={"size": "outer_fit_count"})
            )
        elif spec.algorithm == "random_forest":
            tuning_counts = (
                selected.groupby(["model", "mtry", "nodesize", "ntree"], as_index=False)
                .size()
                .rename(columns={"size": "outer_fit_count"})
            )
        if not tuning_counts.empty:
            tuning_counts.insert(0, "analysis", spec.analysis)

    validation_rows = [
        ("held_out_prediction_rows", len(predictions), expected_rows),
        ("outer_repeats", predictions["outer_repeat"].nunique(), outer_repeats),
        ("outer_split_combinations", predictions[["outer_repeat", "outer_fold"]].drop_duplicates().shape[0], outer_repeats * outer_folds),
        ("participants_per_repeat", int(predictions.groupby("outer_repeat").size().min()), spec.expected_participants),
        ("participant_predictions_per_repeat", int(predictions.groupby(["outer_repeat", "participant_id"]).size().max()), 1),
        ("outcome_consistency", int(outcome_consistency.max()), 1),
        ("repeat_metric_rows", len(repeat_metrics), outer_repeats),
    ]
    validation = pd.DataFrame(validation_rows, columns=["check", "observed", "expected"])
    validation["result"] = np.where(validation["observed"] == validation["expected"], "PASS", "FAIL")

    return {
        "predictions": predictions,
        "repeat_metrics": repeat_metrics,
        "performance_summary": performance_summary,
        "participant_mean_predictions": participant_mean,
        "tuning": tuning,
        "tuning_selection_counts": tuning_counts,
        "coefficients": raw["coefficients"],
        "warnings": raw["warnings"],
        "aggregation_validation": validation,
        "missing_splits": raw["missing_splits"],
    }
