#!/usr/bin/env python3
"""Analyse formal calibration and PRS-coefficient stability from fixed models.

Purpose
-------
Reproduce the final formal repeat-wise calibration analysis and describe the
stability of the eight PRS coefficients across the 100 primary combined-model
outer fits. This command reads fixed held-out predictions and fitted
coefficients; it never refits or retunes a predictive model.

Inputs
------
The primary aggregated prediction, coefficient and tuning files produced by
the ``gendep aggregate-models`` command.

Outputs
-------
Repeat calibration estimates, NumPy-seeded participant-bootstrap intervals,
PRS coefficient long data, stability summaries and validation records.

Usage
-----
gendep model-diagnostics \
  --predictions work/modelling/aggregated/primary/all_held_out_predictions.tsv.gz \
  --coefficients work/modelling/aggregated/primary/all_coefficients.tsv.gz \
  --tuning work/modelling/aggregated/primary/all_tuning_results.tsv.gz \
  --output-dir work/modelling/diagnostics/primary
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from numba import njit, prange
except ImportError:  # pragma: no cover - slow but functional fallback
    def njit(*_args, **_kwargs):
        """Return the original function when Numba is unavailable."""
        def decorator(function):
            """Return the undecorated function in the no-Numba fallback."""
            return function
        return decorator

    def prange(value):
        """Use the built-in range in the no-Numba fallback."""
        return range(value)

ROOT = Path(__file__).resolve().parents[3]
PRS_ORDER = ["PRS_MDD", "PRS_ANX", "PRS_BIP", "PRS_SCZ", "PRS_NEUR", "PRS_INSOM", "PRS_SWB", "PRS_EA"]


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the model-diagnostics command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--coefficients", type=Path, required=True)
    parser.add_argument("--tuning", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--bootstrap-replicates", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--comparator-label", default="clinical")
    parser.add_argument("--combined-label", default="combined")
    parser.add_argument("--expected-participants", type=int, default=430)
    parser.add_argument("--expected-events", type=int, default=166)
    parser.add_argument("--validate-fixed", action="store_true")
    return parser.parse_args()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    """Write a deterministic TSV workflow output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        sep="\t",
        index=False,
        compression="gzip" if path.suffix == ".gz" else None,
        lineterminator="\n",
    )


@njit(cache=False)
def expit_scalar(value: float) -> float:
    """Evaluate the logistic function while avoiding numerical overflow."""
    if value >= 0.0:
        z = math.exp(-value)
        return 1.0 / (1.0 + z)
    z = math.exp(value)
    return z / (1.0 + z)


@njit(cache=False)
def fit_calibration_in_large(y: np.ndarray, x: np.ndarray) -> tuple[float, int]:
    """Estimate calibration-in-the-large with the slope fixed at one."""
    intercept = 0.0
    converged = 0
    for _ in range(60):
        score = 0.0
        information = 0.0
        for index in range(y.shape[0]):
            probability = expit_scalar(intercept + x[index])
            score += y[index] - probability
            information += probability * (1.0 - probability)
        if information <= 1e-14:
            break
        step = score / information
        intercept += step
        if abs(step) < 1e-11:
            converged = 1
            break
    return intercept, converged


@njit(cache=False)
def fit_calibration_slope(y: np.ndarray, x: np.ndarray) -> tuple[float, float, int]:
    """Estimate calibration intercept and slope from predicted logits."""
    intercept = 0.0
    slope = 1.0
    converged = 0
    for _ in range(80):
        score_0 = 0.0
        score_1 = 0.0
        info_00 = 0.0
        info_01 = 0.0
        info_11 = 0.0
        for index in range(y.shape[0]):
            linear_predictor = intercept + slope * x[index]
            probability = expit_scalar(linear_predictor)
            residual = y[index] - probability
            weight = probability * (1.0 - probability)
            score_0 += residual
            score_1 += x[index] * residual
            info_00 += weight
            info_01 += weight * x[index]
            info_11 += weight * x[index] * x[index]
        determinant = info_00 * info_11 - info_01 * info_01
        if abs(determinant) <= 1e-14:
            break
        step_0 = (info_11 * score_0 - info_01 * score_1) / determinant
        step_1 = (-info_01 * score_0 + info_00 * score_1) / determinant
        intercept += step_0
        slope += step_1
        if max(abs(step_0), abs(step_1)) < 1e-11:
            converged = 1
            break
    return intercept, slope, converged


@njit(cache=False)
def repeat_calibration(y: np.ndarray, logits: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Calculate formal calibration estimates for each outer repeat."""
    repeats = logits.shape[1]
    in_large = np.empty(repeats, dtype=np.float64)
    joint_intercept = np.empty(repeats, dtype=np.float64)
    slopes = np.empty(repeats, dtype=np.float64)
    failures = 0
    for repeat_index in range(repeats):
        intercept, ok_1 = fit_calibration_in_large(y, logits[:, repeat_index])
        joint, slope, ok_2 = fit_calibration_slope(y, logits[:, repeat_index])
        in_large[repeat_index] = intercept
        joint_intercept[repeat_index] = joint
        slopes[repeat_index] = slope
        if ok_1 == 0 or ok_2 == 0:
            failures += 1
    return in_large, joint_intercept, slopes, failures


@njit(parallel=True, cache=False)
def bootstrap_calibration(
    y: np.ndarray,
    comparator_logits: np.ndarray,
    combined_logits: np.ndarray,
    sampled_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Bootstrap participant-level calibration contrasts using fixed predictions."""
    replicates = sampled_indices.shape[0]
    repeats = comparator_logits.shape[1]
    output = np.empty((replicates, 6), dtype=np.float64)
    failure_counts = np.zeros(replicates, dtype=np.int64)
    for bootstrap_index in prange(replicates):
        indices = sampled_indices[bootstrap_index]
        sampled_y = y[indices]
        comparator_intercept_sum = 0.0
        combined_intercept_sum = 0.0
        comparator_slope_sum = 0.0
        combined_slope_sum = 0.0
        failures = 0
        for repeat_index in range(repeats):
            comparator_x = comparator_logits[indices, repeat_index]
            combined_x = combined_logits[indices, repeat_index]
            comparator_intercept, ok_1 = fit_calibration_in_large(sampled_y, comparator_x)
            _, comparator_slope, ok_2 = fit_calibration_slope(sampled_y, comparator_x)
            combined_intercept, ok_3 = fit_calibration_in_large(sampled_y, combined_x)
            _, combined_slope, ok_4 = fit_calibration_slope(sampled_y, combined_x)
            comparator_intercept_sum += comparator_intercept
            combined_intercept_sum += combined_intercept
            comparator_slope_sum += comparator_slope
            combined_slope_sum += combined_slope
            if ok_1 == 0 or ok_2 == 0 or ok_3 == 0 or ok_4 == 0:
                failures += 1
        comparator_intercept_mean = comparator_intercept_sum / repeats
        combined_intercept_mean = combined_intercept_sum / repeats
        comparator_slope_mean = comparator_slope_sum / repeats
        combined_slope_mean = combined_slope_sum / repeats
        output[bootstrap_index, 0] = comparator_intercept_mean
        output[bootstrap_index, 1] = combined_intercept_mean
        output[bootstrap_index, 2] = combined_intercept_mean - comparator_intercept_mean
        output[bootstrap_index, 3] = comparator_slope_mean
        output[bootstrap_index, 4] = combined_slope_mean
        output[bootstrap_index, 5] = combined_slope_mean - comparator_slope_mean
        failure_counts[bootstrap_index] = failures
    return output, failure_counts


def reshape_predictions(
    path: Path,
    expected_participants: int,
    expected_events: int,
) -> tuple[pd.DataFrame, np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Reshape paired model predictions into a calibration-ready long table."""
    frame = pd.read_csv(path, sep="\t")
    required = {
        "participant_id", "outcome", "comparator_probability",
        "combined_probability", "outer_repeat", "outer_fold",
    }
    if not required.issubset(frame.columns):
        raise ValueError(f"Prediction columns missing: {sorted(required - set(frame.columns))}")
    expected_rows = expected_participants * 10
    if len(frame) != expected_rows:
        raise ValueError(f"Expected {expected_rows} predictions; observed {len(frame)}")
    if frame[["participant_id", "outer_repeat"]].duplicated().any():
        raise ValueError("Duplicate participant-repeat predictions")
    if set(frame["outer_repeat"].unique()) != set(range(1, 11)):
        raise ValueError("Outer repeats are not exactly 1 through 10")
    for column in ["comparator_probability", "combined_probability"]:
        values = frame[column].to_numpy(float)
        if not np.isfinite(values).all() or np.any(values <= 0.0) or np.any(values >= 1.0):
            raise ValueError(f"{column} contains invalid probabilities")

    participants = sorted(frame["participant_id"].astype(str).unique())
    if len(participants) != expected_participants:
        raise ValueError("Unexpected participant count")
    index = {participant: position for position, participant in enumerate(participants)}
    outcome = np.full(expected_participants, -1, dtype=np.int64)
    comparator = np.empty((expected_participants, 10), dtype=np.float64)
    combined = np.empty((expected_participants, 10), dtype=np.float64)
    for row in frame.itertuples(index=False):
        position = index[str(row.participant_id)]
        repeat_index = int(row.outer_repeat) - 1
        if outcome[position] == -1:
            outcome[position] = int(row.outcome)
        elif outcome[position] != int(row.outcome):
            raise ValueError(f"Inconsistent outcome for {row.participant_id}")
        comparator[position, repeat_index] = float(row.comparator_probability)
        combined[position, repeat_index] = float(row.combined_probability)
    if int((outcome == 1).sum()) != expected_events:
        raise ValueError("Unexpected event count")
    epsilon = 1e-15
    comparator_bounded = np.clip(comparator, epsilon, 1 - epsilon)
    combined_bounded = np.clip(combined, epsilon, 1 - epsilon)
    comparator_logits = np.log(comparator_bounded / (1 - comparator_bounded))
    combined_logits = np.log(combined_bounded / (1 - combined_bounded))
    return frame, outcome, comparator_logits, combined_logits, participants


def analyse_coefficients(
    coefficient_path: Path,
    tuning_path: Path,
    combined_label: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Summarise PRS coefficient selection, sign and magnitude across fits."""
    coefficients = pd.read_csv(coefficient_path, sep="\t")
    required = {"outer_repeat", "outer_fold", "model", "term", "coefficient", "PRS_term", "nonzero"}
    if not required.issubset(coefficients.columns):
        raise ValueError(f"Coefficient columns missing: {sorted(required - set(coefficients.columns))}")
    prs = coefficients.loc[
        (coefficients["model"] == combined_label)
        & coefficients["term"].isin(PRS_ORDER)
    ].copy()
    if len(prs) != 800:
        raise ValueError(f"Expected 800 PRS coefficient rows; observed {len(prs)}")
    if prs[["outer_repeat", "outer_fold", "term"]].duplicated().any():
        raise ValueError("Duplicate PRS coefficient within an outer fit")

    tuning = pd.read_csv(tuning_path, sep="\t")
    selected = tuning.loc[
        (tuning["model"] == combined_label)
        & tuning["selected"].astype(str).str.lower().isin({"true", "1", "t", "yes"})
    ].copy()
    required_tuning = {"outer_repeat", "outer_fold", "alpha", "selected_lambda"}
    if not required_tuning.issubset(selected.columns):
        raise ValueError("Selected tuning rows have an invalid schema")
    if len(selected) != 100:
        raise ValueError(f"Expected 100 selected combined-model fits; observed {len(selected)}")
    selected = selected.rename(columns={"alpha": "selected_alpha"})[
        ["outer_repeat", "outer_fold", "selected_alpha", "selected_lambda"]
    ]
    nonzero_count = (
        prs.assign(nonzero_bool=prs["nonzero"].astype(str).str.lower().isin({"true", "1", "t", "yes"}))
        .groupby(["outer_repeat", "outer_fold"], as_index=False)["nonzero_bool"]
        .sum()
        .rename(columns={"nonzero_bool": "nonzero_PRS_predictors"})
    )
    fit_summary = selected.merge(nonzero_count, on=["outer_repeat", "outer_fold"], validate="one_to_one")
    fit_summary["fit_id"] = (
        "R" + fit_summary["outer_repeat"].astype(int).astype(str).str.zfill(2)
        + "F" + fit_summary["outer_fold"].astype(int).astype(str).str.zfill(2)
    )
    fit_summary["sparse_fit"] = fit_summary["selected_alpha"] > 0.0
    fit_summary = fit_summary.sort_values(["outer_repeat", "outer_fold"]).reset_index(drop=True)

    prs = prs.merge(
        fit_summary[["outer_repeat", "outer_fold", "selected_alpha", "selected_lambda", "fit_id", "sparse_fit"]],
        on=["outer_repeat", "outer_fold"],
        validate="many_to_one",
    )
    prs["nonzero"] = prs["nonzero"].astype(str).str.lower().isin({"true", "1", "t", "yes"})
    prs["coefficient_sign"] = np.select(
        [prs["coefficient"] > 0.0, prs["coefficient"] < 0.0],
        ["positive", "negative"],
        default="zero",
    )

    rows: list[dict[str, object]] = []
    for term in PRS_ORDER:
        trait = prs.loc[prs["term"] == term].copy()
        sparse = trait.loc[trait["sparse_fit"]]
        nonzero = trait.loc[trait["nonzero"]]
        sparse_nonzero = sparse.loc[sparse["nonzero"]]
        rows.append(
            {
                "PRS": term.replace("PRS_", ""),
                "outer_fits": len(trait),
                "ridge_fits_alpha_0": int((trait["selected_alpha"] == 0.0).sum()),
                "sparse_fits_alpha_gt_0": len(sparse),
                "overall_nonzero_fits": int(trait["nonzero"].sum()),
                "overall_nonzero_percent": 100.0 * trait["nonzero"].mean(),
                "sparse_nonzero_fits": int(sparse["nonzero"].sum()),
                "sparse_nonzero_percent": 100.0 * sparse["nonzero"].mean(),
                "positive_nonzero_fits": int((nonzero["coefficient"] > 0.0).sum()),
                "negative_nonzero_fits": int((nonzero["coefficient"] < 0.0).sum()),
                "positive_percent_among_nonzero": 100.0 * (nonzero["coefficient"] > 0.0).mean(),
                "median_coefficient_all_fits": trait["coefficient"].median(),
                "median_absolute_coefficient_all_fits": trait["coefficient"].abs().median(),
                "median_coefficient_nonzero_fits": nonzero["coefficient"].median(),
                "q1_coefficient_nonzero_fits": nonzero["coefficient"].quantile(0.25),
                "q3_coefficient_nonzero_fits": nonzero["coefficient"].quantile(0.75),
                "median_coefficient_sparse_nonzero_fits": (
                    sparse_nonzero["coefficient"].median() if len(sparse_nonzero) else np.nan
                ),
            }
        )
    stability = pd.DataFrame(rows)
    long = prs[
        [
            "outer_repeat", "outer_fold", "fit_id", "term", "coefficient", "nonzero",
            "coefficient_sign", "selected_alpha", "selected_lambda", "sparse_fit",
        ]
    ].copy()
    long["term"] = long["term"].str.replace("PRS_", "", regex=False)
    long = long.sort_values(["outer_repeat", "outer_fold", "term"]).reset_index(drop=True)
    return long, stability, fit_summary


def percentile_summary(values: np.ndarray) -> tuple[float, float, float, float]:
    """Return an estimate with percentile bootstrap limits."""
    return (
        float(np.mean(values)),
        float(np.std(values, ddof=1)),
        float(np.quantile(values, 0.025)),
        float(np.quantile(values, 0.975)),
    )


def main() -> None:
    """Generate formal calibration and PRS-coefficient diagnostics."""
    args = parse_args()
    if args.bootstrap_replicates < 1_000:
        raise SystemExit("At least 1,000 bootstrap replicates are required")
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    _, outcome, comparator_logits, combined_logits, participants = reshape_predictions(
        args.predictions, args.expected_participants, args.expected_events
    )
    comparator_intercepts, comparator_joint, comparator_slopes, comparator_failures = repeat_calibration(
        outcome, comparator_logits
    )
    combined_intercepts, combined_joint, combined_slopes, combined_failures = repeat_calibration(
        outcome, combined_logits
    )
    if comparator_failures or combined_failures:
        raise RuntimeError("Point-estimate calibration fitting did not converge")

    repeat_rows: list[dict[str, object]] = []
    for repeat in range(10):
        for model, intercepts, joint, slopes in [
            (args.comparator_label, comparator_intercepts, comparator_joint, comparator_slopes),
            (args.combined_label, combined_intercepts, combined_joint, combined_slopes),
        ]:
            repeat_rows.append(
                {
                    "outer_repeat": repeat + 1,
                    "model": model,
                    "calibration_in_the_large": intercepts[repeat],
                    "joint_recalibration_intercept": joint[repeat],
                    "calibration_slope": slopes[repeat],
                }
            )
    repeat_estimates = pd.DataFrame(repeat_rows)

    rng = np.random.default_rng(args.seed)
    case_indices = np.flatnonzero(outcome == 1)
    control_indices = np.flatnonzero(outcome == 0)
    sampled_cases = rng.choice(case_indices, size=(args.bootstrap_replicates, len(case_indices)), replace=True)
    sampled_controls = rng.choice(control_indices, size=(args.bootstrap_replicates, len(control_indices)), replace=True)
    sampled_indices = np.concatenate([sampled_cases, sampled_controls], axis=1)
    bootstrap_values, failure_counts = bootstrap_calibration(
        outcome, comparator_logits, combined_logits, sampled_indices
    )
    if np.any(failure_counts != 0):
        raise RuntimeError(f"Calibration failed in {int(np.sum(failure_counts != 0))} bootstrap replicates")

    columns = [
        "comparator_calibration_in_the_large",
        "combined_calibration_in_the_large",
        "combined_minus_comparator_calibration_in_the_large",
        "comparator_calibration_slope",
        "combined_calibration_slope",
        "combined_minus_comparator_calibration_slope",
    ]
    bootstrap_distribution = pd.DataFrame(bootstrap_values, columns=columns)
    bootstrap_distribution.insert(0, "bootstrap_replicate", np.arange(1, args.bootstrap_replicates + 1))
    point = {
        "comparator_calibration_in_the_large": float(np.mean(comparator_intercepts)),
        "combined_calibration_in_the_large": float(np.mean(combined_intercepts)),
        "combined_minus_comparator_calibration_in_the_large": float(np.mean(combined_intercepts - comparator_intercepts)),
        "comparator_calibration_slope": float(np.mean(comparator_slopes)),
        "combined_calibration_slope": float(np.mean(combined_slopes)),
        "combined_minus_comparator_calibration_slope": float(np.mean(combined_slopes - comparator_slopes)),
    }
    labels = {
        "comparator_calibration_in_the_large": ("Comparator", "Calibration-in-the-large"),
        "combined_calibration_in_the_large": ("Combined", "Calibration-in-the-large"),
        "combined_minus_comparator_calibration_in_the_large": ("Combined minus comparator", "Calibration-in-the-large difference"),
        "comparator_calibration_slope": ("Comparator", "Calibration slope"),
        "combined_calibration_slope": ("Combined", "Calibration slope"),
        "combined_minus_comparator_calibration_slope": ("Combined minus comparator", "Calibration slope difference"),
    }
    summary_rows = []
    for column in columns:
        mean, standard_error, lower, upper = percentile_summary(bootstrap_distribution[column].to_numpy(float))
        model_or_contrast, metric = labels[column]
        summary_rows.append(
            {
                "model_or_contrast": model_or_contrast,
                "metric": metric,
                "point_estimate_mean_across_repeats": point[column],
                "bootstrap_mean": mean,
                "bootstrap_standard_error": standard_error,
                "lower_95_percentile": lower,
                "upper_95_percentile": upper,
                "bootstrap_replicates": args.bootstrap_replicates,
                "seed": args.seed,
            }
        )
    calibration_summary = pd.DataFrame(summary_rows)

    coefficient_long, stability_summary, tuning_summary = analyse_coefficients(
        args.coefficients, args.tuning, args.combined_label
    )

    checks: list[tuple[str, object, object]] = [
        ("participants", len(participants), args.expected_participants),
        ("events", int(np.sum(outcome == 1)), args.expected_events),
        ("outer_repeats", comparator_logits.shape[1], 10),
        ("point_calibration_converged", comparator_failures + combined_failures, 0),
        ("bootstrap_calibration_converged", int(np.sum(failure_counts)), 0),
        ("bootstrap_replicates", len(bootstrap_distribution), args.bootstrap_replicates),
        ("PRS_coefficient_rows", len(coefficient_long), 800),
        ("combined_outer_fits", len(tuning_summary), 100),
        ("PRS_stability_rows", len(stability_summary), 8),
        ("ridge_combined_fits", int((tuning_summary["selected_alpha"] == 0.0).sum()), 53),
        ("sparse_combined_fits", int((tuning_summary["selected_alpha"] > 0.0).sum()), 47),
    ]

    if args.validate_fixed:
        expected_calibration = pd.read_csv(
            ROOT / "results" / "validation" / "primary_calibration_expected.tsv", sep="\t"
        )
        lookup = {
            ("Clinical + ancestry PCs", "Calibration-in-the-large"): ("Comparator", "Calibration-in-the-large"),
            ("Clinical + ancestry PCs + 8 PRS", "Calibration-in-the-large"): ("Combined", "Calibration-in-the-large"),
            ("Combined minus clinical", "Calibration-in-the-large difference"): ("Combined minus comparator", "Calibration-in-the-large difference"),
            ("Clinical + ancestry PCs", "Calibration slope"): ("Comparator", "Calibration slope"),
            ("Clinical + ancestry PCs + 8 PRS", "Calibration slope"): ("Combined", "Calibration slope"),
            ("Combined minus clinical", "Calibration slope difference"): ("Combined minus comparator", "Calibration slope difference"),
        }
        for expected in expected_calibration.itertuples(index=False):
            model, metric = lookup[(expected.model_or_contrast, expected.metric)]
            observed = calibration_summary.loc[
                (calibration_summary["model_or_contrast"] == model)
                & (calibration_summary["metric"] == metric)
            ].iloc[0]
            checks.extend(
                [
                    (f"fixed_{model}_{metric}_estimate", abs(observed.point_estimate_mean_across_repeats - expected.point_estimate_mean_across_repeats) <= 1e-10, True),
                    (f"fixed_{model}_{metric}_lower", abs(observed.lower_95_percentile - expected.lower_95_percentile) <= 1e-10, True),
                    (f"fixed_{model}_{metric}_upper", abs(observed.upper_95_percentile - expected.upper_95_percentile) <= 1e-10, True),
                ]
            )
        fixed_stability = pd.read_csv(
            ROOT / "results" / "figure_data" / "supplementary_figure_s8_prs_coefficient_stability_summary.tsv",
            sep="\t",
        )
        for observed_row, expected_row in zip(stability_summary.itertuples(index=False), fixed_stability.itertuples(index=False), strict=True):
            checks.append((f"fixed_{observed_row.PRS}_overall_nonzero", observed_row.overall_nonzero_fits, expected_row.overall_nonzero_fits))
            checks.append((f"fixed_{observed_row.PRS}_sparse_nonzero", observed_row.sparse_nonzero_fits, expected_row.sparse_nonzero_fits))

    validation = pd.DataFrame(checks, columns=["check", "observed", "expected"])
    validation["result"] = np.where(validation["observed"].astype(str) == validation["expected"].astype(str), "PASS", "FAIL")

    write_tsv(repeat_estimates, output_dir / "calibration_repeat_estimates.tsv")
    write_tsv(calibration_summary, output_dir / "calibration_summary.tsv")
    write_tsv(bootstrap_distribution, output_dir / "calibration_bootstrap_distribution.tsv.gz")
    write_tsv(coefficient_long, output_dir / "prs_coefficient_outer_fit_long.tsv.gz")
    write_tsv(stability_summary, output_dir / "prs_coefficient_stability_summary.tsv")
    write_tsv(tuning_summary, output_dir / "outer_fit_tuning_summary.tsv")
    write_tsv(validation, output_dir / "diagnostics_validation.tsv")

    failed = validation.loc[validation["result"] != "PASS"]
    if not failed.empty:
        raise SystemExit(f"Model diagnostics validation failed: {failed['check'].tolist()}")
    (output_dir / "_SUCCESS").write_text(
        "\n".join(
            [
                "status=PASS",
                f"participants={args.expected_participants}",
                f"bootstrap_replicates={args.bootstrap_replicates}",
                f"bootstrap_seed={args.seed}",
                "models_refitted=FALSE",
                "hyperparameters_retuned=FALSE",
                "outer_training_sets_overlap=TRUE",
                "coefficient_frequencies_are_descriptive=TRUE",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"MODEL_DIAGNOSTICS=PASS output={output_dir}")


if __name__ == "__main__":
    main()
