#!/usr/bin/env python3
"""Generate Figure 3: primary repeated nested Elastic Net performance.

Purpose
-------
Visualise the primary paired repeated-cross-validation comparison, descriptive
participant-mean discrimination and calibration, and conditional paired-
bootstrap increments after addition of eight Bayesian polygenic scores.

Inputs
------
- results/figure_data/figure_03_repeat_auc.tsv
- results/figure_data/figure_03_roc_coordinates.tsv
- results/figure_data/figure_03_calibration_bins.tsv
- results/figure_data/figure_03_increment_intervals.tsv
- results/figure_data/figure_03_model_summary.tsv
- results/figure_data/figure_03_formal_calibration_summary.tsv

Outputs
-------
- results/figures/Figure_3_primary_performance.png
- results/figures/Figure_3_primary_performance.svg
- results/figures/Figure_3_primary_performance.pdf
- results/manifests/figure_03_sha256.txt

Scope
-----
The script reads deidentified aggregate plotting tables only. It does not
access participant identifiers, fit or tune models, calculate polygenic
scores, generate resampling folds or repeat the paired bootstrap.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final = Path(__file__).resolve().parents[3]
DEFAULT_DATA_DIR: Final = PROJECT_ROOT / "results" / "figure_data"
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

# Marker shape, line style and panel position also encode model identity, so
# the figure does not depend on colour alone.
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 10.2,
        "axes.labelsize": 9.3,
        "xtick.labelsize": 8.3,
        "ytick.labelsize": 8.3,
        "legend.fontsize": 7.8,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

MODEL_ORDER: Final = ("clinical", "combined")
MODEL_SHORT_LABELS: Final = {
    "clinical": "Clinical + ancestry PCs",
    "combined": "Clinical + ancestry PCs + 8 PRS",
}
MODEL_MARKERS: Final = {
    "clinical": "o",
    "combined": "s",
}
METRIC_ORDER: Final = ("AUC", "Brier score", "Log loss")

EXPECTED_REPEAT_MEANS: Final = {
    "clinical": 0.693930266520628,
    "combined": 0.688378514056225,
}
EXPECTED_PARTICIPANT_MEAN_AUC: Final = {
    "clinical": 0.699981745162468,
    "combined": 0.694596568090544,
}


@dataclass(frozen=True)
class RepeatAUC:
    """One paired repetition-level AUC comparison."""

    outer_repeat: int
    clinical_auc: float
    combined_auc: float
    difference: float


@dataclass(frozen=True)
class ROCCoordinate:
    """One point on a deidentified participant-mean ROC curve."""

    model_key: str
    point_order: int
    false_positive_rate: float
    true_positive_rate: float


@dataclass(frozen=True)
class CalibrationBin:
    """One quantile-bin calibration point."""

    model_key: str
    bin_number: int
    mean_predicted_probability: float
    observed_remission_proportion: float
    participants: int


@dataclass(frozen=True)
class IncrementInterval:
    """One paired-bootstrap increment and percentile interval."""

    metric_order: int
    metric: str
    estimate: float
    lower_95: float
    upper_95: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the primary-performance figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory containing the release-safe Figure 3 plotting tables.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for PNG, SVG and PDF outputs.",
    )
    parser.add_argument(
        "--provenance-dir",
        type=Path,
        default=DEFAULT_PROVENANCE_DIR,
        help="Directory for SHA-256 provenance records.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a tab-separated plotting table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_repeat_auc(path: Path) -> list[RepeatAUC]:
    """Read and validate the ten paired repetition-level AUC values."""

    rows = read_tsv(path)
    required_columns = {
        "outer_repeat",
        "clinical_auc",
        "combined_auc",
        "combined_minus_clinical_auc",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected repetition-level Figure 3 schema.")
    if len(rows) != 10:
        raise ValueError(f"Expected 10 repetition rows; observed {len(rows)}.")

    values = [
        RepeatAUC(
            outer_repeat=int(row["outer_repeat"]),
            clinical_auc=float(row["clinical_auc"]),
            combined_auc=float(row["combined_auc"]),
            difference=float(row["combined_minus_clinical_auc"]),
        )
        for row in rows
    ]
    values.sort(key=lambda row: row.outer_repeat)

    if [row.outer_repeat for row in values] != list(range(1, 11)):
        raise ValueError("Repetition identifiers must be the integers 1 through 10.")

    for row in values:
        if not np.isclose(
            row.combined_auc - row.clinical_auc,
            row.difference,
            rtol=0.0,
            atol=1e-14,
        ):
            raise ValueError(f"Inconsistent AUC difference in repetition {row.outer_repeat}.")

    for model_key, observed_mean in (
        ("clinical", float(np.mean([row.clinical_auc for row in values]))),
        ("combined", float(np.mean([row.combined_auc for row in values]))),
    ):
        if not np.isclose(
            observed_mean,
            EXPECTED_REPEAT_MEANS[model_key],
            rtol=0.0,
            atol=1e-14,
        ):
            raise ValueError(f"Unexpected repeat-mean AUC for {model_key}.")

    if sum(row.difference > 0 for row in values) != 1:
        raise ValueError("The combined model must have higher AUC in 1/10 repetitions.")
    return values


def read_roc(path: Path) -> dict[str, list[ROCCoordinate]]:
    """Read and validate the two deidentified ROC curves."""

    rows = read_tsv(path)
    required_columns = {
        "model_key",
        "model_label",
        "point_order",
        "false_positive_rate",
        "true_positive_rate",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected ROC-coordinate Figure 3 schema.")

    curves: dict[str, list[ROCCoordinate]] = {model_key: [] for model_key in MODEL_ORDER}
    for row in rows:
        model_key = row["model_key"]
        if model_key not in curves:
            raise ValueError(f"Unexpected ROC model key: {model_key}")
        curves[model_key].append(
            ROCCoordinate(
                model_key=model_key,
                point_order=int(row["point_order"]),
                false_positive_rate=float(row["false_positive_rate"]),
                true_positive_rate=float(row["true_positive_rate"]),
            )
        )

    for model_key, coordinates in curves.items():
        coordinates.sort(key=lambda row: row.point_order)
        if len(coordinates) < 2:
            raise ValueError(f"Insufficient ROC coordinates for {model_key}.")
        false_positive = np.array(
            [row.false_positive_rate for row in coordinates],
            dtype=float,
        )
        true_positive = np.array(
            [row.true_positive_rate for row in coordinates],
            dtype=float,
        )
        if np.any(np.diff(false_positive) < 0):
            raise ValueError(f"ROC false-positive rates decrease for {model_key}.")
        if not (
            np.isclose(false_positive[0], 0.0)
            and np.isclose(true_positive[0], 0.0)
            and np.isclose(false_positive[-1], 1.0)
            and np.isclose(true_positive[-1], 1.0)
        ):
            raise ValueError(f"ROC endpoints are invalid for {model_key}.")

        auc_value = float(np.trapezoid(true_positive, false_positive))
        if not np.isclose(
            auc_value,
            EXPECTED_PARTICIPANT_MEAN_AUC[model_key],
            rtol=0.0,
            atol=1e-14,
        ):
            raise ValueError(f"Unexpected participant-mean ROC AUC for {model_key}.")
    return curves


def read_calibration(path: Path) -> dict[str, list[CalibrationBin]]:
    """Read and validate ten quantile-bin calibration points per model."""

    rows = read_tsv(path)
    required_columns = {
        "model_key",
        "model_label",
        "bin",
        "mean_predicted_probability",
        "observed_remission_proportion",
        "participants",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected calibration-bin Figure 3 schema.")

    curves: dict[str, list[CalibrationBin]] = {model_key: [] for model_key in MODEL_ORDER}
    for row in rows:
        model_key = row["model_key"]
        if model_key not in curves:
            raise ValueError(f"Unexpected calibration model key: {model_key}")
        curves[model_key].append(
            CalibrationBin(
                model_key=model_key,
                bin_number=int(row["bin"]),
                mean_predicted_probability=float(row["mean_predicted_probability"]),
                observed_remission_proportion=float(
                    row["observed_remission_proportion"]
                ),
                participants=int(row["participants"]),
            )
        )

    for model_key, bins in curves.items():
        bins.sort(key=lambda row: row.bin_number)
        if [row.bin_number for row in bins] != list(range(1, 11)):
            raise ValueError(f"Calibration bins must be numbered 1 through 10 for {model_key}.")
        if sum(row.participants for row in bins) != 430:
            raise ValueError(f"Calibration-bin counts do not sum to 430 for {model_key}.")
        if any(
            not 0.0 <= value <= 1.0
            for row in bins
            for value in (
                row.mean_predicted_probability,
                row.observed_remission_proportion,
            )
        ):
            raise ValueError(f"Calibration coordinates fall outside [0, 1] for {model_key}.")
    return curves


def read_increments(path: Path) -> list[IncrementInterval]:
    """Read and validate the three primary paired-bootstrap increments."""

    rows = read_tsv(path)
    required_columns = {
        "metric_order",
        "metric",
        "estimate",
        "lower_95",
        "upper_95",
        "positive_value_interpretation",
        "uncertainty_scope",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected increment-interval Figure 3 schema.")

    intervals = [
        IncrementInterval(
            metric_order=int(row["metric_order"]),
            metric=row["metric"],
            estimate=float(row["estimate"]),
            lower_95=float(row["lower_95"]),
            upper_95=float(row["upper_95"]),
        )
        for row in rows
    ]
    intervals.sort(key=lambda row: row.metric_order)

    if [row.metric for row in intervals] != list(METRIC_ORDER):
        raise ValueError("Increment metrics are missing or in the wrong order.")
    for row in intervals:
        if not row.lower_95 <= row.estimate <= row.upper_95:
            raise ValueError(f"The {row.metric} estimate is outside its interval.")
    return intervals


def read_model_summary(path: Path) -> dict[str, dict[str, float | int | str]]:
    """Read the participant-mean model labels and AUC values."""

    rows = read_tsv(path)
    required_columns = {
        "model_key",
        "model_label",
        "participants",
        "events",
        "participant_mean_auc",
        "mean_predicted_probability",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected model-summary Figure 3 schema.")

    summary: dict[str, dict[str, float | int | str]] = {}
    for row in rows:
        model_key = row["model_key"]
        if model_key not in MODEL_ORDER or model_key in summary:
            raise ValueError(f"Unexpected or duplicate model key: {model_key}")
        summary[model_key] = {
            "model_label": row["model_label"],
            "participants": int(row["participants"]),
            "events": int(row["events"]),
            "participant_mean_auc": float(row["participant_mean_auc"]),
            "mean_predicted_probability": float(row["mean_predicted_probability"]),
        }

    if set(summary) != set(MODEL_ORDER):
        raise ValueError("The model summary must contain both primary comparator models.")
    for model_key in MODEL_ORDER:
        if summary[model_key]["participants"] != 430 or summary[model_key]["events"] != 166:
            raise ValueError(f"Unexpected cohort counts for {model_key}.")
    return summary


def clean_axes(ax: plt.Axes, *, grid_axis: str | None = None) -> None:
    """Apply a consistent publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis is not None:
        ax.grid(axis=grid_axis, alpha=0.22, linewidth=0.7)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel identifier outside the plotting region."""

    ax.text(
        -0.12,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=11.5,
        fontweight="bold",
        va="top",
    )


def add_repeat_auc_panel(ax: plt.Axes, values: list[RepeatAUC]) -> None:
    """Draw the primary paired repetition-level AUC comparison."""

    clinical_auc = np.array([row.clinical_auc for row in values], dtype=float)
    combined_auc = np.array([row.combined_auc for row in values], dtype=float)

    # One low-emphasis segment per repetition preserves pairing without
    # obscuring the two model-specific point distributions.
    segments = [
        [(0.0, clinical), (1.0, combined)]
        for clinical, combined in zip(clinical_auc, combined_auc, strict=True)
    ]
    ax.add_collection(LineCollection(segments, linewidths=0.9, alpha=0.25))

    ax.scatter(
        np.zeros(len(values)),
        clinical_auc,
        marker=MODEL_MARKERS["clinical"],
        s=30,
        label=MODEL_SHORT_LABELS["clinical"],
        zorder=3,
    )
    ax.scatter(
        np.ones(len(values)),
        combined_auc,
        marker=MODEL_MARKERS["combined"],
        s=30,
        label=MODEL_SHORT_LABELS["combined"],
        zorder=3,
    )

    # Short horizontal bars show the two primary repetition-mean estimands.
    clinical_mean = float(clinical_auc.mean())
    combined_mean = float(combined_auc.mean())
    ax.hlines(clinical_mean, -0.12, 0.12, linewidth=2.5)
    ax.hlines(combined_mean, 0.88, 1.12, linewidth=2.5)

    # A compact in-panel summary avoids placing mean labels over individual
    # repetition points while preserving the primary numerical comparison.
    ax.text(
        0.98,
        0.055,
        f"Mean AUC {clinical_mean:.3f} vs {combined_mean:.3f}\n"
        "PRS model higher in 1/10 repetitions",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.7,
        linespacing=1.15,
    )

    ax.set_xticks(
        (0, 1),
        (
            "Clinical +\nancestry PCs",
            "Clinical + ancestry PCs\n+ 8 PRS",
        ),
    )
    ax.set_xlim(-0.35, 1.35)
    ax.set_ylabel("ROC AUC")
    ax.set_title("Primary repeat-wise held-out AUC", loc="left")
    clean_axes(ax, grid_axis="y")
    add_panel_label(ax, "A")


def add_roc_panel(
    ax: plt.Axes,
    curves: dict[str, list[ROCCoordinate]],
    summary: dict[str, dict[str, float | int | str]],
) -> None:
    """Draw descriptive participant-mean out-of-fold ROC curves."""

    for model_key in MODEL_ORDER:
        coordinates = curves[model_key]
        ax.plot(
            [row.false_positive_rate for row in coordinates],
            [row.true_positive_rate for row in coordinates],
            linewidth=1.7,
            label=(
                (
                    "Comparator"
                    if model_key == "clinical"
                    else "+ 8 PRS"
                )
                + f" (AUC {float(summary[model_key]['participant_mean_auc']):.3f})"
            ),
        )

    # The diagonal reference indicates chance discrimination and is not a
    # third fitted model.
    ax.plot((0, 1), (0, 1), linestyle="--", linewidth=0.9, alpha=0.45)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xlabel("False-positive rate")
    ax.set_ylabel("True-positive rate")
    ax.set_title("Descriptive participant-mean OOF ROC", loc="left")
    ax.legend(loc="lower right", frameon=False)
    clean_axes(ax)
    add_panel_label(ax, "B")


def read_formal_calibration(path: Path) -> dict[str, dict[str, float]]:
    """Read formal repeat-wise calibration point estimates and intervals."""

    rows = read_tsv(path)
    required_columns = {
        "model_or_contrast",
        "metric",
        "point_estimate_mean_across_repeats",
        "lower_95_percentile",
        "upper_95_percentile",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError("Unexpected formal-calibration Figure 3 schema.")
    label_to_key = {
        "Clinical + ancestry PCs": "clinical",
        "Clinical + ancestry PCs + 8 PRS": "combined",
    }
    output: dict[str, dict[str, float]] = {key: {} for key in MODEL_ORDER}
    for row in rows:
        if row["model_or_contrast"] not in label_to_key:
            raise ValueError("Unexpected formal-calibration model label.")
        model_key = label_to_key[row["model_or_contrast"]]
        metric_key = (
            "intercept"
            if row["metric"] == "Calibration-in-the-large"
            else "slope"
            if row["metric"] == "Calibration slope"
            else None
        )
        if metric_key is None:
            raise ValueError("Unexpected formal-calibration metric.")
        output[model_key][metric_key] = float(
            row["point_estimate_mean_across_repeats"]
        )
        output[model_key][f"{metric_key}_lower"] = float(row["lower_95_percentile"])
        output[model_key][f"{metric_key}_upper"] = float(row["upper_95_percentile"])
    for model_key in MODEL_ORDER:
        if set(output[model_key]) != {
            "intercept",
            "intercept_lower",
            "intercept_upper",
            "slope",
            "slope_lower",
            "slope_upper",
        }:
            raise ValueError(f"Incomplete formal calibration for {model_key}.")
    return output


def add_calibration_panel(
    ax: plt.Axes,
    curves: dict[str, list[CalibrationBin]],
    formal_calibration: dict[str, dict[str, float]],
) -> None:
    """Draw descriptive ten-quantile-bin calibration curves."""

    for model_key in MODEL_ORDER:
        bins = curves[model_key]
        ax.plot(
            [row.mean_predicted_probability for row in bins],
            [row.observed_remission_proportion for row in bins],
            marker=MODEL_MARKERS[model_key],
            markersize=4.5,
            linewidth=1.3,
            label=(
                "Clinical + ancestry PCs"
                if model_key == "clinical"
                else "+ 8 PRS"
            ),
        )

    # The 45-degree line represents perfect agreement between mean predicted
    # probability and observed remission proportion.
    ax.plot((0, 1), (0, 1), linestyle="--", linewidth=0.9, alpha=0.45)
    ax.set_xlim(0, 0.85)
    ax.set_ylim(0, 0.85)
    ax.set_xlabel("Mean predicted probability")
    ax.set_ylabel("Observed remission proportion")
    ax.set_title("Descriptive calibration by prediction decile", loc="left")
    ax.legend(loc="upper left", frameon=False)
    annotation = (
        "Formal repeat-wise calibration\n"
        "(clinical | + 8 PRS)\n"
        f"Intercept: {formal_calibration['clinical']['intercept']:.3f} | "
        f"{formal_calibration['combined']['intercept']:.3f}\n"
        f"Slope: {formal_calibration['clinical']['slope']:.3f} | "
        f"{formal_calibration['combined']['slope']:.3f}"
    )
    ax.text(
        0.97,
        0.05,
        annotation,
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7.2,
        linespacing=1.15,
        bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "0.75", "linewidth": 0.6, "alpha": 0.92},
    )
    clean_axes(ax)
    add_panel_label(ax, "C")


def add_increment_panel(
    ax: plt.Axes,
    intervals: list[IncrementInterval],
) -> None:
    """Draw conditional paired-bootstrap increments and 95% intervals."""

    y_positions = np.arange(len(intervals))
    estimates = np.array([row.estimate for row in intervals], dtype=float)
    lower_errors = np.array(
        [row.estimate - row.lower_95 for row in intervals],
        dtype=float,
    )
    upper_errors = np.array(
        [row.upper_95 - row.estimate for row in intervals],
        dtype=float,
    )

    ax.errorbar(
        estimates,
        y_positions,
        xerr=np.vstack([lower_errors, upper_errors]),
        fmt="o",
        elinewidth=1.2,
        capsize=3.0,
        markersize=5.5,
    )
    ax.axvline(0, linewidth=0.9, linestyle="--")
    ax.set_yticks(y_positions, [row.metric for row in intervals])
    ax.invert_yaxis()
    ax.set_xlabel("Paired increment (positive favours PRS)")
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.set_title("Conditional paired-bootstrap increments", loc="left")
    clean_axes(ax, grid_axis="x")
    add_panel_label(ax, "D")

    # Compact numerical labels make the plotted estimates independently
    # readable without reproducing a full results table inside the panel.
    x_span = ax.get_xlim()[1] - ax.get_xlim()[0]
    for y_position, row in zip(y_positions, intervals, strict=True):
        ax.text(
            row.upper_95 + x_span * 0.02,
            y_position,
            f"{row.estimate:.4f}",
            va="center",
            ha="left",
            fontsize=7.6,
        )


def render(
    repeat_values: list[RepeatAUC],
    roc_curves: dict[str, list[ROCCoordinate]],
    calibration_curves: dict[str, list[CalibrationBin]],
    intervals: list[IncrementInterval],
    model_summary: dict[str, dict[str, float | int | str]],
    formal_calibration: dict[str, dict[str, float]],
    output_dir: Path,
) -> list[Path]:
    """Render the complete four-panel primary performance figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # The balanced two-by-two layout separates the primary repeated-CV
    # estimand from descriptive participant-mean plots and bootstrap intervals.
    fig, axes = plt.subplots(
        2,
        2,
        figsize=(7.55, 7.15),
        constrained_layout=True,
    )

    add_repeat_auc_panel(axes[0, 0], repeat_values)
    add_roc_panel(axes[0, 1], roc_curves, model_summary)
    add_calibration_panel(axes[1, 0], calibration_curves, formal_calibration)
    add_increment_panel(axes[1, 1], intervals)

    stem = output_dir / "Figure_3_primary_performance"
    outputs = [stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf")]

    # PNG supports direct dissertation insertion; SVG and PDF retain vector
    # geometry and editable text for archival and publication use.
    fig.savefig(outputs[0], dpi=600, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[2], bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return outputs


def write_provenance(
    input_paths: list[Path],
    outputs: list[Path],
    provenance_dir: Path,
) -> Path:
    """Record the plotting tables, script and rendered output checksums."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "figure_03_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]
    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render primary model performance and paired uncertainty."""
    args = parse_args()
    data_dir = args.data_dir.resolve()

    input_paths = [
        data_dir / "figure_03_repeat_auc.tsv",
        data_dir / "figure_03_roc_coordinates.tsv",
        data_dir / "figure_03_calibration_bins.tsv",
        data_dir / "figure_03_increment_intervals.tsv",
        data_dir / "figure_03_model_summary.tsv",
        data_dir / "figure_03_formal_calibration_summary.tsv",
    ]

    # Every deidentified source table is validated before the first plotting
    # object is created.
    repeat_values = read_repeat_auc(input_paths[0])
    roc_curves = read_roc(input_paths[1])
    calibration_curves = read_calibration(input_paths[2])
    intervals = read_increments(input_paths[3])
    model_summary = read_model_summary(input_paths[4])
    formal_calibration = read_formal_calibration(input_paths[5])

    outputs = render(
        repeat_values,
        roc_curves,
        calibration_curves,
        intervals,
        model_summary,
        formal_calibration,
        args.output_dir.resolve(),
    )
    provenance_record = write_provenance(
        input_paths,
        outputs,
        args.provenance_dir.resolve(),
    )

    print("FIGURE_03_GENERATION=PASS")
    print("participant_identifiers_accessed=FALSE")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
