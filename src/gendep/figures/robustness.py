#!/usr/bin/env python3
"""Generate Figure 4: integrated robustness of predictive increments.

Purpose
-------
Compare paired changes in discrimination and probabilistic accuracy after
addition of eight Bayesian polygenic scores across the predefined robustness
analyses and the completed orientation-concordance sensitivities.

Inputs
------
- results/figure_data/figure_04_robustness_increments_integrated.tsv

Outputs
-------
- results/figures/Figure_4_robustness_increments_integrated.png
- results/figures/Figure_4_robustness_increments_integrated.svg
- results/figures/Figure_4_robustness_increments_integrated.pdf
- results/manifests/figure_04_integrated_sha256.txt

Scope
-----
The script reads a fixed 21-row aggregate source table only. It does not
access participant-level predictions, alter PRS, fit models, perform
cross-validation or repeat the paired bootstrap.
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
DEFAULT_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "figure_04_robustness_increments_integrated.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

# Analysis order, interval geometry, marker shape and marker fill jointly
# encode the results, so interpretation does not depend on colour alone.
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass(frozen=True)
class Increment:
    """One paired increment and 95% percentile interval."""

    analysis_order: int
    analysis_key: str
    analysis_label: str
    analysis_detail: str
    participants: int
    metric_order: int
    metric: str
    estimate: float
    lower_95: float
    upper_95: float
    interval_excludes_zero: bool
    analysis_class: str
    bootstrap_replicates: int


ANALYSIS_ORDER: Final = (
    "primary_elastic_net",
    "orientation_primary",
    "strict_eur",
    "prs_focused",
    "summary_predictor",
    "orientation_summary",
    "random_forest",
)
METRIC_ORDER: Final = ("AUC", "Brier score", "Log loss")

ANALYSIS_MARKERS: Final = {
    "primary_analysis": "o",
    "orientation_sensitivity": "D",
}

# Metric-specific domains use the full interval range plus a small
# margin, preserving resolution despite different numerical scales.
PANEL_SPECIFICATION: Final = {
    "AUC": {
        "left": 0.43,
        "right": 0.600,
        "minimum": -0.040,
        "maximum": 0.010,
        "ticks": (-0.04, -0.02, 0.00),
        "title": "AUC increment",
    },
    "Brier score": {
        "left": 0.645,
        "right": 0.795,
        "minimum": -0.007,
        "maximum": 0.003,
        "ticks": (-0.006, -0.003, 0.000, 0.003),
        "title": "Brier improvement",
    },
    "Log loss": {
        "left": 0.840,
        "right": 0.975,
        "minimum": -0.016,
        "maximum": 0.010,
        "ticks": (-0.015, 0.000, 0.010),
        "title": "Log-loss improvement",
    },
}

ROW_Y: Final = {
    "primary_elastic_net": 0.755,
    "orientation_primary": 0.665,
    "strict_eur": 0.555,
    "prs_focused": 0.465,
    "summary_predictor": 0.355,
    "orientation_summary": 0.265,
    "random_forest": 0.155,
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the robustness figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--provenance-dir",
        type=Path,
        default=DEFAULT_PROVENANCE_DIR,
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_increments(path: Path) -> list[Increment]:
    """Read and validate the release-safe 21-row integrated source table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")

    required = {
        "analysis_order",
        "analysis_key",
        "analysis_label",
        "analysis_detail",
        "participants",
        "metric_order",
        "metric",
        "estimate",
        "lower_95",
        "upper_95",
        "interval_excludes_zero",
        "analysis_class",
        "positive_value_interpretation",
        "bootstrap_replicates",
        "source_record",
        "source_stage",
    }

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError("Unexpected integrated robustness-table schema.")

        increments: list[Increment] = []
        observed_pairs: set[tuple[str, str]] = set()

        # Every analysis-metric pair must occur exactly once, preventing silent
        # omission or duplication of any robustness result.
        for row_number, row in enumerate(reader, start=2):
            analysis_key = row["analysis_key"]
            metric = row["metric"]
            pair = (analysis_key, metric)

            if analysis_key not in ANALYSIS_ORDER:
                raise ValueError(f"Unexpected analysis on row {row_number}: {analysis_key}")
            if metric not in METRIC_ORDER:
                raise ValueError(f"Unexpected metric on row {row_number}: {metric}")
            if pair in observed_pairs:
                raise ValueError(f"Duplicate analysis-metric pair: {pair}")
            observed_pairs.add(pair)

            analysis_class = row["analysis_class"]
            if analysis_class not in ANALYSIS_MARKERS:
                raise ValueError(f"Unexpected analysis class on row {row_number}.")
            if row["positive_value_interpretation"] != "Favours addition of PRS":
                raise ValueError(f"Unexpected metric direction on row {row_number}.")
            if int(row["bootstrap_replicates"]) != 10000:
                raise ValueError(f"Unexpected bootstrap count on row {row_number}.")

            estimate = float(row["estimate"])
            lower = float(row["lower_95"])
            upper = float(row["upper_95"])
            interval_excludes_zero = upper < 0.0 or lower > 0.0
            if not lower <= estimate <= upper:
                raise ValueError(f"Estimate outside interval on row {row_number}.")
            if (
                row["interval_excludes_zero"].upper() == "TRUE"
            ) != interval_excludes_zero:
                raise ValueError(f"Incorrect zero-exclusion flag on row {row_number}.")

            increments.append(
                Increment(
                    analysis_order=int(row["analysis_order"]),
                    analysis_key=analysis_key,
                    analysis_label=row["analysis_label"],
                    analysis_detail=row["analysis_detail"],
                    participants=int(row["participants"]),
                    metric_order=int(row["metric_order"]),
                    metric=metric,
                    estimate=estimate,
                    lower_95=lower,
                    upper_95=upper,
                    interval_excludes_zero=interval_excludes_zero,
                    analysis_class=analysis_class,
                    bootstrap_replicates=int(row["bootstrap_replicates"]),
                )
            )

    expected_pairs = {
        (analysis_key, metric)
        for analysis_key in ANALYSIS_ORDER
        for metric in METRIC_ORDER
    }
    if observed_pairs != expected_pairs or len(increments) != 21:
        raise ValueError("The integrated robustness result set is incomplete.")
    if any(row.estimate >= 0.0 for row in increments):
        raise ValueError("All integrated point estimates must remain negative.")
    return increments


def map_value(metric: str, value: float) -> float:
    """Map one metric value into its assigned horizontal panel."""

    panel = PANEL_SPECIFICATION[metric]
    minimum = float(panel["minimum"])
    maximum = float(panel["maximum"])
    if not minimum <= value <= maximum:
        raise ValueError(
            f"{metric} value {value} falls outside the configured panel domain."
        )
    proportion = (value - minimum) / (maximum - minimum)
    return float(panel["left"]) + proportion * (
        float(panel["right"]) - float(panel["left"])
    )


def plot_with_default_style(
    ax: plt.Axes,
    x_values: object,
    y_values: object,
    **kwargs: object,
) -> list[plt.Line2D]:
    """Draw one element using the first default style-cycle entry.

    Resetting the property cycle keeps the forest plot visually consistent
    without defining a project-specific colour in the script.
    """

    ax.set_prop_cycle(None)
    return ax.plot(x_values, y_values, **kwargs)


def add_analysis_labels(
    ax: plt.Axes,
    increments: list[Increment],
) -> None:
    """Add comparator labels, cohort sizes and sensitivity indentation."""

    metadata: dict[str, tuple[str, str, int, str]] = {}
    for row in increments:
        metadata.setdefault(
            row.analysis_key,
            (
                row.analysis_label,
                row.analysis_detail,
                row.participants,
                row.analysis_class,
            ),
        )

    for analysis_key in ANALYSIS_ORDER:
        label, detail, participants, analysis_class = metadata[analysis_key]
        y = ROW_Y[analysis_key]
        x = 0.070 if analysis_class == "primary_analysis" else 0.095

        ax.text(
            x,
            y + 0.013,
            label,
            ha="left",
            va="center",
            fontsize=9.0,
            fontweight=(
                "semibold"
                if analysis_class == "primary_analysis"
                else "normal"
            ),
        )
        ax.text(
            x,
            y - 0.020,
            f"{detail} · n={participants}",
            ha="left",
            va="center",
            fontsize=7.2,
        )

    # Thin separators preserve the pairing of the primary and summary
    # analyses with their orientation-sensitivity counterparts.
    for y in (0.610, 0.210):
        plot_with_default_style(
            ax,
            [0.055, 0.975],
            [y, y],
            linewidth=0.65,
            alpha=0.22,
        )


def add_metric_panel(
    ax: plt.Axes,
    metric: str,
    increments: list[Increment],
) -> None:
    """Draw one metric-specific forest-plot column."""

    panel = PANEL_SPECIFICATION[metric]
    left = float(panel["left"])
    right = float(panel["right"])
    top = 0.805
    bottom = 0.115

    ax.text(
        (left + right) / 2,
        0.865,
        str(panel["title"]),
        ha="center",
        va="center",
        fontsize=9.5,
        fontweight="bold",
    )

    # Each metric keeps its own numerical scale. The aligned rows permit
    # comparison across analyses without implying that the metrics share units.
    for tick in panel["ticks"]:
        x = map_value(metric, float(tick))
        plot_with_default_style(
            ax,
            [x, x],
            [bottom, top],
            linewidth=0.55,
            alpha=0.16,
        )
        ax.text(
            x,
            0.092,
            f"{float(tick):.3f}",
            ha="center",
            va="top",
            fontsize=7.0,
        )

    zero_x = map_value(metric, 0.0)
    plot_with_default_style(
        ax,
        [zero_x, zero_x],
        [bottom, top],
        linestyle="--",
        linewidth=0.9,
    )

    metric_rows = {
        row.analysis_key: row for row in increments if row.metric == metric
    }

    for analysis_key in ANALYSIS_ORDER:
        row = metric_rows[analysis_key]
        y = ROW_Y[analysis_key]
        lower_x = map_value(metric, row.lower_95)
        upper_x = map_value(metric, row.upper_95)
        estimate_x = map_value(metric, row.estimate)

        # The horizontal segment and caps represent the 95% percentile interval.
        plot_with_default_style(
            ax,
            [lower_x, upper_x],
            [y, y],
            linewidth=1.2,
        )
        plot_with_default_style(
            ax,
            [lower_x, lower_x],
            [y - 0.009, y + 0.009],
            linewidth=0.85,
        )
        plot_with_default_style(
            ax,
            [upper_x, upper_x],
            [y - 0.009, y + 0.009],
            linewidth=0.85,
        )

        # Marker shape identifies primary versus orientation-sensitivity
        # analyses. Marker fill independently records whether the interval
        # excludes zero.
        marker = ANALYSIS_MARKERS[row.analysis_class]
        marker_kwargs: dict[str, object] = {
            "marker": marker,
            "markersize": 5.5,
        }
        if not row.interval_excludes_zero:
            marker_kwargs["fillstyle"] = "none"

        plot_with_default_style(
            ax,
            estimate_x,
            y,
            **marker_kwargs,
        )

    plot_with_default_style(
        ax,
        [left, right],
        [0.105, 0.105],
        linewidth=0.75,
    )


def add_legend(ax: plt.Axes) -> None:
    """Add a compact shape and fill legend below the forest plot."""

    legend_y = 0.040
    positions = (0.47, 0.61, 0.76, 0.90)
    entries = (
        ("o", None, "Primary analysis"),
        ("D", None, "Orientation sensitivity"),
        ("o", "none", "Interval includes zero"),
        ("o", "full", "Interval excludes zero"),
    )

    for x, (marker, fillstyle, label) in zip(positions, entries, strict=True):
        kwargs: dict[str, object] = {
            "marker": marker,
            "markersize": 5.4,
        }
        if fillstyle == "none":
            kwargs["fillstyle"] = "none"
        plot_with_default_style(ax, x, legend_y, **kwargs)
        ax.text(
            x + 0.012,
            legend_y,
            label,
            ha="left",
            va="center",
            fontsize=7.0,
        )


def render(
    increments: list[Increment],
    output_dir: Path,
) -> list[Path]:
    """Render the integrated seven-analysis robustness forest plot."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # A single normalised canvas aligns all seven analyses while allowing each
    # metric to retain a scientifically appropriate scale.
    fig, ax = plt.subplots(figsize=(11.3, 8.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    fig.text(
        0.055,
        0.955,
        "Robustness of incremental predictive performance",
        ha="left",
        va="top",
        fontsize=14,
        fontweight="bold",
    )
    fig.text(
        0.055,
        0.920,
        "Primary analyses and orientation-concordance sensitivities",
        ha="left",
        va="top",
        fontsize=9.2,
    )

    add_analysis_labels(ax, increments)
    for metric in METRIC_ORDER:
        add_metric_panel(ax, metric, increments)
    add_legend(ax)

    # The uncertainty intervals condition on fixed saved out-of-fold
    # predictions; no model was refitted during the bootstrap.
    fig.text(
        0.055,
        0.010,
        "Positive values favour addition of PRS. Points show paired estimates; "
        "horizontal bars show 95% percentile intervals from 10,000 "
        "participant-level paired bootstrap replicates conditional on fixed "
        "saved predictions.",
        ha="left",
        va="bottom",
        fontsize=7.4,
    )

    stem = output_dir / "Figure_4_robustness_increments_integrated"
    outputs = [stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf")]

    # PNG supports dissertation insertion; SVG and PDF retain vector geometry
    # and editable text for archival and publication use.
    fig.savefig(outputs[0], dpi=600, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[2], bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return outputs


def write_provenance(
    input_path: Path,
    outputs: list[Path],
    provenance_dir: Path,
) -> Path:
    """Record checksums for the integrated table, script and outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "figure_04_integrated_sha256.txt"
    paths = [input_path, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render the predefined robustness comparison figure."""
    args = parse_args()
    input_path = args.input.resolve()

    # The complete integrated source is validated before any plotting object is
    # created.
    increments = read_increments(input_path)
    outputs = render(increments, args.output_dir.resolve())
    provenance_record = write_provenance(
        input_path,
        outputs,
        args.provenance_dir.resolve(),
    )

    print("FIGURE_04_INTEGRATED_GENERATION=PASS")
    print(f"rows={len(increments)}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
