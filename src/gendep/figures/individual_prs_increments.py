#!/usr/bin/env python3
"""Generate Supplementary Figure S7: individual-PRS predictive increments.

Purpose
-------
Display AUC, Brier-score and log-loss increments for each of the eight
individual Bayesian polygenic scores, with participant-level paired-bootstrap
intervals.

Inputs
------
- results/figure_data/supplementary_figure_s7_individual_prs_increments.tsv
- results/figure_data/supplementary_figure_s7_provenance.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S7_individual_PRS_increments.png
- results/figures/Supplementary_Figure_S7_individual_PRS_increments.svg
- results/figures/Supplementary_Figure_S7_individual_PRS_increments.pdf
- results/manifests/supplementary_figure_s7_sha256.txt

Scope
-----
The script reads release-safe aggregate trait-level results only. It does not access
participant-level predictions, fit models, perform cross-validation or repeat
the paired bootstrap.
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
    / "supplementary_figure_s7_individual_prs_increments.tsv"
)
DEFAULT_PROVENANCE: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s7_provenance.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# The highlighted row, interval geometry and exact panel scales communicate the
# result independently of colour alone.
STANDARD_POINT_COLOUR: Final = "#0072B2"
HIGHLIGHT_POINT_COLOUR: Final = "#D55E00"
INTERVAL_COLOUR: Final = "#666666"
HIGHLIGHT_BAND_COLOUR: Final = "#FFF1C7"
ZERO_LINE_COLOUR: Final = "#222222"
GRID_COLOUR: Final = "#D9D9D9"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 11.0,
        "axes.labelsize": 9.8,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.8,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

TRAIT_ORDER: Final = (
    "MDD",
    "ANX",
    "BIP",
    "SCZ",
    "NEUR",
    "INSOM",
    "SWB",
    "EA",
)
METRIC_ORDER: Final = ("AUC", "Brier", "Log loss")
PANEL_LABELS: Final = {
    "AUC": "A",
    "Brier": "B",
    "Log loss": "C",
}
PANEL_TITLES: Final = {
    "AUC": "AUC increment",
    "Brier": "Brier improvement",
    "Log loss": "Log-loss improvement",
}
PANEL_LIMITS: Final = {
    "AUC": (-0.030, 0.050),
    "Brier": (-0.0045, 0.0055),
    "Log loss": (-0.0095, 0.0115),
}
PANEL_TICKS: Final = {
    "AUC": (-0.02, 0.00, 0.02, 0.04),
    "Brier": (-0.004, -0.002, 0.000, 0.002, 0.004),
    "Log loss": (-0.005, 0.000, 0.005, 0.010),
}


@dataclass(frozen=True)
class TraitIncrement:
    """One individual-PRS point estimate and uncertainty interval."""

    trait_order: int
    trait: str
    metric_order: int
    metric: str
    point_estimate: float
    lower_95: float
    upper_95: float
    interval_excludes_zero: bool
    highlight_trait: bool
    bootstrap_replicates: int


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the individual-PRS increment figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument(
        "--provenance",
        type=Path,
        default=DEFAULT_PROVENANCE,
    )
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


def read_increments(path: Path) -> list[TraitIncrement]:
    """Read and validate the fixed 24-row plotting table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")

    required_columns = {
        "trait_order",
        "trait",
        "metric_order",
        "metric",
        "source_metric",
        "participants",
        "outer_repetitions",
        "point_estimate",
        "bootstrap_standard_error",
        "lower_95",
        "upper_95",
        "interval_excludes_zero",
        "highlight_trait",
        "positive_value_interpretation",
        "bootstrap_positive_proportion",
        "bootstrap_negative_proportion",
        "bootstrap_zero_proportion",
        "bootstrap_replicates",
    }

    increments: list[TraitIncrement] = []
    observed_pairs: set[tuple[str, str]] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected Supplementary Figure S7 plotting-table schema."
            )

        # Every trait-metric pair must occur exactly once.
        for row_number, row in enumerate(reader, start=2):
            trait = row["trait"]
            metric = row["metric"]
            pair = (trait, metric)

            if trait not in TRAIT_ORDER:
                raise ValueError(
                    f"Unexpected trait on row {row_number}: {trait}"
                )
            if metric not in METRIC_ORDER:
                raise ValueError(
                    f"Unexpected metric on row {row_number}: {metric}"
                )
            if pair in observed_pairs:
                raise ValueError(f"Duplicate trait-metric pair: {pair}")
            observed_pairs.add(pair)

            if int(row["participants"]) != 430:
                raise ValueError(
                    f"Unexpected participant count on row {row_number}."
                )
            if int(row["outer_repetitions"]) != 10:
                raise ValueError(
                    f"Unexpected outer-repetition count on row {row_number}."
                )
            if int(row["bootstrap_replicates"]) != 10000:
                raise ValueError(
                    f"Unexpected bootstrap count on row {row_number}."
                )
            if row["positive_value_interpretation"] != "Favours individual PRS":
                raise ValueError(
                    f"Unexpected direction definition on row {row_number}."
                )

            estimate = float(row["point_estimate"])
            lower = float(row["lower_95"])
            upper = float(row["upper_95"])
            if not lower <= estimate <= upper:
                raise ValueError(
                    f"Point estimate outside interval on row {row_number}."
                )

            excludes_zero = lower > 0.0 or upper < 0.0
            if (
                row["interval_excludes_zero"].upper() == "TRUE"
            ) != excludes_zero:
                raise ValueError(
                    f"Incorrect zero-exclusion flag on row {row_number}."
                )

            highlight_trait = row["highlight_trait"].upper() == "TRUE"
            if highlight_trait != (trait == "SWB"):
                raise ValueError(
                    f"Incorrect highlighted-trait flag on row {row_number}."
                )

            increments.append(
                TraitIncrement(
                    trait_order=int(row["trait_order"]),
                    trait=trait,
                    metric_order=int(row["metric_order"]),
                    metric=metric,
                    point_estimate=estimate,
                    lower_95=lower,
                    upper_95=upper,
                    interval_excludes_zero=excludes_zero,
                    highlight_trait=highlight_trait,
                    bootstrap_replicates=int(row["bootstrap_replicates"]),
                )
            )

    expected_pairs = {
        (trait, metric)
        for trait in TRAIT_ORDER
        for metric in METRIC_ORDER
    }
    if observed_pairs != expected_pairs or len(increments) != 24:
        raise ValueError("The Supplementary Figure S7 result set is incomplete.")

    excluding_zero = {
        (row.trait, row.metric)
        for row in increments
        if row.interval_excludes_zero
    }
    expected_excluding_zero = {
        ("SWB", "AUC"),
        ("SWB", "Brier"),
        ("SWB", "Log loss"),
    }
    if excluding_zero != expected_excluding_zero:
        raise ValueError(
            "Expected only the three SWB intervals to exclude zero."
        )
    return increments


def validate_provenance(path: Path) -> None:
    """Validate the two-row provenance manifest."""

    if not path.is_file():
        raise FileNotFoundError(f"Source-run table not found: {path}")

    required_columns = {
        "source_order",
        "source_role",
        "source_file",
        "source_sha256",
        "source_stage",
        "included_in_public_package",
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows or set(rows[0]) != required_columns:
        raise ValueError(
            "Unexpected Supplementary Figure S7 source-run schema."
        )
    if len(rows) != 2:
        raise ValueError(f"Expected two source rows; observed {len(rows)}.")
    if any(row["included_in_public_package"].upper() != "TRUE" for row in rows):
        raise ValueError("Both aggregate source tables must be packaged.")


def clean_axes(ax: plt.Axes) -> None:
    """Apply the common publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", color=GRID_COLOUR, linewidth=0.65, alpha=0.75)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel label outside the plotting area."""

    ax.text(
        -0.13,
        1.055,
        label,
        transform=ax.transAxes,
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
    )


def add_metric_panel(
    ax: plt.Axes,
    metric: str,
    increments: list[TraitIncrement],
    *,
    show_trait_labels: bool,
) -> None:
    """Draw one metric-specific individual-PRS interval panel."""

    metric_rows = sorted(
        (row for row in increments if row.metric == metric),
        key=lambda row: row.trait_order,
    )
    y_positions = np.arange(len(TRAIT_ORDER), dtype=float)

    # The highlighted band spans the SWB row across the full plotting width.
    swb_index = TRAIT_ORDER.index("SWB")
    ax.axhspan(
        swb_index - 0.36,
        swb_index + 0.36,
        facecolor=HIGHLIGHT_BAND_COLOUR,
        edgecolor="none",
        zorder=0,
    )

    for y_position, row in zip(y_positions, metric_rows, strict=True):
        lower_error = row.point_estimate - row.lower_95
        upper_error = row.upper_95 - row.point_estimate

        ax.errorbar(
            row.point_estimate,
            y_position,
            xerr=np.array([[lower_error], [upper_error]]),
            fmt="none",
            ecolor=INTERVAL_COLOUR,
            elinewidth=1.15,
            capsize=2.8,
            capthick=1.0,
            zorder=2,
        )
        ax.plot(
            row.point_estimate,
            y_position,
            marker="o",
            markersize=5.7,
            color=(
                HIGHLIGHT_POINT_COLOUR
                if row.highlight_trait
                else STANDARD_POINT_COLOUR
            ),
            linestyle="none",
            zorder=3,
        )

    ax.axvline(
        0.0,
        color=ZERO_LINE_COLOUR,
        linestyle="--",
        linewidth=0.9,
        zorder=1,
    )
    ax.set_xlim(*PANEL_LIMITS[metric])
    ax.set_xticks(PANEL_TICKS[metric])
    ax.set_ylim(len(TRAIT_ORDER) - 0.5, -0.5)
    ax.set_yticks(
        y_positions,
        TRAIT_ORDER if show_trait_labels else [""] * len(TRAIT_ORDER),
    )
    ax.set_xlabel("Increment")
    ax.set_title(PANEL_TITLES[metric], pad=8)
    clean_axes(ax)
    add_panel_label(ax, PANEL_LABELS[metric])


def render(
    increments: list[TraitIncrement],
    output_dir: Path,
) -> list[Path]:
    """Render the complete three-panel individual-PRS uncertainty figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # Manual axes preserve three aligned panels while allowing metric-specific
    # x scales and one shared trait order.
    fig = plt.figure(figsize=(10.6, 5.8))
    axes = {
        "AUC": fig.add_axes((0.070, 0.17, 0.285, 0.72)),
        "Brier": fig.add_axes((0.385, 0.17, 0.285, 0.72)),
        "Log loss": fig.add_axes((0.700, 0.17, 0.285, 0.72)),
    }

    for metric in METRIC_ORDER:
        add_metric_panel(
            axes[metric],
            metric,
            increments,
            show_trait_labels=(metric == "AUC"),
        )

    fig.text(
        0.53,
        0.055,
        "Positive values favour the individual PRS. SWB is highlighted; "
        "the 24 exploratory predictive contrasts were not jointly "
        "multiplicity-corrected.",
        ha="center",
        va="bottom",
        fontsize=7.7,
        color="#666666",
    )

    stem = output_dir / "Supplementary_Figure_S7_individual_PRS_increments"
    outputs = [stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf")]

    # PNG supports dissertation insertion; SVG and PDF retain vector geometry
    # and editable text for archival and publication use.
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
    """Record checksums for inputs, script and rendered outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "supplementary_figure_s7_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n"
            )
    return record


def main() -> None:
    """Render the individual-PRS predictive-increment figure."""
    args = parse_args()
    input_path = args.input.resolve()
    provenance_path = args.provenance.resolve()

    # Both release-safe aggregate inputs are validated before the first plotting
    # object is created.
    increments = read_increments(input_path)
    validate_provenance(provenance_path)

    outputs = render(increments, args.output_dir.resolve())
    provenance_record = write_provenance(
        [input_path, provenance_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S7_GENERATION=PASS")
    print(f"trait_metric_rows={len(increments)}")
    print("intervals_excluding_zero=3")
    print("intervals_excluding_zero_traits=SWB")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
