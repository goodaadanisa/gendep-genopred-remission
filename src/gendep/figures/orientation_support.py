#!/usr/bin/env python3
"""Generate Supplementary Figure S1: frequency-orientation evidence.

Purpose
-------
Visualise the relationship between the supplied code-2 sample frequency and
the GenoPred European major-allele frequency, together with predefined
frequency-edge categories used in the orientation audit.

Inputs
------
- results/figure_data/supplementary_figure_s1_frequency_points.tsv.gz
- results/figure_data/supplementary_figure_s1_summary.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S1_frequency_orientation.png
- results/figures/Supplementary_Figure_S1_frequency_orientation.svg
- results/figures/Supplementary_Figure_S1_frequency_orientation.pdf
- results/manifests/supplementary_figure_s1_sha256.txt

Scope
-----
The script reads release-safe identifier-free frequency records only. It does not
access participant-level genotypes, alter the selected mapping, rebuild the
target or execute polygenic-score analyses.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import os
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Final

PROJECT_ROOT: Final = Path(__file__).resolve().parents[3]
DEFAULT_POINTS: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s1_frequency_points.tsv.gz"
)
DEFAULT_SUMMARY: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s1_summary.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Density, marker shape and panel position jointly encode meaning, so the
# figure remains interpretable without relying on colour alone.
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 10.2,
        "axes.labelsize": 9.4,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.2,
        "legend.fontsize": 7.8,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

CATEGORY_ORDER: Final = (
    "retained_sample_frequency_flip",
    "retained_sample_frequency_tie",
    "excluded_eur_frequency_tie",
)

CATEGORY_LABELS: Final = {
    "retained_sample_frequency_flip": "Retained sample-frequency flips",
    "retained_sample_frequency_tie": "Retained sample-frequency ties",
    "excluded_eur_frequency_tie": "Excluded EUR-frequency ties",
}

CATEGORY_MARKERS: Final = {
    "retained_sample_frequency_flip": "o",
    "retained_sample_frequency_tie": "s",
    "excluded_eur_frequency_tie": "^",
}

EXPECTED_VARIANTS: Final = 147438
EXPECTED_CORRELATION: Final = 0.9869353564843381
EXPECTED_COUNTS: Final = {
    "retained_standard": 146366,
    "retained_sample_frequency_flip": 899,
    "retained_sample_frequency_tie": 105,
    "excluded_eur_frequency_tie": 68,
}


@dataclass(frozen=True)
class FrequencyPoint:
    """One identifier-free variant frequency record."""

    code2_frequency: float
    eur_major_frequency: float
    category: str


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the orientation-support figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=Path, default=DEFAULT_POINTS)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
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


def read_points(path: Path) -> list[FrequencyPoint]:
    """Read and validate the release-safe identifier-free frequency table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure point table not found: {path}")

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_columns = {
            "gendep_code2_frequency",
            "eur_major_frequency",
            "category",
        }
        if reader.fieldnames is None or set(reader.fieldnames) != expected_columns:
            raise ValueError("Unexpected Supplementary Figure S1 point-table schema.")

        points: list[FrequencyPoint] = []
        category_counts = {category: 0 for category in EXPECTED_COUNTS}

        # Every coordinate is checked before plotting because a single invalid
        # probability could distort the density scale or zoomed panel.
        for row_number, row in enumerate(reader, start=2):
            code2_frequency = float(row["gendep_code2_frequency"])
            eur_major_frequency = float(row["eur_major_frequency"])
            category = row["category"].strip()

            if category not in EXPECTED_COUNTS:
                raise ValueError(f"Unexpected category on row {row_number}: {category}")
            for value in (code2_frequency, eur_major_frequency):
                if not np.isfinite(value) or not 0.0 <= value <= 1.0:
                    raise ValueError(f"Invalid frequency on row {row_number}: {value}")

            category_counts[category] += 1
            points.append(
                FrequencyPoint(
                    code2_frequency=code2_frequency,
                    eur_major_frequency=eur_major_frequency,
                    category=category,
                )
            )

    if len(points) != EXPECTED_VARIANTS:
        raise ValueError(f"Expected {EXPECTED_VARIANTS} points; observed {len(points)}.")
    if category_counts != EXPECTED_COUNTS:
        raise ValueError(
            f"Unexpected Figure S1 category counts: {category_counts}"
        )

    observed_correlation = float(
        np.corrcoef(
            [point.code2_frequency for point in points],
            [point.eur_major_frequency for point in points],
        )[0, 1]
    )
    if not np.isclose(
        observed_correlation,
        EXPECTED_CORRELATION,
        rtol=0.0,
        atol=1e-14,
    ):
        raise ValueError(
            f"Unexpected Pearson correlation: {observed_correlation:.17g}"
        )
    return points


def read_summary(path: Path) -> dict[str, float]:
    """Read and validate the annotation summary."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure summary table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows or set(rows[0]) != {"metric", "value", "display_label"}:
        raise ValueError("Unexpected Supplementary Figure S1 summary schema.")

    summary = {row["metric"]: float(row["value"]) for row in rows}
    required_metrics = {
        "matched_variants",
        "pearson_correlation",
        *EXPECTED_COUNTS,
    }
    if set(summary) != required_metrics:
        raise ValueError("The Supplementary Figure S1 summary is incomplete.")

    if int(summary["matched_variants"]) != EXPECTED_VARIANTS:
        raise ValueError("Unexpected matched-variant count in summary.")
    if not np.isclose(
        summary["pearson_correlation"],
        EXPECTED_CORRELATION,
        rtol=0.0,
        atol=1e-14,
    ):
        raise ValueError("Unexpected Pearson correlation in summary.")
    for category, expected_count in EXPECTED_COUNTS.items():
        if int(summary[category]) != expected_count:
            raise ValueError(f"Unexpected summary count for {category}.")
    return summary


def clean_axes(ax: plt.Axes) -> None:
    """Apply the shared publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(alpha=0.14, linewidth=0.6)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel identifier outside the plotting region."""

    ax.text(
        -0.12,
        1.055,
        label,
        transform=ax.transAxes,
        fontsize=11.5,
        fontweight="bold",
        va="top",
    )


def add_overall_panel(
    ax: plt.Axes,
    points: list[FrequencyPoint],
    summary: dict[str, float],
) -> None:
    """Draw the complete 147,438-variant density relationship."""

    x_values = np.array(
        [point.eur_major_frequency for point in points],
        dtype=float,
    )
    y_values = np.array(
        [point.code2_frequency for point in points],
        dtype=float,
    )

    # Log-scaled hexagon counts preserve the dense diagonal structure while
    # retaining visibility of lower-density departures from identity.
    density = ax.hexbin(
        x_values,
        y_values,
        gridsize=72,
        extent=(0.5, 1.0, 0.0, 1.0),
        mincnt=1,
        bins="log",
        linewidths=0.0,
    )

    # The identity line represents perfect agreement between sample code-2
    # frequency and the selected European major-allele frequency.
    ax.plot((0.5, 1.0), (0.5, 1.0), linestyle="--", linewidth=1.0)
    ax.set_xlim(0.5, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.set_xlabel("GenoPred EUR major-allele frequency")
    ax.set_ylabel("GENDEP code-2 sample frequency")
    ax.set_title("All matched variants", loc="left")
    clean_axes(ax)
    add_panel_label(ax, "A")

    ax.text(
        0.04,
        0.95,
        f"n = {int(summary['matched_variants']):,}\n"
        f"Pearson r = {summary['pearson_correlation']:.4f}",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=8.3,
    )

    # A compact colour bar reports the log-density scale without requiring raw
    # point overplotting across the complete matched set.
    colour_bar = ax.figure.colorbar(density, ax=ax, fraction=0.046, pad=0.03)
    colour_bar.set_label("log10 hexagon count", fontsize=8.2)
    colour_bar.ax.tick_params(labelsize=7.3)


def add_edge_panel(
    ax: plt.Axes,
    points: list[FrequencyPoint],
    summary: dict[str, float],
) -> None:
    """Draw the predefined edge categories near frequency one half."""

    for category in CATEGORY_ORDER:
        selected = [point for point in points if point.category == category]
        x_values = [point.eur_major_frequency for point in selected]
        y_values = [point.code2_frequency for point in selected]

        # Marker shape distinguishes the three audit categories independently
        # of the default plotting colours.
        ax.scatter(
            x_values,
            y_values,
            marker=CATEGORY_MARKERS[category],
            s=24,
            alpha=0.72,
            label=(
                f"{CATEGORY_LABELS[category]} "
                f"(n={int(summary[category]):,})"
            ),
        )

    ax.plot((0.495, 0.61), (0.495, 0.61), linestyle="--", linewidth=1.0)
    ax.axhline(0.5, linestyle=":", linewidth=0.8)
    ax.axvline(0.5, linestyle=":", linewidth=0.8)
    ax.set_xlim(0.495, 0.61)
    ax.set_ylim(0.45, 0.61)
    ax.set_xlabel("GenoPred EUR major-allele frequency")
    ax.set_ylabel("GENDEP code-2 sample frequency")
    ax.set_title("Frequency-edge categories", loc="left")
    ax.legend(loc="upper right", frameon=False)
    clean_axes(ax)
    add_panel_label(ax, "B")

    ax.text(
        0.03,
        0.04,
        "Dashed: identity\nDotted: frequency 0.5",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=7.5,
    )


def render(
    points: list[FrequencyPoint],
    summary: dict[str, float],
    output_dir: Path,
) -> list[Path]:
    """Render the complete two-panel frequency-orientation figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # The wide two-panel layout separates the global density relationship from
    # the near-half-frequency categories that would otherwise be obscured.
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.4, 4.6),
        gridspec_kw={"width_ratios": (1.05, 1.0)},
        constrained_layout=True,
    )

    add_overall_panel(axes[0], points, summary)
    add_edge_panel(axes[1], points, summary)

    stem = output_dir / "Supplementary_Figure_S1_frequency_orientation"
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
    """Record checksums for the plotting tables, script and outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "supplementary_figure_s1_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]
    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render the allele-orientation support figure."""
    args = parse_args()
    points_path = args.points.resolve()
    summary_path = args.summary.resolve()

    # Both release-safe inputs are validated before any plotting object is created.
    points = read_points(points_path)
    summary = read_summary(summary_path)
    outputs = render(points, summary, args.output_dir.resolve())
    provenance_record = write_provenance(
        [points_path, summary_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S1_GENERATION=PASS")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
