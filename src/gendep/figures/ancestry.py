#!/usr/bin/env python3
"""Generate Supplementary Figure S4: ancestry probability and EUR membership.

Purpose
-------
Compare maximum European-ancestry assignment probabilities between
participants retained and not retained by the strict model-based European
keep file.

Inputs
------
- a controlled 430-row deidentified probability table supplied at runtime
- results/figure_data/supplementary_figure_s4_summary.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S4_ancestry_probability.png
- results/figures/Supplementary_Figure_S4_ancestry_probability.svg
- results/figures/Supplementary_Figure_S4_ancestry_probability.pdf
- results/manifests/supplementary_figure_s4_sha256.txt

Scope
-----
The script reads a controlled deidentified probability table and an aggregate
summary. The participant-level plotting values are intentionally not committed
to the public repository. It does not access source identifiers, execute ancestry
assignment or alter strict keep-file membership.
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
DEFAULT_PROBABILITIES: Final = (
    PROJECT_ROOT / "work" / "reporting" / "supplementary_figure_s4_probabilities.tsv"
)
DEFAULT_SUMMARY: Final = (
    PROJECT_ROOT / "results" / "figure_data" / "supplementary_figure_s4_summary.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Colours reproduce the supplied reference. Line form, group position and
# boxplot geometry also encode membership, so interpretation does not rely on
# colour alone.
STRICT_COLOUR: Final = "#0072B2"
OUTSIDE_COLOUR: Final = "#D55E00"
BOX_EDGE: Final = "#8C8C8C"
GRID_COLOUR: Final = "#D9D9D9"
PLOT_SEED: Final = 20260715

GROUP_ORDER: Final = ("strict_eur_keep", "not_strict_eur_keep")
GROUP_LABELS: Final = {
    "strict_eur_keep": "Strict EUR keep",
    "not_strict_eur_keep": "Not in strict EUR keep",
}
GROUP_COLOURS: Final = {
    "strict_eur_keep": STRICT_COLOUR,
    "not_strict_eur_keep": OUTSIDE_COLOUR,
}
EXPECTED_COUNTS: Final = {
    "strict_eur_keep": 418,
    "not_strict_eur_keep": 12,
}

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 11.0,
        "axes.labelsize": 9.8,
        "xtick.labelsize": 8.4,
        "ytick.labelsize": 8.6,
        "legend.fontsize": 8.3,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass(frozen=True)
class ProbabilityRow:
    """One deidentified maximum-EUR-probability observation."""

    display_index: int
    group_key: str
    probability: float


@dataclass(frozen=True)
class GroupSummary:
    """One committed descriptive membership summary."""

    group_order: int
    group_key: str
    group_label: str
    participants: int
    minimum: float
    q1: float
    median: float
    q3: float
    maximum: float
    mean: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the ancestry figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--probabilities", type=Path, default=DEFAULT_PROBABILITIES)
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


def read_probabilities(path: Path) -> list[ProbabilityRow]:
    """Read and validate the deidentified 430-row plotting table."""

    expected_columns = {
        "display_index",
        "group_key",
        "group_label",
        "maximum_eur_assignment_probability",
    }
    if not path.is_file():
        raise FileNotFoundError(f"Probability source table not found: {path}")

    rows: list[ProbabilityRow] = []
    seen_indices: set[int] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != expected_columns:
            raise ValueError("Unexpected Supplementary Figure S4 probability schema.")

        # The display index is a deidentified plotting key, not the source
        # participant identifier.
        for row_number, row in enumerate(reader, start=2):
            display_index = int(row["display_index"])
            if display_index in seen_indices:
                raise ValueError(f"Duplicate display index on row {row_number}.")
            seen_indices.add(display_index)

            group_key = row["group_key"].strip()
            if group_key not in EXPECTED_COUNTS:
                raise ValueError(f"Unexpected group on row {row_number}: {group_key}")
            if row["group_label"].strip() != GROUP_LABELS[group_key]:
                raise ValueError(f"Unexpected group label on row {row_number}.")

            probability = float(row["maximum_eur_assignment_probability"])
            if not np.isfinite(probability) or not 0.0 <= probability <= 1.0:
                raise ValueError(f"Invalid probability on row {row_number}.")

            rows.append(
                ProbabilityRow(
                    display_index=display_index,
                    group_key=group_key,
                    probability=probability,
                )
            )

    rows.sort(key=lambda row: row.display_index)
    if [row.display_index for row in rows] != list(range(1, 431)):
        raise ValueError("Display indices must be the integers 1 through 430.")

    counts = {
        group_key: sum(row.group_key == group_key for row in rows)
        for group_key in GROUP_ORDER
    }
    if counts != EXPECTED_COUNTS:
        raise ValueError(f"Unexpected strict-EUR membership counts: {counts}")
    return rows


def read_summary(path: Path) -> dict[str, GroupSummary]:
    """Read and validate the committed two-row descriptive summary."""

    expected_columns = {
        "group_order",
        "group_key",
        "group_label",
        "participants",
        "minimum",
        "q1",
        "median",
        "q3",
        "maximum",
        "mean",
    }
    if not path.is_file():
        raise FileNotFoundError(f"Summary source table not found: {path}")

    summaries: dict[str, GroupSummary] = {}

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != expected_columns:
            raise ValueError("Unexpected Supplementary Figure S4 summary schema.")

        for row in reader:
            group_key = row["group_key"].strip()
            if group_key not in EXPECTED_COUNTS or group_key in summaries:
                raise ValueError(f"Unexpected or duplicate summary group: {group_key}")

            summary = GroupSummary(
                group_order=int(row["group_order"]),
                group_key=group_key,
                group_label=row["group_label"].strip(),
                participants=int(row["participants"]),
                minimum=float(row["minimum"]),
                q1=float(row["q1"]),
                median=float(row["median"]),
                q3=float(row["q3"]),
                maximum=float(row["maximum"]),
                mean=float(row["mean"]),
            )
            if summary.group_label != GROUP_LABELS[group_key]:
                raise ValueError(f"Unexpected summary label for {group_key}.")
            if summary.participants != EXPECTED_COUNTS[group_key]:
                raise ValueError(f"Unexpected participant count for {group_key}.")
            if not (
                0.0 <= summary.minimum <= summary.q1 <= summary.median
                <= summary.q3 <= summary.maximum <= 1.0
            ):
                raise ValueError(f"Invalid summary quantile order for {group_key}.")
            summaries[group_key] = summary

    if set(summaries) != set(GROUP_ORDER):
        raise ValueError("The Supplementary Figure S4 summary is incomplete.")
    return summaries


def clean_axes(ax: plt.Axes, *, grid_axis: str) -> None:
    """Apply the common publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.grid(axis=grid_axis, color=GRID_COLOUR, linewidth=0.65, alpha=0.75)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel label outside the plotting area."""

    ax.text(
        -0.12,
        1.06,
        label,
        transform=ax.transAxes,
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
    )


def empirical_cdf(values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return sorted observations and empirical cumulative proportions."""

    sorted_values = np.sort(values)
    cumulative = np.arange(1, sorted_values.size + 1, dtype=float)
    cumulative /= sorted_values.size
    return sorted_values, cumulative


def add_ecdf_panel(ax: plt.Axes, rows: list[ProbabilityRow]) -> None:
    """Draw membership-specific empirical cumulative distributions."""

    for group_key in GROUP_ORDER:
        values = np.array(
            [row.probability for row in rows if row.group_key == group_key],
            dtype=float,
        )
        x_values, y_values = empirical_cdf(values)

        # The leading zero anchors each curve at its minimum observation
        # without extrapolating below the observed range.
        x_values = np.r_[x_values[0], x_values]
        y_values = np.r_[0.0, y_values]

        ax.step(
            x_values,
            y_values,
            where="post",
            linewidth=1.8,
            color=GROUP_COLOURS[group_key],
            label=f"{GROUP_LABELS[group_key]} (n={EXPECTED_COUNTS[group_key]})",
        )

    ax.set_title("Empirical cumulative distributions", pad=8)
    ax.set_xlabel("Maximum EUR-assignment probability")
    ax.set_ylabel("Cumulative proportion")
    ax.set_xlim(0.74, 1.00)
    ax.set_ylim(0.00, 1.01)
    ax.set_xticks((0.75, 0.80, 0.85, 0.90, 0.95, 1.00))
    ax.set_yticks(np.linspace(0.0, 1.0, 6))
    ax.legend(loc="upper left", frameon=False)
    clean_axes(ax, grid_axis="y")
    add_panel_label(ax, "A")


def boxplot_statistics(summary: GroupSummary) -> dict[str, object]:
    """Convert a committed summary into Matplotlib boxplot statistics."""

    return {
        "label": summary.group_label,
        "whislo": summary.minimum,
        "q1": summary.q1,
        "med": summary.median,
        "q3": summary.q3,
        "whishi": summary.maximum,
        "fliers": [],
    }


def add_distribution_panel(
    ax: plt.Axes,
    rows: list[ProbabilityRow],
    summaries: dict[str, GroupSummary],
) -> None:
    """Draw horizontal boxplots and deterministic deidentified points."""

    positions = {
        "strict_eur_keep": 0.0,
        "not_strict_eur_keep": 1.0,
    }
    statistics = [boxplot_statistics(summaries[key]) for key in GROUP_ORDER]

    ax.bxp(
        statistics,
        positions=[positions[key] for key in GROUP_ORDER],
        vert=False,
        widths=0.42,
        showfliers=False,
        patch_artist=True,
        boxprops={
            "facecolor": "white",
            "edgecolor": BOX_EDGE,
            "linewidth": 1.0,
        },
        whiskerprops={"color": BOX_EDGE, "linewidth": 1.0},
        capprops={"color": BOX_EDGE, "linewidth": 1.0},
        medianprops={"color": BOX_EDGE, "linewidth": 1.1},
    )

    rng = np.random.default_rng(PLOT_SEED)

    for group_key in GROUP_ORDER:
        values = np.array(
            [row.probability for row in rows if row.group_key == group_key],
            dtype=float,
        )

        # Jitter is deterministic and purely graphical. It reveals point
        # density without encoding identity or source row order.
        width = 0.15 if group_key == "strict_eur_keep" else 0.16
        jitter = rng.uniform(-width, width, size=values.size)

        ax.scatter(
            values,
            positions[group_key] + jitter,
            s=8.0 if group_key == "strict_eur_keep" else 15.0,
            color=GROUP_COLOURS[group_key],
            alpha=0.48 if group_key == "strict_eur_keep" else 0.64,
            linewidths=0.0,
            zorder=3,
        )

    ax.set_title("Participant-level distributions", pad=8)
    ax.set_xlabel("Maximum EUR-assignment probability")
    ax.set_xlim(0.74, 1.00)
    ax.set_xticks((0.75, 0.80, 0.85, 0.90, 0.95, 1.00))
    ax.set_ylim(-0.52, 1.52)
    ax.set_yticks(
        (0.0, 1.0),
        ("Strict EUR keep\n(n=418)", "Not in strict EUR keep\n(n=12)"),
    )
    clean_axes(ax, grid_axis="x")
    add_panel_label(ax, "B")


def render(
    rows: list[ProbabilityRow],
    summaries: dict[str, GroupSummary],
    output_dir: Path,
) -> list[Path]:
    """Render the complete two-panel ancestry-probability figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # Separate axes preserve the old reference composition while allowing the
    # ECDF and horizontal distribution panels to use different y scales.
    fig = plt.figure(figsize=(9.2, 5.25))
    ecdf_ax = fig.add_axes((0.075, 0.17, 0.42, 0.70))
    distribution_ax = fig.add_axes((0.705, 0.17, 0.265, 0.70))

    add_ecdf_panel(ecdf_ax, rows)
    add_distribution_panel(distribution_ax, rows, summaries)

    # The note keeps maximum-probability assignment distinct from the stricter
    # model-based keep-file definition.
    fig.text(
        0.50,
        0.055,
        "All 430 participants had EUR as the maximum-probability class; "
        "strict model-based EUR membership retained 418.",
        ha="center",
        va="bottom",
        fontsize=8.0,
        color="#666666",
    )

    stem = output_dir / "Supplementary_Figure_S4_ancestry_probability"
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
    record = provenance_dir / "supplementary_figure_s4_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render the controlled ancestry-probability figure."""
    args = parse_args()
    probability_path = args.probabilities.resolve()
    summary_path = args.summary.resolve()

    # Both deidentified inputs are validated before any plotting object is
    # created.
    rows = read_probabilities(probability_path)
    summaries = read_summary(summary_path)

    outputs = render(rows, summaries, args.output_dir.resolve())
    provenance_record = write_provenance(
        [probability_path, summary_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S4_GENERATION=PASS")
    print("participant_identifiers_accessed=FALSE")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
