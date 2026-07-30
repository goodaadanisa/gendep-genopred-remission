#!/usr/bin/env python3
"""Generate Figure 2: target genotype reconstruction and marker retention.

Purpose
-------
Render a publication figure summarising target genotype reconstruction,
scoring-marker retention, orientation diagnostics and PLINK2 round-trip
validation.

Inputs
------
- results/figure_data/figure_02_target_reconstruction_retention.tsv

Outputs
-------
- results/figures/Figure_2_target_reconstruction_and_marker_retention.png
- results/figures/Figure_2_target_reconstruction_and_marker_retention.svg
- results/figures/Figure_2_target_reconstruction_and_marker_retention.pdf
- results/manifests/figure_02_validation.txt

Scope
-----
The script visualises release-safe aggregate analysis outputs only. It does not access controlled genotype inputs, generate PRS, or build PLINK2 target files.
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
DEFAULT_INPUT: Final = PROJECT_ROOT / "results" / "figure_data" / "figure_02_target_reconstruction_retention.tsv"
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402
from matplotlib.ticker import FuncFormatter, MaxNLocator  # noqa: E402

INK = "#1B263B"
MUTED_TEXT = "#4B5563"
LINE = "#6B7280"
WHITE = "#FFFFFF"
NAVY = "#0B4F94"
BLUE = "#2C7DD1"
LIGHT_BLUE = "#DCEAF8"
MID_BLUE = "#7EB3E6"
ORANGE = "#F28E2B"
PANEL_LINE = "#CBD5E1"
GREY_FILL = "#F8FAFC"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.unicode_minus": True,
        "figure.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "axes.facecolor": WHITE,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

PROCESSING_IDS: Final = [
    "supplied_variants",
    "autosomal_snvs",
    "official_genopred_reference_overlap",
    "orientation_resolved_target",
    "prscs_compatible_scoring_variants",
]
RETAINED_BAR_IDS: Final = [
    "official_genopred_reference_matches",
    "prscs_compatible_subset",
    "conservative_orientation_sensitivity_retained",
    "three_resource_diagnostic_subset",
]
DIAGNOSTIC_BAR_IDS: Final = [
    "same_order_alleles",
    "swapped_order_alleles_corrected",
    "absent_from_prscs_after_orientation",
]
VALIDATION_IDS: Final = [
    "validated_chromosomes",
    "validated_participants",
    "validated_variants",
    "roundtrip_genotype_comparisons",
    "validation_mismatches",
]
BAR_COLOURS: Final = {
    "official_genopred_reference_matches": NAVY,
    "prscs_compatible_subset": BLUE,
    "conservative_orientation_sensitivity_retained": MID_BLUE,
    "three_resource_diagnostic_subset": LIGHT_BLUE,
    "same_order_alleles": LIGHT_BLUE,
    "swapped_order_alleles_corrected": ORANGE,
    "absent_from_prscs_after_orientation": "#FFB55A",
}
FILL_GRADIENT: Final = [NAVY, "#1566B8", "#2886E0", "#61A7E8", "#D9EAF8"]


@dataclass(frozen=True)
class Metric:
    """Represent Metric as an immutable workflow record."""
    section: str
    order_index: int
    metric_id: str
    label: str
    value: int
    percent_of_start: float | None
    drop_from_previous: int | None
    drop_percent_from_previous: float | None
    note: str


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the target-reconstruction figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Release-safe Figure 2 metric table.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Directory for figure files.")
    parser.add_argument(
        "--provenance-dir",
        type=Path,
        default=DEFAULT_PROVENANCE_DIR,
        help="Directory for validation output.",
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_optional_float(value: str) -> float | None:
    """Parse optional float from command-line or configuration input."""
    value = value.strip()
    return None if value == "" else float(value)


def parse_optional_int(value: str) -> int | None:
    """Parse optional int from command-line or configuration input."""
    value = value.strip()
    return None if value == "" else int(value)


def read_metrics(path: Path) -> dict[str, Metric]:
    """Read metrics from disk and validate its basic structure."""
    if not path.is_file():
        raise FileNotFoundError(f"Figure 2 metric table not found: {path}")

    expected_columns = {
        "section",
        "order_index",
        "metric_id",
        "label",
        "value",
        "percent_of_start",
        "drop_from_previous",
        "drop_percent_from_previous",
        "note",
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != expected_columns:
            raise ValueError(
                "Unexpected Figure 2 metric-table schema. "
                f"Expected {sorted(expected_columns)}; observed {reader.fieldnames}."
            )

        metrics: dict[str, Metric] = {}
        for row_number, row in enumerate(reader, start=2):
            metric_id = row["metric_id"].strip()
            if metric_id in metrics:
                raise ValueError(f"Duplicate metric_id at row {row_number}: {metric_id}")
            metrics[metric_id] = Metric(
                section=row["section"].strip(),
                order_index=int(row["order_index"]),
                metric_id=metric_id,
                label=row["label"].strip(),
                value=int(row["value"]),
                percent_of_start=parse_optional_float(row["percent_of_start"]),
                drop_from_previous=parse_optional_int(row["drop_from_previous"]),
                drop_percent_from_previous=parse_optional_float(row["drop_percent_from_previous"]),
                note=row["note"].strip(),
            )

    required_ids = set(PROCESSING_IDS) | set(RETAINED_BAR_IDS) | set(DIAGNOSTIC_BAR_IDS) | set(VALIDATION_IDS)
    missing_ids = sorted(required_ids - set(metrics))
    if missing_ids:
        raise ValueError(f"Missing Figure 2 metrics: {missing_ids}")
    return metrics


def human_count(value: int) -> str:
    """Format an integer count with thousands separators."""
    return f"{value:,}"


def negative_change_text(count: int | None, percent: float | None) -> str:
    """Format a retained-count reduction for figure annotation."""
    if count is None or percent is None:
        return ""
    return f"−{human_count(count)}\n(−{percent:.2f}%)"


def draw_processing_panel(ax: plt.Axes, metrics: dict[str, Metric]) -> None:
    """Draw the target-reconstruction attrition panel."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.005, 0.985, "A", fontsize=15, fontweight="bold", color=INK, va="top")
    ax.text(
        0.165,
        0.985,
        "Variant processing and retention",
        fontsize=14.2,
        fontweight="bold",
        color=INK,
        va="top",
    )

    # A dashed guide and stage numerals encode process order without non-data artwork.
    ax.plot([0.085, 0.085], [0.155, 0.855], color=PANEL_LINE, linewidth=1.0, linestyle=(0, (3, 3)))
    ax.text(0.025, 0.50, "Processing stage", rotation=90, va="center", ha="center", fontsize=11.8, color=INK)

    y_positions = [0.805, 0.655, 0.505, 0.355, 0.205]
    x_box = 0.450
    box_widths = [0.380, 0.370, 0.350, 0.330, 0.310]
    box_height = 0.096

    label_replacements = {
        "Supplied rsID-only genotype matrix": "Supplied rsID-only\ngenotype matrix",
        "Reconstructed autosomal SNVs": "Reconstructed\nautosomal SNVs",
        "Official GenoPred reference overlap": "Official GenoPred\nreference overlap",
        "Orientation-resolved target markers": "Orientation-resolved\ntarget markers",
        "PRS-CS-compatible scoring variants": "PRS-CS-compatible\nscoring variants",
    }

    for index, (metric_id, y_pos, box_width, fill_colour) in enumerate(
        zip(PROCESSING_IDS, y_positions, box_widths, FILL_GRADIENT, strict=True),
        start=1,
    ):
        metric = metrics[metric_id]
        label = label_replacements.get(metric.label, metric.label)
        centre_y = y_pos + box_height / 2

        ax.text(0.055, centre_y, str(index), ha="center", va="center", fontsize=16, fontweight="bold", color=NAVY)
        ax.text(0.12, centre_y, label, ha="left", va="center", fontsize=10.7, color=INK, linespacing=1.18)

        patch = FancyBboxPatch(
            (x_box, y_pos),
            box_width,
            box_height,
            boxstyle="round,pad=0.012,rounding_size=0.012",
            linewidth=0.0,
            edgecolor=fill_colour,
            facecolor=fill_colour,
        )
        ax.add_patch(patch)

        text_colour = WHITE if index < 5 else INK
        ax.text(
            x_box + box_width / 2,
            y_pos + box_height * 0.62,
            human_count(metric.value),
            ha="center",
            va="center",
            fontsize=18.2,
            fontweight="bold",
            color=text_colour,
        )
        ax.text(
            x_box + box_width / 2,
            y_pos + box_height * 0.28,
            f"({metric.percent_of_start:.2f}% of start)",
            ha="center",
            va="center",
            fontsize=10.5,
            color=text_colour,
        )

        if index < len(PROCESSING_IDS):
            next_y = y_positions[index]
            arrow_x = x_box + box_width / 2
            arrow = FancyArrowPatch(
                (arrow_x, y_pos - 0.004),
                (arrow_x, next_y + box_height + 0.004),
                arrowstyle="simple",
                mutation_scale=13,
                linewidth=0.0,
                color=LINE,
                alpha=0.85,
            )
            ax.add_patch(arrow)

            # Difference labels sit in dedicated whitespace and no longer collide
            # with the processing boxes or the right-hand panel.
            ax.text(
                0.855,
                (y_pos + next_y + box_height) / 2,
                negative_change_text(metric.drop_from_previous, metric.drop_percent_from_previous),
                ha="left",
                va="center",
                fontsize=9.9,
                color=INK,
                linespacing=1.15,
            )

    first_value = metrics["supplied_variants"].value
    last_value = metrics["prscs_compatible_scoring_variants"].value
    total_reduction = first_value - last_value
    retention_from_final_target = 100 * last_value / metrics["orientation_resolved_target"].value

    # The summary box concentrates the end-point retention metrics used in the text.
    summary_box = FancyBboxPatch(
        (0.105, 0.025),
        0.79,
        0.105,
        boxstyle="round,pad=0.012,rounding_size=0.008",
        linewidth=1.0,
        edgecolor=LINE,
        facecolor=GREY_FILL,
        linestyle=(0, (3, 3)),
    )
    ax.add_patch(summary_box)
    ax.text(
        0.50,
        0.088,
        f"Total reduction: -{human_count(total_reduction)} variants (-{100 * total_reduction / first_value:.2f}% of start)",
        ha="center",
        va="center",
        fontsize=10.6,
        color=INK,
    )
    ax.text(
        0.50,
        0.050,
        f"Retention from orientation-resolved to PRS-CS: {retention_from_final_target:.2f}%",
        ha="center",
        va="center",
        fontsize=10.6,
        color=INK,
    )


def draw_diagnostic_panel(ax: plt.Axes, metrics: dict[str, Metric]) -> None:
    """Draw the orientation and compatibility diagnostics panel."""
    ax.set_xlim(0, 170000)
    ax.set_ylim(-0.65, 9.55)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(INK)
    ax.spines["bottom"].set_color(INK)
    ax.tick_params(axis="both", colors=INK, labelsize=8.9)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=6))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{int(value):,}"))
    ax.grid(False)

    ax.text(-0.29, 1.018, "B", transform=ax.transAxes, fontsize=15, fontweight="bold", color=INK, va="bottom")
    ax.text(
        0.47,
        1.018,
        "Compatibility and validation diagnostics",
        transform=ax.transAxes,
        fontsize=14.2,
        fontweight="bold",
        color=INK,
        va="bottom",
        ha="center",
    )

    display_ids = list(RETAINED_BAR_IDS) + list(DIAGNOSTIC_BAR_IDS)
    y_positions = [8.15, 6.98, 5.81, 4.64, 2.66, 1.49, 0.32]
    y_labels = [
        "Official GenoPred\nreference matches",
        "PRS-CS-compatible\nscoring variants",
        "Conservative sensitivity\nretained",
        "Three-resource\ndiagnostic subset",
        "Same-order alleles",
        "Swapped-order\nalleles corrected",
        "Absent from PRS-CS\nafter orientation",
    ]

    for metric_id, y_pos in zip(display_ids, y_positions, strict=True):
        metric = metrics[metric_id]
        ax.barh(y_pos, metric.value, height=0.68, color=BAR_COLOURS[metric_id], edgecolor="none")
        ax.text(
            metric.value + 5000,
            y_pos,
            human_count(metric.value),
            va="center",
            ha="left",
            fontsize=9.8,
            color=INK,
            fontweight="bold",
            clip_on=False,
        )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(y_labels, fontsize=9.2, linespacing=1.15)
    ax.set_xlabel("Markers (count)", color=INK, labelpad=8)
    ax.axhline(3.72, color=LINE, linewidth=0.9, linestyle=(0, (4, 4)))
    ax.text(1300, 9.05, "Retained marker sets", fontsize=10.1, color=INK, fontweight="semibold", ha="left")
    ax.text(1300, 3.12, "Orientation diagnostics", fontsize=10.1, color=INK, fontweight="semibold", ha="left")


def draw_validation_box(ax: plt.Axes, metrics: dict[str, Metric]) -> None:
    """Draw one exact-validation summary box."""
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    box = FancyBboxPatch(
        (0.10, 0.11),
        0.80,
        0.77,
        boxstyle="round,pad=0.018,rounding_size=0.028",
        linewidth=1.4,
        edgecolor="#1D4ED8",
        facecolor=GREY_FILL,
    )
    ax.add_patch(box)

    ax.text(0.50, 0.76, "PLINK2 validation passed", ha="center", va="center", fontsize=15.4, fontweight="bold", color=INK)
    summary_lines = [
        f"{metrics['validated_chromosomes'].value} chromosomes",
        f"{metrics['validated_participants'].value} participants",
        f"{human_count(metrics['validated_variants'].value)} variants",
        f"{human_count(metrics['roundtrip_genotype_comparisons'].value)} genotype comparisons",
        f"{metrics['validation_mismatches'].value} mismatches",
    ]
    for offset, line in enumerate(summary_lines):
        ax.text(0.50, 0.60 - offset * 0.095, line, ha="center", va="center", fontsize=11.7, color=INK)


def create_figure(metrics: dict[str, Metric]) -> plt.Figure:
    """Assemble the complete target-reconstruction figure."""
    fig = plt.figure(figsize=(14.0, 10.5), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.12, 1.0],
        height_ratios=[1.0, 0.39],
        left=0.055,
        right=0.965,
        top=0.855,
        bottom=0.075,
        wspace=0.20,
        hspace=0.10,
    )

    ax_left = fig.add_subplot(grid[:, 0])
    ax_right_top = fig.add_subplot(grid[0, 1])
    ax_right_bottom = fig.add_subplot(grid[1, 1])

    draw_processing_panel(ax_left, metrics)
    draw_diagnostic_panel(ax_right_top, metrics)
    draw_validation_box(ax_right_bottom, metrics)

    # A figure-level title keeps the export interpretable as a standalone item.
    fig.suptitle(
        "Target genotype reconstruction and scoring-marker retention",
        fontsize=18.5,
        fontweight="bold",
        color=INK,
        y=0.975,
    )
    fig.add_artist(Line2D([0.545, 0.545], [0.10, 0.845], transform=fig.transFigure, color=PANEL_LINE, linewidth=1.0))
    return fig


def write_validation(path: Path, metrics: dict[str, Metric], outputs: list[Path]) -> None:
    """Write validation to a deterministic workflow output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        handle.write("figure_02_target_reconstruction_retention=PASS\n")
        handle.write(f"input_table={DEFAULT_INPUT}\n")
        handle.write(f"processing_stage_count={len(PROCESSING_IDS)}\n")
        handle.write(f"retained_bar_count={len(RETAINED_BAR_IDS)}\n")
        handle.write(f"diagnostic_bar_count={len(DIAGNOSTIC_BAR_IDS)}\n")
        handle.write(f"validation_box_values={len(VALIDATION_IDS)}\n")
        handle.write(f"total_start_variants={metrics['supplied_variants'].value}\n")
        handle.write(f"prscs_compatible_variants={metrics['prscs_compatible_scoring_variants'].value}\n")
        handle.write(f"conservative_sensitivity_retained={metrics['conservative_orientation_sensitivity_retained'].value}\n")
        handle.write(f"roundtrip_genotype_comparisons={metrics['roundtrip_genotype_comparisons'].value}\n")
        handle.write(f"validation_mismatches={metrics['validation_mismatches'].value}\n")
        for output_path in outputs:
            handle.write(f"sha256_{output_path.name}={sha256(output_path)}\n")


def main() -> None:
    """Render the target-reconstruction and validation figure."""
    args = parse_args()
    metrics = read_metrics(args.input)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.provenance_dir.mkdir(parents=True, exist_ok=True)

    figure = create_figure(metrics)
    output_base = args.output_dir / "Figure_2_target_reconstruction_and_marker_retention"
    png_path = output_base.with_suffix(".png")
    svg_path = output_base.with_suffix(".svg")
    pdf_path = output_base.with_suffix(".pdf")

    figure.savefig(png_path, dpi=300)
    figure.savefig(svg_path)
    figure.savefig(pdf_path)
    plt.close(figure)

    validation_path = args.provenance_dir / "figure_02_validation.txt"
    write_validation(validation_path, metrics, [png_path, svg_path, pdf_path])

    print(f"Wrote: {png_path}")
    print(f"Wrote: {svg_path}")
    print(f"Wrote: {pdf_path}")
    print(f"Wrote: {validation_path}")


if __name__ == "__main__":
    main()
