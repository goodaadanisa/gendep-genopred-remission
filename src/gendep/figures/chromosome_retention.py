#!/usr/bin/env python3
"""Generate Supplementary Figure S2: chromosome-specific target retention.

Purpose
-------
Present chromosome-level matched and retained target-variant counts, expand
the small European-reference frequency-tie exclusions, and summarise the exact
PLINK2 round-trip validation.

Inputs
------
- results/figure_data/supplementary_figure_s2_chromosome_results.tsv
- results/figure_data/supplementary_figure_s2_genomewide_summary.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S2_chromosome_retention.png
- results/figures/Supplementary_Figure_S2_chromosome_retention.svg
- results/figures/Supplementary_Figure_S2_chromosome_retention.pdf
- results/manifests/supplementary_figure_s2_sha256.txt

Scope
-----
The script reads release-safe aggregate chromosome-level tables only. It does not
access participant-level genotypes, alter allele assignments, rebuild PLINK2
files or repeat the round-trip comparison.
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
DEFAULT_CHROMOSOME_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s2_chromosome_results.tsv"
)
DEFAULT_SUMMARY_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s2_genomewide_summary.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# The palette follows the approved visual reference. Bar position, outline and
# panel structure also encode the series, so the figure remains readable in
# greyscale and does not rely on colour alone.
MATCHED_COLOUR: Final = "#D9D9D9"
RETAINED_COLOUR: Final = "#0072B2"
EXCLUDED_COLOUR: Final = "#D55E00"
EDGE_COLOUR: Final = "#4D4D4D"
GRID_COLOUR: Final = "#D9D9D9"
ANNOTATION_EDGE: Final = "#4C72B0"
ANNOTATION_FACE: Final = "#F4F7FB"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 11.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 8.4,
        "ytick.labelsize": 8.4,
        "legend.fontsize": 8.4,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "axes.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

EXPECTED_CHROMOSOMES: Final = tuple(range(1, 23))
EXPECTED_TOTALS: Final = {
    "matched_variants": 147438,
    "retained_variants": 147370,
    "excluded_eur_frequency_ties": 68,
    "chromosomes": 22,
    "participants": 430,
    "participant_genotype_comparisons": 63369100,
    "participant_genotype_mismatches": 0,
}


@dataclass(frozen=True)
class ChromosomeResult:
    """One chromosome-level retention and round-trip result."""

    chromosome: int
    matched_variants: int
    retained_variants: int
    excluded_ties: int
    retention_percentage: float
    participants: int
    genotype_comparisons: int
    genotype_mismatches: int
    roundtrip_pass: bool


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the chromosome-retention figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--chromosomes",
        type=Path,
        default=DEFAULT_CHROMOSOME_INPUT,
        help="Release-safe 22-row chromosome-level audit table.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY_INPUT,
        help="Fixed genome-wide audit summary.",
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


def read_chromosome_results(path: Path) -> list[ChromosomeResult]:
    """Read and validate the fixed 22-row chromosome table."""

    if not path.is_file():
        raise FileNotFoundError(f"Chromosome source table not found: {path}")

    required_columns = {
        "chromosome",
        "matched_variants",
        "retained_variants",
        "excluded_eur_frequency_ties",
        "retention_percentage",
        "participants",
        "genotype_comparisons",
        "genotype_mismatches",
        "roundtrip_pass",
    }

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected Supplementary Figure S2 chromosome schema."
            )

        results: list[ChromosomeResult] = []
        observed_chromosomes: set[int] = set()

        # Every chromosome is validated before plotting so a malformed count
        # cannot be hidden by the small graphical difference between the bars.
        for row_number, row in enumerate(reader, start=2):
            chromosome = int(row["chromosome"])
            if chromosome in observed_chromosomes:
                raise ValueError(f"Duplicate chromosome on row {row_number}: {chromosome}")
            observed_chromosomes.add(chromosome)

            result = ChromosomeResult(
                chromosome=chromosome,
                matched_variants=int(row["matched_variants"]),
                retained_variants=int(row["retained_variants"]),
                excluded_ties=int(row["excluded_eur_frequency_ties"]),
                retention_percentage=float(row["retention_percentage"]),
                participants=int(row["participants"]),
                genotype_comparisons=int(row["genotype_comparisons"]),
                genotype_mismatches=int(row["genotype_mismatches"]),
                roundtrip_pass=row["roundtrip_pass"].upper() == "TRUE",
            )

            # The matched-minus-retained difference must exactly equal the
            # chromosome-specific European-reference frequency-tie exclusion.
            if (
                result.matched_variants - result.retained_variants
                != result.excluded_ties
            ):
                raise ValueError(
                    f"Retention arithmetic failed on source row {row_number}."
                )

            # Each chromosome was validated across the same 430 participants.
            if result.participants != EXPECTED_TOTALS["participants"]:
                raise ValueError(
                    f"Unexpected participant count on source row {row_number}."
                )
            if (
                result.retained_variants * result.participants
                != result.genotype_comparisons
            ):
                raise ValueError(
                    f"Comparison arithmetic failed on source row {row_number}."
                )
            if result.genotype_mismatches != 0 or not result.roundtrip_pass:
                raise ValueError(
                    f"Round-trip validation failed on source row {row_number}."
                )

            expected_retention = (
                100.0 * result.retained_variants / result.matched_variants
            )
            if not np.isclose(
                result.retention_percentage,
                expected_retention,
                rtol=0.0,
                atol=1e-12,
            ):
                raise ValueError(
                    f"Retention percentage is inconsistent on row {row_number}."
                )

            results.append(result)

    results.sort(key=lambda row: row.chromosome)
    if tuple(row.chromosome for row in results) != EXPECTED_CHROMOSOMES:
        raise ValueError("Chromosome rows must cover autosomes 1 through 22.")

    # Chromosome-level totals are checked independently of the genome-wide
    # summary so both aggregate inputs must agree with the target summary.
    chromosome_totals = {
        "matched_variants": sum(row.matched_variants for row in results),
        "retained_variants": sum(row.retained_variants for row in results),
        "excluded_eur_frequency_ties": sum(row.excluded_ties for row in results),
        "participant_genotype_comparisons": sum(
            row.genotype_comparisons for row in results
        ),
        "participant_genotype_mismatches": sum(
            row.genotype_mismatches for row in results
        ),
    }
    for metric, observed in chromosome_totals.items():
        if observed != EXPECTED_TOTALS[metric]:
            raise ValueError(
                f"Unexpected chromosome-summed {metric}: {observed}"
            )

    return results


def read_summary(path: Path) -> dict[str, str]:
    """Read and validate the fixed genome-wide summary."""

    if not path.is_file():
        raise FileNotFoundError(f"Genome-wide summary not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows or set(rows[0]) != {"metric", "value"}:
        raise ValueError(
            "Unexpected Supplementary Figure S2 summary schema."
        )

    summary = {row["metric"]: row["value"] for row in rows}
    required_metrics = {*EXPECTED_TOTALS, "all_chromosomes_pass"}
    if set(summary) != required_metrics:
        missing = sorted(required_metrics - set(summary))
        unexpected = sorted(set(summary) - required_metrics)
        raise ValueError(
            f"Genome-wide summary mismatch. Missing={missing}; "
            f"unexpected={unexpected}"
        )

    for metric, expected in EXPECTED_TOTALS.items():
        if int(summary[metric]) != expected:
            raise ValueError(f"Unexpected genome-wide value for {metric}.")
    if summary["all_chromosomes_pass"].upper() != "TRUE":
        raise ValueError("Not all chromosome-level round-trip checks passed.")
    return summary


def clean_axes(ax: plt.Axes) -> None:
    """Apply the common publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.8)
    ax.spines["bottom"].set_linewidth(0.8)
    ax.tick_params(width=0.8, length=3.5)
    ax.grid(
        axis="y",
        linestyle="-",
        linewidth=0.65,
        color=GRID_COLOUR,
        alpha=0.75,
    )
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel label outside the plotting area."""

    ax.text(
        -0.072,
        1.055,
        label,
        transform=ax.transAxes,
        fontsize=12.0,
        fontweight="bold",
        ha="left",
        va="top",
    )


def add_matched_retained_panel(
    ax: plt.Axes,
    results: list[ChromosomeResult],
) -> None:
    """Draw side-by-side matched and retained chromosome counts."""

    chromosomes = np.array([row.chromosome for row in results], dtype=int)
    matched = np.array([row.matched_variants for row in results], dtype=int)
    retained = np.array([row.retained_variants for row in results], dtype=int)
    x_positions = np.arange(len(results), dtype=float)
    bar_width = 0.36

    # Side-by-side bars make the source and final target independently visible;
    # their near equality communicates the very high retention rate.
    ax.bar(
        x_positions - bar_width / 2,
        matched,
        width=bar_width,
        color=MATCHED_COLOUR,
        edgecolor=EDGE_COLOUR,
        linewidth=0.55,
        label="Matched to GenoPred",
        zorder=2,
    )
    ax.bar(
        x_positions + bar_width / 2,
        retained,
        width=bar_width,
        color=RETAINED_COLOUR,
        edgecolor=EDGE_COLOUR,
        linewidth=0.45,
        label="Final retained",
        zorder=3,
    )

    ax.set_title(
        "Matched and retained target variants",
        loc="center",
        fontweight="semibold",
        pad=10,
    )
    ax.set_ylabel("Variants")
    ax.set_xticks(x_positions, chromosomes)
    ax.set_xlim(-0.8, len(results) - 0.2)
    ax.set_ylim(0, 12800)
    ax.set_yticks(np.arange(0, 12001, 2000))
    ax.legend(
        loc="upper right",
        frameon=False,
        ncol=2,
        columnspacing=1.6,
        handlelength=1.7,
        borderaxespad=0.0,
    )
    clean_axes(ax)
    add_panel_label(ax, "A")


def add_exclusion_panel(
    ax: plt.Axes,
    results: list[ChromosomeResult],
    summary: dict[str, str],
) -> None:
    """Draw chromosome-specific excluded tie counts and validation summary."""

    chromosomes = np.array([row.chromosome for row in results], dtype=int)
    excluded = np.array([row.excluded_ties for row in results], dtype=int)
    x_positions = np.arange(len(results), dtype=float)

    bars = ax.bar(
        x_positions,
        excluded,
        width=0.50,
        color=EXCLUDED_COLOUR,
        edgecolor=EXCLUDED_COLOUR,
        linewidth=0.55,
        zorder=3,
    )

    # Integer labels are retained above every chromosome, including chromosome
    # 22 where the zero count would otherwise be invisible.
    for bar, count in zip(bars, excluded, strict=True):
        y_position = max(float(count) + 0.15, 0.16)
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            y_position,
            str(int(count)),
            ha="center",
            va="bottom",
            fontsize=8.0,
        )

    ax.set_title(
        "Small retention differences made visible",
        loc="center",
        fontweight="semibold",
        pad=7,
    )
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("EUR ties excluded")
    ax.set_xticks(x_positions, chromosomes)
    ax.set_xlim(-0.8, len(results) - 0.2)
    ax.set_ylim(0, 8.3)
    ax.set_yticks(np.arange(0, 9, 1))
    clean_axes(ax)
    add_panel_label(ax, "B")

    # The audit outcome is placed within unused upper-right panel space rather
    # than below the chart, keeping the validation result visually connected
    # to the chromosome-level exclusions.
    ax.text(
        0.985,
        0.895,
        "All 22 chromosomes passed\n"
        f"{int(summary['participant_genotype_comparisons']):,} comparisons; "
        f"{int(summary['participant_genotype_mismatches']):,} mismatches",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.2,
        linespacing=1.20,
        bbox={
            "boxstyle": "round,pad=0.38",
            "facecolor": ANNOTATION_FACE,
            "edgecolor": ANNOTATION_EDGE,
            "linewidth": 0.8,
        },
    )


def render(
    results: list[ChromosomeResult],
    summary: dict[str, str],
    output_dir: Path,
) -> list[Path]:
    """Render the complete two-panel chromosome-retention figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # The stacked layout mirrors the approved reference: Panel A preserves the
    # chromosome-scale counts, while Panel B expands the otherwise invisible
    # differences onto an integer exclusion scale.
    fig, axes = plt.subplots(
        2,
        1,
        figsize=(10.6, 7.6),
        gridspec_kw={"height_ratios": (1.25, 0.78)},
    )
    fig.subplots_adjust(
        left=0.085,
        right=0.985,
        top=0.955,
        bottom=0.085,
        hspace=0.34,
    )

    add_matched_retained_panel(axes[0], results)
    add_exclusion_panel(axes[1], results, summary)

    stem = output_dir / "Supplementary_Figure_S2_chromosome_retention"
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
    """Record checksums for the plotting tables, script and outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "supplementary_figure_s2_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    # The checksum record links the fixed inputs, generator and three rendered
    # formats without publishing participant-level data.
    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n"
            )
    return record


def main() -> None:
    """Render chromosome-level target retention and validation results."""
    args = parse_args()
    chromosome_path = args.chromosomes.resolve()
    summary_path = args.summary.resolve()

    # Both fixed inputs are validated before the first plotting object is
    # created.
    results = read_chromosome_results(chromosome_path)
    summary = read_summary(summary_path)

    outputs = render(
        results,
        summary,
        args.output_dir.resolve(),
    )
    provenance_record = write_provenance(
        [chromosome_path, summary_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S2_GENERATION=PASS")
    print(f"source_chromosomes={chromosome_path}")
    print(f"source_summary={summary_path}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
