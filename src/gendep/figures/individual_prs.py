#!/usr/bin/env python3
"""Generate Figure 5: individual PRS association and predictive results.

Purpose
-------
Visualise adjusted remission associations and cross-validated discrimination
changes for each of the eight Bayesian polygenic scores.

Inputs
------
- results/figure_data/figure_05_individual_prs_results.tsv

Outputs
-------
- results/figures/Figure_5_individual_PRS_results.png
- results/figures/Figure_5_individual_PRS_results.svg
- results/figures/Figure_5_individual_PRS_results.pdf
- results/manifests/figure_05_sha256.txt

Scope
-----
The script reads a fixed aggregate source table only. It does not access
participant-level data, fit logistic models, perform cross-validation or repeat
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
    PROJECT_ROOT / "results" / "figure_data" / "figure_05_individual_prs_results.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.ticker import FixedFormatter, FixedLocator  # noqa: E402

# Trait position, interval geometry and marker fill jointly encode the results,
# so interpretation does not depend on colour alone.
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.3,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.6,
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


@dataclass(frozen=True)
class TraitResult:
    """One trait-level association and predictive-increment result."""

    trait_order: int
    trait: str
    trait_label: str
    participants: int
    events: int
    odds_ratio: float
    or_lower_95: float
    or_upper_95: float
    wald_p_value: float
    bh_fdr_p_value: float
    association_survives_fdr: bool
    auc_increment: float
    auc_lower_95: float
    auc_upper_95: float
    auc_interval_excludes_zero: bool
    auc_positive_repetitions: int
    bootstrap_replicates: int


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the individual-PRS figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Fixed aggregate individual-PRS table.",
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


def read_results(path: Path) -> list[TraitResult]:
    """Read and validate the fixed eight-row Figure 5 source table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")

    required_columns = {
        "trait_order",
        "trait",
        "trait_label",
        "participants",
        "events",
        "adjusted_odds_ratio_per_SD",
        "or_lower_95",
        "or_upper_95",
        "wald_p_value",
        "bh_fdr_p_value",
        "association_survives_fdr_0_05",
        "mean_auc_increment",
        "auc_lower_95",
        "auc_upper_95",
        "auc_interval_excludes_zero",
        "auc_positive_repetitions",
        "bootstrap_replicates",
    }

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected Figure 5 source schema. "
                f"Expected {sorted(required_columns)}; observed {reader.fieldnames}."
            )

        results: list[TraitResult] = []
        observed_traits: set[str] = set()

        # Every configured trait must occur exactly once so no result can be
        # silently omitted or duplicated in the forest plots.
        for row_number, row in enumerate(reader, start=2):
            trait = row["trait"].strip()
            if trait not in TRAIT_ORDER:
                raise ValueError(f"Unexpected trait on row {row_number}: {trait}")
            if trait in observed_traits:
                raise ValueError(f"Duplicate trait on row {row_number}: {trait}")
            observed_traits.add(trait)

            odds_ratio = float(row["adjusted_odds_ratio_per_SD"])
            or_lower = float(row["or_lower_95"])
            or_upper = float(row["or_upper_95"])
            auc_increment = float(row["mean_auc_increment"])
            auc_lower = float(row["auc_lower_95"])
            auc_upper = float(row["auc_upper_95"])
            association_survives_fdr = (
                row["association_survives_fdr_0_05"].upper() == "TRUE"
            )
            auc_excludes_zero = row["auc_interval_excludes_zero"].upper() == "TRUE"

            if not or_lower <= odds_ratio <= or_upper:
                raise ValueError(f"Odds ratio outside interval for {trait}.")
            if not auc_lower <= auc_increment <= auc_upper:
                raise ValueError(f"AUC increment outside interval for {trait}.")
            if association_survives_fdr != (float(row["bh_fdr_p_value"]) < 0.05):
                raise ValueError(f"Incorrect FDR flag for {trait}.")
            if auc_excludes_zero != (auc_upper < 0.0 or auc_lower > 0.0):
                raise ValueError(f"Incorrect AUC zero-exclusion flag for {trait}.")

            results.append(
                TraitResult(
                    trait_order=int(row["trait_order"]),
                    trait=trait,
                    trait_label=row["trait_label"].strip(),
                    participants=int(row["participants"]),
                    events=int(row["events"]),
                    odds_ratio=odds_ratio,
                    or_lower_95=or_lower,
                    or_upper_95=or_upper,
                    wald_p_value=float(row["wald_p_value"]),
                    bh_fdr_p_value=float(row["bh_fdr_p_value"]),
                    association_survives_fdr=association_survives_fdr,
                    auc_increment=auc_increment,
                    auc_lower_95=auc_lower,
                    auc_upper_95=auc_upper,
                    auc_interval_excludes_zero=auc_excludes_zero,
                    auc_positive_repetitions=int(row["auc_positive_repetitions"]),
                    bootstrap_replicates=int(row["bootstrap_replicates"]),
                )
            )

    results.sort(key=lambda row: row.trait_order)
    if [row.trait for row in results] != list(TRAIT_ORDER):
        raise ValueError("Figure 5 traits are missing or in the wrong order.")
    if any(row.participants != 430 or row.events != 166 for row in results):
        raise ValueError("Unexpected analysis-cohort counts in Figure 5.")
    if any(row.bootstrap_replicates != 10000 for row in results):
        raise ValueError("All AUC intervals must use 10,000 bootstrap replicates.")
    if any(row.association_survives_fdr for row in results):
        raise ValueError("No individual association should survive BH-FDR correction.")
    if [row.trait for row in results if row.auc_interval_excludes_zero] != ["SWB"]:
        raise ValueError("Only the SWB AUC interval should exclude zero.")
    return results


def plot_with_default_style(
    ax: plt.Axes,
    x_value: float,
    y_value: float,
    **kwargs: object,
) -> list[plt.Line2D]:
    """Draw one marker using the first default style-cycle entry.

    Resetting the property cycle keeps all trait markers visually consistent
    without defining a project-specific colour in the script.
    """

    ax.set_prop_cycle(None)
    return ax.plot(x_value, y_value, **kwargs)


def clean_axes(ax: plt.Axes) -> None:
    """Apply the shared publication treatment to one forest-plot panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", alpha=0.18, linewidth=0.7)


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


def add_association_panel(
    ax: plt.Axes,
    results: list[TraitResult],
    y_positions: np.ndarray,
) -> None:
    """Draw adjusted odds ratios and 95% Wald intervals."""

    odds_ratios = np.array([row.odds_ratio for row in results], dtype=float)
    lower_errors = np.array(
        [row.odds_ratio - row.or_lower_95 for row in results],
        dtype=float,
    )
    upper_errors = np.array(
        [row.or_upper_95 - row.odds_ratio for row in results],
        dtype=float,
    )

    # The odds-ratio scale is logarithmic so equal proportional changes around
    # the null value of one occupy equal horizontal distances.
    ax.errorbar(
        odds_ratios,
        y_positions,
        xerr=np.vstack([lower_errors, upper_errors]),
        fmt="o",
        markersize=5.5,
        elinewidth=1.2,
        capsize=3.0,
    )
    ax.axvline(1.0, linestyle="--", linewidth=0.9)
    ax.set_xscale("log")
    ax.set_xlim(0.64, 1.80)
    ax.xaxis.set_major_locator(FixedLocator((0.7, 1.0, 1.4)))
    ax.xaxis.set_major_formatter(FixedFormatter(("0.7", "1.0", "1.4")))
    ax.set_xlabel("Adjusted odds ratio per source-scale SD")
    ax.set_title("Adjusted remission association", loc="left")
    ax.set_yticks(y_positions, [row.trait for row in results])
    ax.invert_yaxis()

    # Reserve a dedicated header band above the first trait row so the compact
    # BH-FDR column heading cannot collide with the MDD value.
    ax.set_ylim(len(results) - 0.5, -0.75)
    clean_axes(ax)
    add_panel_label(ax, "A")

    # BH-FDR values are retained as a compact audit column rather than using a
    # significance colour scale; none of the eight associations survives FDR.
    ax.text(
        0.985,
        0.955,
        "BH-FDR q",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=7.4,
        fontweight="semibold",
    )
    for y_position, row in zip(y_positions, results, strict=True):
        ax.text(
            1.77,
            y_position,
            f"{row.bh_fdr_p_value:.3f}",
            ha="right",
            va="center",
            fontsize=7.4,
        )


def add_auc_panel(
    ax: plt.Axes,
    results: list[TraitResult],
    y_positions: np.ndarray,
) -> None:
    """Draw cross-validated AUC increments and paired-bootstrap intervals."""

    increments = np.array([row.auc_increment for row in results], dtype=float)
    lower_errors = np.array(
        [row.auc_increment - row.auc_lower_95 for row in results],
        dtype=float,
    )
    upper_errors = np.array(
        [row.auc_upper_95 - row.auc_increment for row in results],
        dtype=float,
    )

    # Interval segments are drawn first, allowing marker fill to separately
    # encode whether the 95% interval excludes zero.
    ax.errorbar(
        increments,
        y_positions,
        xerr=np.vstack([lower_errors, upper_errors]),
        fmt="none",
        elinewidth=1.2,
        capsize=3.0,
    )

    for y_position, row in zip(y_positions, results, strict=True):
        if row.auc_interval_excludes_zero:
            plot_with_default_style(
                ax,
                row.auc_increment,
                y_position,
                marker="o",
                markersize=5.8,
            )
        else:
            plot_with_default_style(
                ax,
                row.auc_increment,
                y_position,
                marker="o",
                markersize=5.8,
                fillstyle="none",
            )

    ax.axvline(0.0, linestyle="--", linewidth=0.9)
    ax.set_xlim(-0.055, 0.065)
    ax.set_xticks((-0.05, -0.025, 0.0, 0.025, 0.05))
    ax.set_xlabel("Mean AUC increment (positive favours individual PRS)")
    ax.set_title("Cross-validated AUC increment", loc="left")
    ax.set_yticks(y_positions, [])
    ax.invert_yaxis()

    # The same header band separates the repeat-count heading from the first
    # 0/10 annotation while preserving row alignment with Panel A.
    ax.set_ylim(len(results) - 0.5, -0.75)
    clean_axes(ax)
    add_panel_label(ax, "B")

    # The repeat count supplements the participant-level interval without
    # replacing it or treating repeated folds as independent observations.
    ax.text(
        0.985,
        0.955,
        "Positive repeats",
        transform=ax.transAxes,
        ha="right",
        va="center",
        fontsize=7.4,
        fontweight="semibold",
    )
    for y_position, row in zip(y_positions, results, strict=True):
        ax.text(
            0.062,
            y_position,
            f"{row.auc_positive_repetitions}/10",
            ha="right",
            va="center",
            fontsize=7.4,
        )


def render(
    results: list[TraitResult],
    output_dir: Path,
) -> list[Path]:
    """Render the complete two-panel individual-PRS figure."""

    output_dir.mkdir(parents=True, exist_ok=True)
    y_positions = np.arange(len(results), dtype=float)

    # The shared row order lets association and predictive results be compared
    # trait by trait without implying that the two estimands are equivalent.
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.7, 5.9),
        gridspec_kw={"width_ratios": (1.08, 1.15)},
    )
    fig.subplots_adjust(
        left=0.10,
        right=0.97,
        top=0.84,
        bottom=0.19,
        wspace=0.28,
    )

    add_association_panel(axes[0], results, y_positions)
    add_auc_panel(axes[1], results, y_positions)

    fig.suptitle(
        "Individual PRS association and predictive results",
        fontsize=13,
        fontweight="semibold",
    )

    # Compact panel notes preserve the distinction between association
    # multiplicity control and predictive interval interpretation.
    axes[0].text(
        0.0,
        -0.17,
        "No association survived BH-FDR correction.",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
    )
    axes[1].text(
        0.0,
        -0.17,
        "Open: interval includes zero; filled: interval excludes zero.",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=7.5,
    )

    stem = output_dir / "Figure_5_individual_PRS_results"
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
    """Record checksums for the source table, script and rendered outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "figure_05_sha256.txt"
    paths = [input_path, Path(__file__).resolve(), *outputs]
    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render the individual-PRS association figure."""
    args = parse_args()
    input_path = args.input.resolve()

    # Source validation is completed before any plotting object is created.
    results = read_results(input_path)
    outputs = render(results, args.output_dir.resolve())
    provenance_record = write_provenance(
        input_path,
        outputs,
        args.provenance_dir.resolve(),
    )

    print("FIGURE_05_GENERATION=PASS")
    print(f"source={input_path}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
