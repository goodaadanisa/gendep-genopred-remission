#!/usr/bin/env python3
"""Generate Supplementary Figure S3: PRS distributions and correlations.

Purpose
-------
Visualise within-score-standardised participant PRS distributions and the
lower-triangular Pearson correlation matrix for the eight Bayesian
polygenic scores.

Inputs
------
- results/figure_data/supplementary_figure_s3_violin_density.tsv
- results/figure_data/supplementary_figure_s3_correlation_matrix.tsv
- results/figure_data/supplementary_figure_s3_summary.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S3_PRS_distributions_correlations.png
- results/figures/Supplementary_Figure_S3_PRS_distributions_correlations.svg
- results/figures/Supplementary_Figure_S3_PRS_distributions_correlations.pdf
- results/manifests/supplementary_figure_s3_sha256.txt

Scope
-----
The script reads release-safe aggregate density and correlation tables only. It does
not access participant identifiers or participant-level PRS values, calculate
new scores or fit predictive models.
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
DEFAULT_DENSITY_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s3_violin_density.tsv"
)
DEFAULT_CORRELATION_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s3_correlation_matrix.tsv"
)
DEFAULT_SUMMARY_INPUT: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s3_summary.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Position, symmetric violin geometry, quartile lines and annotated heatmap
# cells jointly encode the figure, so interpretation does not rely on colour
# alone.
VIOLIN_FACE: Final = "#56B4E9"
VIOLIN_EDGE: Final = "#4D4D4D"
QUARTILE_COLOUR: Final = "#666666"
MEDIAN_COLOUR: Final = "#FFFFFF"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 10.8,
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

TRAITS: Final = (
    "MDD",
    "ANX",
    "BIP",
    "SCZ",
    "NEUR",
    "INSOM",
    "SWB",
    "EA",
)
EXPECTED_PARTICIPANTS: Final = 430
EXPECTED_GRID_POINTS: Final = 401
EXPECTED_MAXIMUM_ABSOLUTE_OFF_DIAGONAL: Final = 0.444443895393558


@dataclass(frozen=True)
class DensityCurve:
    """One aggregate PRS density curve and quartile summary."""

    trait_order: int
    trait: str
    x_values: np.ndarray
    density_values: np.ndarray
    q1: float
    median: float
    q3: float
    minimum: float
    maximum: float
    participants: int


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the PRS-distribution figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--density", type=Path, default=DEFAULT_DENSITY_INPUT)
    parser.add_argument(
        "--correlations",
        type=Path,
        default=DEFAULT_CORRELATION_INPUT,
    )
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY_INPUT)
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


def read_density_curves(path: Path) -> list[DensityCurve]:
    """Read and validate the aggregate violin-density table."""

    if not path.is_file():
        raise FileNotFoundError(f"Density source table not found: {path}")

    required_columns = {
        "trait_order",
        "trait",
        "grid_order",
        "within_score_sd",
        "density",
        "q1",
        "median",
        "q3",
        "minimum",
        "maximum",
        "participants",
    }

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError("Unexpected Supplementary Figure S3 density schema.")

        grouped: dict[str, list[dict[str, str]]] = {trait: [] for trait in TRAITS}
        for row_number, row in enumerate(reader, start=2):
            trait = row["trait"].strip()
            if trait not in grouped:
                raise ValueError(f"Unexpected trait on row {row_number}: {trait}")
            grouped[trait].append(row)

    curves: list[DensityCurve] = []
    for expected_order, trait in enumerate(TRAITS, start=1):
        rows = grouped[trait]
        rows.sort(key=lambda row: int(row["grid_order"]))
        if len(rows) != EXPECTED_GRID_POINTS:
            raise ValueError(
                f"Expected {EXPECTED_GRID_POINTS} density points for {trait}; "
                f"observed {len(rows)}."
            )
        if [int(row["grid_order"]) for row in rows] != list(
            range(1, EXPECTED_GRID_POINTS + 1)
        ):
            raise ValueError(f"Density grid order is incomplete for {trait}.")

        trait_orders = {int(row["trait_order"]) for row in rows}
        participants = {int(row["participants"]) for row in rows}
        if trait_orders != {expected_order}:
            raise ValueError(f"Unexpected trait order for {trait}.")
        if participants != {EXPECTED_PARTICIPANTS}:
            raise ValueError(f"Unexpected participant count for {trait}.")

        x_values = np.array(
            [float(row["within_score_sd"]) for row in rows],
            dtype=float,
        )
        density_values = np.array(
            [float(row["density"]) for row in rows],
            dtype=float,
        )
        if np.any(np.diff(x_values) <= 0.0):
            raise ValueError(f"Density grid is not increasing for {trait}.")
        if np.any(~np.isfinite(density_values)) or np.any(density_values < 0.0):
            raise ValueError(f"Invalid density values for {trait}.")

        # Summary fields are repeated on every grid row so one missing or
        # conflicting value is detected before the curve is rendered.
        q1_values = {float(row["q1"]) for row in rows}
        median_values = {float(row["median"]) for row in rows}
        q3_values = {float(row["q3"]) for row in rows}
        minimum_values = {float(row["minimum"]) for row in rows}
        maximum_values = {float(row["maximum"]) for row in rows}

        if not all(
            len(values) == 1
            for values in (
                q1_values,
                median_values,
                q3_values,
                minimum_values,
                maximum_values,
            )
        ):
            raise ValueError(f"Repeated summary values disagree for {trait}.")

        q1 = q1_values.pop()
        median = median_values.pop()
        q3 = q3_values.pop()
        minimum = minimum_values.pop()
        maximum = maximum_values.pop()
        if not minimum <= q1 <= median <= q3 <= maximum:
            raise ValueError(f"Invalid quartile order for {trait}.")

        curves.append(
            DensityCurve(
                trait_order=expected_order,
                trait=trait,
                x_values=x_values,
                density_values=density_values,
                q1=q1,
                median=median,
                q3=q3,
                minimum=minimum,
                maximum=maximum,
                participants=EXPECTED_PARTICIPANTS,
            )
        )

    return curves


def read_correlation_matrix(path: Path) -> np.ndarray:
    """Read and validate the release-safe long-format correlation table."""

    if not path.is_file():
        raise FileNotFoundError(f"Correlation source table not found: {path}")

    required_columns = {
        "row_order",
        "row_trait",
        "column_order",
        "column_trait",
        "pearson_r",
    }
    matrix = np.full((len(TRAITS), len(TRAITS)), np.nan, dtype=float)
    observed_pairs: set[tuple[int, int]] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError("Unexpected Supplementary Figure S3 correlation schema.")

        for row_number, row in enumerate(reader, start=2):
            row_order = int(row["row_order"])
            column_order = int(row["column_order"])
            pair = (row_order, column_order)
            if pair in observed_pairs:
                raise ValueError(f"Duplicate correlation cell: {pair}")
            observed_pairs.add(pair)

            if not 1 <= row_order <= len(TRAITS):
                raise ValueError(f"Invalid row order on source row {row_number}.")
            if not 1 <= column_order <= len(TRAITS):
                raise ValueError(f"Invalid column order on source row {row_number}.")
            if row["row_trait"] != TRAITS[row_order - 1]:
                raise ValueError(f"Row trait mismatch on source row {row_number}.")
            if row["column_trait"] != TRAITS[column_order - 1]:
                raise ValueError(f"Column trait mismatch on source row {row_number}.")

            value = float(row["pearson_r"])
            if not np.isfinite(value) or not -1.0 <= value <= 1.0:
                raise ValueError(f"Invalid Pearson correlation on row {row_number}.")
            matrix[row_order - 1, column_order - 1] = value

    if len(observed_pairs) != len(TRAITS) ** 2:
        raise ValueError("The correlation matrix is incomplete.")
    if not np.allclose(matrix, matrix.T, rtol=0.0, atol=1e-15):
        raise ValueError("The correlation matrix is not symmetric.")
    if not np.allclose(np.diag(matrix), 1.0, rtol=0.0, atol=1e-15):
        raise ValueError("The correlation diagonal is not exactly one.")

    off_diagonal = matrix.copy()
    np.fill_diagonal(off_diagonal, np.nan)
    maximum_absolute = float(np.nanmax(np.abs(off_diagonal)))
    if not np.isclose(
        maximum_absolute,
        EXPECTED_MAXIMUM_ABSOLUTE_OFF_DIAGONAL,
        rtol=0.0,
        atol=1e-15,
    ):
        raise ValueError(
            f"Unexpected maximum absolute off-diagonal correlation: "
            f"{maximum_absolute:.17g}"
        )
    return matrix


def read_summary(path: Path) -> dict[str, str]:
    """Read and validate the aggregate cohort summary."""

    if not path.is_file():
        raise FileNotFoundError(f"Summary source table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows or set(rows[0]) != {"metric", "value"}:
        raise ValueError("Unexpected Supplementary Figure S3 summary schema.")

    summary = {row["metric"]: row["value"] for row in rows}
    required_metrics = {
        "participants",
        "traits",
        "maximum_absolute_off_diagonal_r",
        "maximum_absolute_pair",
        "minimum_standardised_value",
        "maximum_standardised_value",
    }
    if set(summary) != required_metrics:
        raise ValueError("The Supplementary Figure S3 summary is incomplete.")
    if int(summary["participants"]) != EXPECTED_PARTICIPANTS:
        raise ValueError("Unexpected participant count in Figure S3 summary.")
    if int(summary["traits"]) != len(TRAITS):
        raise ValueError("Unexpected trait count in Figure S3 summary.")
    return summary


def clean_axes(ax: plt.Axes, *, grid_axis: str | None = None) -> None:
    """Apply the shared publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if grid_axis is not None:
        ax.grid(axis=grid_axis, alpha=0.20, linewidth=0.7)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel label outside the plotting area."""

    ax.text(
        -0.12,
        1.055,
        label,
        transform=ax.transAxes,
        fontsize=11.8,
        fontweight="bold",
        ha="left",
        va="top",
    )


def add_distribution_panel(
    ax: plt.Axes,
    curves: list[DensityCurve],
) -> None:
    """Draw eight horizontal aggregate violin distributions."""

    y_positions = np.arange(len(curves), dtype=float)

    for y_position, curve in zip(y_positions, curves, strict=True):
        # Each density is scaled to the same maximum half-width, making shape
        # comparable across traits without using density magnitude as a second
        # quantitative variable.
        maximum_density = float(curve.density_values.max())
        if maximum_density <= 0.0:
            raise ValueError(f"Non-positive maximum density for {curve.trait}.")
        half_width = 0.40 * curve.density_values / maximum_density

        ax.fill_between(
            curve.x_values,
            y_position - half_width,
            y_position + half_width,
            facecolor=VIOLIN_FACE,
            edgecolor=VIOLIN_EDGE,
            linewidth=0.8,
            alpha=0.88,
        )

        # Quartile and median lines are clipped to the local violin width,
        # matching the original presentation while remaining code-generated.
        for quantile, line_style, colour, line_width in (
            (curve.q1, ":", QUARTILE_COLOUR, 0.9),
            (curve.median, "--", MEDIAN_COLOUR, 1.1),
            (curve.q3, ":", QUARTILE_COLOUR, 0.9),
        ):
            local_width = float(
                np.interp(
                    quantile,
                    curve.x_values,
                    half_width,
                )
            )
            ax.plot(
                [quantile, quantile],
                [y_position - local_width, y_position + local_width],
                linestyle=line_style,
                color=colour,
                linewidth=line_width,
                zorder=4,
            )

    ax.axvline(0.0, linestyle="--", linewidth=0.8, alpha=0.55)
    ax.set_yticks(y_positions, [curve.trait for curve in curves])
    ax.invert_yaxis()
    ax.set_xlim(-4.25, 4.25)
    ax.set_xticks((-4, -2, 0, 2, 4))
    ax.set_xlabel("Within-score standard deviation")
    ax.set_ylabel("")
    ax.set_title(
        "Participant PRS distributions (n=430)",
        loc="left",
        pad=8,
    )
    clean_axes(ax, grid_axis="x")
    add_panel_label(ax, "A")


def add_correlation_panel(
    ax: plt.Axes,
    matrix: np.ndarray,
) -> None:
    """Draw the lower-triangular annotated Pearson correlation matrix."""

    # Upper-triangle cells are masked because the matrix is symmetric. The
    # diagonal remains visible to establish the 1.00 self-correlation scale.
    mask = np.triu(np.ones_like(matrix, dtype=bool), k=1)
    masked_matrix = np.ma.array(matrix, mask=mask)

    image = ax.imshow(
        masked_matrix,
        cmap="RdBu_r",
        vmin=-0.5,
        vmax=0.5,
        interpolation="none",
        aspect="equal",
    )

    # White cell boundaries retain the table-like reading order and make the
    # lower-triangular structure explicit.
    ax.set_xticks(np.arange(-0.5, len(TRAITS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(TRAITS), 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=1.0)
    ax.tick_params(which="minor", bottom=False, left=False)

    ax.set_xticks(np.arange(len(TRAITS)), TRAITS, rotation=45, ha="right")
    ax.set_yticks(np.arange(len(TRAITS)), TRAITS)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Pearson correlation matrix", loc="center", pad=8)

    # Numerical labels make the exact correlations accessible independently
    # of the diverging colour scale.
    for row_index in range(len(TRAITS)):
        for column_index in range(row_index + 1):
            value = float(matrix[row_index, column_index])
            text_colour = "white" if abs(value) >= 0.33 else "#222222"
            ax.text(
                column_index,
                row_index,
                f"{value:.2f}",
                ha="center",
                va="center",
                fontsize=7.7,
                color=text_colour,
            )

    for spine in ax.spines.values():
        spine.set_visible(False)
    add_panel_label(ax, "B")

    colour_bar = ax.figure.colorbar(
        image,
        ax=ax,
        fraction=0.046,
        pad=0.05,
        shrink=0.86,
    )
    colour_bar.set_label("Pearson r")
    colour_bar.outline.set_visible(False)
    colour_bar.ax.tick_params(labelsize=7.8)


def render(
    curves: list[DensityCurve],
    matrix: np.ndarray,
    summary: dict[str, str],
    output_dir: Path,
) -> list[Path]:
    """Render the complete two-panel PRS distribution figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # The asymmetric panel widths preserve readable violin distributions while
    # allowing square heatmap cells and a dedicated colour bar.
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(8.45, 4.95),
        gridspec_kw={"width_ratios": (0.90, 1.22)},
    )
    fig.subplots_adjust(
        left=0.085,
        right=0.955,
        top=0.90,
        bottom=0.16,
        wspace=0.30,
    )

    add_distribution_panel(axes[0], curves)
    add_correlation_panel(axes[1], matrix)

    # This statement reinforces the descriptive scope without placing
    # inferential claims inside the figure.
    fig.text(
        0.085,
        0.035,
        "Within-score standardisation was used for display only; "
        f"largest absolute off-diagonal r = "
        f"{float(summary['maximum_absolute_off_diagonal_r']):.3f} "
        f"({summary['maximum_absolute_pair']}).",
        ha="left",
        va="bottom",
        fontsize=7.6,
    )

    stem = output_dir / "Supplementary_Figure_S3_PRS_distributions_correlations"
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
    """Record checksums for aggregate inputs, script and rendered outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "supplementary_figure_s3_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n"
            )
    return record


def main() -> None:
    """Render the final PRS distribution figure."""
    args = parse_args()
    density_path = args.density.resolve()
    correlation_path = args.correlations.resolve()
    summary_path = args.summary.resolve()

    # All aggregate inputs are validated before the first plotting object is
    # created.
    curves = read_density_curves(density_path)
    matrix = read_correlation_matrix(correlation_path)
    summary = read_summary(summary_path)

    outputs = render(
        curves,
        matrix,
        summary,
        args.output_dir.resolve(),
    )
    provenance_record = write_provenance(
        [density_path, correlation_path, summary_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S3_GENERATION=PASS")
    print("participant_identifiers_accessed=FALSE")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
