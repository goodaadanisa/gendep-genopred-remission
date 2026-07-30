#!/usr/bin/env python3
"""Generate Supplementary Figure S6: Elastic Net alpha selection.

Purpose
-------
Visualise the frequency with which each Elastic Net mixing parameter was
selected across 100 outer fits for seven predefined model specifications.

Inputs
------
- results/figure_data/supplementary_figure_s6_alpha_selection_counts.tsv
- results/figure_data/supplementary_figure_s6_provenance.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S6_alpha_selection.png
- results/figures/Supplementary_Figure_S6_alpha_selection.svg
- results/figures/Supplementary_Figure_S6_alpha_selection.pdf
- results/manifests/supplementary_figure_s6_sha256.txt

Scope
-----
The script reads release-safe aggregate model-selection counts only. It does not
access participant-level data, fit models, perform cross-validation or change
any recorded hyperparameter selection.
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
    / "supplementary_figure_s6_alpha_selection_counts.tsv"
)
DEFAULT_PROVENANCE: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s6_provenance.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

# The blue count scale and orange ridge outline follow the supplied visual
# reference. Cell values and the explicit alpha axis also encode the complete
# result independently of colour.
RIDGE_COLOUR: Final = "#E66101"
CELL_EDGE_COLOUR: Final = "#FFFFFF"
LOW_COUNT_COLOUR: Final = "#F2F4F6"
HIGH_COUNT_COLOUR: Final = "#0072B2"
TEXT_DARK: Final = "#222222"
TEXT_LIGHT: Final = "#FFFFFF"

COUNT_CMAP: Final = LinearSegmentedColormap.from_list(
    "selection_frequency_blue",
    (LOW_COUNT_COLOUR, HIGH_COUNT_COLOUR),
)

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 12.0,
        "axes.labelsize": 10.0,
        "xtick.labelsize": 8.6,
        "ytick.labelsize": 8.8,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

MODEL_ORDER: Final = (
    "primary_clinical",
    "primary_combined",
    "strict_eur_clinical",
    "strict_eur_combined",
    "prs_focused",
    "summary_clinical",
    "summary_combined",
)
ALPHA_VALUES: Final = tuple(round(value / 10, 1) for value in range(11))


@dataclass(frozen=True)
class SelectionCell:
    """One model-alpha selection-frequency cell."""

    model_order: int
    model_key: str
    model_label: str
    alpha_order: int
    selected_alpha: float
    outer_fit_count: int
    is_ridge: bool


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the alpha-selection figure command."""
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


def read_cells(path: Path) -> list[SelectionCell]:
    """Read and validate the release-safe 77-cell plotting table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")

    required_columns = {
        "model_order",
        "model_key",
        "model_label",
        "alpha_order",
        "selected_alpha",
        "outer_fit_count",
        "outer_fit_proportion",
        "is_ridge",
        "source_key",
        "source_model",
    }

    cells: list[SelectionCell] = []
    observed_pairs: set[tuple[str, float]] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected Supplementary Figure S6 plotting-table schema."
            )

        # Every model-alpha combination must occur exactly once.
        for row_number, row in enumerate(reader, start=2):
            model_key = row["model_key"]
            alpha = float(row["selected_alpha"])
            pair = (model_key, alpha)

            if model_key not in MODEL_ORDER:
                raise ValueError(
                    f"Unexpected model on row {row_number}: {model_key}"
                )
            if alpha not in ALPHA_VALUES:
                raise ValueError(
                    f"Unexpected alpha on row {row_number}: {alpha}"
                )
            if pair in observed_pairs:
                raise ValueError(f"Duplicate model-alpha pair: {pair}")
            observed_pairs.add(pair)

            count = int(row["outer_fit_count"])
            proportion = float(row["outer_fit_proportion"])
            if count < 0 or count > 100:
                raise ValueError(
                    f"Invalid selection count on row {row_number}: {count}"
                )
            if abs(proportion - count / 100.0) > 1e-12:
                raise ValueError(
                    f"Count/proportion mismatch on row {row_number}."
                )

            is_ridge = row["is_ridge"].upper() == "TRUE"
            if is_ridge != (alpha == 0.0):
                raise ValueError(
                    f"Incorrect ridge flag on row {row_number}."
                )

            cells.append(
                SelectionCell(
                    model_order=int(row["model_order"]),
                    model_key=model_key,
                    model_label=row["model_label"],
                    alpha_order=int(row["alpha_order"]),
                    selected_alpha=alpha,
                    outer_fit_count=count,
                    is_ridge=is_ridge,
                )
            )

    expected_pairs = {
        (model_key, alpha)
        for model_key in MODEL_ORDER
        for alpha in ALPHA_VALUES
    }
    if observed_pairs != expected_pairs or len(cells) != 77:
        raise ValueError("The Supplementary Figure S6 result set is incomplete.")

    for model_key in MODEL_ORDER:
        model_total = sum(
            cell.outer_fit_count
            for cell in cells
            if cell.model_key == model_key
        )
        if model_total != 100:
            raise ValueError(
                f"Selection counts for {model_key} sum to {model_total}, not 100."
            )
    return cells


def validate_provenance(path: Path) -> None:
    """Validate the four-row provenance manifest."""

    if not path.is_file():
        raise FileNotFoundError(f"Source-run table not found: {path}")

    required_columns = {
        "source_order",
        "source_key",
        "source_description",
        "source_file",
        "source_sha256",
    }
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    if not rows or set(rows[0]) != required_columns:
        raise ValueError(
            "Unexpected Supplementary Figure S6 source-run schema."
        )
    if len(rows) != 4:
        raise ValueError(f"Expected four source rows; observed {len(rows)}.")


def build_matrix(
    cells: list[SelectionCell],
) -> tuple[np.ndarray, list[str]]:
    """Convert the validated long-format cells into a seven-by-eleven matrix."""

    matrix = np.zeros((len(MODEL_ORDER), len(ALPHA_VALUES)), dtype=int)
    model_labels = [""] * len(MODEL_ORDER)

    for cell in cells:
        row_index = cell.model_order - 1
        column_index = cell.alpha_order - 1
        matrix[row_index, column_index] = cell.outer_fit_count
        model_labels[row_index] = cell.model_label

    if np.any(matrix.sum(axis=1) != 100):
        raise ValueError("At least one model row does not sum to 100 outer fits.")
    if any(not label for label in model_labels):
        raise ValueError("At least one model label is missing.")
    return matrix, model_labels


def annotation_colour(value: int, maximum: int) -> str:
    """Choose legible cell text based on the normalised count intensity."""

    return TEXT_LIGHT if value / maximum >= 0.56 else TEXT_DARK


def render(
    matrix: np.ndarray,
    model_labels: list[str],
    output_dir: Path,
) -> list[Path]:
    """Render the annotated alpha-selection heatmap."""

    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9.5, 6.2))
    maximum_count = int(matrix.max())

    # pcolormesh provides explicit cell edges without requiring Seaborn and
    # keeps the matrix fully vectorised in SVG and PDF outputs.
    mesh = ax.pcolormesh(
        np.arange(matrix.shape[1] + 1),
        np.arange(matrix.shape[0] + 1),
        matrix,
        cmap=COUNT_CMAP,
        vmin=0,
        vmax=maximum_count,
        edgecolors=CELL_EDGE_COLOUR,
        linewidth=0.75,
        shading="flat",
    )

    # Place exact counts at cell centres so the data remain accessible
    # independently of the colour scale.
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = int(matrix[row_index, column_index])
            ax.text(
                column_index + 0.5,
                row_index + 0.5,
                str(value),
                ha="center",
                va="center",
                fontsize=8.3,
                color=annotation_colour(value, maximum_count),
            )

    # The first column is the alpha=0 ridge solution. Its orange outline is a
    # visual annotation only and does not alter the heatmap scale.
    ax.add_patch(
        Rectangle(
            (0, 0),
            1,
            matrix.shape[0],
            fill=False,
            edgecolor=RIDGE_COLOUR,
            linewidth=1.8,
            clip_on=False,
        )
    )
    ax.text(
        0.5,
        -0.34,
        "ridge",
        ha="center",
        va="bottom",
        fontsize=8.3,
        color=RIDGE_COLOUR,
    )

    ax.set_xlim(0, matrix.shape[1])
    ax.set_ylim(matrix.shape[0], 0)
    ax.set_xticks(
        np.arange(len(ALPHA_VALUES)) + 0.5,
        [f"{alpha:.1f}" for alpha in ALPHA_VALUES],
    )
    ax.set_yticks(
        np.arange(len(model_labels)) + 0.5,
        model_labels,
    )
    ax.set_xlabel(r"Selected Elastic Net $\alpha$")
    ax.set_ylabel("")
    ax.set_title(
        "Elastic Net mixing-parameter selection across outer fits",
        pad=12,
        fontweight="semibold",
    )

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.tick_params(axis="both", length=0)

    colour_bar = fig.colorbar(
        mesh,
        ax=ax,
        fraction=0.035,
        pad=0.045,
        shrink=0.86,
    )
    colour_bar.set_label("Outer fits (of 100)")
    colour_bar.outline.set_visible(False)
    colour_bar.set_ticks(np.arange(0, 71, 10))

    fig.tight_layout()

    stem = output_dir / "Supplementary_Figure_S6_alpha_selection"
    outputs = [stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf")]

    # PNG supports direct dissertation insertion; SVG and PDF retain vector
    # geometry and editable text for archival and publication use.
    fig.savefig(outputs[0], dpi=600, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.08)
    fig.savefig(outputs[2], bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)
    return outputs


def write_provenance(
    input_paths: list[Path],
    outputs: list[Path],
    provenance_dir: Path,
) -> Path:
    """Record checksums for inputs, script and rendered outputs."""

    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "supplementary_figure_s6_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n"
            )
    return record


def main() -> None:
    """Render the supplementary Elastic Net alpha-selection figure."""
    args = parse_args()
    input_path = args.input.resolve()
    provenance_path = args.provenance.resolve()

    # Both fixed aggregate inputs are validated before the first plotting
    # object is created.
    cells = read_cells(input_path)
    validate_provenance(provenance_path)
    matrix, model_labels = build_matrix(cells)

    outputs = render(matrix, model_labels, args.output_dir.resolve())
    provenance_record = write_provenance(
        [input_path, provenance_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S6_GENERATION=PASS")
    print(f"model_alpha_cells={len(cells)}")
    print(f"maximum_selection_count={int(matrix.max())}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
