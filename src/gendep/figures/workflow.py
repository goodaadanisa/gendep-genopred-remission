#!/usr/bin/env python3
"""Generate Figure 1: study design and analysis workflow.

Purpose
-------
Create a deterministic schematic of the clinical, genetic and modelling
workflow used to evaluate the incremental predictive value of eight Bayesian
polygenic risk scores.

Inputs
------
- results/figure_data/figure_01_study_design_nodes.tsv

Outputs
-------
- results/figures/Figure_1_study_design_workflow.png
- results/figures/Figure_1_study_design_workflow.svg
- results/figures/Figure_1_study_design_workflow.pdf
- results/manifests/figure_01_sha256.txt

Scope
-----
The script reads aggregate design metadata only. It does not access
participant-level data, regenerate polygenic scores or fit statistical models.
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
DEFAULT_INPUT: Final = PROJECT_ROOT / "results" / "figure_data" / "figure_01_study_design_nodes.tsv"
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

# Colour-blind-safe, print-compatible palette. Meaning is also encoded through
# row position and border treatment so the figure remains interpretable in grey scale.
INK = "#1B263B"
MUTED_TEXT = "#4B5563"
LINE = "#6B7280"
INPUT_FACE = "#F3F6FA"
PROCESS_FACE = "#EDF3F8"
INTEGRATION_FACE = "#E8F4F1"
ANALYSIS_FACE = "#F4F4F3"
EVALUATION_FACE = "#EAF0F7"
NAVY = "#234B72"
TEAL = "#197A6A"
WHITE = "#FFFFFF"

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.unicode_minus": True,
        "figure.facecolor": WHITE,
        "savefig.facecolor": WHITE,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


@dataclass(frozen=True)
class Node:
    """One workflow node loaded from the fixed source table."""

    node_id: str
    stage: int
    category: str
    title: str
    detail: str
    parent_ids: tuple[str, ...]


REQUIRED_NODE_IDS: Final = {
    "clinical_input",
    "genotype_input",
    "clinical_processing",
    "genetic_processing",
    "analysis_base",
    "primary_analysis",
    "secondary_analyses",
    "evaluation",
}

ALLOWED_CATEGORIES: Final = {"input", "processing", "integration", "analysis", "evaluation"}

# Fixed positions preserve a stable publication layout while all displayed text
# and counts are supplied by the source table.
LAYOUT: Final = {
    "clinical_input": (0.17, 0.855, 0.34, 0.105),
    "genotype_input": (0.57, 0.855, 0.34, 0.105),
    "clinical_processing": (0.17, 0.675, 0.34, 0.105),
    "genetic_processing": (0.57, 0.675, 0.34, 0.105),
    "analysis_base": (0.30, 0.495, 0.48, 0.105),
    "primary_analysis": (0.17, 0.275, 0.34, 0.135),
    "secondary_analyses": (0.57, 0.275, 0.34, 0.135),
    "evaluation": (0.30, 0.065, 0.48, 0.115),
}

STAGE_LABELS: Final = {
    1: "INPUTS",
    2: "PROCESSING",
    3: "INTEGRATION",
    4: "ANALYSIS",
    5: "EVALUATION",
}

STAGE_Y: Final = {1: 0.907, 2: 0.727, 3: 0.547, 4: 0.342, 5: 0.122}

CATEGORY_STYLE: Final = {
    "input": (INPUT_FACE, NAVY, 1.35),
    "processing": (PROCESS_FACE, NAVY, 1.35),
    "integration": (INTEGRATION_FACE, TEAL, 1.7),
    "analysis": (ANALYSIS_FACE, NAVY, 1.35),
    "evaluation": (EVALUATION_FACE, NAVY, 1.7),
}


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the study-workflow figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Workflow-node TSV.")
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Directory for PNG, SVG and PDF outputs."
    )
    parser.add_argument(
        "--provenance-dir", type=Path, default=DEFAULT_PROVENANCE_DIR, help="Directory for checksum records."
    )
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_nodes(path: Path) -> dict[str, Node]:
    """Read nodes from disk and validate its basic structure."""
    if not path.is_file():
        raise FileNotFoundError(f"Workflow source table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"node_id", "stage", "category", "title", "detail", "parent_ids"}
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected workflow-table schema. "
                f"Expected {sorted(required_columns)}; observed {reader.fieldnames}."
            )

        nodes: dict[str, Node] = {}
        for row_number, row in enumerate(reader, start=2):
            node_id = row["node_id"].strip()
            if not node_id:
                raise ValueError(f"Missing node_id on row {row_number}.")
            if node_id in nodes:
                raise ValueError(f"Duplicate node_id: {node_id}")

            try:
                stage = int(row["stage"])
            except ValueError as exc:
                raise ValueError(f"Non-integer stage for {node_id}: {row['stage']}") from exc

            category = row["category"].strip()
            if category not in ALLOWED_CATEGORIES:
                raise ValueError(f"Unsupported category for {node_id}: {category}")

            title = row["title"].strip()
            detail = row["detail"].strip()
            if not title or not detail:
                raise ValueError(f"Title and detail are required for {node_id}.")

            parent_ids = tuple(part.strip() for part in row["parent_ids"].split(";") if part.strip())
            nodes[node_id] = Node(node_id, stage, category, title, detail, parent_ids)

    observed_ids = set(nodes)
    if observed_ids != REQUIRED_NODE_IDS:
        missing = sorted(REQUIRED_NODE_IDS - observed_ids)
        unexpected = sorted(observed_ids - REQUIRED_NODE_IDS)
        raise ValueError(f"Workflow-node mismatch. Missing={missing}; unexpected={unexpected}")

    for node in nodes.values():
        for parent_id in node.parent_ids:
            if parent_id not in nodes:
                raise ValueError(f"Unknown parent {parent_id!r} for node {node.node_id!r}.")
            if nodes[parent_id].stage >= node.stage:
                raise ValueError(
                    f"Parent {parent_id!r} must precede child {node.node_id!r} in the workflow."
                )

    return nodes


def add_node(ax: plt.Axes, node: Node) -> None:
    """Draw one labelled workflow node."""
    x, y, width, height = LAYOUT[node.node_id]
    face, edge, line_width = CATEGORY_STYLE[node.category]

    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.016",
        linewidth=line_width,
        edgecolor=edge,
        facecolor=face,
        zorder=3,
    )
    ax.add_patch(patch)

    # Titles and details use separate vertical anchors to preserve hierarchy and
    # prevent long explanatory text from competing with the node heading.
    ax.text(
        x + width / 2,
        y + height * 0.66,
        node.title.replace("|", "\n"),
        ha="center",
        va="center",
        fontsize=9.55,
        fontweight="semibold",
        color=INK,
        linespacing=1.05,
        zorder=4,
    )
    ax.text(
        x + width / 2,
        y + height * 0.28,
        node.detail.replace("|", "\n"),
        ha="center",
        va="center",
        fontsize=7.75,
        color=MUTED_TEXT,
        linespacing=1.12,
        zorder=4,
    )


def anchor(node_id: str, position: str) -> tuple[float, float]:
    """Return a named connection point for a workflow node."""
    x, y, width, height = LAYOUT[node_id]
    if position == "top":
        return x + width / 2, y + height
    if position == "bottom":
        return x + width / 2, y
    if position == "bottom_left":
        return x + width * 0.33, y
    if position == "bottom_right":
        return x + width * 0.67, y
    if position == "top_left":
        return x + width * 0.33, y + height
    if position == "top_right":
        return x + width * 0.67, y + height
    raise ValueError(f"Unsupported anchor position: {position}")


def add_arrow(
    ax: plt.Axes,
    start: tuple[float, float],
    end: tuple[float, float],
    *,
    connectionstyle: str = "arc3,rad=0",
) -> None:
    """Draw a directional connection between workflow stages."""
    arrow = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=11.5,
        linewidth=1.25,
        color=LINE,
        connectionstyle=connectionstyle,
        shrinkA=4,
        shrinkB=4,
        zorder=2,
    )
    ax.add_patch(arrow)


def draw_connections(ax: plt.Axes) -> None:
    """Draw the fixed workflow connections in dependency order."""
    add_arrow(ax, anchor("clinical_input", "bottom"), anchor("clinical_processing", "top"))
    add_arrow(ax, anchor("genotype_input", "bottom"), anchor("genetic_processing", "top"))

    add_arrow(
        ax,
        anchor("clinical_processing", "bottom_right"),
        anchor("analysis_base", "top_left"),
        connectionstyle="arc3,rad=-0.04",
    )
    add_arrow(
        ax,
        anchor("genetic_processing", "bottom_left"),
        anchor("analysis_base", "top_right"),
        connectionstyle="arc3,rad=0.04",
    )

    add_arrow(
        ax,
        anchor("analysis_base", "bottom_left"),
        anchor("primary_analysis", "top_right"),
        connectionstyle="arc3,rad=0.04",
    )
    add_arrow(
        ax,
        anchor("analysis_base", "bottom_right"),
        anchor("secondary_analyses", "top_left"),
        connectionstyle="arc3,rad=-0.04",
    )

    add_arrow(
        ax,
        anchor("primary_analysis", "bottom_right"),
        anchor("evaluation", "top_left"),
        connectionstyle="arc3,rad=-0.03",
    )
    add_arrow(
        ax,
        anchor("secondary_analyses", "bottom_left"),
        anchor("evaluation", "top_right"),
        connectionstyle="arc3,rad=0.03",
    )


def add_stage_labels(ax: plt.Axes) -> None:
    """Add stage labels and execution-boundary annotations."""
    for stage, label in STAGE_LABELS.items():
        ax.text(
            0.055,
            STAGE_Y[stage],
            label,
            ha="left",
            va="center",
            fontsize=7.4,
            fontweight="bold",
            color=NAVY,
            rotation=90,
            rotation_mode="anchor",
        )


def render(nodes: dict[str, Node], output_dir: Path) -> list[Path]:
    """Render the complete workflow diagram to the requested outputs."""
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(11.4, 8.0))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_stage_labels(ax)
    draw_connections(ax)
    for node_id in sorted(nodes, key=lambda value: (nodes[value].stage, value)):
        add_node(ax, nodes[node_id])

    stem = output_dir / "Figure_1_study_design_workflow"
    outputs = [stem.with_suffix(suffix) for suffix in (".png", ".svg", ".pdf")]
    fig.savefig(outputs[0], dpi=600, bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[1], bbox_inches="tight", pad_inches=0.06)
    fig.savefig(outputs[2], bbox_inches="tight", pad_inches=0.06)
    plt.close(fig)
    return outputs


def write_provenance(input_path: Path, outputs: list[Path], provenance_dir: Path) -> Path:
    """Write provenance to a deterministic workflow output."""
    provenance_dir.mkdir(parents=True, exist_ok=True)
    record = provenance_dir / "figure_01_sha256.txt"
    paths = [input_path, Path(__file__).resolve(), *outputs]
    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n")
    return record


def main() -> None:
    """Render the dissertation workflow figure."""
    args = parse_args()
    input_path = args.input.resolve()
    nodes = read_nodes(input_path)
    outputs = render(nodes, args.output_dir.resolve())
    provenance_record = write_provenance(input_path, outputs, args.provenance_dir.resolve())

    print("FIGURE_01_GENERATION=PASS")
    print(f"source={input_path}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
