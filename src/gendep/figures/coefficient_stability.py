#!/usr/bin/env python3
"""Generate Supplementary Figure S8: PRS coefficient stability.

The figure is descriptive. It summarises the eight PRS coefficients from the
100 primary clinical-plus-PRS outer fits. Ridge fits (alpha = 0) do not
perform variable selection, so Panel B distinguishes overall non-zero frequency
from frequency within the 47 sparse fits with alpha > 0.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import tempfile
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATA = ROOT / "results" / "figure_data"
DEFAULT_OUTPUT = ROOT / "results" / "figures"
DEFAULT_PROVENANCE = ROOT / "results" / "manifests"
PRS_ORDER = ["MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA"]

MPL_CACHE = Path(tempfile.gettempdir()) / "gendep_s8_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.2,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "xtick.labelsize": 8.2,
        "ytick.labelsize": 8.5,
        "legend.fontsize": 8.0,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the coefficient-stability figure command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--provenance-dir", type=Path, default=DEFAULT_PROVENANCE)
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inputs(data_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Load inputs as typed workflow records."""
    long = pd.read_csv(
        data_dir / "supplementary_figure_s8_prs_coefficient_outer_fit_long.tsv",
        sep="\t",
    )
    summary = pd.read_csv(
        data_dir / "supplementary_figure_s8_prs_coefficient_stability_summary.tsv",
        sep="\t",
    )
    tuning = pd.read_csv(
        data_dir / "supplementary_figure_s8_outer_fit_tuning_summary.tsv",
        sep="\t",
    )
    required_long = {
        "outer_repeat",
        "outer_fold",
        "fit_id",
        "term",
        "coefficient",
        "nonzero",
        "selected_alpha",
        "sparse_fit",
    }
    if not required_long.issubset(long.columns):
        raise ValueError("Unexpected long-form coefficient schema.")
    if len(long) != 800 or long["fit_id"].nunique() != 100:
        raise ValueError("Expected eight PRS across 100 outer fits.")
    if set(long["term"]) != set(PRS_ORDER):
        raise ValueError("Unexpected PRS names.")
    if len(summary) != 8 or set(summary["PRS"]) != set(PRS_ORDER):
        raise ValueError("Unexpected coefficient summary.")
    if len(tuning) != 100 or int(tuning["sparse_fit"].sum()) != 47:
        raise ValueError("Expected 100 fits, including 47 alpha > 0 fits.")
    return long, summary, tuning


def panel_label(axis: plt.Axes, label: str) -> None:
    """Add a consistent panel label to an axes object."""
    axis.text(
        -0.12,
        1.04,
        label,
        transform=axis.transAxes,
        fontsize=15,
        fontweight="bold",
        va="top",
        ha="left",
    )


def main() -> None:
    """Render the PRS coefficient-stability figure."""
    args = parse_args()
    long, summary, tuning = load_inputs(args.data_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.provenance_dir.mkdir(parents=True, exist_ok=True)

    fit_order = (
        tuning.sort_values(["outer_repeat", "outer_fold"])["fit_id"].tolist()
    )
    matrix = (
        long.pivot(index="term", columns="fit_id", values="coefficient")
        .reindex(index=PRS_ORDER, columns=fit_order)
        .to_numpy(float)
    )
    absolute_nonzero = np.abs(matrix[np.abs(matrix) > 0])
    colour_limit = float(np.quantile(absolute_nonzero, 0.98))
    colour_limit = max(colour_limit, 0.06)

    fig = plt.figure(figsize=(10.8, 7.1), constrained_layout=False)
    grid = fig.add_gridspec(
        2,
        1,
        height_ratios=[1.55, 1.0],
        left=0.09,
        right=0.96,
        bottom=0.10,
        top=0.94,
        hspace=0.42,
    )
    ax_heat = fig.add_subplot(grid[0, 0])
    ax_freq = fig.add_subplot(grid[1, 0])

    norm = TwoSlopeNorm(vmin=-colour_limit, vcenter=0.0, vmax=colour_limit)
    image = ax_heat.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap="coolwarm",
        norm=norm,
    )
    ax_heat.set_yticks(np.arange(len(PRS_ORDER)), PRS_ORDER)
    ax_heat.set_xticks(np.arange(4.5, 100, 10), [str(index) for index in range(1, 11)])
    ax_heat.set_xlabel("Outer repetition (ten folds per repetition)")
    ax_heat.set_ylabel("Polygenic score")
    ax_heat.set_title("PRS coefficients across 100 outer fits", loc="left", pad=10)
    for boundary in np.arange(9.5, 100, 10):
        ax_heat.axvline(boundary, color="black", linewidth=0.45, alpha=0.65)

    sparse_by_fit = tuning.set_index("fit_id").reindex(fit_order)["sparse_fit"].to_numpy(bool)
    for column, is_sparse in enumerate(sparse_by_fit):
        if is_sparse:
            ax_heat.plot(column, -0.72, marker="|", markersize=8, color="black", clip_on=False)
    cbar = fig.colorbar(image, ax=ax_heat, fraction=0.022, pad=0.015)
    cbar.set_label("Fitted coefficient")
    panel_label(ax_heat, "A")

    summary = summary.set_index("PRS").reindex(PRS_ORDER).reset_index()
    y_positions = np.arange(len(PRS_ORDER))
    overall = summary["overall_nonzero_percent"].to_numpy(float)
    sparse = summary["sparse_nonzero_percent"].to_numpy(float)

    for y, left, right in zip(y_positions, overall, sparse):
        ax_freq.plot([left, right], [y, y], color="0.70", linewidth=1.3, zorder=1)
    ax_freq.scatter(
        overall,
        y_positions,
        marker="s",
        s=46,
        facecolors="white",
        edgecolors="black",
        linewidths=1.0,
        label="All 100 fits",
        zorder=3,
    )
    ax_freq.scatter(
        sparse,
        y_positions,
        marker="o",
        s=48,
        facecolors="black",
        edgecolors="black",
        linewidths=0.8,
        label="47 sparse fits (alpha > 0)",
        zorder=4,
    )
    ax_freq.set_yticks(y_positions, PRS_ORDER)
    ax_freq.invert_yaxis()
    ax_freq.set_xlim(-2, 102)
    ax_freq.set_xticks(np.arange(0, 101, 20))
    ax_freq.set_xlabel("Outer fits with a non-zero PRS coefficient (%)")
    ax_freq.set_ylabel("Polygenic score")
    ax_freq.set_title("Non-zero frequency overall and when variable selection was possible", loc="left", pad=10)
    ax_freq.grid(axis="x", color="0.88", linewidth=0.7)
    ax_freq.set_axisbelow(True)
    ax_freq.spines["top"].set_visible(False)
    ax_freq.spines["right"].set_visible(False)

    for y, row in summary.iterrows():
        predominant_sign = "+" if row["positive_percent_among_nonzero"] >= 50 else "-"
        ax_freq.text(
            101.2,
            y,
            f"{predominant_sign} sign: {max(row['positive_percent_among_nonzero'], 100-row['positive_percent_among_nonzero']):.0f}%",
            ha="left",
            va="center",
            fontsize=7.4,
            clip_on=False,
        )
    ax_freq.legend(
        handles=[
            Line2D([0], [0], marker="s", markerfacecolor="white", markeredgecolor="black", linestyle="none", label="All 100 fits"),
            Line2D([0], [0], marker="o", markerfacecolor="black", markeredgecolor="black", linestyle="none", label="47 sparse fits (alpha > 0)"),
        ],
        loc="lower right",
        frameon=False,
    )
    ax_freq.text(
        0.0,
        -0.27,
        "Ridge fits (alpha = 0; 53/100) keep all penalised coefficients non-zero, so the sparse-fit frequency is the more informative selection-stability summary.",
        transform=ax_freq.transAxes,
        ha="left",
        va="top",
        fontsize=7.8,
    )
    panel_label(ax_freq, "B")

    output_paths = [
        args.output_dir / "Supplementary_Figure_S8_PRS_coefficient_stability.png",
        args.output_dir / "Supplementary_Figure_S8_PRS_coefficient_stability.svg",
        args.output_dir / "Supplementary_Figure_S8_PRS_coefficient_stability.pdf",
    ]
    fig.savefig(output_paths[0], dpi=400, bbox_inches="tight")
    fig.savefig(output_paths[1], bbox_inches="tight")
    fig.savefig(output_paths[2], bbox_inches="tight")
    plt.close(fig)

    inputs = sorted(args.data_dir.glob("*.tsv"))
    with (args.provenance_dir / "supplementary_figure_s8_sha256.txt").open("w", encoding="utf-8") as handle:
        for path in inputs + output_paths + [Path(__file__)]:
            handle.write(f"{sha256(path)}  {path.relative_to(ROOT)}\n")
    print("Supplementary Figure S8 generation: PASS")


if __name__ == "__main__":
    main()
