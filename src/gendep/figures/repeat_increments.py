#!/usr/bin/env python3
"""Generate Supplementary Figure S5: repeat-level performance increments.

Purpose
-------
Display repetition-level AUC, Brier-score and log-loss increments across the
five predefined robustness analyses, with the ten outer-repeat values and their
scenario means.

Inputs
------
- results/figure_data/supplementary_figure_s5_repeat_level_increments.tsv
- results/figure_data/supplementary_figure_s5_provenance.tsv

Outputs
-------
- results/figures/Supplementary_Figure_S5_repeat_level_increments.png
- results/figures/Supplementary_Figure_S5_repeat_level_increments.svg
- results/figures/Supplementary_Figure_S5_repeat_level_increments.pdf
- results/manifests/supplementary_figure_s5_sha256.txt

Scope
-----
The script reads release-safe aggregate repetition-level metrics only. It does not
access participant-level predictions, fit models, perform cross-validation or
repeat uncertainty analyses.
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
    / "supplementary_figure_s5_repeat_level_increments.tsv"
)
DEFAULT_PROVENANCE: Final = (
    PROJECT_ROOT
    / "results"
    / "figure_data"
    / "supplementary_figure_s5_provenance.tsv"
)
DEFAULT_OUTPUT_DIR: Final = PROJECT_ROOT / "results" / "figures"
DEFAULT_PROVENANCE_DIR: Final = PROJECT_ROOT / "results" / "manifests"

MPL_CACHE: Final = Path(tempfile.gettempdir()) / "gendep_figure_matplotlib_cache"
MPL_CACHE.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE))

import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

# Scenario position, point jitter and the mean diamond jointly encode the
# result, so interpretation remains possible without relying on colour alone.
mpl.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.4,
        "axes.titlesize": 11.0,
        "axes.labelsize": 9.8,
        "xtick.labelsize": 8.4,
        "ytick.labelsize": 8.5,
        "figure.facecolor": "white",
        "savefig.facecolor": "white",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

SCENARIO_ORDER: Final = (
    "primary",
    "strict_eur",
    "prs_focused",
    "summary",
    "rf",
)
SCENARIO_LABELS: Final = {
    "primary": "Primary",
    "strict_eur": "Strict\nEUR",
    "prs_focused": "PRS-\nfocused",
    "summary": "Summary",
    "rf": "RF",
}
METRIC_ORDER: Final = ("AUC", "Brier", "Log loss")
METRIC_PANEL_LABELS: Final = {
    "AUC": "A",
    "Brier": "B",
    "Log loss": "C",
}
METRIC_LIMITS: Final = {
    "AUC": (-0.042, 0.010),
    "Brier": (-0.0071, 0.0050),
    "Log loss": (-0.0255, 0.0255),
}
METRIC_TICKS: Final = {
    "AUC": (-0.04, -0.03, -0.02, -0.01, 0.00),
    "Brier": (-0.006, -0.004, -0.002, 0.000, 0.002, 0.004),
    "Log loss": (-0.02, -0.01, 0.00, 0.01, 0.02),
}
JITTER_SEED: Final = 20260715


@dataclass(frozen=True)
class RepeatIncrement:
    """One repetition-level metric increment."""

    scenario_order: int
    scenario_key: str
    scenario_label: str
    outer_repeat: int
    metric_order: int
    metric: str
    value: float


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the repeat-level increment figure command."""
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


def read_increments(path: Path) -> list[RepeatIncrement]:
    """Read and validate the fixed 150-row plotting table."""

    if not path.is_file():
        raise FileNotFoundError(f"Figure source table not found: {path}")

    required_columns = {
        "scenario_order",
        "scenario_key",
        "scenario_label",
        "outer_repeat",
        "participants",
        "events",
        "non_events",
        "metric_order",
        "metric",
        "increment_value",
        "positive_value_interpretation",
    }

    increments: list[RepeatIncrement] = []
    observed_keys: set[tuple[str, int, str]] = set()

    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or set(reader.fieldnames) != required_columns:
            raise ValueError(
                "Unexpected Supplementary Figure S5 plotting-table schema."
            )

        # Every scenario-repeat-metric combination must occur exactly once.
        for row_number, row in enumerate(reader, start=2):
            scenario_key = row["scenario_key"]
            metric = row["metric"]
            outer_repeat = int(row["outer_repeat"])
            combination = (scenario_key, outer_repeat, metric)

            if scenario_key not in SCENARIO_ORDER:
                raise ValueError(
                    f"Unexpected scenario on row {row_number}: {scenario_key}"
                )
            if metric not in METRIC_ORDER:
                raise ValueError(
                    f"Unexpected metric on row {row_number}: {metric}"
                )
            if combination in observed_keys:
                raise ValueError(
                    f"Duplicate scenario-repeat-metric combination: {combination}"
                )
            observed_keys.add(combination)

            if not 1 <= outer_repeat <= 10:
                raise ValueError(
                    f"Invalid outer repetition on row {row_number}."
                )
            if row["positive_value_interpretation"] != "Favours addition of PRS":
                raise ValueError(
                    f"Unexpected metric direction on row {row_number}."
                )

            value = float(row["increment_value"])
            if not np.isfinite(value):
                raise ValueError(
                    f"Non-finite increment on row {row_number}."
                )

            increments.append(
                RepeatIncrement(
                    scenario_order=int(row["scenario_order"]),
                    scenario_key=scenario_key,
                    scenario_label=row["scenario_label"],
                    outer_repeat=outer_repeat,
                    metric_order=int(row["metric_order"]),
                    metric=metric,
                    value=value,
                )
            )

    expected_keys = {
        (scenario_key, outer_repeat, metric)
        for scenario_key in SCENARIO_ORDER
        for outer_repeat in range(1, 11)
        for metric in METRIC_ORDER
    }
    if observed_keys != expected_keys or len(increments) != 150:
        raise ValueError("The Supplementary Figure S5 result set is incomplete.")
    return increments


def validate_provenance(path: Path) -> None:
    """Validate the five-row source-run manifest."""

    if not path.is_file():
        raise FileNotFoundError(f"Source-run table not found: {path}")

    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required_columns = {
        "scenario_order",
        "scenario_key",
        "scenario_label",
        "source_file",
        "source_sha256",
        "participants",
        "events",
        "non_events",
        "outer_repetitions",
        "comparison",
    }
    if not rows or set(rows[0]) != required_columns:
        raise ValueError(
            "Unexpected Supplementary Figure S5 source-run schema."
        )
    if [row["scenario_key"] for row in rows] != list(SCENARIO_ORDER):
        raise ValueError("Source-run scenarios are missing or in the wrong order.")
    if any(int(row["outer_repetitions"]) != 10 for row in rows):
        raise ValueError("Every source run must contain ten outer repetitions.")


def default_scenario_styles() -> dict[str, str]:
    """Return the first five colours from Matplotlib's default style cycle."""

    colours = mpl.rcParams["axes.prop_cycle"].by_key()["color"]
    if len(colours) < len(SCENARIO_ORDER):
        raise ValueError("The active Matplotlib style has fewer than five colours.")
    return {
        scenario_key: colours[index]
        for index, scenario_key in enumerate(SCENARIO_ORDER)
    }


def clean_axes(ax: plt.Axes) -> None:
    """Apply the common publication treatment to one panel."""

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", alpha=0.20, linewidth=0.7)
    ax.set_axisbelow(True)


def add_panel_label(ax: plt.Axes, label: str) -> None:
    """Place a stable panel label outside the plotting area."""

    ax.text(
        -0.11,
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
    increments: list[RepeatIncrement],
    scenario_colours: dict[str, str],
    jitter_by_repeat: dict[int, float],
) -> None:
    """Draw one metric-specific repeat-level scatter panel."""

    metric_rows = [row for row in increments if row.metric == metric]

    for scenario_index, scenario_key in enumerate(SCENARIO_ORDER):
        values_by_repeat = {
            row.outer_repeat: row.value
            for row in metric_rows
            if row.scenario_key == scenario_key
        }
        if set(values_by_repeat) != set(range(1, 11)):
            raise ValueError(
                f"Missing repeat values for {scenario_key}, {metric}."
            )

        values = np.array(
            [values_by_repeat[repeat] for repeat in range(1, 11)],
            dtype=float,
        )
        x_values = np.array(
            [
                scenario_index + jitter_by_repeat[repeat]
                for repeat in range(1, 11)
            ],
            dtype=float,
        )

        # Identical jitter is reused across the three panels, allowing the same
        # outer repetition to occupy the same relative horizontal position.
        ax.scatter(
            x_values,
            values,
            s=24,
            color=scenario_colours[scenario_key],
            alpha=0.85,
            linewidths=0.0,
            zorder=3,
        )

        # The black diamond is the arithmetic mean across the ten outer
        # repetitions, not an independently fitted or bootstrapped estimate.
        ax.plot(
            scenario_index,
            float(values.mean()),
            marker="D",
            markersize=7.0,
            color="#222222",
            linestyle="none",
            zorder=5,
        )

    ax.axhline(0.0, linestyle="--", linewidth=0.9, color="#222222")
    ax.set_xlim(-0.40, len(SCENARIO_ORDER) - 0.60)
    ax.set_ylim(*METRIC_LIMITS[metric])
    ax.set_yticks(METRIC_TICKS[metric])
    ax.set_xticks(
        np.arange(len(SCENARIO_ORDER)),
        [SCENARIO_LABELS[key] for key in SCENARIO_ORDER],
    )
    ax.set_title(metric, pad=8)
    clean_axes(ax)
    add_panel_label(ax, METRIC_PANEL_LABELS[metric])


def render(
    increments: list[RepeatIncrement],
    output_dir: Path,
) -> list[Path]:
    """Render the complete three-panel repeat-level performance figure."""

    output_dir.mkdir(parents=True, exist_ok=True)

    # Manual axes preserve three aligned panels while allowing each metric to
    # retain its own scientifically appropriate y scale.
    fig = plt.figure(figsize=(10.9, 5.8))
    axes = {
        "AUC": fig.add_axes((0.070, 0.17, 0.285, 0.73)),
        "Brier": fig.add_axes((0.385, 0.17, 0.285, 0.73)),
        "Log loss": fig.add_axes((0.700, 0.17, 0.285, 0.73)),
    }

    scenario_colours = default_scenario_styles()
    rng = np.random.default_rng(JITTER_SEED)
    jitter_by_repeat = {
        repeat: float(jitter)
        for repeat, jitter in zip(
            range(1, 11),
            rng.uniform(-0.14, 0.14, size=10),
            strict=True,
        )
    }

    for metric in METRIC_ORDER:
        add_metric_panel(
            axes[metric],
            metric,
            increments,
            scenario_colours,
            jitter_by_repeat,
        )

    axes["AUC"].set_ylabel("Increment")

    fig.text(
        0.53,
        0.055,
        "Positive values favour addition of PRS; black diamonds show "
        "the mean across ten outer repetitions.",
        ha="center",
        va="bottom",
        fontsize=7.8,
    )

    stem = output_dir / "Supplementary_Figure_S5_repeat_level_increments"
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
    record = provenance_dir / "supplementary_figure_s5_sha256.txt"
    paths = [*input_paths, Path(__file__).resolve(), *outputs]

    with record.open("w", encoding="utf-8") as handle:
        for path in paths:
            handle.write(
                f"{sha256(path)}  {path.relative_to(PROJECT_ROOT)}\n"
            )
    return record


def main() -> None:
    """Render repeat-level performance increments."""
    args = parse_args()
    input_path = args.input.resolve()
    provenance_path = args.provenance.resolve()

    # Both fixed aggregate inputs are validated before the first plotting
    # object is created.
    increments = read_increments(input_path)
    validate_provenance(provenance_path)

    outputs = render(increments, args.output_dir.resolve())
    provenance_record = write_provenance(
        [input_path, provenance_path],
        outputs,
        args.provenance_dir.resolve(),
    )

    print("SUPPLEMENTARY_FIGURE_S5_GENERATION=PASS")
    print(f"plotting_rows={len(increments)}")
    for output in outputs:
        print(f"output={output} sha256={sha256(output)}")
    print(f"provenance={provenance_record}")


if __name__ == "__main__":
    main()
