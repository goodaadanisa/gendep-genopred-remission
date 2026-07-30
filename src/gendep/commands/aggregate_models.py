#!/usr/bin/env python3
"""Aggregate all completed outer splits for one configured modelling analysis.

Purpose
-------
Combine split-level held-out predictions, recalculate pooled metrics within each
outer repeat, summarise tuning and coefficients, and enforce split completeness.

Inputs
------
The directory written by ``scripts/run_models.R`` for one analysis.

Outputs
-------
Repeat-level metrics, mean performance, all held-out predictions, participant-
mean descriptive predictions, tuning summaries, coefficients and validation.

Usage
-----
gendep aggregate-models \
  --analysis primary \
  --split-root work/modelling \
  --output-dir work/modelling/aggregated/primary
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.modeling.aggregation import (  # noqa: E402
    aggregate_outputs,
    load_split_outputs,
    read_analysis_spec,
    read_parameters,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the model-aggregation command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", required=True)
    parser.add_argument("--split-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--analysis-config", type=Path, default=ROOT / "config" / "model_analyses.tsv"
    )
    parser.add_argument(
        "--parameter-config", type=Path, default=ROOT / "config" / "model_parameters.tsv"
    )
    parser.add_argument("--allow-incomplete", action="store_true")
    return parser.parse_args()


def write_frame(frame: pd.DataFrame, path: Path) -> None:
    """Write frame to a deterministic workflow output."""
    path.parent.mkdir(parents=True, exist_ok=True)
    compression = "gzip" if path.suffix == ".gz" else None
    frame.to_csv(path, sep="\t", index=False, compression=compression, lineterminator="\n")


def main() -> None:
    """Aggregate and validate every completed split for one configured analysis."""
    args = parse_args()
    spec = read_analysis_spec(args.analysis_config, args.analysis)
    parameters = read_parameters(args.parameter_config)
    repeats = int(parameters["outer_repeats"])
    folds = int(parameters["outer_folds"])
    epsilon = float(parameters["probability_epsilon"])

    raw = load_split_outputs(
        args.split_root.resolve(),
        spec,
        repeats,
        folds,
        allow_incomplete=args.allow_incomplete,
    )
    if args.allow_incomplete and not raw["missing_splits"].empty:
        write_frame(raw["missing_splits"], args.output_dir / "missing_splits.tsv")
        raise SystemExit(
            "Incomplete analyses cannot be aggregated into final performance estimates."
        )

    outputs = aggregate_outputs(raw, spec, repeats, folds, epsilon)
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    names = {
        "predictions": "all_held_out_predictions.tsv.gz",
        "repeat_metrics": "repeat_metrics.tsv",
        "performance_summary": "performance_summary.tsv",
        "participant_mean_predictions": "participant_mean_predictions.tsv.gz",
        "tuning": "all_tuning_results.tsv.gz",
        "tuning_selection_counts": "tuning_selection_counts.tsv",
        "coefficients": "all_coefficients.tsv.gz",
        "warnings": "all_warnings.tsv",
        "aggregation_validation": "aggregation_validation.tsv",
        "missing_splits": "missing_splits.tsv",
    }
    for key, filename in names.items():
        write_frame(outputs[key], output_dir / filename)

    validation = outputs["aggregation_validation"]
    if not validation.empty and (validation["result"] != "PASS").any():
        failed = validation.loc[validation["result"] != "PASS", "check"].tolist()
        raise SystemExit(f"Aggregation validation failed: {', '.join(failed)}")

    (output_dir / "_SUCCESS").write_text(
        "\n".join(
            [
                "status=PASS",
                f"analysis={spec.analysis}",
                f"algorithm={spec.algorithm}",
                f"participants={spec.expected_participants}",
                f"outer_repeats={repeats}",
                f"outer_folds={folds}",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    print(
        f"MODEL_AGGREGATION=PASS analysis={spec.analysis} "
        f"rows={len(outputs['predictions'])} output={output_dir}"
    )


if __name__ == "__main__":
    main()
