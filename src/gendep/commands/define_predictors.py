#!/usr/bin/env python3
"""Create the final outcome-blind predictor and encoding specification.

Purpose
-------
Apply the fixed source-order duplicate rule, produce comprehensive and summary
predictor sets, and document semantic encodings for Elastic Net and Random
Forest modelling.

Inputs
------
The authorised clinical workbook. Remission is used only to identify and exclude
the outcome column; it is never used for predictor selection.

Outputs
-------
Release-safe predictor inventories, model lists, binary mappings and validation
records. No participant-level rows are written.

Usage
-----
gendep define-predictors \
  --clinical /authorised/path/data-rem.xlsx \
  --output-dir work/analysis/predictor_policy
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.clinical.io import read_clinical_workbook, write_records  # noqa: E402
from gendep.clinical.predictors import build_predictor_policy  # noqa: E402
from gendep.clinical.schema import MODEL_PREDICTOR_SETS  # noqa: E402


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the predictor-definition command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--output-dir", required=True)
    return parser.parse_args()


def main() -> None:
    """Generate and validate the fixed outcome-blind predictor policy."""
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    _, frame = read_clinical_workbook(args.clinical)
    result = build_predictor_policy(frame)
    write_records(result["inventory"], output_dir / "clinical_predictor_inventory.tsv")
    write_records(result["model_rows"], output_dir / "model_predictor_sets.tsv")
    write_records(result["semantic_rows"], output_dir / "predictor_semantic_type_counts.tsv")
    write_records(result["encoding_rows"], output_dir / "predictor_encoding_policy.tsv")
    write_records(result["binary_rows"], output_dir / "binary_source_code_mappings.tsv")
    write_records(result["checks"], output_dir / "predictor_policy_validation.tsv")
    machine = {name: list(values) for name, values in MODEL_PREDICTOR_SETS.items()}
    (output_dir / "model_predictor_sets.json").write_text(
        json.dumps(machine, indent=2) + "\n", encoding="utf-8"
    )
    summary = [
        "stage=outcome_blind_predictor_policy",
        "baseline_clinical_candidates=135",
        "removed_later_exact_duplicates=7",
        "comprehensive_clinical_predictors=128",
        "summary_clinical_predictors=46",
        "primary_clinical_model_predictors=134",
        "primary_combined_model_predictors=142",
        "secondary_prs_focused_predictors=15",
        "summary_clinical_model_predictors=52",
        "summary_combined_model_predictors=60",
        "outcome_used_for_selection=FALSE",
        f"validation={'PASS' if result['passed'] else 'FAIL'}",
    ]
    (output_dir / "predictor_policy_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"PREDICTOR_POLICY={'PASS' if result['passed'] else 'FAIL'}")
    print(f"output_dir={output_dir}")
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
