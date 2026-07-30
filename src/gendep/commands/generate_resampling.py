#!/usr/bin/env python3
"""Generate deterministic repeated nested-CV assignments for one or more cohorts.

Purpose
-------
Create jointly outcome-by-treatment-stratified outer and inner fold manifests
using the fixed repeated nested-CV design.

Inputs
------
One or more controlled analysis-base files containing participant identifiers,
remission and treatment-group codes.

Outputs
-------
Release-excluded outer- and inner-fold assignment tables plus balance and
integrity summaries under the requested output directory.

Examples
--------
gendep resampling \
  --cohort primary_430=data/processed/analysis_primary.tsv.gz:0:430 \
  --cohort EUR_418=data/processed/analysis_eur.tsv.gz:50000000:418 \
  --output-dir data/processed/resampling

Each ``--cohort`` value has the form ``NAME=PATH:SEED_OFFSET[:EXPECTED_N]``.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.resampling import (  # noqa: E402
    CohortSpec,
    FoldDesign,
    generate_cohort_folds,
    read_participants,
    write_tsv,
)


def parse_cohort(value: str) -> CohortSpec:
    """Parse cohort from command-line or configuration input."""
    if "=" not in value:
        raise argparse.ArgumentTypeError("Cohort must be NAME=PATH:SEED_OFFSET[:EXPECTED_N]")
    name, payload = value.split("=", 1)
    parts = payload.rsplit(":", 2)
    if len(parts) not in {2, 3}:
        raise argparse.ArgumentTypeError("Cohort must be NAME=PATH:SEED_OFFSET[:EXPECTED_N]")
    path = Path(parts[0])
    try:
        seed_offset = int(parts[1])
        expected_n = int(parts[2]) if len(parts) == 3 else None
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return CohortSpec(name=name, path=path, seed_offset=seed_offset, expected_n=expected_n)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the resampling command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cohort", action="append", type=parse_cohort, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=ROOT / "config" / "analyses.yml")
    parser.add_argument("--id-column", default="Row.names")
    parser.add_argument("--outcome-column", default="hdremit.all")
    parser.add_argument("--treatment-column", default="drug")
    return parser.parse_args()


def main() -> None:
    """Generate deterministic nested-CV assignments for the requested cohorts."""
    args = parse_args()
    with args.config.open("r", encoding="utf-8") as handle:
        config = yaml.safe_load(handle)
    design = FoldDesign(
        outer_folds=int(config["outer_folds"]),
        outer_repeats=int(config["outer_repeats"]),
        inner_folds=int(config["inner_folds"]),
        master_seed=int(config["master_seed"]),
    )

    all_balance: list[dict[str, object]] = []
    summary: list[dict[str, object]] = []
    cohort_ids: dict[str, set[str]] = {}

    for cohort in args.cohort:
        participants = read_participants(
            cohort.path,
            id_column=args.id_column,
            outcome_column=args.outcome_column,
            treatment_column=args.treatment_column,
        )
        if cohort.expected_n is not None and len(participants) != cohort.expected_n:
            raise SystemExit(f"{cohort.name}: expected {cohort.expected_n} rows; found {len(participants)}")
        cohort_ids[cohort.name] = {p.participant_id for p in participants}
        outer, inner, balance = generate_cohort_folds(
            participants,
            cohort_name=cohort.name,
            seed_offset=cohort.seed_offset,
            design=design,
        )
        write_tsv(args.output_dir / f"{cohort.name}_outer_folds.tsv.gz", outer)
        write_tsv(args.output_dir / f"{cohort.name}_inner_folds.tsv.gz", inner)
        all_balance.extend(balance)
        summary.append(
            {
                "cohort": cohort.name,
                "participants": len(participants),
                "outer_manifest_rows": len(outer),
                "inner_manifest_rows": len(inner),
                "seed_offset": cohort.seed_offset,
                "status": "PASS",
            }
        )

    # A named EUR cohort, when supplied, must be a subset of the primary cohort.
    if "primary_430" in cohort_ids and "EUR_418" in cohort_ids:
        if not cohort_ids["EUR_418"].issubset(cohort_ids["primary_430"]):
            raise SystemExit("EUR_418 is not a subset of primary_430")

    write_tsv(args.output_dir / "joint_stratified_fold_balance.tsv.gz", all_balance)
    write_tsv(args.output_dir / "resampling_summary.tsv", summary)
    with (args.output_dir / "resampling_design.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["parameter", "value"])
        writer.writerows(
            [
                ("master_seed", design.master_seed),
                ("outer_folds", design.outer_folds),
                ("outer_repeats", design.outer_repeats),
                ("inner_folds", design.inner_folds),
                ("stratification", f"{args.outcome_column}+{args.treatment_column}"),
            ]
        )

    print("NESTED_RESAMPLING_GENERATION=PASS")
    for row in summary:
        print(
            f"{row['cohort']}: participants={row['participants']} "
            f"outer_rows={row['outer_manifest_rows']} inner_rows={row['inner_manifest_rows']}"
        )


if __name__ == "__main__":
    main()
