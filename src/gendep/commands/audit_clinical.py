#!/usr/bin/env python3
"""Audit the controlled GENDEP clinical workbook without exporting participants.

Purpose
-------
Verify the fixed 430 x 139 clinical source, identifiers, remission and treatment
counts, missingness, storage classes and exact duplicate groups.

Inputs
------
The authorised data-rem.xlsx workbook.

Outputs
-------
Release-safe aggregate audit tables only. No participant-level clinical values
or identifiers are written.

Usage
-----
gendep audit-clinical \
  --clinical /authorised/path/data-rem.xlsx \
  --output-dir work/analysis/clinical_audit
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.clinical.audit import audit_clinical, sha256_file  # noqa: E402
from gendep.clinical.io import read_clinical_workbook, write_records  # noqa: E402


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the clinical-audit command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--expected-rows", type=int, default=430)
    return parser.parse_args()


def main() -> None:
    """Audit the authorised clinical workbook and write aggregate evidence only."""
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    sheet, frame = read_clinical_workbook(args.clinical)
    result = audit_clinical(frame, expected_rows=args.expected_rows)
    write_records(result["inventory"], output_dir / "clinical_variable_inventory.tsv")
    write_records(result["duplicate_rows"], output_dir / "clinical_exact_duplicate_groups.tsv")
    write_records(result["cohort"], output_dir / "clinical_cohort_summary.tsv")
    write_records(result["checks"], output_dir / "clinical_audit_validation.tsv")
    summary = [
        "stage=clinical_source_audit",
        f"clinical_workbook={Path(args.clinical).resolve()}",
        f"clinical_workbook_sha256={sha256_file(args.clinical)}",
        f"worksheet={sheet}",
        f"participant_rows={len(frame)}",
        f"clinical_columns={len(frame.columns)}",
        f"validation={'PASS' if result['passed'] else 'FAIL'}",
        "participant_level_values_exported=FALSE",
    ]
    (output_dir / "clinical_audit_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print(f"CLINICAL_AUDIT={'PASS' if result['passed'] else 'FAIL'}")
    print(f"output_dir={output_dir}")
    if not result["passed"]:
        failed = [row for row in result["checks"] if row["result"] == "FAIL"]
        for row in failed:
            print(f"FAIL {row['metric']}: observed={row['observed']} expected={row['expected']}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()
