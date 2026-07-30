#!/usr/bin/env python3
"""Build controlled primary and strict-EUR modelling datasets by participant ID.

Purpose
-------
Use the clinical workbook as the master table, align the eight PRS and six
projected ancestry PCs by IID, and add strict model-based EUR membership.

Inputs
------
Authorised clinical workbook, final eight-PRS matrix, projected PCs, ancestry
probabilities and the model-based EUR keep file.

Outputs
-------
Controlled participant-level primary and strict-EUR TSV files beneath work/,
plus a release-safe aggregate merge summary. The participant-level files are
excluded from Git.

Usage
-----
gendep build-analysis \
  --clinical /authorised/path/data-rem.xlsx \
  --prs work/prs/final/gendep_prscs_auto_8trait_prs.tsv.gz \
  --pcs /path/to/projected_pcs.eigenvec \
  --ancestry /path/to/GENDEP.Ancestry.model_pred \
  --eur-keep /path/to/EUR.keep \
  --output-dir work/analysis/datasets \
  --audit-dir work/analysis/integration_audit
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.clinical.integration import integrate_analysis_data, read_keep_file  # noqa: E402
from gendep.clinical.io import read_clinical_workbook, read_table, write_records, write_tsv  # noqa: E402


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the analysis-dataset build command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True)
    parser.add_argument("--prs", required=True)
    parser.add_argument("--pcs", required=True)
    parser.add_argument("--ancestry", required=True)
    parser.add_argument("--eur-keep", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--audit-dir", required=True)
    return parser.parse_args()


def main() -> None:
    """Build the controlled primary and strict-EUR analysis matrices."""
    args = parse_args()
    output_dir = Path(args.output_dir).resolve()
    audit_dir = Path(args.audit_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    _, clinical = read_clinical_workbook(args.clinical)
    prs = read_table(args.prs)
    pcs = read_table(args.pcs)
    ancestry = read_table(args.ancestry)
    eur_keep = read_keep_file(args.eur_keep)
    result = integrate_analysis_data(clinical, prs, pcs, eur_keep, ancestry)
    primary_path = output_dir / "gendep_clinical_prs_pc_analysis_base.tsv.gz"
    eur_path = output_dir / "gendep_clinical_prs_pc_analysis_base_eur418.tsv.gz"
    write_tsv(result["primary"], primary_path)
    write_tsv(result["strict_eur"], eur_path)
    write_records(
        [{"metric": key, "value": value} for key, value in result["summary"].items()],
        audit_dir / "analysis_dataset_merge_summary.tsv",
    )
    summary = [
        "stage=clinical_prs_pc_integration",
        "merge_key=clinical.Row.names_to_genetic.IID",
        "clinical_master_order_preserved=TRUE",
        f"primary_analysis_base={primary_path}",
        f"strict_eur_analysis_base={eur_path}",
        *(f"{key}={value}" for key, value in result["summary"].items()),
        "participant_level_outputs_release_safe=FALSE",
        "validation=PASS",
    ]
    (audit_dir / "analysis_dataset_merge_summary.txt").write_text("\n".join(summary) + "\n", encoding="utf-8")
    print("ANALYSIS_DATASET_BUILD=PASS")
    print(f"primary={primary_path}")
    print(f"strict_eur={eur_path}")


if __name__ == "__main__":
    main()
