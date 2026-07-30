#!/usr/bin/env python3
"""Audit and standardise the eight fixed discovery GWAS inputs.

Purpose
-------
Convert heterogeneous source files into the column schemas consumed by
GenoUtils within GenoPred. Processing is streaming and outcome-blind.

Inputs
------
``config/project.yml`` for authorised paths and ``config/gwas.yml`` for the
fixed trait-specific mappings and sample-size rules.

Outputs
-------
One deterministic gzip-compressed table per trait, an aggregate manifest,
source headers, sample-size strategies and a validation record.

Usage
-----
gendep prepare-gwas --config config/project.yml
gendep prepare-gwas --config config/project.yml --no-strict-counts
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

from gendep.prs.gwas import (  # noqa: E402
    load_trait_specs,
    resolve_source,
    standardise_trait,
    write_manifest,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the GWAS-preparation command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="Project configuration YAML.")
    parser.add_argument(
        "--gwas-config",
        default=str(ROOT / "config" / "gwas.yml"),
        help="Trait-specific GWAS mapping YAML.",
    )
    parser.add_argument("--gwas-directory", help="Override the configured raw GWAS directory.")
    parser.add_argument("--output-dir", help="Override the configured standardised output directory.")
    parser.add_argument("--audit-dir", help="Override the configured audit directory.")
    parser.add_argument(
        "--traits",
        nargs="+",
        help="Optional subset of trait codes. The primary analysis uses all eight.",
    )
    parser.add_argument(
        "--strict-counts",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Require exact configured source and autosomal row counts (default: true).",
    )
    return parser.parse_args()


def resolve(value: str, project_root: Path) -> Path:
    """Resolve a project-relative path against the configured project root."""
    path = Path(value).expanduser()
    return path if path.is_absolute() else project_root / path


def main() -> None:
    """Standardise and validate the configured GWAS summary statistics."""
    args = parse_args()
    config_path = Path(args.config).expanduser().resolve()
    project = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    project_root = Path(project.get("project_root", ROOT)).expanduser()
    if not project_root.is_absolute():
        project_root = (config_path.parent / project_root).resolve()

    external = project.get("external_resources", {})
    outputs = project.get("outputs", {})
    gwas_directory = resolve(
        args.gwas_directory or external.get("gwas_directory", ""), project_root
    )
    output_dir = resolve(
        args.output_dir or outputs.get("gwas_standardised", "work/prs/gwas_standardised"),
        project_root,
    )
    audit_dir = resolve(
        args.audit_dir or outputs.get("gwas_audit", "work/prs/gwas_audit"),
        project_root,
    )
    if not gwas_directory.is_dir():
        raise SystemExit(f"Raw GWAS directory does not exist: {gwas_directory}")

    selected = set(args.traits or [])
    specs = load_trait_specs(args.gwas_config)
    if selected:
        unknown = selected.difference(spec.trait for spec in specs)
        if unknown:
            raise SystemExit(f"Unknown trait code(s): {', '.join(sorted(unknown))}")
        specs = [spec for spec in specs if spec.trait in selected]

    output_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)
    records = []
    header_rows = []

    for spec in specs:
        source = resolve_source(gwas_directory, spec)
        output = output_dir / f"{spec.trait}.standardised.tsv.gz"
        print(f"[{spec.trait}] {source.name}", flush=True)
        record, header = standardise_trait(
            spec,
            source,
            output,
            strict_counts=args.strict_counts,
        )
        records.append(record)
        header_rows.append(
            {
                "trait": spec.trait,
                "source": str(source.resolve()),
                "delimiter": spec.delimiter,
                "header_columns": len(header),
                "header": "|".join(header),
            }
        )

    write_manifest(audit_dir / "gwas_standardisation_manifest.tsv", records)

    with (audit_dir / "gwas_source_headers.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(header_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(header_rows)

    with (audit_dir / "sample_size_strategy.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["trait", "strategy"])
        for spec in specs:
            writer.writerow([spec.trait, spec.sample_size_strategy])

    all_passed = all(record.validation == "PASS" for record in records)
    validation_lines = [
        f"trait_count={len(records)}",
        f"traits={','.join(record.trait for record in records)}",
        f"strict_counts={str(args.strict_counts).upper()}",
        f"output_directory={output_dir.resolve()}",
        f"manifest={(audit_dir / 'gwas_standardisation_manifest.tsv').resolve()}",
        f"validation={'PASS' if all_passed else 'COUNT_DIFFERENCE'}",
    ]
    (audit_dir / "gwas_standardisation_validation.txt").write_text(
        "\n".join(validation_lines) + "\n", encoding="utf-8"
    )

    print(f"traits={len(records)}")
    print(f"output={output_dir}")
    print(f"audit={audit_dir}")
    print(f"GWAS_STANDARDISATION={'PASS' if all_passed else 'COUNT_DIFFERENCE'}")
    if args.strict_counts and not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
