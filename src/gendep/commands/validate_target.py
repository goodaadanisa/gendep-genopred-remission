#!/usr/bin/env python3
"""Build PLINK2 pfiles and validate the complete target round trip.

Purpose
-------
Import each chromosome VCF into PLINK2, export it again and compare every
participant-level genotype call, sample, variant and ALT count with the
matrix-derived input VCF and manifest.

Inputs
------
Chromosome-specific target VCFs and manifests plus a compatible PLINK2 executable.

Outputs
-------
PLINK2 pfiles, exported round-trip VCFs, chromosome-level comparison tables and
a genome-wide exact mismatch summary.

Usage
-----
gendep validate-target \
  --target-dir work/genotype/target \
  --output-dir work/genotype/target_validation \
  --plink2 plink2 --expected config/expected_results.yml
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
SRC = REPOSITORY_ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.genotype.common import sha256, write_tsv  # noqa: E402
from gendep.genotype.validation import compare_vcfs, run_command  # noqa: E402


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the target-validation command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--plink2", default="plink2")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument(
        "--expected",
        default=str(REPOSITORY_ROOT / "config" / "expected_results.yml"),
    )
    parser.add_argument(
        "--reuse-existing",
        action="store_true",
        help="Reuse existing PLINK and round-trip files instead of executing PLINK2 again.",
    )
    parser.add_argument("--keep-roundtrip-vcf", action="store_true")
    return parser.parse_args()


def read_whitespace_table(path: Path, skip_double_hash: bool = False) -> list[dict[str, str]]:
    """Read whitespace table from disk and validate its basic structure."""
    rows: list[dict[str, str]] = []
    header: list[str] | None = None
    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            if skip_double_hash and raw_line.startswith("##"):
                continue
            if not raw_line.strip():
                continue
            fields = raw_line.strip().split()
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed table row in {path}")
            rows.append(dict(zip(header, fields)))
    if header is None:
        raise ValueError(f"No header found in {path}")
    return rows


def read_manifest(path: Path) -> list[dict[str, str]]:
    """Read manifest from disk and validate its basic structure."""
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sample_matches(observed: str, expected: str) -> bool:
    """Compare PSAM identifiers with the expected participant order."""
    return observed in {expected, f"{expected}_{expected}"}


def append_output_suffix(prefix: Path, suffix: str) -> Path:
    """Append an output extension while preserving names such as GENDEP.chr1."""
    if not suffix.startswith("."):
        raise ValueError("Output suffix must begin with a dot")
    return Path(f"{prefix}{suffix}")


def validate_plink_tables(
    chromosome: int,
    manifest: list[dict[str, str]],
    samples: list[str],
    pvar_path: Path,
    psam_path: Path,
    acount_path: Path,
    vmiss_path: Path,
    smiss_path: Path,
) -> list[str]:
    """Validate plink tables against the fixed project invariants."""
    reasons: list[str] = []
    pvar = read_whitespace_table(pvar_path, skip_double_hash=True)
    psam = read_whitespace_table(psam_path)
    acount = read_whitespace_table(acount_path)
    vmiss = read_whitespace_table(vmiss_path)
    smiss = read_whitespace_table(smiss_path)

    if not all(len(table) == len(manifest) for table in (pvar, acount, vmiss)):
        reasons.append("variant_table_row_count_mismatch")
    if not all(len(table) == len(samples) for table in (psam, smiss)):
        reasons.append("sample_table_row_count_mismatch")

    for index, expected in enumerate(manifest):
        if index >= len(pvar) or index >= len(acount) or index >= len(vmiss):
            break
        expected_metadata = (
            str(chromosome),
            expected["POS"],
            expected["ID"],
            expected["REF"],
            expected["ALT"],
        )
        observed_pvar = (
            pvar[index].get("#CHROM", ""),
            pvar[index].get("POS", ""),
            pvar[index].get("ID", ""),
            pvar[index].get("REF", ""),
            pvar[index].get("ALT", ""),
        )
        if observed_pvar != expected_metadata:
            reasons.append(f"pvar_metadata_mismatch:{index + 1}")

        observed_acount = acount[index]
        if (
            observed_acount.get("#CHROM", "") != str(chromosome)
            or observed_acount.get("ID", "") != expected["ID"]
            or observed_acount.get("REF", "") != expected["REF"]
            or observed_acount.get("ALT", "") != expected["ALT"]
        ):
            reasons.append(f"acount_metadata_mismatch:{index + 1}")
        try:
            if int(float(observed_acount["ALT_CTS"])) != int(float(expected["EXPECTED_ALT_COUNT"])):
                reasons.append(f"acount_alt_count_mismatch:{index + 1}")
            if int(float(observed_acount["OBS_CT"])) != int(float(expected["EXPECTED_OBS_CT"])):
                reasons.append(f"acount_observation_count_mismatch:{index + 1}")
        except (KeyError, ValueError):
            reasons.append(f"invalid_acount_record:{index + 1}")

        observed_vmiss = vmiss[index]
        if observed_vmiss.get("ID", "") != expected["ID"]:
            reasons.append(f"vmiss_order_mismatch:{index + 1}")
        try:
            if int(float(observed_vmiss["MISSING_CT"])) != 0:
                reasons.append(f"variant_missing_genotypes:{index + 1}")
            if int(float(observed_vmiss["OBS_CT"])) != len(samples):
                reasons.append(f"unexpected_variant_sample_count:{index + 1}")
        except (KeyError, ValueError):
            reasons.append(f"invalid_vmiss_record:{index + 1}")

    for index, expected_sample in enumerate(samples):
        if index >= len(psam) or index >= len(smiss):
            break
        observed_psam = psam[index]
        if observed_psam.get("#FID", "") != expected_sample or observed_psam.get("IID", "") != expected_sample:
            reasons.append(f"psam_sample_mismatch:{index + 1}")
        observed_smiss = smiss[index]
        if observed_smiss.get("#FID", "") != expected_sample or observed_smiss.get("IID", "") != expected_sample:
            reasons.append(f"smiss_sample_mismatch:{index + 1}")
        try:
            if int(float(observed_smiss["MISSING_CT"])) != 0:
                reasons.append(f"sample_missing_genotypes:{index + 1}")
            if int(float(observed_smiss["OBS_CT"])) != len(manifest):
                reasons.append(f"unexpected_sample_variant_count:{index + 1}")
        except (KeyError, ValueError):
            reasons.append(f"invalid_smiss_record:{index + 1}")

    return reasons


def validate_expected(summary: dict[str, int | str], path: Path) -> list[dict[str, object]]:
    """Validate expected against the fixed project invariants."""
    if not path.is_file():
        return []
    expected = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    variants = expected.get("variants", {})
    checks = {
        "chromosomes": 22,
        "participants": expected.get("cohort", {}).get("participants"),
        "variants": variants.get("final_oriented_target"),
        "participant_genotype_comparisons": variants.get("roundtrip_comparisons"),
        "participant_genotype_mismatches": variants.get("roundtrip_mismatches"),
    }
    rows = []
    failures = []
    for metric, expected_value in checks.items():
        if expected_value is None:
            continue
        observed = summary[metric]
        passed = int(observed) == int(expected_value)
        rows.append(
            {
                "metric": metric,
                "expected": expected_value,
                "observed": observed,
                "pass": "TRUE" if passed else "FALSE",
            }
        )
        if not passed:
            failures.append(f"{metric}: expected {expected_value}, observed {observed}")
    if failures:
        raise RuntimeError("Target round-trip regression checks failed:\n" + "\n".join(failures))
    return rows


def main() -> None:
    """Validate PLINK2 target metadata and exact genotype round trips."""
    args = parse_args()
    if args.threads < 1:
        raise ValueError("--threads must be positive")

    target_dir = Path(args.target_dir).resolve()
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    pfile_dir = output_dir / "pfiles"
    audit_dir = output_dir / "audit"
    roundtrip_dir = output_dir / "roundtrip"
    log_dir = output_dir / "logs"
    for directory in (pfile_dir, audit_dir, roundtrip_dir, log_dir):
        directory.mkdir(parents=True, exist_ok=True)

    samples_path = target_dir / "validated_sample_ids.txt"
    if not samples_path.is_file():
        raise FileNotFoundError(f"Validated sample list not found: {samples_path}")
    samples = [line.strip() for line in samples_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not samples or len(samples) != len(set(samples)):
        raise ValueError("Validated sample list is empty or duplicated")

    plink2 = shutil.which(args.plink2) or (args.plink2 if Path(args.plink2).is_file() else None)
    if not args.reuse_existing and plink2 is None:
        raise FileNotFoundError(f"PLINK2 executable not found: {args.plink2}")

    chromosome_rows = []
    mismatch_examples: list[dict[str, object]] = []

    for chromosome in range(1, 23):
        input_vcf = target_dir / "vcf" / f"GENDEP.chr{chromosome}.vcf"
        manifest_path = target_dir / "manifests" / f"GENDEP.chr{chromosome}.variant_manifest.tsv"
        if not input_vcf.is_file() or not manifest_path.is_file():
            raise FileNotFoundError(f"Chromosome {chromosome} VCF or manifest is missing")

        pfile_prefix = pfile_dir / f"GENDEP.chr{chromosome}"
        audit_prefix = audit_dir / f"GENDEP.chr{chromosome}"
        roundtrip_prefix = roundtrip_dir / f"GENDEP.chr{chromosome}"
        roundtrip_vcf = append_output_suffix(roundtrip_prefix, ".vcf")

        if not args.reuse_existing:
            run_command(
                [
                    str(plink2),
                    "--vcf",
                    str(input_vcf),
                    "--double-id",
                    "--make-pgen",
                    "--threads",
                    str(args.threads),
                    "--out",
                    str(pfile_prefix),
                ],
                log_dir / f"chr{chromosome}.import.log",
            )
            run_command(
                [
                    str(plink2),
                    "--pfile",
                    str(pfile_prefix),
                    "--freq",
                    "counts",
                    "--threads",
                    str(args.threads),
                    "--out",
                    str(audit_prefix),
                ],
                log_dir / f"chr{chromosome}.freq.log",
            )
            run_command(
                [
                    str(plink2),
                    "--pfile",
                    str(pfile_prefix),
                    "--missing",
                    "--threads",
                    str(args.threads),
                    "--out",
                    str(audit_prefix),
                ],
                log_dir / f"chr{chromosome}.missing.log",
            )
            run_command(
                [
                    str(plink2),
                    "--pfile",
                    str(pfile_prefix),
                    "--export",
                    "vcf",
                    "--threads",
                    str(args.threads),
                    "--out",
                    str(roundtrip_prefix),
                ],
                log_dir / f"chr{chromosome}.export.log",
            )

        required = [
            append_output_suffix(pfile_prefix, ".pgen"),
            append_output_suffix(pfile_prefix, ".pvar"),
            append_output_suffix(pfile_prefix, ".psam"),
            append_output_suffix(audit_prefix, ".acount"),
            append_output_suffix(audit_prefix, ".vmiss"),
            append_output_suffix(audit_prefix, ".smiss"),
            roundtrip_vcf,
        ]
        missing = [str(path) for path in required if not path.is_file()]
        if missing:
            raise FileNotFoundError("Missing PLINK2 validation outputs:\n" + "\n".join(missing))

        comparison, examples = compare_vcfs(input_vcf, roundtrip_vcf, samples)
        manifest = read_manifest(manifest_path)
        table_reasons = validate_plink_tables(
            chromosome,
            manifest,
            samples,
            append_output_suffix(pfile_prefix, ".pvar"),
            append_output_suffix(pfile_prefix, ".psam"),
            append_output_suffix(audit_prefix, ".acount"),
            append_output_suffix(audit_prefix, ".vmiss"),
            append_output_suffix(audit_prefix, ".smiss"),
        )
        mismatch_examples.extend({"chromosome": chromosome, **example} for example in examples)
        mismatch_examples.extend(
            {
                "chromosome": chromosome,
                "variant_row": "",
                "variant_id": "",
                "sample_id": "",
                "input_gt": "",
                "roundtrip_gt": "",
                "reason": reason,
            }
            for reason in table_reasons[:50]
        )

        passed = (
            comparison.variants == len(manifest)
            and comparison.samples == len(samples)
            and comparison.genotype_mismatches == 0
            and comparison.metadata_mismatches == 0
            and not table_reasons
        )
        chromosome_rows.append(
            {
                "chromosome": chromosome,
                "participants": comparison.samples,
                "variants": comparison.variants,
                "genotype_comparisons": comparison.genotype_comparisons,
                "genotype_mismatches": comparison.genotype_mismatches,
                "metadata_mismatches": comparison.metadata_mismatches,
                "plink_table_validation_errors": len(table_reasons),
                "roundtrip_pass": "TRUE" if passed else "FALSE",
            }
        )
        if not passed:
            raise RuntimeError(f"Chromosome {chromosome} failed target round-trip validation")
        if not args.keep_roundtrip_vcf:
            roundtrip_vcf.unlink()
        print(f"Validated chromosome {chromosome}: {comparison.variants:,} variants")

    summary = {
        "chromosomes": len(chromosome_rows),
        "participants": len(samples),
        "variants": sum(int(row["variants"]) for row in chromosome_rows),
        "participant_genotype_comparisons": sum(
            int(row["genotype_comparisons"]) for row in chromosome_rows
        ),
        "participant_genotype_mismatches": sum(
            int(row["genotype_mismatches"]) for row in chromosome_rows
        ),
        "metadata_mismatches": sum(int(row["metadata_mismatches"]) for row in chromosome_rows),
        "all_chromosomes_pass": "TRUE",
    }
    write_tsv(
        output_dir / "target_roundtrip_by_chromosome.tsv",
        list(chromosome_rows[0].keys()),
        chromosome_rows,
    )
    write_tsv(
        output_dir / "target_roundtrip_summary.tsv",
        ["metric", "value"],
        ({"metric": key, "value": value} for key, value in summary.items()),
    )
    write_tsv(
        output_dir / "target_roundtrip_mismatch_examples.tsv",
        [
            "chromosome",
            "variant_row",
            "variant_id",
            "sample_id",
            "input_gt",
            "roundtrip_gt",
            "reason",
        ],
        mismatch_examples,
    )

    regression_rows = validate_expected(summary, Path(args.expected)) if args.expected else []
    if regression_rows:
        write_tsv(
            output_dir / "target_roundtrip_regression_validation.tsv",
            ["metric", "expected", "observed", "pass"],
            regression_rows,
        )

    files = sorted(path for path in output_dir.rglob("*") if path.is_file() and path.suffix != ".pgen")
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "target_dir": str(target_dir),
        "plink2": str(plink2 or args.plink2),
        "threads": args.threads,
        "summary": summary,
        "checksums": {str(path.relative_to(output_dir)): sha256(path) for path in files},
    }
    (output_dir / "target_validation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print("Target round-trip validation completed successfully.")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
