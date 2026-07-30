"""VCF and PLINK2 round-trip validation helpers."""

from __future__ import annotations

import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO


@dataclass(frozen=True)
class VcfComparison:
    """Represent VcfComparison as an immutable workflow record."""
    samples: int
    variants: int
    genotype_comparisons: int
    genotype_mismatches: int
    metadata_mismatches: int


def read_vcf_header(handle: TextIO, source: str | Path) -> list[str]:
    """Read vcf header from disk and validate its basic structure."""
    for raw_line in handle:
        if raw_line.startswith("#CHROM"):
            fields = raw_line.rstrip("\r\n").split("\t")
            if len(fields) < 10:
                raise ValueError(f"VCF header has too few fields in {source}")
            return fields[9:]
        if not raw_line.startswith("#"):
            raise ValueError(f"Variant data appeared before the VCF header in {source}")
    raise ValueError(f"VCF header not found in {source}")


def next_variant(handle: TextIO) -> list[str] | None:
    """Advance a VCF iterator and return the next parsed variant."""
    for raw_line in handle:
        if not raw_line.startswith("#"):
            return raw_line.rstrip("\r\n").split("\t")
    return None


def normalise_gt(sample_field: str, format_field: str) -> str | None:
    """Normalise gt to the canonical project representation."""
    format_parts = format_field.split(":")
    if "GT" not in format_parts:
        return None
    gt_index = format_parts.index("GT")
    sample_parts = sample_field.split(":")
    if gt_index >= len(sample_parts):
        return None
    genotype = sample_parts[gt_index].replace("|", "/")
    if genotype in {".", "./."}:
        return None
    alleles = genotype.split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return "/".join(sorted(alleles))


def compare_vcfs(
    input_vcf: str | Path,
    roundtrip_vcf: str | Path,
    expected_samples: Iterable[str] | None = None,
    max_examples: int = 50,
) -> tuple[VcfComparison, list[dict[str, object]]]:
    """Compare source and round-tripped VCFs across sample order, variant order and every genotype call."""
    examples: list[dict[str, object]] = []
    genotype_mismatches = 0
    metadata_mismatches = 0
    variants = 0
    genotype_comparisons = 0

    with Path(input_vcf).open("r", encoding="utf-8", newline="") as input_handle, Path(
        roundtrip_vcf
    ).open("r", encoding="utf-8", newline="") as roundtrip_handle:
        input_samples = read_vcf_header(input_handle, input_vcf)
        roundtrip_samples = read_vcf_header(roundtrip_handle, roundtrip_vcf)
        expected = list(expected_samples) if expected_samples is not None else input_samples
        # PLINK2 can export FID_IID when FID and IID are identical.
        if len(roundtrip_samples) != len(expected):
            raise ValueError("Round-trip VCF sample count differs from the expected count")
        for observed, wanted in zip(roundtrip_samples, expected):
            if observed not in {wanted, f"{wanted}_{wanted}"}:
                raise ValueError(f"Round-trip sample order mismatch: expected {wanted}, found {observed}")
        if input_samples != expected:
            raise ValueError("Input VCF samples do not match the expected order")

        while True:
            left = next_variant(input_handle)
            right = next_variant(roundtrip_handle)
            if left is None and right is None:
                break
            if left is None or right is None:
                raise ValueError("Input and round-trip VCFs contain different numbers of variants")
            variants += 1
            if tuple(left[:5]) != tuple(right[:5]):
                metadata_mismatches += 1
                if len(examples) < max_examples:
                    examples.append(
                        {
                            "variant_row": variants,
                            "variant_id": left[2] if len(left) > 2 else "",
                            "sample_id": "",
                            "input_gt": "",
                            "roundtrip_gt": "",
                            "reason": "metadata_mismatch",
                        }
                    )
            if len(left) != 9 + len(expected) or len(right) != 9 + len(expected):
                raise ValueError(f"Unexpected VCF field count at variant row {variants}")
            for index, sample_id in enumerate(expected):
                left_gt = normalise_gt(left[9 + index], left[8])
                right_gt = normalise_gt(right[9 + index], right[8])
                genotype_comparisons += 1
                if left_gt != right_gt or left_gt is None:
                    genotype_mismatches += 1
                    if len(examples) < max_examples:
                        examples.append(
                            {
                                "variant_row": variants,
                                "variant_id": left[2],
                                "sample_id": sample_id,
                                "input_gt": left_gt or "INVALID",
                                "roundtrip_gt": right_gt or "INVALID",
                                "reason": "genotype_mismatch",
                            }
                        )

    return (
        VcfComparison(
            samples=len(expected),
            variants=variants,
            genotype_comparisons=genotype_comparisons,
            genotype_mismatches=genotype_mismatches,
            metadata_mismatches=metadata_mismatches,
        ),
        examples,
    )


def run_command(command: list[str], log_path: str | Path) -> None:
    """Run command and record its validation status."""
    log = Path(log_path)
    log.parent.mkdir(parents=True, exist_ok=True)
    with log.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, text=True)
    if completed.returncode != 0:
        raise RuntimeError(f"Command failed ({completed.returncode}): {' '.join(command)}; see {log}")
