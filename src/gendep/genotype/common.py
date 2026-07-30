"""Shared I/O and validation helpers for genotype reconstruction."""

from __future__ import annotations

import csv
import gzip
import hashlib
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence, TextIO

VALID_BASES = frozenset({"A", "C", "G", "T"})
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def open_text(path: str | Path, mode: str = "rt") -> TextIO:
    """Open plain or gzip-compressed text using UTF-8."""
    file_path = Path(path)
    if "b" in mode:
        raise ValueError("open_text only supports text modes")
    if file_path.suffix == ".gz":
        return gzip.open(file_path, mode, encoding="utf-8", newline="")
    return file_path.open(mode, encoding="utf-8", newline="")


def normalise_chromosome(value: object) -> str:
    """Normalise chromosome to the canonical project representation."""
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return text


def normalise_variant_id(value: object) -> str:
    """Normalise variant id to the canonical project representation."""
    return str(value).strip().lower()


def complement(allele: str) -> str:
    """Return the Watson-Crick complement of one nucleotide allele."""
    return COMPLEMENT.get(allele.upper(), "")


def is_simple_snv(ref: str, alt: str) -> bool:
    """Return whether both alleles define a simple biallelic SNV."""
    return (
        len(ref) == 1
        and len(alt) == 1
        and ref.upper() in VALID_BASES
        and alt.upper() in VALID_BASES
        and ref.upper() != alt.upper()
    )


def is_snv_record(ref: str, alt_field: str) -> bool:
    """Return True for single-base REF and one or more single-base ALT alleles."""
    ref = ref.upper()
    alts = [value.upper() for value in alt_field.split(",") if value]
    return (
        len(ref) == 1
        and ref in VALID_BASES
        and bool(alts)
        and all(len(alt) == 1 and alt in VALID_BASES and alt != ref for alt in alts)
    )


def allele_relation(ref: str, alt: str, a1: str, a2: str) -> str | None:
    """Classify two biallelic SNP representations, including strand complements."""
    ref, alt, a1, a2 = (value.upper() for value in (ref, alt, a1, a2))
    if (ref, alt) == (a1, a2):
        return "same_order"
    if (ref, alt) == (a2, a1):
        return "swapped_order"
    comp_ref, comp_alt = complement(ref), complement(alt)
    if not comp_ref or not comp_alt:
        return None
    if (comp_ref, comp_alt) == (a1, a2):
        return "strand_complement"
    if (comp_ref, comp_alt) == (a2, a1):
        return "strand_complement_swapped"
    return None


def sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_variant_ids(path: str | Path) -> list[str]:
    """Read variant ids from disk and validate its basic structure."""
    values = [line.strip() for line in Path(path).read_text(encoding="utf-8").splitlines()]
    values = [value for value in values if value]
    if not values:
        raise ValueError(f"Variant list is empty: {path}")
    if len(set(values)) != len(values):
        raise ValueError(f"Variant list contains duplicate identifiers: {path}")
    return values


def require_columns(fieldnames: Sequence[str] | None, required: Iterable[str], source: str) -> None:
    """Require columns and raise a clear error when absent."""
    if fieldnames is None:
        raise ValueError(f"No header was found in {source}")
    missing = sorted(set(required).difference(fieldnames))
    if missing:
        raise ValueError(f"{source} is missing required columns: {', '.join(missing)}")


def read_tsv(path: str | Path) -> Iterator[dict[str, str]]:
    """Read tsv from disk and validate its basic structure."""
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"No header was found in {path}")
        for row in reader:
            yield dict(row)


def write_tsv(
    path: str | Path,
    fieldnames: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    """Write a deterministic TSV workflow output."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with open_text(output, "wt") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
