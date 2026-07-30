"""Streaming audit and standardisation of the eight discovery GWAS files."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable, Mapping, TextIO

import yaml


@dataclass(frozen=True)
class TraitSpec:
    """Configuration for one source GWAS."""

    trait: str
    filename: str
    transform: str
    delimiter: str
    expected_source_rows: int
    expected_autosomal_rows: int
    sample_size_strategy: str
    fixed_n: float | None = None
    cases: float | None = None
    controls: float | None = None


@dataclass(frozen=True)
class StandardisationResult:
    """Aggregate record for one standardised GWAS."""

    trait: str
    source: str
    output: str
    source_compression: str
    delimiter: str
    header_columns: int
    source_rows: int
    autosomal_rows: int
    non_autosomal_rows: int
    invalid_sample_size_rows: int
    minimum_n: float
    mean_n: float
    maximum_n: float
    source_sha256: str
    output_sha256: str
    expected_source_rows: int
    expected_autosomal_rows: int
    validation: str


REQUIRED_COLUMNS: dict[str, tuple[str, ...]] = {
    "mdd": ("#CHROM", "POS", "ID", "EA", "NEA", "BETA", "SE", "PVAL", "NEFF", "FCAS", "FCON", "NCAS", "NCON", "IMPINFO"),
    "anx": ("CHR", "BP", "SNP", "A1", "A2", "OR", "SE", "P", "Neff_half", "FRQ_A_122083", "FRQ_U_729602", "Nca", "Nco", "INFO"),
    "bip": ("CHR", "BP", "SNP", "A1", "A2", "OR", "SE", "P", "Neff_half", "HRC_FRQ_A1", "INFO"),
    "scz": ("CHROM", "POS", "ID", "A1", "A2", "BETA", "SE", "PVAL", "NEFF", "FCAS", "FCON", "NCAS", "NCON", "IMPINFO"),
    "neur": ("CHR", "POS", "RSID", "A1", "A2", "Z", "P", "N", "EAF_UKB", "INFO_UKB"),
    "insom": ("CHR", "BP", "RSID_UKB", "A1", "A2", "OR", "SE", "P", "INFO_UKB"),
    "swb": ("CHR", "POS", "MarkerName", "A1", "A2", "Beta", "SE", "Pval", "EAF"),
    "ea": ("Chr", "BP", "rsID", "Effect_allele", "Other_allele", "Beta", "SE", "P", "EAF_HRC"),
}


def sha256(path: str | Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _is_gzip(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def open_source_text(path: str | Path) -> tuple[TextIO, str]:
    """Open an input by file signature rather than extension.

    This intentionally supports the SWB source whose filename ends
    in ``.gz`` although the downloaded bytes were plain text.
    """

    source = Path(path)
    if not source.is_file():
        raise FileNotFoundError(source)
    if _is_gzip(source):
        return gzip.open(source, "rt", encoding="utf-8", errors="replace", newline=""), "gzip"
    return source.open("rt", encoding="utf-8", errors="replace", newline=""), "plain"


def open_deterministic_gzip(path: str | Path) -> TextIO:
    """Open deterministic gzip while handling plain-text and compressed files."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    raw = output.open("wb")
    zipped = gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0)
    return io.TextIOWrapper(zipped, encoding="utf-8", newline="")


def split_line(line: str, delimiter: str) -> list[str]:
    """Split one summary-statistics row using the detected delimiter."""
    text = line.rstrip("\r\n")
    if delimiter == "tab":
        return text.split("\t")
    if delimiter == "comma":
        return next(csv.reader([text]))
    if delimiter == "whitespace":
        return text.split()
    raise ValueError(f"Unsupported delimiter: {delimiter}")


def clean_float(value: object, field: str, trait: str) -> float:
    """Parse a finite floating-point value or return missing."""
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{trait}: invalid numeric value in {field}: {value!r}") from error
    if not math.isfinite(result):
        raise ValueError(f"{trait}: non-finite value in {field}: {value!r}")
    return result


def format_number(value: float) -> str:
    """Format number for deterministic output."""
    return f"{value:.12g}"


def weighted_frequency(freq_cases: float, freq_controls: float, cases: float, controls: float) -> float:
    """Combine case and control allele frequencies by sample size."""
    denominator = cases + controls
    if denominator <= 0:
        raise ValueError("Case/control denominator must be positive")
    result = (freq_cases * cases + freq_controls * controls) / denominator
    if result < 0 or result > 1:
        raise ValueError(f"Weighted allele frequency outside [0, 1]: {result}")
    return result


def effective_sample_size(cases: float, controls: float) -> float:
    """Calculate the effective sample size for a binary GWAS."""
    if cases <= 0 or controls <= 0:
        raise ValueError("Case and control counts must be positive")
    return 4 * cases * controls / (cases + controls)


def normalise_chromosome(value: object) -> int | None:
    """Normalise chromosome to the canonical project representation."""
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    try:
        chromosome = int(text)
    except ValueError:
        return None
    return chromosome if 1 <= chromosome <= 22 else None


def _mdd(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    cases = clean_float(row["NCAS"], "NCAS", spec.trait)
    controls = clean_float(row["NCON"], "NCON", spec.trait)
    freq = weighted_frequency(
        clean_float(row["FCAS"], "FCAS", spec.trait),
        clean_float(row["FCON"], "FCON", spec.trait),
        cases,
        controls,
    )
    return {"CHR": row["#CHROM"], "BP": row["POS"], "SNP": row["ID"], "A1": row["EA"], "A2": row["NEA"], "BETA": row["BETA"], "SE": row["SE"], "P": row["PVAL"], "N": row["NEFF"], "FREQ": format_number(freq), "INFO": row["IMPINFO"]}


def _anx(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    cases = clean_float(row["Nca"], "Nca", spec.trait)
    controls = clean_float(row["Nco"], "Nco", spec.trait)
    freq = weighted_frequency(
        clean_float(row["FRQ_A_122083"], "FRQ_A_122083", spec.trait),
        clean_float(row["FRQ_U_729602"], "FRQ_U_729602", spec.trait),
        cases,
        controls,
    )
    n_value = 2 * clean_float(row["Neff_half"], "Neff_half", spec.trait)
    return {"CHR": row["CHR"], "BP": row["BP"], "SNP": row["SNP"], "A1": row["A1"], "A2": row["A2"], "OR": row["OR"], "SE": row["SE"], "P": row["P"], "N": format_number(n_value), "FREQ": format_number(freq), "INFO": row["INFO"]}


def _bip(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    n_value = 2 * clean_float(row["Neff_half"], "Neff_half", spec.trait)
    return {"CHR": row["CHR"], "BP": row["BP"], "SNP": row["SNP"], "A1": row["A1"], "A2": row["A2"], "OR": row["OR"], "SE": row["SE"], "P": row["P"], "N": format_number(n_value), "FREQ": row["HRC_FRQ_A1"], "INFO": row["INFO"]}


def _scz(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    cases = clean_float(row["NCAS"], "NCAS", spec.trait)
    controls = clean_float(row["NCON"], "NCON", spec.trait)
    freq = weighted_frequency(
        clean_float(row["FCAS"], "FCAS", spec.trait),
        clean_float(row["FCON"], "FCON", spec.trait),
        cases,
        controls,
    )
    return {"CHR": row["CHROM"], "BP": row["POS"], "SNP": row["ID"], "A1": row["A1"], "A2": row["A2"], "BETA": row["BETA"], "SE": row["SE"], "P": row["PVAL"], "N": row["NEFF"], "FREQ": format_number(freq), "INFO": row["IMPINFO"]}


def _neur(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    return {"CHR": row["CHR"], "BP": row["POS"], "SNP": row["RSID"], "A1": row["A1"], "A2": row["A2"], "Z": row["Z"], "P": row["P"], "N": row["N"], "FREQ": row["EAF_UKB"], "INFO": row["INFO_UKB"]}


def _insom(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    if spec.cases is None or spec.controls is None:
        raise ValueError("INSOM requires cases and controls in config/gwas.yml")
    n_value = effective_sample_size(spec.cases, spec.controls)
    return {"CHR": row["CHR"], "BP": row["BP"], "SNP": row["RSID_UKB"], "A1": row["A1"], "A2": row["A2"], "OR": row["OR"], "SE": row["SE"], "P": row["P"], "N": format_number(n_value), "INFO": row["INFO_UKB"]}


def _swb(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    if spec.fixed_n is None:
        raise ValueError("SWB requires fixed_n in config/gwas.yml")
    return {"CHR": row["CHR"], "BP": row["POS"], "SNP": row["MarkerName"], "A1": row["A1"], "A2": row["A2"], "BETA": row["Beta"], "SE": row["SE"], "P": row["Pval"], "N": format_number(spec.fixed_n), "FREQ": row["EAF"]}


def _ea(row: Mapping[str, str], spec: TraitSpec) -> dict[str, str]:
    if spec.fixed_n is None:
        raise ValueError("EA requires fixed_n in config/gwas.yml")
    return {"CHR": row["Chr"], "BP": row["BP"], "SNP": row["rsID"], "A1": row["Effect_allele"], "A2": row["Other_allele"], "BETA": row["Beta"], "SE": row["SE"], "P": row["P"], "N": format_number(spec.fixed_n), "FREQ": row["EAF_HRC"]}


TRANSFORMS: dict[str, Callable[[Mapping[str, str], TraitSpec], dict[str, str]]] = {
    "mdd": _mdd,
    "anx": _anx,
    "bip": _bip,
    "scz": _scz,
    "neur": _neur,
    "insom": _insom,
    "swb": _swb,
    "ea": _ea,
}


def load_trait_specs(path: str | Path) -> list[TraitSpec]:
    """Load trait specs as typed workflow records."""
    payload = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or not isinstance(payload.get("traits"), dict):
        raise ValueError("GWAS configuration must contain a traits mapping")
    specs: list[TraitSpec] = []
    for trait, values in payload["traits"].items():
        if not isinstance(values, dict):
            raise ValueError(f"Invalid configuration for {trait}")
        specs.append(
            TraitSpec(
                trait=str(trait),
                filename=str(values["filename"]),
                transform=str(values["transform"]).lower(),
                delimiter=str(values.get("delimiter", "tab")).lower(),
                expected_source_rows=int(values["expected_source_rows"]),
                expected_autosomal_rows=int(values["expected_autosomal_rows"]),
                sample_size_strategy=str(values["sample_size_strategy"]),
                fixed_n=float(values["fixed_n"]) if values.get("fixed_n") is not None else None,
                cases=float(values["cases"]) if values.get("cases") is not None else None,
                controls=float(values["controls"]) if values.get("controls") is not None else None,
            )
        )
    return specs


def resolve_source(gwas_directory: str | Path, spec: TraitSpec) -> Path:
    """Resolve source from configuration and validated paths."""
    candidate = Path(spec.filename)
    return candidate if candidate.is_absolute() else Path(gwas_directory) / candidate


def _read_header(handle: TextIO, delimiter: str) -> list[str]:
    for line in handle:
        if not line.strip() or line.startswith("##"):
            continue
        return split_line(line, delimiter)
    raise ValueError("GWAS header was not found")


def standardise_trait(
    spec: TraitSpec,
    source: str | Path,
    output: str | Path,
    *,
    strict_counts: bool = True,
) -> tuple[StandardisationResult, list[str]]:
    """Stream one source GWAS into the GenoUtils-compatible schema."""

    transform = TRANSFORMS.get(spec.transform)
    if transform is None:
        raise ValueError(f"Unknown transform for {spec.trait}: {spec.transform}")

    input_handle, compression = open_source_text(source)
    source_rows = autosomal_rows = non_autosomal_rows = invalid_n_rows = 0
    n_values_count = 0
    n_sum = 0.0
    n_min = math.inf
    n_max = -math.inf
    writer: csv.DictWriter[str] | None = None

    with input_handle, open_deterministic_gzip(output) as output_handle:
        header = _read_header(input_handle, spec.delimiter)
        missing = sorted(set(REQUIRED_COLUMNS[spec.transform]).difference(header))
        if missing:
            raise ValueError(f"{spec.trait}: missing required source columns: {', '.join(missing)}")

        for line_number, line in enumerate(input_handle, start=2):
            if not line.strip() or line.startswith("#"):
                continue
            fields = split_line(line, spec.delimiter)
            if len(fields) != len(header):
                raise ValueError(
                    f"{spec.trait}: line {line_number} has {len(fields)} fields; expected {len(header)}"
                )
            source_rows += 1
            result = transform(dict(zip(header, fields)), spec)
            chromosome = normalise_chromosome(result["CHR"])
            if chromosome is None:
                non_autosomal_rows += 1
                continue
            result["CHR"] = str(chromosome)

            n_value = clean_float(result["N"], "N", spec.trait)
            if n_value <= 0:
                invalid_n_rows += 1
                continue
            result["N"] = format_number(n_value)

            if "FREQ" in result:
                frequency = clean_float(result["FREQ"], "FREQ", spec.trait)
                if not 0 <= frequency <= 1:
                    raise ValueError(f"{spec.trait}: FREQ outside [0, 1] on line {line_number}")
                result["FREQ"] = format_number(frequency)

            if writer is None:
                writer = csv.DictWriter(
                    output_handle,
                    fieldnames=list(result),
                    delimiter="\t",
                    lineterminator="\n",
                    extrasaction="raise",
                )
                writer.writeheader()
            writer.writerow(result)
            autosomal_rows += 1
            n_values_count += 1
            n_sum += n_value
            n_min = min(n_min, n_value)
            n_max = max(n_max, n_value)

    if writer is None or n_values_count == 0:
        raise ValueError(f"{spec.trait}: no valid autosomal rows were written")

    counts_pass = (
        source_rows == spec.expected_source_rows
        and autosomal_rows == spec.expected_autosomal_rows
        and invalid_n_rows == 0
    )
    if strict_counts and not counts_pass:
        raise ValueError(
            f"{spec.trait}: row counts differ from configured targets: "
            f"source {source_rows}/{spec.expected_source_rows}, "
            f"autosomal {autosomal_rows}/{spec.expected_autosomal_rows}, "
            f"invalid N {invalid_n_rows}/0"
        )

    record = StandardisationResult(
        trait=spec.trait,
        source=str(Path(source).resolve()),
        output=str(Path(output).resolve()),
        source_compression=compression,
        delimiter=spec.delimiter,
        header_columns=len(header),
        source_rows=source_rows,
        autosomal_rows=autosomal_rows,
        non_autosomal_rows=non_autosomal_rows,
        invalid_sample_size_rows=invalid_n_rows,
        minimum_n=n_min,
        mean_n=n_sum / n_values_count,
        maximum_n=n_max,
        source_sha256=sha256(source),
        output_sha256=sha256(output),
        expected_source_rows=spec.expected_source_rows,
        expected_autosomal_rows=spec.expected_autosomal_rows,
        validation="PASS" if counts_pass else "COUNT_DIFFERENCE",
    )
    return record, header


def write_manifest(path: str | Path, records: Iterable[StandardisationResult]) -> None:
    """Write manifest to a deterministic workflow output."""
    rows = list(records)
    if not rows:
        raise ValueError("No GWAS standardisation records were supplied")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].__dict__), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            values = dict(row.__dict__)
            for field in ("minimum_n", "mean_n", "maximum_n"):
                values[field] = format_number(float(values[field]))
            writer.writerow(values)
