"""Reference reconstruction and GenoPred target matching."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from .common import (
    is_simple_snv,
    is_snv_record,
    normalise_chromosome,
    normalise_variant_id,
    open_text,
)


@dataclass(frozen=True)
class EnsemblRecord:
    """Represent EnsemblRecord as an immutable workflow record."""
    chromosome: str
    position: int
    variant_id: str
    ref: str
    alt: str
    variant_class: str
    is_multiallelic: bool
    is_ambiguous: bool


@dataclass(frozen=True)
class GenoPredRecord:
    """Represent GenoPredRecord as an immutable workflow record."""
    chromosome: str
    position: int
    variant_id: str
    ref: str
    alt: str
    eur_alt_freq: float
    eur_obs_ct: int


def _ambiguous_pair(ref: str, alt_field: str) -> bool:
    pairs = {(ref.upper(), alt.upper()) for alt in alt_field.split(",")}
    return any(pair in {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")} for pair in pairs)


def scan_ensembl_vcfs(
    requested_ids: Iterable[str],
    vcf_pattern: str,
    chromosomes: Iterable[int] = range(1, 23),
) -> tuple[list[EnsemblRecord], dict[str, str], Counter[str]]:
    """Scan autosomal Ensembl VCFs and classify every requested rsID."""
    requested = {normalise_variant_id(value): value for value in requested_ids}
    snv_ids: set[str] = set()
    non_snv_ids: set[str] = set()
    records: list[EnsemblRecord] = []
    counts: Counter[str] = Counter()

    for chromosome in chromosomes:
        path = Path(vcf_pattern.format(chromosome=chromosome, chr=chromosome))
        if not path.is_file():
            raise FileNotFoundError(f"Missing Ensembl VCF: {path}")

        with open_text(path, "rt") as handle:
            for line in handle:
                if not line or line.startswith("#"):
                    continue
                fields = line.rstrip("\r\n").split("\t")
                if len(fields) < 8:
                    raise ValueError(f"Malformed Ensembl VCF record in {path}")
                chrom, pos_text, variant_id, ref, alt = fields[:5]
                key = normalise_variant_id(variant_id)
                if key not in requested:
                    continue

                counts["matched_reference_rows"] += 1
                if is_snv_record(ref, alt):
                    snv_ids.add(key)
                    records.append(
                        EnsemblRecord(
                            chromosome=normalise_chromosome(chrom),
                            position=int(pos_text),
                            variant_id=requested[key],
                            ref=ref.upper(),
                            alt=alt.upper(),
                            variant_class="SNV",
                            is_multiallelic="," in alt,
                            is_ambiguous=_ambiguous_pair(ref, alt),
                        )
                    )
                else:
                    non_snv_ids.add(key)

    classification: dict[str, str] = {}
    for key, original in requested.items():
        if key in snv_ids:
            classification[original] = "reconstructed_autosomal_snv"
        elif key in non_snv_ids:
            classification[original] = "excluded_non_snv"
        else:
            classification[original] = "unresolved_rsid"

    counts["requested_unique_rsids"] = len(requested)
    counts["reconstructed_autosomal_snvs"] = sum(
        value == "reconstructed_autosomal_snv" for value in classification.values()
    )
    counts["excluded_non_snv"] = sum(value == "excluded_non_snv" for value in classification.values())
    counts["unresolved_rsids"] = sum(value == "unresolved_rsid" for value in classification.values())
    counts["multiallelic_snv_rows"] = sum(record.is_multiallelic for record in records)
    counts["ambiguous_snv_rows"] = sum(record.is_ambiguous for record in records)
    return records, classification, counts


def _read_pvar_matches(path: Path, requested: set[str]) -> dict[str, list[dict[str, object]]]:
    matches: dict[str, list[dict[str, object]]] = defaultdict(list)
    header: list[str] | None = None
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if raw_line.startswith("##") or not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split()
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed PVAR row in {path}")
            row = dict(zip(header, fields))
            variant_id = row.get("ID", "")
            key = normalise_variant_id(variant_id)
            if key not in requested:
                continue
            ref, alt = row.get("REF", "").upper(), row.get("ALT", "").upper()
            if not is_simple_snv(ref, alt):
                continue
            matches[key].append(
                {
                    "chromosome": normalise_chromosome(row.get("#CHROM", row.get("CHROM", ""))),
                    "position": int(row["POS"]),
                    "variant_id": variant_id,
                    "ref": ref,
                    "alt": alt,
                }
            )
    return matches


def _read_afreq(path: Path, selected_ids: set[str]) -> dict[str, dict[str, object]]:
    result: dict[str, dict[str, object]] = {}
    header: list[str] | None = None
    with open_text(path, "rt") as handle:
        for raw_line in handle:
            if not raw_line.strip():
                continue
            fields = raw_line.rstrip("\r\n").split()
            if header is None:
                header = fields
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed AFREQ row in {path}")
            row = dict(zip(header, fields))
            key = normalise_variant_id(row.get("ID", ""))
            if key not in selected_ids:
                continue
            if key in result:
                raise ValueError(f"Duplicate AFREQ identifier {row.get('ID')} in {path}")
            alt_freq_text = row.get("ALT_FREQS", "")
            if "," in alt_freq_text:
                raise ValueError(f"Multiple ALT frequencies for {row.get('ID')} in {path}")
            result[key] = {
                "chromosome": normalise_chromosome(row.get("#CHROM", row.get("CHROM", ""))),
                "variant_id": row.get("ID", ""),
                "ref": row.get("REF", "").upper(),
                "alt": row.get("ALT", "").upper(),
                "eur_alt_freq": float(alt_freq_text),
                "eur_obs_ct": int(float(row.get("OBS_CT", 0))),
            }
    return result


def match_genopred_reference(
    requested_ids: Iterable[str],
    pvar_pattern: str,
    afreq_pattern: str,
    chromosomes: Iterable[int] = range(1, 23),
) -> tuple[list[GenoPredRecord], dict[str, str], Counter[str]]:
    """Match original rsIDs directly to the GenoPred reference and EUR frequencies."""
    requested_map = {normalise_variant_id(value): value for value in requested_ids}
    remaining = set(requested_map)
    output: list[GenoPredRecord] = []
    exclusions: dict[str, str] = {}
    counts: Counter[str] = Counter()

    for chromosome in chromosomes:
        pvar_path = Path(pvar_pattern.format(chromosome=chromosome, chr=chromosome))
        afreq_path = Path(afreq_pattern.format(chromosome=chromosome, chr=chromosome))
        if not pvar_path.is_file():
            raise FileNotFoundError(f"Missing GenoPred PVAR: {pvar_path}")
        if not afreq_path.is_file():
            raise FileNotFoundError(f"Missing GenoPred EUR AFREQ: {afreq_path}")

        pvar_matches = _read_pvar_matches(pvar_path, remaining)
        unambiguous = {key: rows[0] for key, rows in pvar_matches.items() if len(rows) == 1}
        for key, rows in pvar_matches.items():
            if len(rows) > 1:
                exclusions[requested_map[key]] = "ambiguous_genopred_rsid"
                counts["ambiguous_genopred_rsid"] += 1

        frequencies = _read_afreq(afreq_path, set(unambiguous))
        for key, pvar in unambiguous.items():
            frequency = frequencies.get(key)
            original_id = requested_map[key]
            if frequency is None:
                exclusions[original_id] = "missing_eur_frequency"
                counts["missing_eur_frequency"] += 1
                continue
            if (
                frequency["chromosome"] != pvar["chromosome"]
                or frequency["ref"] != pvar["ref"]
                or frequency["alt"] != pvar["alt"]
            ):
                exclusions[original_id] = "pvar_afreq_metadata_mismatch"
                counts["pvar_afreq_metadata_mismatch"] += 1
                continue
            alt_freq = float(frequency["eur_alt_freq"])
            if not 0.0 <= alt_freq <= 1.0:
                exclusions[original_id] = "invalid_eur_alt_frequency"
                counts["invalid_eur_alt_frequency"] += 1
                continue
            output.append(
                GenoPredRecord(
                    chromosome=str(pvar["chromosome"]),
                    position=int(pvar["position"]),
                    variant_id=original_id,
                    ref=str(pvar["ref"]),
                    alt=str(pvar["alt"]),
                    eur_alt_freq=alt_freq,
                    eur_obs_ct=int(frequency["eur_obs_ct"]),
                )
            )
            remaining.discard(key)
            counts["matched"] += 1

    for key in remaining:
        original_id = requested_map[key]
        exclusions.setdefault(original_id, "absent_from_genopred_reference")
    counts["absent_from_genopred_reference"] = sum(
        value == "absent_from_genopred_reference" for value in exclusions.values()
    )
    output.sort(key=lambda record: (int(record.chromosome), record.position, record.variant_id))
    return output, exclusions, counts
