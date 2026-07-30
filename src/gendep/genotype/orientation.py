"""Frequency-based target orientation and PRS-CS compatibility."""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .common import allele_relation, normalise_chromosome, normalise_variant_id, open_text
from .reference import GenoPredRecord


@dataclass(frozen=True)
class GenotypeFrequency:
    """Represent GenotypeFrequency as an immutable workflow record."""
    variant_id: str
    n_code_0: int
    n_code_1: int
    n_code_2: int
    n_missing: int
    code2_frequency: float


@dataclass(frozen=True)
class OrientedMarker:
    """Represent OrientedMarker as an immutable workflow record."""
    SNP: str
    ORIGINAL_CATEGORY: str
    GENOPRED_CHR: str
    GENOPRED_BP: int
    GENOPRED_REF: str
    GENOPRED_ALT: str
    EUR_ALT_FREQ: float
    EUR_REF_FREQ: float
    EUR_MAJOR_FREQ: float
    EUR_OBS_CT: int
    N_CODE_0: int
    N_CODE_1: int
    N_CODE_2: int
    N_MISSING: int
    GENDEP_CODE2_FREQ: float
    DIFFERENCE_TO_EUR_MAJOR_FREQ: float
    CODE2_SIDE: str
    CODE0_SIDE: str
    CODE2_ALLELE: str
    CODE0_ALLELE: str
    GENDEP_CODE2_SAMPLE_STATUS: str
    ORIENTATION_DECISION: str
    KEEP_FOR_TARGET: bool


def read_genotype_frequencies(path: str | Path) -> dict[str, GenotypeFrequency]:
    """Read genotype frequencies from disk and validate its basic structure."""
    import csv

    required = {
        "SNP",
        "N_CODE_0",
        "N_CODE_1",
        "N_CODE_2",
        "N_MISSING",
        "GENDEP_CODE2_FREQ",
    }
    output: dict[str, GenotypeFrequency] = {}
    expected_sample_count: int | None = None
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Frequency table is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            variant_id = row["SNP"].strip()
            if variant_id in output:
                raise ValueError(f"Duplicate frequency record for {variant_id}")
            record = GenotypeFrequency(
                variant_id=variant_id,
                n_code_0=int(row["N_CODE_0"]),
                n_code_1=int(row["N_CODE_1"]),
                n_code_2=int(row["N_CODE_2"]),
                n_missing=int(row["N_MISSING"]),
                code2_frequency=float(row["GENDEP_CODE2_FREQ"]),
            )
            if not 0.0 <= record.code2_frequency <= 1.0:
                raise ValueError(f"Invalid code-2 frequency for {variant_id}")
            sample_count = record.n_code_0 + record.n_code_1 + record.n_code_2 + record.n_missing
            if sample_count <= 0:
                raise ValueError(f"Invalid aggregate genotype count for {variant_id}")
            if expected_sample_count is None:
                expected_sample_count = sample_count
            elif sample_count != expected_sample_count:
                raise ValueError(
                    f"Aggregate genotype counts do not use one sample count: {variant_id}"
                )
            nonmissing = sample_count - record.n_missing
            if nonmissing <= 0:
                raise ValueError(f"No nonmissing genotypes for {variant_id}")
            expected_frequency = (record.n_code_1 + 2 * record.n_code_2) / (2 * nonmissing)
            if abs(expected_frequency - record.code2_frequency) > 1e-12:
                raise ValueError(
                    f"Code-2 frequency does not agree with aggregate counts for {variant_id}"
                )
            output[variant_id] = record
    return output


def orient_markers(
    records: Iterable[GenoPredRecord],
    frequencies: Mapping[str, GenotypeFrequency],
    source_classification: Mapping[str, str],
    tie_tolerance: float = 1e-12,
) -> tuple[list[OrientedMarker], list[OrientedMarker], Counter[str]]:
    """Infer the nucleotide counted by source code 2 using the fixed European major-allele rule."""
    retained: list[OrientedMarker] = []
    excluded: list[OrientedMarker] = []
    counts: Counter[str] = Counter()

    for record in records:
        frequency = frequencies.get(record.variant_id)
        if frequency is None:
            raise ValueError(f"Missing GENDEP genotype frequency for {record.variant_id}")
        eur_ref_freq = 1.0 - record.eur_alt_freq
        eur_major_freq = max(eur_ref_freq, record.eur_alt_freq)
        eur_tie = abs(record.eur_alt_freq - 0.5) <= tie_tolerance
        if eur_tie:
            code2_side, code0_side, code2_allele, code0_allele = "", "", "", ""
            decision = "exclude_eur_frequency_tie"
            keep = False
        else:
            code2_side = "ALT" if record.eur_alt_freq > 0.5 else "REF"
            code0_side = "REF" if code2_side == "ALT" else "ALT"
            code2_allele = record.alt if code2_side == "ALT" else record.ref
            code0_allele = record.ref if code0_side == "REF" else record.alt
            keep = True
            if abs(frequency.code2_frequency - 0.5) <= tie_tolerance:
                decision = "retain_eur_major_sample_tie"
            elif frequency.code2_frequency > 0.5:
                decision = "retain_sample_and_eur_major_agree"
            else:
                decision = "retain_eur_major_sample_flip_near_half"

        sample_status = (
            "equal_frequency"
            if abs(frequency.code2_frequency - 0.5) <= tie_tolerance
            else "code2_more_frequent"
            if frequency.code2_frequency > 0.5
            else "code2_less_frequent"
        )
        marker = OrientedMarker(
            SNP=record.variant_id,
            ORIGINAL_CATEGORY=source_classification.get(record.variant_id, "unknown"),
            GENOPRED_CHR=record.chromosome,
            GENOPRED_BP=record.position,
            GENOPRED_REF=record.ref,
            GENOPRED_ALT=record.alt,
            EUR_ALT_FREQ=record.eur_alt_freq,
            EUR_REF_FREQ=eur_ref_freq,
            EUR_MAJOR_FREQ=eur_major_freq,
            EUR_OBS_CT=record.eur_obs_ct,
            N_CODE_0=frequency.n_code_0,
            N_CODE_1=frequency.n_code_1,
            N_CODE_2=frequency.n_code_2,
            N_MISSING=frequency.n_missing,
            GENDEP_CODE2_FREQ=frequency.code2_frequency,
            DIFFERENCE_TO_EUR_MAJOR_FREQ=abs(frequency.code2_frequency - eur_major_freq),
            CODE2_SIDE=code2_side,
            CODE0_SIDE=code0_side,
            CODE2_ALLELE=code2_allele,
            CODE0_ALLELE=code0_allele,
            GENDEP_CODE2_SAMPLE_STATUS=sample_status,
            ORIENTATION_DECISION=decision,
            KEEP_FOR_TARGET=keep,
        )
        counts[decision] += 1
        counts[sample_status] += 1
        (retained if keep else excluded).append(marker)

    counts["matched_genopred_variants"] = len(retained) + len(excluded)
    counts["retained_oriented_variants"] = len(retained)
    counts["excluded_eur_frequency_ties"] = len(excluded)
    return retained, excluded, counts


def marker_to_dict(marker: OrientedMarker) -> dict[str, object]:
    """Serialise an oriented marker record for tabular output."""
    output = asdict(marker)
    output["KEEP_FOR_TARGET"] = "TRUE" if marker.KEEP_FOR_TARGET else "FALSE"
    return output


def pearson_code2_vs_eur_major(markers: Iterable[OrientedMarker]) -> float:
    """Calculate the concordance between GENDEP code-2 frequency and European major-allele frequency."""
    import math

    pairs = [(marker.GENDEP_CODE2_FREQ, marker.EUR_MAJOR_FREQ) for marker in markers]
    if len(pairs) < 2:
        return float("nan")
    x = [pair[0] for pair in pairs]
    y = [pair[1] for pair in pairs]
    mean_x = sum(x) / len(x)
    mean_y = sum(y) / len(y)
    numerator = sum((a - mean_x) * (b - mean_y) for a, b in zip(x, y))
    denominator = math.sqrt(
        sum((a - mean_x) ** 2 for a in x) * sum((b - mean_y) ** 2 for b in y)
    )
    return numerator / denominator if denominator else float("nan")


def load_prscs(
    path: str | Path,
    requested_ids: set[str] | None = None,
    requested_coordinates: set[tuple[str, int]] | None = None,
):
    """Load prscs as typed workflow records."""
    by_id: dict[str, list[dict[str, object]]] = defaultdict(list)
    by_coordinate: dict[tuple[str, int], list[dict[str, object]]] = defaultdict(list)
    rows = 0
    with open_text(path, "rt") as handle:
        header_fields = handle.readline().strip().split()
        header = {name: index for index, name in enumerate(header_fields)}
        required = {"CHR", "SNP", "BP", "A1", "A2", "MAF"}
        missing = required.difference(header)
        if missing:
            raise ValueError(f"PRS-CS SNP info is missing columns: {', '.join(sorted(missing))}")
        for line in handle:
            fields = line.strip().split()
            if not fields:
                continue
            record = {
                "CHR": normalise_chromosome(fields[header["CHR"]]),
                "SNP": fields[header["SNP"]],
                "BP": int(fields[header["BP"]]),
                "A1": fields[header["A1"]].upper(),
                "A2": fields[header["A2"]].upper(),
                "MAF": fields[header["MAF"]],
            }
            rows += 1
            record_id = normalise_variant_id(record["SNP"])
            coordinate = (record["CHR"], record["BP"])
            if (
                requested_ids is not None
                and requested_coordinates is not None
                and record_id not in requested_ids
                and coordinate not in requested_coordinates
            ):
                continue
            by_id[record_id].append(record)
            by_coordinate[coordinate].append(record)
    return by_id, by_coordinate, rows


def intersect_prscs(
    markers: Iterable[OrientedMarker],
    snpinfo_path: str | Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]], Counter[str]]:
    """Match oriented target variants to PRS-CS by rsID, coordinate and allele relation."""
    marker_list = list(markers)
    requested_ids = {normalise_variant_id(marker.SNP) for marker in marker_list}
    requested_coordinates = {
        (marker.GENOPRED_CHR, marker.GENOPRED_BP) for marker in marker_list
    }
    by_id, by_coordinate, reference_rows = load_prscs(
        snpinfo_path,
        requested_ids=requested_ids,
        requested_coordinates=requested_coordinates,
    )
    compatible: list[dict[str, object]] = []
    excluded: list[dict[str, object]] = []
    counts: Counter[str] = Counter(prscs_reference_rows=reference_rows)

    def candidates(records, marker):
        """Return valid orientation candidates for the current marker."""
        output = []
        for record in records:
            if record["CHR"] != marker.GENOPRED_CHR or record["BP"] != marker.GENOPRED_BP:
                continue
            relation = allele_relation(
                marker.GENOPRED_REF,
                marker.GENOPRED_ALT,
                str(record["A1"]),
                str(record["A2"]),
            )
            if relation is not None:
                output.append((record, relation))
        return output

    for marker in marker_list:
        counts["oriented_input_variants"] += 1
        exact_records = by_id.get(normalise_variant_id(marker.SNP), [])
        matches = candidates(exact_records, marker)
        match_method = ""
        status = ""
        selected = None
        if len(matches) == 1:
            selected = matches[0]
            match_method = "exact_rsid"
            status = "compatible"
            counts["exact_rsid_compatible"] += 1
        elif len(matches) > 1:
            status = "ambiguous_prscs_match"
        else:
            coordinate_matches = candidates(
                by_coordinate.get((marker.GENOPRED_CHR, marker.GENOPRED_BP), []), marker
            )
            if len(coordinate_matches) == 1:
                selected = coordinate_matches[0]
                match_method = "coordinate_and_alleles"
                status = "compatible"
                counts["coordinate_allele_recovered"] += 1
            elif len(coordinate_matches) > 1:
                status = "ambiguous_prscs_match"
            elif exact_records:
                status = "exact_rsid_metadata_mismatch"
            else:
                status = "absent_from_prscs"

        base = marker_to_dict(marker)
        if selected is None:
            counts[status] += 1
            base["PRSCS_INTERSECTION_STATUS"] = status
            excluded.append(base)
            continue

        record, relation = selected
        base.update(
            {
                "PRSCS_SNP": record["SNP"],
                "PRSCS_CHR": record["CHR"],
                "PRSCS_BP": record["BP"],
                "PRSCS_A1": record["A1"],
                "PRSCS_A2": record["A2"],
                "PRSCS_MAF": record["MAF"],
                "PRSCS_MATCH_METHOD": match_method,
                "PRSCS_ALLELE_RELATION": relation,
                "PRSCS_INTERSECTION_STATUS": status,
            }
        )
        compatible.append(base)
        counts["final_prscs_compatible"] += 1
        counts[relation] += 1

    return compatible, excluded, counts
