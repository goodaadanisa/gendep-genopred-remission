"""Conservative target-orientation sensitivity utilities.

The sensitivity retains markers whose GENDEP code-2 direction agrees with the
European-reference major-allele direction and zero-weights the tied or
near-half-frequency discordant markers. It never modifies the primary files.
"""

from __future__ import annotations

import csv
import gzip
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TextIO

import numpy as np
import pandas as pd

AGREEMENT_DECISION = "retain_sample_and_eur_major_agree"
EXPECTED_TRAITS = ("MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA")
PRS_COLUMNS = tuple(f"PRS_{trait}" for trait in EXPECTED_TRAITS)


@dataclass(frozen=True)
class OrientationPartition:
    """Represent OrientationPartition as an immutable workflow record."""
    total: int
    retained_ids: frozenset[str]
    excluded_ids: frozenset[str]
    decision_counts: dict[str, int]


def open_text(path: Path, mode: str = "rt") -> TextIO:
    """Open text while handling plain-text and compressed files."""
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def partition_orientation_markers(path: Path) -> OrientationPartition:
    """Split a PRS-CS-compatible orientation table into retained and flagged IDs."""

    retained: set[str] = set()
    excluded: set[str] = set()
    counts: dict[str, int] = {}
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"SNP", "ORIENTATION_DECISION"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Orientation table is missing: {', '.join(sorted(missing))}")
        for row_number, row in enumerate(reader, start=2):
            snp = row["SNP"].strip()
            decision = row["ORIENTATION_DECISION"].strip()
            if not snp or not decision:
                raise ValueError(f"Blank SNP or decision at {path}:{row_number}")
            if snp in retained or snp in excluded:
                raise ValueError(f"Duplicate SNP in orientation table: {snp}")
            counts[decision] = counts.get(decision, 0) + 1
            if decision == AGREEMENT_DECISION:
                retained.add(snp)
            else:
                excluded.add(snp)
    if not retained or not excluded:
        raise ValueError("Orientation partition must contain retained and excluded variants")
    return OrientationPartition(
        total=len(retained) + len(excluded),
        retained_ids=frozenset(retained),
        excluded_ids=frozenset(excluded),
        decision_counts=counts,
    )


def _split_whitespace(line: str) -> list[str]:
    return line.rstrip("\n").split()


def zero_score_file(
    source: Path,
    destination: Path,
    excluded_ids: set[str] | frozenset[str],
    *,
    weight_column: str = "SCORE_phi_auto",
    expected_rows: int = 147370,
) -> dict[str, int | float | str]:
    """Copy a GenoPred score file while setting flagged posterior weights to zero."""

    destination.parent.mkdir(parents=True, exist_ok=True)
    rows = source_nonzero = excluded_present = excluded_nonzero = 0
    absolute_removed = 0.0
    seen_excluded: set[str] = set()

    with open_text(source, "rt") as input_handle, open_text(destination, "wt") as output_handle:
        header_line = input_handle.readline()
        if not header_line:
            raise ValueError(f"Empty score file: {source}")
        header = _split_whitespace(header_line)
        required = {"SNP", "A1", "A2", weight_column}
        missing = required.difference(header)
        if missing:
            raise ValueError(f"Score file is missing: {', '.join(sorted(missing))}")
        snp_index = header.index("SNP")
        weight_index = header.index(weight_column)
        output_handle.write(" ".join(header) + "\n")

        for row_number, line in enumerate(input_handle, start=2):
            fields = _split_whitespace(line)
            if not fields:
                continue
            if len(fields) != len(header):
                raise ValueError(f"Ragged score row at {source}:{row_number}")
            rows += 1
            snp = fields[snp_index]
            weight = float(fields[weight_index])
            if not math.isfinite(weight):
                raise ValueError(f"Non-finite score weight at {source}:{row_number}")
            if weight != 0.0:
                source_nonzero += 1
            if snp in excluded_ids:
                seen_excluded.add(snp)
                excluded_present += 1
                if weight != 0.0:
                    excluded_nonzero += 1
                    absolute_removed += abs(weight)
                fields[weight_index] = "0"
            output_handle.write(" ".join(fields) + "\n")

    if rows != expected_rows:
        raise ValueError(f"Expected {expected_rows:,} score rows in {source}; found {rows:,}")
    missing_excluded = set(excluded_ids).difference(seen_excluded)
    if missing_excluded:
        preview = ", ".join(sorted(missing_excluded)[:5])
        raise ValueError(
            f"{len(missing_excluded)} flagged variants were absent from {source}; examples: {preview}"
        )
    return {
        "score_file": str(source),
        "score_rows": rows,
        "source_nonzero_weights": source_nonzero,
        "excluded_variants_present": excluded_present,
        "excluded_nonzero_weights": excluded_nonzero,
        "sum_absolute_excluded_weights": absolute_removed,
    }


def _canonical_id(values: Iterable[object]) -> pd.Series:
    series = pd.Series(list(values), dtype="string").str.strip()
    if series.isna().any() or (series == "").any():
        raise ValueError("Participant identifiers contain missing or blank values")
    return series


def build_conservative_analysis_base(
    accepted_base_path: Path,
    conservative_prs_path: Path,
    output_path: Path,
    *,
    id_column: str = "Row.names",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Replace only the eight PRS columns in the primary analysis base.

    Returns the conservative base and a release-safe aggregate comparison table.
    """

    primary = pd.read_csv(accepted_base_path, sep="\t", compression="infer", dtype={id_column: "string"})
    conservative = pd.read_csv(conservative_prs_path, sep="\t", compression="infer", dtype="string")
    if id_column not in primary.columns:
        raise ValueError(f"Primary analysis base is missing {id_column}")
    if "IID" not in conservative.columns:
        raise ValueError("Conservative PRS matrix is missing IID")
    missing_prs = set(PRS_COLUMNS).difference(conservative.columns)
    if missing_prs:
        raise ValueError(f"Conservative PRS matrix is missing: {', '.join(sorted(missing_prs))}")
    if primary[id_column].duplicated().any() or conservative["IID"].duplicated().any():
        raise ValueError("Participant identifiers must be unique")

    accepted_ids = _canonical_id(primary[id_column])
    conservative_ids = _canonical_id(conservative["IID"])
    if set(accepted_ids) != set(conservative_ids):
        raise ValueError("Primary analysis-base IDs and conservative PRS IDs differ")
    conservative = conservative.assign(_participant_id=conservative_ids).set_index("_participant_id")

    output = primary.copy()
    summaries: list[dict[str, object]] = []
    for column in PRS_COLUMNS:
        accepted_values = pd.to_numeric(primary[column], errors="raise").to_numpy(dtype=float)
        new_values = pd.to_numeric(
            conservative.loc[accepted_ids.tolist(), column], errors="raise"
        ).to_numpy(dtype=float)
        if not np.isfinite(accepted_values).all() or not np.isfinite(new_values).all():
            raise ValueError(f"Non-finite values found for {column}")
        output[column] = new_values
        difference = new_values - accepted_values
        summaries.append(
            {
                "trait": column.removeprefix("PRS_"),
                "participants": len(new_values),
                "identifiers_match": "TRUE",
                "pearson_correlation": float(np.corrcoef(accepted_values, new_values)[0, 1]),
                "spearman_correlation": float(pd.Series(accepted_values).corr(pd.Series(new_values), method="spearman")),
                "maximum_absolute_difference": float(np.max(np.abs(difference))),
                "mean_absolute_difference": float(np.mean(np.abs(difference))),
                "participants_with_changed_score": int(np.sum(difference != 0)),
                "accepted_standard_deviation": float(np.std(accepted_values, ddof=1)),
                "conservative_standard_deviation": float(np.std(new_values, ddof=1)),
            }
        )

    non_prs = [column for column in primary.columns if column not in PRS_COLUMNS]
    if not primary[non_prs].equals(output[non_prs]):
        raise RuntimeError("Non-PRS cells changed while building the conservative analysis base")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output.to_csv(output_path, sep="\t", index=False, compression="infer", lineterminator="\n")
    return output, pd.DataFrame(summaries)


def replace_outdir(source_config: Path, destination_config: Path, new_outdir: Path) -> None:
    """Copy a GenoPred YAML configuration while changing only its outdir field."""

    import re

    text = source_config.read_text(encoding="utf-8")
    updated, replacements = re.subn(
        r"(?m)^(\s*outdir\s*:\s*).*$",
        lambda match: f"{match.group(1)}{new_outdir}",
        text,
        count=1,
    )
    if replacements != 1:
        raise ValueError(
            f"Expected one outdir setting in {source_config}; found {replacements}."
        )
    destination_config.parent.mkdir(parents=True, exist_ok=True)
    destination_config.write_text(updated, encoding="utf-8")


def make_symlink(source: Path, destination: Path) -> None:
    """Create a relative-independent symlink after strict source/destination checks."""

    source = source.resolve()
    if not source.exists():
        raise FileNotFoundError(f"Required source was not found: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.is_symlink():
        if destination.resolve() == source:
            return
        raise FileExistsError(f"Destination links to a different source: {destination}")
    if destination.exists():
        raise FileExistsError(f"Destination already exists: {destination}")
    destination.symlink_to(source, target_is_directory=source.is_dir())


def read_genopred_profile(path: Path) -> tuple[pd.DataFrame, str]:
    """Read a one-score GenoPred profile and return its score-column name."""

    frame = pd.read_csv(path, sep=r"\s+", dtype={"FID": "string", "IID": "string"})
    frame = frame.rename(columns={column: column.removeprefix("#") for column in frame.columns})
    required = {"FID", "IID"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Profile is missing {', '.join(sorted(missing))}: {path}")
    score_columns = [column for column in frame.columns if column not in required]
    if len(score_columns) != 1:
        raise ValueError(f"Expected one score column in {path}; observed {len(score_columns)}")
    score_column = score_columns[0]
    frame[score_column] = pd.to_numeric(frame[score_column], errors="raise")
    if len(frame) != 430 or frame[["FID", "IID"]].duplicated().any():
        raise ValueError(f"Profile does not contain 430 unique participants: {path}")
    if not np.isfinite(frame[score_column].to_numpy(dtype=float)).all():
        raise ValueError(f"Profile contains non-finite scores: {path}")
    return frame, score_column


def prepare_genopred_scoring_bundle(
    *,
    accepted_output: Path,
    accepted_config: Path,
    output_dir: Path,
    conservative_score_root: Path,
    target_name: str = "GENDEP",
) -> dict[str, Path]:
    """Prepare exact baseline and conservative GenoPred target-scoring directories.

    The baseline mode links the primary score/model files. The conservative mode
    links the same model files but uses copied score files in which only the
    flagged posterior weights were set to zero. Both modes link the primary
    genotype, projected-PC and completed ancestry resources.
    """

    accepted_output = accepted_output.resolve()
    accepted_config = accepted_config.resolve()
    output_dir = output_dir.resolve()
    conservative_score_root = conservative_score_root.resolve()
    ancestry_marker = (
        accepted_output / "reference" / "target_checks" / target_name / "ancestry_reporter.done"
    )
    resources = {
        "genotypes": accepted_output / target_name / "geno",
        "projected_pcs": accepted_output / target_name / "pcs",
        "ancestry_marker": ancestry_marker,
    }
    for label, path in resources.items():
        if not path.exists():
            raise FileNotFoundError(f"Primary GenoPred {label} resource was not found: {path}")

    for mode in ("baseline", "conservative"):
        mode_root = output_dir / mode
        replace_outdir(accepted_config, mode_root / "config.yaml", mode_root)
        make_symlink(resources["genotypes"], mode_root / target_name / "geno")
        make_symlink(resources["projected_pcs"], mode_root / target_name / "pcs")
        make_symlink(
            ancestry_marker,
            mode_root / "reference" / "target_checks" / target_name / "ancestry_reporter.done",
        )

    for trait in EXPECTED_TRAITS:
        accepted_trait = accepted_output / "reference" / "pgs_score_files" / "prscs" / trait
        source_score = accepted_trait / f"ref-{trait}.score.gz"
        source_model = accepted_trait / f"ref-{trait}-TRANS.model.rds"
        conservative_score = conservative_score_root / trait / f"ref-{trait}.score.gz"
        for required in (source_score, source_model, conservative_score):
            if not required.exists():
                raise FileNotFoundError(f"Required scoring resource was not found: {required}")

        baseline_trait = output_dir / "baseline" / "reference" / "pgs_score_files" / "prscs" / trait
        conservative_trait = (
            output_dir / "conservative" / "reference" / "pgs_score_files" / "prscs" / trait
        )
        make_symlink(source_score, baseline_trait / source_score.name)
        make_symlink(source_model, baseline_trait / source_model.name)
        make_symlink(source_model, conservative_trait / source_model.name)
        make_symlink(conservative_score, conservative_trait / conservative_score.name)

    return {
        "baseline_config": output_dir / "baseline" / "config.yaml",
        "conservative_config": output_dir / "conservative" / "config.yaml",
        "ancestry_marker": ancestry_marker,
    }


def collect_orientation_profiles(
    *,
    accepted_output: Path,
    scoring_bundle: Path,
    output_matrix: Path,
    target_name: str = "GENDEP",
    population: str = "TRANS",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Validate baseline reproduction and collect conservative participant PRS."""

    accepted_output = accepted_output.resolve()
    scoring_bundle = scoring_bundle.resolve()
    baseline_rows: list[dict[str, object]] = []
    comparison_rows: list[dict[str, object]] = []
    conservative_matrix: pd.DataFrame | None = None

    for trait in EXPECTED_TRAITS:
        relative = Path(target_name) / "pgs" / population / "prscs" / trait / f"{target_name}-{trait}-{population}.profiles"
        primary, accepted_column = read_genopred_profile(accepted_output / relative)
        baseline, baseline_column = read_genopred_profile(scoring_bundle / "baseline" / relative)
        conservative, conservative_column = read_genopred_profile(scoring_bundle / "conservative" / relative)

        accepted_ids = primary[["FID", "IID"]].astype("string").reset_index(drop=True)
        baseline_ids = baseline[["FID", "IID"]].astype("string").reset_index(drop=True)
        conservative_ids = conservative[["FID", "IID"]].astype("string").reset_index(drop=True)
        baseline_ids_match = accepted_ids.equals(baseline_ids)
        conservative_ids_match = accepted_ids.equals(conservative_ids)
        if not baseline_ids_match or not conservative_ids_match:
            raise ValueError(f"Participant order differs for {trait}")

        accepted_values = primary[accepted_column].to_numpy(dtype=float)
        baseline_values = baseline[baseline_column].to_numpy(dtype=float)
        conservative_values = conservative[conservative_column].to_numpy(dtype=float)
        baseline_difference = np.abs(accepted_values - baseline_values)
        conservative_difference = conservative_values - accepted_values
        baseline_exact = bool(np.array_equal(accepted_values, baseline_values))
        baseline_rows.append(
            {
                "trait": trait,
                "participants": len(accepted_values),
                "identifiers_match": baseline_ids_match,
                "exact_three_decimal_reproduction": baseline_exact,
                "maximum_absolute_difference": float(np.max(baseline_difference)),
                "differing_participants": int(np.sum(baseline_difference != 0)),
            }
        )
        comparison_rows.append(
            {
                "trait": trait,
                "participants": len(accepted_values),
                "identifiers_match": conservative_ids_match,
                "pearson_correlation": float(np.corrcoef(accepted_values, conservative_values)[0, 1]),
                "spearman_correlation": float(
                    pd.Series(accepted_values).corr(pd.Series(conservative_values), method="spearman")
                ),
                "maximum_absolute_difference": float(np.max(np.abs(conservative_difference))),
                "mean_absolute_difference": float(np.mean(np.abs(conservative_difference))),
                "participants_with_changed_score": int(np.sum(conservative_difference != 0)),
                "accepted_standard_deviation": float(np.std(accepted_values, ddof=1)),
                "conservative_standard_deviation": float(np.std(conservative_values, ddof=1)),
            }
        )
        if conservative_matrix is None:
            conservative_matrix = accepted_ids.copy()
        conservative_matrix[f"PRS_{trait}"] = conservative_values

    baseline_table = pd.DataFrame(baseline_rows)
    comparison_table = pd.DataFrame(comparison_rows)
    if not (
        baseline_table["identifiers_match"].all()
        and baseline_table["exact_three_decimal_reproduction"].all()
        and (baseline_table["maximum_absolute_difference"] == 0).all()
        and (baseline_table["differing_participants"] == 0).all()
    ):
        raise RuntimeError("Baseline target scoring did not exactly reproduce primary profiles")
    assert conservative_matrix is not None
    output_matrix.parent.mkdir(parents=True, exist_ok=True)
    conservative_matrix.to_csv(
        output_matrix, sep="\t", index=False, compression="infer", lineterminator="\n"
    )
    return conservative_matrix, baseline_table, comparison_table
