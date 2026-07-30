"""Deterministic joint outcome-by-treatment nested cross-validation.

The implementation preserves the fixed GENDEP fold-allocation rule while
making cohort names, columns, fold counts and seed offsets explicit inputs.
"""

from __future__ import annotations

import csv
import gzip
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping, Sequence


@dataclass(frozen=True)
class Participant:
    """Represent Participant as an immutable workflow record."""
    participant_id: str
    outcome: str
    treatment: str

    @property
    def stratum(self) -> tuple[str, str]:
        """Return the participant outcome-by-treatment stratum."""
        return self.outcome, self.treatment


@dataclass(frozen=True)
class FoldDesign:
    """Represent FoldDesign as an immutable workflow record."""
    outer_folds: int = 10
    outer_repeats: int = 10
    inner_folds: int = 5
    master_seed: int = 20260715


@dataclass(frozen=True)
class CohortSpec:
    """Represent CohortSpec as an immutable workflow record."""
    name: str
    path: Path
    seed_offset: int = 0
    expected_n: int | None = None


def canonical_value(value: object) -> str:
    """Convert value to a stable canonical representation."""
    text = str(value).strip()
    try:
        number = float(text)
    except ValueError:
        return text
    return str(int(number)) if number.is_integer() else text


def open_text(path: Path, mode: str):
    """Open text while handling plain-text and compressed files."""
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_participants(
    path: Path,
    *,
    id_column: str,
    outcome_column: str,
    treatment_column: str,
    allowed_outcomes: Sequence[str] = ("0", "1"),
    allowed_treatments: Sequence[str] = ("1", "2"),
) -> list[Participant]:
    """Read participants from disk and validate its basic structure."""
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {id_column, outcome_column, treatment_column}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            missing = sorted(required.difference(reader.fieldnames or []))
            raise ValueError(f"Missing required columns in {path}: {', '.join(missing)}")

        participants: list[Participant] = []
        for row_number, row in enumerate(reader, start=2):
            participant_id = canonical_value(row[id_column])
            outcome = canonical_value(row[outcome_column])
            treatment = canonical_value(row[treatment_column])
            if not participant_id:
                raise ValueError(f"Missing participant ID at {path}:{row_number}")
            if outcome not in allowed_outcomes:
                raise ValueError(f"Unexpected outcome {outcome!r} for {participant_id}")
            if treatment not in allowed_treatments:
                raise ValueError(f"Unexpected treatment {treatment!r} for {participant_id}")
            participants.append(Participant(participant_id, outcome, treatment))

    counts = Counter(p.participant_id for p in participants)
    duplicates = sorted(key for key, count in counts.items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate participant IDs in {path}: {', '.join(duplicates[:10])}")
    return participants


def expected_strata(participants: Sequence[Participant]) -> tuple[tuple[str, str], ...]:
    """Return the four required outcome-by-treatment strata."""
    outcomes = sorted({p.outcome for p in participants})
    treatments = sorted({p.treatment for p in participants})
    return tuple((outcome, treatment) for outcome in outcomes for treatment in treatments)


def assign_joint_stratified_folds(
    participants: Sequence[Participant],
    *,
    number_of_folds: int,
    seed: int,
    strata: Sequence[tuple[str, str]] | None = None,
) -> dict[str, int]:
    """Assign deterministic balanced folds within each joint stratum.

    Each stratum is independently shuffled. Remainder observations are assigned
    to folds with the smallest current total, with a seeded random tie order.
    This reproduces the specified balancing rule and keeps final fold sizes within
    one participant.
    """

    if number_of_folds < 2:
        raise ValueError("number_of_folds must be at least 2")
    groups: dict[tuple[str, str], list[str]] = defaultdict(list)
    for participant in participants:
        groups[participant.stratum].append(participant.participant_id)

    strata = tuple(strata or expected_strata(participants))
    if set(groups) != set(strata):
        raise ValueError(f"Observed strata {sorted(groups)} do not match expected {sorted(strata)}")
    if any(len(groups[key]) < number_of_folds for key in strata):
        small = {key: len(groups[key]) for key in strata if len(groups[key]) < number_of_folds}
        raise ValueError(f"Every stratum must contain at least one member per fold: {small}")

    assignments: dict[str, int] = {}
    fold_totals = [0] * number_of_folds

    for group_index, group_key in enumerate(strata):
        identifiers = sorted(groups[group_key])
        generator = random.Random(seed + group_index * 100_003)
        generator.shuffle(identifiers)

        base_count, remainder = divmod(len(identifiers), number_of_folds)
        tie_order = list(range(number_of_folds))
        generator.shuffle(tie_order)
        tie_priority = {fold: priority for priority, fold in enumerate(tie_order)}
        extra_folds = set(
            sorted(
                range(number_of_folds),
                key=lambda fold: (fold_totals[fold], tie_priority[fold]),
            )[:remainder]
        )

        position = 0
        for fold_index in range(number_of_folds):
            count = base_count + int(fold_index in extra_folds)
            for participant_id in identifiers[position : position + count]:
                if participant_id in assignments:
                    raise RuntimeError(f"Duplicate assignment for {participant_id}")
                assignments[participant_id] = fold_index + 1
            position += count
            fold_totals[fold_index] += count
        if position != len(identifiers):
            raise RuntimeError(f"Incomplete assignment for stratum {group_key}")

    if len(assignments) != len(participants):
        raise RuntimeError("Not every participant received a fold")
    if max(fold_totals) - min(fold_totals) > 1:
        raise RuntimeError(f"Final fold sizes differ by more than one: {fold_totals}")
    return assignments


def subset_by_fold(
    participants: Sequence[Participant],
    assignments: Mapping[str, int],
    fold: int,
    *,
    selected: bool,
) -> list[Participant]:
    """Select participants assigned to the requested fold set."""
    return [
        p
        for p in participants
        if (assignments[p.participant_id] == fold) is selected
    ]


def count_summary(participants: Iterable[Participant]) -> dict[str, int]:
    """Count summary without loading unnecessary participant-level data."""
    participants = list(participants)
    outcomes = Counter(p.outcome for p in participants)
    treatments = Counter(p.treatment for p in participants)
    joint = Counter(p.stratum for p in participants)
    return {
        "participants": len(participants),
        "outcome_0": outcomes.get("0", 0),
        "outcome_1": outcomes.get("1", 0),
        "drug_1": treatments.get("1", 0),
        "drug_2": treatments.get("2", 0),
        "outcome_0_drug_1": joint.get(("0", "1"), 0),
        "outcome_0_drug_2": joint.get(("0", "2"), 0),
        "outcome_1_drug_1": joint.get(("1", "1"), 0),
        "outcome_1_drug_2": joint.get(("1", "2"), 0),
        "all_four_joint_strata_present": "TRUE" if all(value > 0 for value in joint.values()) and len(joint) == 4 else "FALSE",
    }


def generate_cohort_folds(
    participants: Sequence[Participant],
    *,
    cohort_name: str,
    seed_offset: int,
    design: FoldDesign,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    """Generate repeated outer folds and fixed inner folds for one configured cohort."""
    by_id = {p.participant_id: p for p in participants}
    ordered_ids = sorted(by_id)
    strata = expected_strata(participants)
    outer_rows: list[dict[str, object]] = []
    inner_rows: list[dict[str, object]] = []
    balance_rows: list[dict[str, object]] = []
    signatures: set[tuple[int, ...]] = set()

    for repeat in range(1, design.outer_repeats + 1):
        outer_seed = design.master_seed + seed_offset + repeat * 10_000
        outer_assignments = assign_joint_stratified_folds(
            participants,
            number_of_folds=design.outer_folds,
            seed=outer_seed,
            strata=strata,
        )
        signature = tuple(outer_assignments[participant_id] for participant_id in ordered_ids)
        if signature in signatures:
            raise RuntimeError(f"Repeated outer assignment in {cohort_name}, repeat {repeat}")
        signatures.add(signature)

        for participant_id in ordered_ids:
            participant = by_id[participant_id]
            outer_rows.append(
                {
                    "cohort": cohort_name,
                    "outer_repeat": repeat,
                    "participant_id": participant_id,
                    "outcome": participant.outcome,
                    "drug": participant.treatment,
                    "outer_fold": outer_assignments[participant_id],
                    "outer_seed": outer_seed,
                    "joint_stratum": f"{participant.outcome}_{participant.treatment}",
                }
            )

        for outer_fold in range(1, design.outer_folds + 1):
            outer_test = subset_by_fold(participants, outer_assignments, outer_fold, selected=True)
            outer_train = subset_by_fold(participants, outer_assignments, outer_fold, selected=False)
            for role, subset in (("validation", outer_test), ("training", outer_train)):
                balance_rows.append(
                    {
                        "cohort": cohort_name,
                        "split_level": "outer",
                        "outer_repeat": repeat,
                        "outer_fold": outer_fold,
                        "inner_fold": "",
                        "role": role,
                        **count_summary(subset),
                    }
                )

            inner_seed = design.master_seed + seed_offset + repeat * 1_000_000 + outer_fold * 1_000
            inner_assignments = assign_joint_stratified_folds(
                outer_train,
                number_of_folds=design.inner_folds,
                seed=inner_seed,
                strata=strata,
            )
            for participant in sorted(outer_train, key=lambda item: item.participant_id):
                inner_rows.append(
                    {
                        "cohort": cohort_name,
                        "outer_repeat": repeat,
                        "held_out_outer_fold": outer_fold,
                        "participant_id": participant.participant_id,
                        "outcome": participant.outcome,
                        "drug": participant.treatment,
                        "inner_fold": inner_assignments[participant.participant_id],
                        "inner_seed": inner_seed,
                        "joint_stratum": f"{participant.outcome}_{participant.treatment}",
                    }
                )

            for inner_fold in range(1, design.inner_folds + 1):
                inner_valid = subset_by_fold(outer_train, inner_assignments, inner_fold, selected=True)
                inner_train = subset_by_fold(outer_train, inner_assignments, inner_fold, selected=False)
                for role, subset in (("validation", inner_valid), ("training", inner_train)):
                    balance_rows.append(
                        {
                            "cohort": cohort_name,
                            "split_level": "inner",
                            "outer_repeat": repeat,
                            "outer_fold": outer_fold,
                            "inner_fold": inner_fold,
                            "role": role,
                            **count_summary(subset),
                        }
                    )

    validate_generated_folds(participants, outer_rows, inner_rows, balance_rows, design)
    return outer_rows, inner_rows, balance_rows


def validate_generated_folds(
    participants: Sequence[Participant],
    outer_rows: Sequence[Mapping[str, object]],
    inner_rows: Sequence[Mapping[str, object]],
    balance_rows: Sequence[Mapping[str, object]],
    design: FoldDesign,
) -> None:
    """Check fold completeness, exclusivity, balance and absence of train-test overlap."""
    identifiers = {p.participant_id for p in participants}
    expected_outer = len(participants) * design.outer_repeats
    expected_inner = len(participants) * design.outer_repeats * (design.outer_folds - 1)
    if len(outer_rows) != expected_outer:
        raise RuntimeError(f"Outer manifest rows {len(outer_rows)} != {expected_outer}")
    if len(inner_rows) != expected_inner:
        raise RuntimeError(f"Inner manifest rows {len(inner_rows)} != {expected_inner}")

    outer_counts = Counter(str(row["participant_id"]) for row in outer_rows)
    inner_counts = Counter(str(row["participant_id"]) for row in inner_rows)
    if set(outer_counts) != identifiers or set(outer_counts.values()) != {design.outer_repeats}:
        raise RuntimeError("Outer participant appearance counts are invalid")
    expected_inner_appearances = design.outer_repeats * (design.outer_folds - 1)
    if set(inner_counts) != identifiers or set(inner_counts.values()) != {expected_inner_appearances}:
        raise RuntimeError("Inner participant appearance counts are invalid")

    for level in ("outer", "inner"):
        validation = [row for row in balance_rows if row["split_level"] == level and row["role"] == "validation"]
        contexts: dict[tuple[object, ...], list[Mapping[str, object]]] = defaultdict(list)
        for row in validation:
            key = (row["outer_repeat"],) if level == "outer" else (row["outer_repeat"], row["outer_fold"])
            contexts[key].append(row)
        expected_context_size = design.outer_folds if level == "outer" else design.inner_folds
        for key, rows in contexts.items():
            if len(rows) != expected_context_size:
                raise RuntimeError(f"Incomplete {level} context {key}")
            for field in (
                "participants",
                "outcome_0_drug_1",
                "outcome_0_drug_2",
                "outcome_1_drug_1",
                "outcome_1_drug_2",
            ):
                values = [int(row[field]) for row in rows]
                if max(values) - min(values) > 1:
                    raise RuntimeError(f"Unbalanced {level} {field} in context {key}: {values}")
            if not all(str(row["all_four_joint_strata_present"]).upper() == "TRUE" for row in rows):
                raise RuntimeError(f"Missing joint stratum in {level} context {key}")


def write_tsv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    """Write a deterministic TSV workflow output."""
    if not rows:
        raise ValueError(f"Cannot write empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
