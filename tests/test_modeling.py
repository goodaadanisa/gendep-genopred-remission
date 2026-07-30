from __future__ import annotations

# ---- test_modeling_metrics.py ----

import numpy as np

from gendep.modeling.metrics import auc, brier, calibration, log_loss, metric_bundle


def test_binary_metrics_match_reference_values() -> None:
    outcome = np.array([0, 0, 1, 1])
    prediction = np.array([0.1, 0.4, 0.35, 0.8])
    assert auc(outcome, prediction) == 0.75
    assert np.isclose(brier(outcome, prediction), 0.158125)
    assert np.isfinite(log_loss(outcome, prediction))


def test_auc_uses_average_ranks_for_ties() -> None:
    outcome = np.array([0, 1, 0, 1])
    prediction = np.array([0.2, 0.2, 0.8, 0.8])
    assert auc(outcome, prediction) == 0.5


def test_calibration_and_bundle_are_finite() -> None:
    outcome = np.array([0, 0, 0, 1, 1, 1])
    prediction = np.array([0.1, 0.2, 0.4, 0.55, 0.7, 0.9])
    result = calibration(outcome, prediction)
    assert np.isfinite(result.intercept)
    assert np.isfinite(result.slope)
    bundle = metric_bundle(outcome, prediction)
    assert set(bundle) == {
        "AUC",
        "Brier",
        "log_loss",
        "calibration_intercept",
        "calibration_slope",
    }

# ---- test_model_aggregation.py ----

import pandas as pd
import pytest

from gendep.modeling.aggregation import AnalysisSpec, aggregate_outputs


def synthetic_predictions() -> pd.DataFrame:
    rows = []
    participants = ["P1", "P2", "P3", "P4"]
    outcomes = [0, 0, 1, 1]
    for repeat in (1, 2):
        for fold, selected in ((1, (0, 2)), (2, (1, 3))):
            for index in selected:
                base = [0.10, 0.30, 0.65, 0.85][index]
                rows.append(
                    {
                        "analysis": "synthetic",
                        "algorithm": "elastic_net",
                        "outer_repeat": repeat,
                        "outer_fold": fold,
                        "participant_id": participants[index],
                        "outcome": outcomes[index],
                        "drug": 1 + index % 2,
                        "comparator_probability": base,
                        "combined_probability": min(0.99, base + (0.01 if index >= 2 else -0.01)),
                    }
                )
    return pd.DataFrame(rows)


def test_aggregate_outputs_checks_repeat_level_predictions() -> None:
    predictions = synthetic_predictions()
    tuning = pd.DataFrame(
        {
            "model": ["clinical", "combined"] * 4,
            "alpha": [0.0, 0.1] * 4,
            "selected": [True, True] * 4,
        }
    )
    raw = {
        "predictions": predictions,
        "metrics": pd.DataFrame(),
        "comparison": pd.DataFrame(),
        "tuning": tuning,
        "coefficients": pd.DataFrame(),
        "warnings": pd.DataFrame(),
        "missing_splits": pd.DataFrame(columns=["split_directory"]),
    }
    spec = AnalysisSpec(
        analysis="synthetic",
        algorithm="elastic_net",
        cohort_name="synthetic",
        expected_participants=4,
        comparator_label="clinical",
        combined_label="combined",
    )
    outputs = aggregate_outputs(raw, spec, outer_repeats=2, outer_folds=2, epsilon=1e-15)
    assert len(outputs["predictions"]) == 8
    assert len(outputs["repeat_metrics"]) == 2
    assert len(outputs["participant_mean_predictions"]) == 4
    assert outputs["aggregation_validation"]["result"].eq("PASS").all()

# ---- test_orientation_sensitivity.py ----

import csv
import gzip
from pathlib import Path

import numpy as np
import pandas as pd

from gendep.sensitivity.orientation import (
    AGREEMENT_DECISION,
    PRS_COLUMNS,
    build_conservative_analysis_base,
    collect_orientation_profiles,
    partition_orientation_markers,
    prepare_genopred_scoring_bundle,
    zero_score_file,
)


def test_partition_and_zero_weight_score_file(tmp_path: Path) -> None:
    orientation = tmp_path / "orientation.tsv.gz"
    with gzip.open(orientation, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SNP", "ORIENTATION_DECISION"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"SNP": "rs1", "ORIENTATION_DECISION": AGREEMENT_DECISION})
        writer.writerow({"SNP": "rs2", "ORIENTATION_DECISION": "retain_eur_major_sample_flip_near_half"})
        writer.writerow({"SNP": "rs3", "ORIENTATION_DECISION": "retain_eur_major_sample_tie"})
    partition = partition_orientation_markers(orientation)
    assert partition.total == 3
    assert partition.retained_ids == frozenset({"rs1"})
    assert partition.excluded_ids == frozenset({"rs2", "rs3"})

    source = tmp_path / "score.gz"
    destination = tmp_path / "conservative.score.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write("SNP A1 A2 SCORE_phi_auto\n")
        handle.write("rs1 A G 0.5\n")
        handle.write("rs2 C T -0.2\n")
        handle.write("rs3 A C 0\n")
    summary = zero_score_file(source, destination, partition.excluded_ids, expected_rows=3)
    assert summary["score_rows"] == 3
    assert summary["source_nonzero_weights"] == 2
    assert summary["excluded_nonzero_weights"] == 1
    assert np.isclose(summary["sum_absolute_excluded_weights"], 0.2)
    with gzip.open(destination, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "rs1 A G 0.5" in text
    assert "rs2 C T 0" in text


def test_analysis_base_replacement_changes_only_prs(tmp_path: Path) -> None:
    identifiers = [f"P{index:03d}" for index in range(1, 11)]
    primary = pd.DataFrame({"Row.names": identifiers, "drug": [1, 2] * 5, "hdremit.all": [0, 1] * 5})
    conservative = pd.DataFrame({"FID": identifiers, "IID": identifiers})
    for position, column in enumerate(PRS_COLUMNS):
        values = np.linspace(-1, 1, len(identifiers)) + position
        primary[column] = values
        conservative[column] = values + 0.01 * (position + 1)
    accepted_path = tmp_path / "primary.tsv.gz"
    conservative_path = tmp_path / "conservative.tsv.gz"
    output_path = tmp_path / "output.tsv.gz"
    primary.to_csv(accepted_path, sep="\t", index=False, compression="gzip")
    conservative = conservative.sample(frac=1, random_state=10)
    conservative.to_csv(conservative_path, sep="\t", index=False, compression="gzip")
    output, summary = build_conservative_analysis_base(accepted_path, conservative_path, output_path)
    pd.testing.assert_frame_equal(
        output[["Row.names", "drug", "hdremit.all"]],
        primary[["Row.names", "drug", "hdremit.all"]],
        check_dtype=False,
    )
    assert len(summary) == 8
    assert (summary["participants_with_changed_score"] == 10).all()



def _write_profile(path: Path, values: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame = pd.DataFrame(
        {
            "#FID": [f"P{index:03d}" for index in range(1, 431)],
            "IID": [f"P{index:03d}" for index in range(1, 431)],
            "SCORE": values,
        }
    )
    frame.to_csv(path, sep="\t", index=False)


def test_prepare_scoring_bundle_and_collect_profiles(tmp_path: Path) -> None:
    # The controlled orientation workflow is executed on Linux HPC. Windows may
    # prohibit symlink creation for ordinary users, so probe capability and skip
    # this platform-specific integration test rather than reporting a code failure.
    probe_source = tmp_path / "symlink_probe_source"
    probe_destination = tmp_path / "symlink_probe_destination"
    probe_source.mkdir()
    try:
        probe_destination.symlink_to(probe_source, target_is_directory=True)
    except OSError as exc:
        pytest.skip(f"Symbolic links are unavailable on this host: {exc}")
    else:
        probe_destination.unlink()
        probe_source.rmdir()
    primary = tmp_path / "primary"
    (primary / "GENDEP" / "geno").mkdir(parents=True)
    (primary / "GENDEP" / "pcs").mkdir(parents=True)
    marker = primary / "reference" / "target_checks" / "GENDEP" / "ancestry_reporter.done"
    marker.parent.mkdir(parents=True)
    marker.write_text("done\n", encoding="utf-8")
    accepted_config = tmp_path / "config.yaml"
    accepted_config.write_text("outdir: /old/output\nother: fixed\n", encoding="utf-8")
    conservative_scores = tmp_path / "score_files"

    baseline_values = np.round(np.linspace(-2, 2, 430), 3)
    for trait_index, trait in enumerate(("MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA")):
        score_dir = primary / "reference" / "pgs_score_files" / "prscs" / trait
        score_dir.mkdir(parents=True, exist_ok=True)
        (score_dir / f"ref-{trait}.score.gz").write_bytes(b"score")
        (score_dir / f"ref-{trait}-TRANS.model.rds").write_bytes(b"model")
        conservative_dir = conservative_scores / trait
        conservative_dir.mkdir(parents=True)
        (conservative_dir / f"ref-{trait}.score.gz").write_bytes(b"conservative")
        relative = Path("GENDEP") / "pgs" / "TRANS" / "prscs" / trait / f"GENDEP-{trait}-TRANS.profiles"
        accepted_trait_values = baseline_values + trait_index
        _write_profile(primary / relative, accepted_trait_values)

    bundle = tmp_path / "bundle"
    paths = prepare_genopred_scoring_bundle(
        accepted_output=primary,
        accepted_config=accepted_config,
        output_dir=bundle,
        conservative_score_root=conservative_scores,
    )
    assert f"outdir: {bundle / 'baseline'}" in paths["baseline_config"].read_text()
    assert (bundle / "baseline" / "GENDEP" / "geno").is_symlink()
    assert (bundle / "conservative" / "reference" / "pgs_score_files" / "prscs" / "MDD" / "ref-MDD.score.gz").is_symlink()

    for trait_index, trait in enumerate(("MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA")):
        relative = Path("GENDEP") / "pgs" / "TRANS" / "prscs" / trait / f"GENDEP-{trait}-TRANS.profiles"
        values = baseline_values + trait_index
        _write_profile(bundle / "baseline" / relative, values)
        _write_profile(bundle / "conservative" / relative, values + 0.001 * (trait_index + 1))

    matrix, baseline, comparison = collect_orientation_profiles(
        accepted_output=primary,
        scoring_bundle=bundle,
        output_matrix=tmp_path / "conservative_prs.tsv.gz",
    )
    assert matrix.shape == (430, 10)
    assert baseline["exact_three_decimal_reproduction"].all()
    assert (comparison["participants_with_changed_score"] == 430).all()
