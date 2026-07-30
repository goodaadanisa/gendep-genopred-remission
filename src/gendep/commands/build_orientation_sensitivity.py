#!/usr/bin/env python3
"""Prepare the conservative target-orientation sensitivity without changing primary outputs.

Purpose
-------
Identify the tied or discordant PRS-CS-compatible markers, zero their primary
posterior weights in copied score files, and optionally replace only the eight
PRS columns in a copied analysis base after conservative target scoring.

Inputs
------
The Stage 2 PRS-CS-compatible orientation table and primary trait-specific
score files. Optional GenoPred paths execute exact baseline and conservative
target scoring and collect the conservative 430 x 8 PRS matrix.

Outputs
-------
Flagged/retained marker lists, copied conservative score files, an exact
baseline-reproduction audit, conservative participant scores, a copied analysis
base and aggregate score comparisons.

Usage
-----
gendep orientation-sensitivity \
  --orientation-table work/genotype/reconstruction/prscs_compatible_markers.tsv.gz \
  --primary-score-root work/prs/genopred/pipeline_output/reference/pgs_score_files/prscs \
  --output-dir work/sensitivity/orientation

For the complete primary target-scoring route, also provide:
  --primary-genopred-output work/prs/genopred/pipeline_output \
  --primary-genopred-config work/prs/genopred/config/config.yaml \
  --genopred-root /path/to/GenoPred --rscript /path/to/Rscript \
  --plink2 /path/to/plink2 --run-scoring \
  --primary-analysis-base work/analysis/datasets/gendep_clinical_prs_pc_analysis_base.tsv.gz
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from gendep.sensitivity.orientation import (  # noqa: E402
    EXPECTED_TRAITS,
    build_conservative_analysis_base,
    collect_orientation_profiles,
    partition_orientation_markers,
    prepare_genopred_scoring_bundle,
    zero_score_file,
)


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the orientation-sensitivity build command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--orientation-table", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--primary-score-root", type=Path)
    parser.add_argument("--score-template", default="{trait}/ref-{trait}.score.gz")
    parser.add_argument("--weight-column", default="SCORE_phi_auto")
    parser.add_argument("--primary-analysis-base", type=Path)
    parser.add_argument("--conservative-prs", type=Path)
    parser.add_argument("--primary-genopred-output", type=Path)
    parser.add_argument("--primary-genopred-config", type=Path)
    parser.add_argument("--genopred-root", type=Path)
    parser.add_argument("--target-scoring-script", type=Path)
    parser.add_argument("--rscript", type=Path)
    parser.add_argument("--plink2", type=Path)
    parser.add_argument("--cores", type=int, default=4)
    parser.add_argument("--run-scoring", action="store_true")
    parser.add_argument("--collect-profiles", action="store_true")
    parser.add_argument("--target-name", default="GENDEP")
    parser.add_argument("--population", default="TRANS")
    parser.add_argument("--expected-score-rows", type=int, default=147370)
    parser.add_argument("--expected-compatible", type=int, default=141038)
    parser.add_argument("--expected-excluded", type=int, default=953)
    parser.add_argument("--expected-retained", type=int, default=140085)
    return parser.parse_args()


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_ids(path: Path, values: frozenset[str]) -> None:
    """Write ids to a deterministic workflow output."""
    path.write_text("\n".join(sorted(values)) + "\n", encoding="utf-8")


def run_target_scoring(
    *,
    mode: str,
    output_dir: Path,
    config_path: Path,
    target_scoring_script: Path,
    rscript: Path,
    plink2: Path,
    target_name: str,
    population: str,
    cores: int,
) -> None:
    """Run target scoring and record its validation status."""
    log_dir = output_dir / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    command = [
        str(rscript),
        str(target_scoring_script),
        "--config",
        str(config_path),
        "--name",
        target_name,
        "--population",
        population,
        "--plink2",
        str(plink2),
        "--test",
        "NA",
        "--n_cores",
        str(cores),
    ]
    with (log_dir / f"{mode}_target_scoring.log").open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            command,
            cwd=target_scoring_script.parents[2] / "pipeline",
            stdout=handle,
            stderr=subprocess.STDOUT,
            check=False,
        )
    if completed.returncode != 0:
        raise SystemExit(f"{mode} target scoring failed with status {completed.returncode}")


def main() -> None:
    """Build the conservative orientation-sensitivity scoring inputs."""
    args = parse_args()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    partition = partition_orientation_markers(args.orientation_table.resolve())
    observed = (partition.total, len(partition.excluded_ids), len(partition.retained_ids))
    expected = (args.expected_compatible, args.expected_excluded, args.expected_retained)
    if observed != expected:
        raise SystemExit(f"Orientation counts differ: observed={observed}, expected={expected}")

    excluded_path = output_dir / "excluded_orientation_variants.txt"
    retained_path = output_dir / "retained_orientation_variants.txt"
    write_ids(excluded_path, partition.excluded_ids)
    write_ids(retained_path, partition.retained_ids)

    decision_rows = [
        {"orientation_decision": key, "variant_count": value}
        for key, value in sorted(partition.decision_counts.items())
    ]
    pd.DataFrame(decision_rows).to_csv(
        output_dir / "orientation_decision_counts.tsv", sep="\t", index=False, lineterminator="\n"
    )

    bundle_requested = bool(args.accepted_genopred_output) or bool(args.accepted_genopred_config)
    if bundle_requested and not (args.accepted_genopred_output and args.accepted_genopred_config):
        raise SystemExit("Provide both --primary-genopred-output and --primary-genopred-config")

    attrition_rows: list[dict[str, object]] = []
    conservative_score_root = output_dir / "score_files"
    if args.accepted_score_root:
        root = args.accepted_score_root.resolve()
        for trait in EXPECTED_TRAITS:
            relative = Path(args.score_template.format(trait=trait))
            source = root / relative
            destination = conservative_score_root / relative
            row = zero_score_file(
                source,
                destination,
                partition.excluded_ids,
                weight_column=args.weight_column,
                expected_rows=args.expected_score_rows,
            )
            row["trait"] = trait
            attrition_rows.append(row)
        pd.DataFrame(attrition_rows)[
            [
                "trait",
                "score_rows",
                "source_nonzero_weights",
                "excluded_variants_present",
                "excluded_nonzero_weights",
                "sum_absolute_excluded_weights",
                "score_file",
            ]
        ].to_csv(output_dir / "trait_weight_attrition.tsv", sep="\t", index=False, lineterminator="\n")

    bundle_paths: dict[str, Path] | None = None
    if bundle_requested:
        if not attrition_rows:
            raise SystemExit("--primary-score-root is required to prepare a scoring bundle")
        bundle_root = output_dir / "scoring_bundle"
        bundle_paths = prepare_genopred_scoring_bundle(
            accepted_output=args.accepted_genopred_output.resolve(),
            accepted_config=args.accepted_genopred_config.resolve(),
            output_dir=bundle_root,
            conservative_score_root=conservative_score_root,
            target_name=args.target_name,
        )

    if args.run_scoring:
        if bundle_paths is None:
            raise SystemExit("--run-scoring requires the primary GenoPred output and config")
        if not args.rscript or not args.plink2:
            raise SystemExit("--run-scoring requires --rscript and --plink2")
        if args.target_scoring_script:
            scoring_script = args.target_scoring_script.resolve()
        elif args.genopred_root:
            scoring_script = args.genopred_root.resolve() / "Scripts" / "target_scoring" / "target_scoring_pipeline.R"
        else:
            raise SystemExit("Provide --target-scoring-script or --genopred-root")
        if not scoring_script.exists():
            raise SystemExit(f"Target-scoring script was not found: {scoring_script}")
        if args.cores < 1:
            raise SystemExit("--cores must be positive")
        for mode in ("baseline", "conservative"):
            run_target_scoring(
                mode=mode,
                output_dir=output_dir,
                config_path=bundle_paths[f"{mode}_config"],
                target_scoring_script=scoring_script,
                rscript=args.rscript.resolve(),
                plink2=args.plink2.resolve(),
                target_name=args.target_name,
                population=args.population,
                cores=args.cores,
            )

    collected_prs_path: Path | None = None
    baseline_validation_path: Path | None = None
    profile_comparison_path: Path | None = None
    if args.collect_profiles or args.run_scoring:
        if not args.accepted_genopred_output or bundle_paths is None:
            raise SystemExit("Profile collection requires the primary GenoPred output and scoring bundle")
        collected_prs_path = output_dir / "conservative_prs.tsv.gz"
        _, baseline_validation, profile_comparison = collect_orientation_profiles(
            accepted_output=args.accepted_genopred_output.resolve(),
            scoring_bundle=output_dir / "scoring_bundle",
            output_matrix=collected_prs_path,
            target_name=args.target_name,
            population=args.population,
        )
        baseline_validation_path = output_dir / "baseline_reproduction_validation.tsv"
        baseline_validation.to_csv(
            baseline_validation_path, sep="\t", index=False, lineterminator="\n"
        )
        profile_comparison_path = output_dir / "conservative_prs_comparison.tsv"
        profile_comparison.to_csv(
            profile_comparison_path, sep="\t", index=False, lineterminator="\n"
        )
        if args.conservative_prs and args.conservative_prs.resolve() != collected_prs_path.resolve():
            raise SystemExit("Do not combine --collect-profiles with a different --conservative-prs file")

    effective_conservative_prs = collected_prs_path or (args.conservative_prs.resolve() if args.conservative_prs else None)
    comparison_path: Path | None = None
    conservative_base_path: Path | None = None
    if bool(args.accepted_analysis_base) != bool(effective_conservative_prs):
        raise SystemExit("Provide both --primary-analysis-base and --conservative-prs, or neither")
    if args.accepted_analysis_base and effective_conservative_prs:
        conservative_base_path = output_dir / "analysis_base" / "gendep_orientation_sensitivity_analysis_base.tsv.gz"
        _, comparison = build_conservative_analysis_base(
            args.accepted_analysis_base.resolve(),
            effective_conservative_prs,
            conservative_base_path,
        )
        comparison_path = output_dir / "analysis_base_prs_replacement_comparison.tsv"
        comparison.to_csv(comparison_path, sep="\t", index=False, lineterminator="\n")
        if profile_comparison_path:
            profile_table = pd.read_csv(profile_comparison_path, sep="\t")
            comparison_columns = [
                "trait",
                "pearson_correlation",
                "spearman_correlation",
                "maximum_absolute_difference",
                "mean_absolute_difference",
                "participants_with_changed_score",
            ]
            if profile_table["trait"].tolist() != comparison["trait"].tolist():
                raise RuntimeError("Profile and analysis-base PRS trait orders differ")
            numeric_columns = [column for column in comparison_columns if column != "trait"]
            profile_values = profile_table[numeric_columns].to_numpy(dtype=float)
            analysis_values = comparison[numeric_columns].to_numpy(dtype=float)
            if not pd.DataFrame(profile_values).round(12).equals(pd.DataFrame(analysis_values).round(12)):
                raise RuntimeError("Profile and analysis-base PRS comparisons differ")

    tracked = [excluded_path, retained_path, output_dir / "orientation_decision_counts.tsv"]
    if attrition_rows:
        tracked.append(output_dir / "trait_weight_attrition.tsv")
    if baseline_validation_path:
        tracked.append(baseline_validation_path)
    if collected_prs_path:
        tracked.append(collected_prs_path)
    if profile_comparison_path:
        tracked.append(profile_comparison_path)
    if comparison_path:
        tracked.extend([comparison_path, conservative_base_path])
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "accepted_outputs_modified": False,
        "prscs_compatible_variants": partition.total,
        "zero_weighted_variants": len(partition.excluded_ids),
        "retained_variants": len(partition.retained_ids),
        "operation": "set_flagged_posterior_weights_to_zero",
        "prs_cs_executed": False,
        "baseline_target_scoring_executed": bool(args.run_scoring),
        "conservative_target_scoring_executed": bool(args.run_scoring),
        "baseline_exact_reproduction_checked": bool(baseline_validation_path),
        "profiles_collected": bool(collected_prs_path),
        "files": {str(path.relative_to(output_dir)): sha256(path) for path in tracked if path},
    }
    (output_dir / "orientation_sensitivity_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "_SUCCESS").write_text("orientation_sensitivity_preparation=PASS\n", encoding="utf-8")
    print("ORIENTATION_SENSITIVITY_PREPARATION=PASS")
    print(f"compatible_variants={partition.total}")
    print(f"zero_weighted_variants={len(partition.excluded_ids)}")
    print(f"retained_variants={len(partition.retained_ids)}")
    print(f"output_dir={output_dir}")


if __name__ == "__main__":
    main()
