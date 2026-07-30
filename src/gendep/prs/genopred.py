"""Prepare the fixed eight-trait GenoPred/PRS-CS-auto execution bundle."""

from __future__ import annotations

import argparse
import csv
import json
import shlex
from pathlib import Path
from typing import Any

import yaml

from .gwas import load_trait_specs


def _resolve(value: str | Path, project_root: Path) -> Path:
    path = Path(value).expanduser()
    return path if path.is_absolute() else project_root / path


def _load_yaml(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a YAML mapping: {path}")
    return payload


def _read_trait_labels(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {
            row["trait"]: row["label"]
            for row in csv.DictReader(handle, delimiter="\t")
        }


def write_gwas_list(
    path: Path,
    standardised_dir: Path,
    gwas_config: Path,
    traits_file: Path,
    *,
    require_inputs: bool,
) -> None:
    """Write gwas list to a deterministic workflow output."""
    labels = _read_trait_labels(traits_file)
    specs = load_trait_specs(gwas_config)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ["name", "path", "population", "n", "sampling", "prevalence", "mean", "sd", "label"]
        )
        for spec in specs:
            source = (standardised_dir / f"{spec.trait}.standardised.tsv.gz").resolve()
            if require_inputs and (not source.is_file() or source.stat().st_size == 0):
                raise FileNotFoundError(f"Standardised GWAS is missing: {source}")
            writer.writerow(
                [
                    spec.trait,
                    str(source),
                    "EUR",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    "NA",
                    '"' + labels.get(spec.trait, spec.trait).replace('"', "") + '"',
                ]
            )


def write_target_list(path: Path, target_prefix: Path, *, require_inputs: bool) -> None:
    """Write target list to a deterministic workflow output."""
    if require_inputs:
        for chromosome in range(1, 23):
            for extension in ("pgen", "pvar", "psam"):
                candidate = Path(f"{target_prefix}.chr{chromosome}.{extension}")
                if not candidate.is_file() or candidate.stat().st_size == 0:
                    raise FileNotFoundError(f"Target pfile component is missing: {candidate}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["name", "path", "type", "indiv_report"])
        writer.writerow(["GENDEP", str(target_prefix.resolve()), "plink2", "F"])


def prepare_run(
    project_config: Path,
    genopred_config: Path,
    gwas_config: Path,
    traits_file: Path,
    *,
    require_inputs: bool = True,
) -> dict[str, str]:
    """Prepare portable GenoPred target, GWAS and configuration files for the eight-trait analysis."""
    project = _load_yaml(project_config)
    fixed = _load_yaml(genopred_config)
    project_root = Path(project.get("project_root", project_config.parent.parent)).expanduser()
    if not project_root.is_absolute():
        project_root = (project_config.parent / project_root).resolve()

    external = project.get("external_resources", {})
    outputs = project.get("outputs", {})
    prs = project.get("prs", {})

    genopred_root = _resolve(external["genopred_root"], project_root)
    base_config = _resolve(external["genopred_base_config"], project_root)
    snakemake = _resolve(external["genopred_snakemake"], project_root)
    rscript = _resolve(external["genopred_runtime_rscript"], project_root)
    runtime_python = _resolve(external["genopred_runtime_python"], project_root)
    plink = _resolve(external["genopred_runtime_plink"], project_root)
    plink2 = _resolve(external["genopred_runtime_plink2"], project_root)
    custom_reference = _resolve(external["genopred_custom_reference"], project_root)
    resource_directory = _resolve(external["genopred_resource_directory"], project_root)

    standardised_dir = _resolve(outputs["gwas_standardised"], project_root)
    run_root = _resolve(outputs["genopred_run"], project_root)
    target_prefix = _resolve(
        prs.get("target_prefix", "work/genotype/target_validation/pfiles/GENDEP"),
        project_root,
    )
    bundle_dir = run_root / "config"
    pipeline_output = run_root / "pipeline_output"
    bundle_dir.mkdir(parents=True, exist_ok=True)
    pipeline_output.mkdir(parents=True, exist_ok=True)

    required_files = [base_config]
    required_directories = [genopred_root / "pipeline", custom_reference, resource_directory]
    required_executables = [snakemake, rscript, runtime_python, plink, plink2]
    if require_inputs:
        for path in required_files:
            if not path.is_file():
                raise FileNotFoundError(path)
        for path in required_directories:
            if not path.is_dir():
                raise FileNotFoundError(path)
        for path in required_executables:
            if not path.is_file():
                raise FileNotFoundError(path)

    gwas_list = bundle_dir / "gwas_list.txt"
    target_list = bundle_dir / "target_list.txt"
    run_config = bundle_dir / "config.yaml"
    write_gwas_list(
        gwas_list,
        standardised_dir,
        gwas_config,
        traits_file,
        require_inputs=require_inputs,
    )
    write_target_list(target_list, target_prefix, require_inputs=require_inputs)

    base = _load_yaml(base_config) if base_config.is_file() else {}
    base.update(
        {
            "outdir": str(pipeline_output.resolve()),
            "resdir": str(resource_directory.resolve()),
            "refdir": str(custom_reference.resolve()),
            "config_file": str(run_config.resolve()),
            "gwas_list": str(gwas_list.resolve()),
            "target_list": str(target_list.resolve()),
            "score_list": "NA",
            "pgs_methods": ["prscs"],
            "prscs_phi": ["auto"],
            "prscs_ldref": "1kg",
            "testing": "NA",
        }
    )
    run_config.write_text(yaml.safe_dump(base, sort_keys=False), encoding="utf-8")

    software = fixed["software"]
    environment = {
        "PROJECT_ROOT": str(project_root.resolve()),
        "GENOPRED_ROOT": str(genopred_root.resolve()),
        "GENOPRED_PIPELINE": str((genopred_root / "pipeline").resolve()),
        "SNAKEMAKE_BIN": str(snakemake.resolve()),
        "RSCRIPT_BIN": str(rscript.resolve()),
        "RUNTIME_PYTHON": str(runtime_python.resolve()),
        "PLINK_BIN": str(plink.resolve()),
        "PLINK2_BIN": str(plink2.resolve()),
        "RUN_CONFIG": str(run_config.resolve()),
        "PIPELINE_OUTPUT": str(pipeline_output.resolve()),
        "RUN_ROOT": str(run_root.resolve()),
        "EXPECTED_GENOPRED_COMMIT": str(software["genopred_commit"]),
        "EXPECTED_GENOUTILS_COMMIT": str(software["genoutils_commit"]),
    }
    env_path = bundle_dir / "run.env"
    env_path.write_text(
        "\n".join(f"{key}={shlex.quote(value)}" for key, value in environment.items()) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "project_config": str(project_config.resolve()),
        "genopred_config": str(genopred_config.resolve()),
        "gwas_config": str(gwas_config.resolve()),
        "traits_file": str(traits_file.resolve()),
        "gwas_list": str(gwas_list.resolve()),
        "target_list": str(target_list.resolve()),
        "run_config": str(run_config.resolve()),
        "pipeline_output": str(pipeline_output.resolve()),
        "genopred_commit": software["genopred_commit"],
        "genoutils_commit": software["genoutils_commit"],
        "pgs_methods": ["prscs"],
        "prscs_phi": ["auto"],
        "prscs_ldref": "1kg",
        "testing": "NA",
    }
    (bundle_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return environment


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments for the GenoPred configuration command."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-config", required=True)
    parser.add_argument("--genopred-config", required=True)
    parser.add_argument("--gwas-config", required=True)
    parser.add_argument("--traits-file", required=True)
    parser.add_argument("--allow-missing-inputs", action="store_true")
    return parser.parse_args()


def main() -> None:
    """Prepare the fixed GenoPred run configuration and manifests."""
    args = parse_args()
    environment = prepare_run(
        Path(args.project_config).resolve(),
        Path(args.genopred_config).resolve(),
        Path(args.gwas_config).resolve(),
        Path(args.traits_file).resolve(),
        require_inputs=not args.allow_missing_inputs,
    )
    print(f"run_config={environment['RUN_CONFIG']}")
    print(f"pipeline_output={environment['PIPELINE_OUTPUT']}")
    print("GENOPRED_BUNDLE=PASS")


if __name__ == "__main__":
    main()
