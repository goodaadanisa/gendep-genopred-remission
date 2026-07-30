from __future__ import annotations

# ---- test_genotype_reconstruction.py ----

import gzip
from pathlib import Path

from gendep.genotype.common import allele_relation, is_snv_record
from gendep.genotype.orientation import (
    GenotypeFrequency,
    intersect_prscs,
    orient_markers,
    pearson_code2_vs_eur_major,
)
from gendep.genotype.reference import match_genopred_reference, scan_ensembl_vcfs


def write_gzip(path: Path, text: str) -> None:
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def test_allele_relation_and_snv_classification() -> None:
    assert is_snv_record("A", "G")
    assert is_snv_record("A", "G,C")
    assert not is_snv_record("A", "AT")
    assert allele_relation("A", "G", "A", "G") == "same_order"
    assert allele_relation("A", "G", "G", "A") == "swapped_order"
    assert allele_relation("A", "G", "T", "C") == "strand_complement"
    assert allele_relation("A", "G", "C", "T") == "strand_complement_swapped"


def test_reference_orientation_and_prscs_intersection(tmp_path: Path) -> None:
    ensembl = tmp_path / "ensembl.chr1.vcf.gz"
    write_gzip(
        ensembl,
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "1\t200\trs2\tA\tAT\t.\tPASS\t.\n"
        "1\t400\trs4\tG\tA\t.\tPASS\t.\n"
        "1\t500\trs5\tA\tC\t.\tPASS\t.\n",
    )
    requested = ["rs1", "rs2", "rs3", "rs4", "rs5"]
    records, classification, counts = scan_ensembl_vcfs(
        requested, str(tmp_path / "ensembl.chr{chromosome}.vcf.gz"), chromosomes=[1]
    )
    assert counts["reconstructed_autosomal_snvs"] == 3
    assert classification["rs2"] == "excluded_non_snv"
    assert classification["rs3"] == "unresolved_rsid"
    assert len(records) == 3

    pvar = tmp_path / "ref.chr1.pvar"
    pvar.write_text(
        "##fileformat=PVARv1.0\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\trs1\tA\tG\n"
        "1\t300\trs3\tC\tT\n"
        "1\t400\trs4\tG\tA\n"
        "1\t500\trs5\tA\tC\n",
        encoding="utf-8",
    )
    afreq = tmp_path / "ref.EUR.chr1.afreq"
    afreq.write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\trs1\tA\tG\t0.7\t100\n"
        "1\trs3\tC\tT\t0.3\t100\n"
        "1\trs4\tG\tA\t0.5\t100\n"
        "1\trs5\tA\tC\t0.8\t100\n",
        encoding="utf-8",
    )
    genopred, exclusions, gp_counts = match_genopred_reference(
        requested,
        str(tmp_path / "ref.chr{chromosome}.pvar"),
        str(tmp_path / "ref.EUR.chr{chromosome}.afreq"),
        chromosomes=[1],
    )
    assert len(genopred) == 4
    assert exclusions["rs2"] == "absent_from_genopred_reference"
    assert gp_counts["matched"] == 4

    frequencies = {
        "rs1": GenotypeFrequency("rs1", 1, 2, 7, 0, 0.8),
        "rs2": GenotypeFrequency("rs2", 5, 0, 5, 0, 0.5),
        "rs3": GenotypeFrequency("rs3", 4, 4, 2, 0, 0.4),
        "rs4": GenotypeFrequency("rs4", 3, 4, 3, 0, 0.5),
        "rs5": GenotypeFrequency("rs5", 1, 2, 7, 0, 0.8),
    }
    retained, ties, orientation_counts = orient_markers(
        genopred, frequencies, classification
    )
    assert len(retained) == 3
    assert len(ties) == 1
    assert orientation_counts["retain_sample_and_eur_major_agree"] == 2
    rs3 = next(marker for marker in retained if marker.SNP == "rs3")
    assert rs3.CODE2_SIDE == "REF"
    assert rs3.ORIENTATION_DECISION == "retain_eur_major_sample_flip_near_half"
    assert pearson_code2_vs_eur_major(retained + ties) > 0.5

    snpinfo = tmp_path / "snpinfo"
    snpinfo.write_text(
        "CHR SNP BP A1 A2 MAF\n"
        "1 rs1 100 A G 0.3\n"
        "1 rs3 300 T C 0.3\n"
        "1 other_id 500 A C 0.2\n",
        encoding="utf-8",
    )
    compatible, excluded, prscs_counts = intersect_prscs(retained, snpinfo)
    assert len(compatible) == 3
    assert not excluded
    assert prscs_counts["same_order"] == 2
    assert prscs_counts["swapped_order"] == 1
    assert prscs_counts["coordinate_allele_recovered"] == 1


def test_reconstruct_target_cli_end_to_end(tmp_path: Path) -> None:
    import subprocess
    import sys
    import yaml

    root = Path(__file__).resolve().parents[1]
    (tmp_path / "ref").mkdir()
    (tmp_path / "variants.txt").write_text("rs1\nrs2\nrs3\nrs4\nrs5\n", encoding="utf-8")
    write_gzip(
        tmp_path / "frequencies.tsv.gz",
        "SNP\tN_CODE_0\tN_CODE_1\tN_CODE_2\tN_MISSING\tGENDEP_CODE2_FREQ\n"
        "rs1\t1\t2\t7\t0\t0.8\n"
        "rs2\t5\t0\t5\t0\t0.5\n"
        "rs3\t4\t4\t2\t0\t0.4\n"
        "rs4\t3\t4\t3\t0\t0.5\n"
        "rs5\t1\t2\t7\t0\t0.8\n",
    )
    write_gzip(
        tmp_path / "ref" / "ensembl.chr1.vcf.gz",
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "1\t100\trs1\tA\tG\t.\tPASS\t.\n"
        "1\t200\trs2\tA\tAT\t.\tPASS\t.\n"
        "1\t400\trs4\tG\tA\t.\tPASS\t.\n"
        "1\t500\trs5\tA\tC\t.\tPASS\t.\n",
    )
    (tmp_path / "ref" / "ref.chr1.pvar").write_text(
        "##fileformat=PVARv1.0\n"
        "#CHROM\tPOS\tID\tREF\tALT\n"
        "1\t100\trs1\tA\tG\n"
        "1\t300\trs3\tC\tT\n"
        "1\t400\trs4\tG\tA\n"
        "1\t500\trs5\tA\tC\n",
        encoding="utf-8",
    )
    (tmp_path / "ref" / "ref.EUR.chr1.afreq").write_text(
        "#CHROM\tID\tREF\tALT\tALT_FREQS\tOBS_CT\n"
        "1\trs1\tA\tG\t0.7\t100\n"
        "1\trs3\tC\tT\t0.3\t100\n"
        "1\trs4\tG\tA\t0.5\t100\n"
        "1\trs5\tA\tC\t0.8\t100\n",
        encoding="utf-8",
    )
    (tmp_path / "ref" / "snpinfo").write_text(
        "CHR SNP BP A1 A2 MAF\n"
        "1 rs1 100 A G 0.3\n"
        "1 rs3 300 T C 0.3\n"
        "1 other 500 A C 0.2\n",
        encoding="utf-8",
    )
    config = {
        "project_root": str(tmp_path),
        "controlled_data": {
            "genotype_variant_list": "variants.txt",
            "genotype_frequency_table": "frequencies.tsv.gz",
        },
        "external_resources": {
            "ensembl_vcf_pattern": "ref/ensembl.chr{chromosome}.vcf.gz",
            "genopred_pvar_pattern": "ref/ref.chr{chromosome}.pvar",
            "genopred_eur_afreq_pattern": "ref/ref.EUR.chr{chromosome}.afreq",
            "prscs_snpinfo": "ref/snpinfo",
        },
        "genotype": {"chromosomes": [1]},
        "outputs": {"genotype_reconstruction": "out"},
    }
    config_path = tmp_path / "config.yml"
    config_path.write_text(yaml.safe_dump(config), encoding="utf-8")

    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "gendep.py"),
                "reconstruct-target",
            "--config",
            str(config_path),
            "--expected",
            "",
        ],
        cwd=root,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, (
        f"reconstruct-target failed with exit code {completed.returncode}\n"
        f"STDOUT:\n{completed.stdout or '<empty>'}\n"
        f"STDERR:\n{completed.stderr or '<empty>'}"
    )
    assert "Target reconstruction completed successfully" in completed.stdout
    summary = (tmp_path / "out" / "reconstruction_summary.tsv").read_text(encoding="utf-8")
    assert "final_oriented_target\t3" in summary
    assert "prscs_compatible\t3" in summary

# ---- test_target_validation.py ----
from pathlib import Path

from gendep.genotype.validation import compare_vcfs


VCF_HEADER = (
    "##fileformat=VCFv4.2\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n"
)


def test_compare_vcfs_exact_and_mismatch(tmp_path: Path) -> None:
    input_vcf = tmp_path / "input.vcf"
    roundtrip_vcf = tmp_path / "roundtrip.vcf"
    body = "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t1/1\n"
    input_vcf.write_text(VCF_HEADER + body, encoding="utf-8")
    roundtrip_vcf.write_text(VCF_HEADER + body, encoding="utf-8")

    summary, examples = compare_vcfs(input_vcf, roundtrip_vcf, ["S1", "S2"])
    assert summary.variants == 1
    assert summary.genotype_comparisons == 2
    assert summary.genotype_mismatches == 0
    assert not examples

    roundtrip_vcf.write_text(
        VCF_HEADER + "1\t100\trs1\tA\tG\t.\tPASS\t.\tGT\t0/0\t0/1\n",
        encoding="utf-8",
    )
    summary, examples = compare_vcfs(input_vcf, roundtrip_vcf, ["S1", "S2"])
    assert summary.genotype_mismatches == 1
    assert examples[0]["sample_id"] == "S2"

# ---- test_prs_stage.py ----

import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from gendep.prs.genopred import prepare_run
from gendep.prs.gwas import TraitSpec, standardise_trait
from gendep.prs.profiles import (
    count_score_weights,
    participant_keys,
    read_table,
    validate_ancestry,
    validate_profile,
    validate_projected_pcs,
)


def test_standardise_plain_file_with_gz_suffix_and_deterministic_output(tmp_path: Path) -> None:
    source = tmp_path / "SWB_Full.txt.gz"
    source.write_text(
        "CHR\tPOS\tMarkerName\tA1\tA2\tBeta\tSE\tPval\tEAF\n"
        "1\t100\trs1\tA\tG\t0.1\t0.01\t0.02\t0.4\n"
        "X\t200\trsX\tC\tT\t0.2\t0.02\t0.03\t0.2\n",
        encoding="utf-8",
    )
    spec = TraitSpec(
        trait="SWB",
        filename=source.name,
        transform="swb",
        delimiter="tab",
        expected_source_rows=2,
        expected_autosomal_rows=1,
        sample_size_strategy="fixed",
        fixed_n=298420,
    )
    first = tmp_path / "first.tsv.gz"
    second = tmp_path / "second.tsv.gz"
    record, header = standardise_trait(spec, source, first, strict_counts=True)
    standardise_trait(spec, source, second, strict_counts=True)
    assert record.source_compression == "plain"
    assert record.autosomal_rows == 1
    assert record.non_autosomal_rows == 1
    assert header[0] == "CHR"
    assert first.read_bytes() == second.read_bytes()
    with gzip.open(first, "rt", encoding="utf-8") as handle:
        text = handle.read()
    assert "N\tFREQ" in text
    assert "298420\t0.4" in text


def test_mdd_weighted_frequency_and_sample_size(tmp_path: Path) -> None:
    source = tmp_path / "mdd.tsv.gz"
    with gzip.open(source, "wt", encoding="utf-8") as handle:
        handle.write(
            "#CHROM\tPOS\tID\tEA\tNEA\tBETA\tSE\tPVAL\tNEFF\tFCAS\tFCON\tNCAS\tNCON\tIMPINFO\n"
            "1\t100\trs1\tA\tG\t0.1\t0.01\t0.02\t1000\t0.2\t0.4\t100\t300\t0.9\n"
        )
    spec = TraitSpec("MDD", source.name, "mdd", "tab", 1, 1, "variant N")
    output = tmp_path / "mdd.standardised.tsv.gz"
    standardise_trait(spec, source, output)
    with gzip.open(output, "rt", encoding="utf-8") as handle:
        frame = pd.read_csv(handle, sep="\t")
    assert frame.loc[0, "FREQ"] == 0.35
    assert frame.loc[0, "N"] == 1000


def test_prepare_genopred_bundle(tmp_path: Path) -> None:
    root = tmp_path
    standardised = root / "work/prs/gwas_standardised"
    target = root / "work/genotype/target_validation/pfiles"
    reference = root / "reference"
    resources = root / "resources"
    genopred = root / "GenoPred"
    tools = root / "tools"
    for directory in (standardised, target, reference, resources, genopred / "pipeline", tools):
        directory.mkdir(parents=True, exist_ok=True)
    traits = ["MDD", "ANX", "BIP", "SCZ", "NEUR", "INSOM", "SWB", "EA"]
    for trait in traits:
        (standardised / f"{trait}.standardised.tsv.gz").write_bytes(b"x")
    for chromosome in range(1, 23):
        for extension in ("pgen", "pvar", "psam"):
            (target / f"GENDEP.chr{chromosome}.{extension}").write_text("x", encoding="utf-8")
    for executable in ("snakemake", "Rscript", "python", "plink", "plink2"):
        path = tools / executable
        path.write_text("#!/bin/sh\n", encoding="utf-8")
    base = root / "base.yml"
    base.write_text("cores_prep_pgs: 4\n", encoding="utf-8")

    repository = Path(__file__).resolve().parents[1]
    project = {
        "project_root": str(root),
        "external_resources": {
            "genopred_root": str(genopred),
            "genopred_base_config": str(base),
            "genopred_snakemake": str(tools / "snakemake"),
            "genopred_runtime_rscript": str(tools / "Rscript"),
            "genopred_runtime_python": str(tools / "python"),
            "genopred_runtime_plink": str(tools / "plink"),
            "genopred_runtime_plink2": str(tools / "plink2"),
            "genopred_custom_reference": str(reference),
            "genopred_resource_directory": str(resources),
        },
        "outputs": {
            "gwas_standardised": "work/prs/gwas_standardised",
            "genopred_run": "work/prs/genopred",
        },
        "prs": {"target_prefix": "work/genotype/target_validation/pfiles/GENDEP"},
    }
    project_path = root / "project.yml"
    project_path.write_text(yaml.safe_dump(project), encoding="utf-8")
    environment = prepare_run(
        project_path,
        repository / "config/genopred.yml",
        repository / "config/gwas.yml",
        repository / "config/traits.tsv",
    )
    config = yaml.safe_load(Path(environment["RUN_CONFIG"]).read_text())
    assert config["pgs_methods"] == ["prscs"]
    assert config["prscs_phi"] == ["auto"]
    assert config["testing"] == "NA"
    gwas_list = Path(config["gwas_list"]).read_text(encoding="utf-8").splitlines()
    assert len(gwas_list) == 9
    assert gwas_list[0].startswith("name\tpath\tpopulation")


def test_profile_score_ancestry_and_pc_validation(tmp_path: Path) -> None:
    output = tmp_path / "output"
    profile_path = output / "GENDEP/pgs/TRANS/prscs/MDD/GENDEP-MDD-TRANS.profiles"
    profile_path.parent.mkdir(parents=True)
    profile_path.write_text("FID\tIID\tSCORE\nF1\tI1\t0.1\nF2\tI2\t0.2\nF3\tI3\t0.4\n", encoding="utf-8")
    psam = tmp_path / "test.psam"
    psam.write_text("#FID\tIID\nF1\tI1\nF2\tI2\nF3\tI3\n", encoding="utf-8")
    expected_keys = participant_keys(read_table(psam))
    record, aligned = validate_profile(profile_path, expected_keys, pd.Series([0.1, 0.2, 0.4]))
    assert record["validation"] == "PASS"
    assert np.allclose(aligned, [0.1, 0.2, 0.4])

    score_path = tmp_path / "score.gz"
    with gzip.open(score_path, "wt", encoding="utf-8") as handle:
        handle.write("CHR SNP BP A1 A2 BETA\n1 rs1 1 A G 0\n1 rs2 2 C T 0.2\n")
    assert count_score_weights(score_path)[:2] == (2, 1)

    ancestry = output / "GENDEP/ancestry/GENDEP.Ancestry.model_pred"
    ancestry.parent.mkdir(parents=True)
    ancestry.write_text(
        "FID IID AFR AMR CSA EAS EUR MID\n"
        "F1 I1 0 0 0 0 1 0\nF2 I2 0 0 0 0 1 0\nF3 I3 0 0 0 0 1 0\n",
        encoding="utf-8",
    )
    keep_dir = output / "GENDEP/ancestry/keep_files/model_based"
    keep_dir.mkdir(parents=True)
    (keep_dir / "EUR.keep").write_text("F1 I1\nF2 I2\n", encoding="utf-8")
    ancestry_record, _ = validate_ancestry(output, expected_keys)
    assert ancestry_record["maximum_population_eur"] == 3
    assert ancestry_record["strict_eur_keep_participants"] == 2
    assert ancestry_record["participants_outside_all_model_keeps"] == 1

    pc = output / "GENDEP/pcs/projected/TRANS/test.eigenvec"
    pc.parent.mkdir(parents=True)
    pc.write_text(
        "#FID IID PC1 PC2 PC3 PC4 PC5 PC6\n"
        "F1 I1 1 2 3 4 5 6\nF2 I2 2 3 4 5 6 7\nF3 I3 3 4 5 6 7 8\n",
        encoding="utf-8",
    )
    pc_record = validate_projected_pcs(output, expected_keys)
    assert pc_record["participant_rows"] == 3
    assert pc_record["pc_columns"] == 6


def test_remaining_six_gwas_schemas(tmp_path: Path) -> None:
    cases = {
        "ANX": (
            "anx", "tab",
            "CHR\tBP\tSNP\tA1\tA2\tOR\tSE\tP\tNeff_half\tFRQ_A_122083\tFRQ_U_729602\tNca\tNco\tINFO\n"
            "1\t100\trs1\tA\tG\t1.1\t0.1\t0.02\t50\t0.2\t0.4\t100\t300\t0.9\n",
            {}, "OR", 100.0,
        ),
        "BIP": (
            "bip", "whitespace",
            "CHR BP SNP A1 A2 OR SE P Neff_half HRC_FRQ_A1 INFO\n"
            "1 100 rs1 A G 1.1 0.1 0.02 50 0.3 0.9\n",
            {}, "OR", 100.0,
        ),
        "SCZ": (
            "scz", "tab",
            "CHROM\tPOS\tID\tA1\tA2\tBETA\tSE\tPVAL\tNEFF\tFCAS\tFCON\tNCAS\tNCON\tIMPINFO\n"
            "1\t100\trs1\tA\tG\t0.1\t0.01\t0.02\t1000\t0.2\t0.4\t100\t300\t0.9\n",
            {}, "BETA", 1000.0,
        ),
        "NEUR": (
            "neur", "tab",
            "CHR\tPOS\tRSID\tA1\tA2\tZ\tP\tN\tEAF_UKB\tINFO_UKB\n"
            "1\t100\trs1\tA\tG\t2.1\t0.02\t1000\t0.3\t0.9\n",
            {}, "Z", 1000.0,
        ),
        "INSOM": (
            "insom", "tab",
            "CHR\tBP\tRSID_UKB\tA1\tA2\tOR\tSE\tP\tINFO_UKB\n"
            "1\t100\trs1\tA\tG\t1.1\t0.1\t0.02\t0.9\n",
            {"cases": 100.0, "controls": 300.0}, "OR", 300.0,
        ),
        "EA": (
            "ea", "tab",
            "Chr\tBP\trsID\tEffect_allele\tOther_allele\tBeta\tSE\tP\tEAF_HRC\n"
            "1\t100\trs1\tA\tG\t0.1\t0.01\t0.02\t0.3\n",
            {"fixed_n": 765283.0}, "BETA", 765283.0,
        ),
    }
    for trait, (transform, delimiter, text, extra, effect_column, expected_n) in cases.items():
        source = tmp_path / f"{trait}.txt.gz"
        with gzip.open(source, "wt", encoding="utf-8") as handle:
            handle.write(text)
        spec = TraitSpec(
            trait=trait,
            filename=source.name,
            transform=transform,
            delimiter=delimiter,
            expected_source_rows=1,
            expected_autosomal_rows=1,
            sample_size_strategy="synthetic",
            fixed_n=extra.get("fixed_n"),
            cases=extra.get("cases"),
            controls=extra.get("controls"),
        )
        output = tmp_path / f"{trait}.standardised.tsv.gz"
        standardise_trait(spec, source, output)
        frame = pd.read_csv(output, sep="\t")
        assert len(frame) == 1
        assert effect_column in frame.columns
        assert np.isclose(frame.loc[0, "N"], expected_n)


def test_plink_output_suffix_preserves_chromosome_token() -> None:
    from gendep.commands.validate_target import append_output_suffix

    prefix = Path("work/pfiles/GENDEP.chr1")
    assert append_output_suffix(prefix, ".pgen").name == "GENDEP.chr1.pgen"
    assert append_output_suffix(prefix, ".pvar").name == "GENDEP.chr1.pvar"
    assert append_output_suffix(prefix, ".psam").name == "GENDEP.chr1.psam"
