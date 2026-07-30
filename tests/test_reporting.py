"""Tests for the release-safe reporting inputs and outputs."""
import csv
import gzip
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]


def test_final_png_set_is_complete() -> None:
    """The repository must contain all five main and eight supplementary figures."""
    pngs = sorted((ROOT / "results" / "figures").glob("*.png"))
    assert len(pngs) == 13


def test_figure_inputs_are_inspectable() -> None:
    """All figure inputs should be visible files rather than an opaque nested archive."""
    data_files = sorted((ROOT / "results" / "figure_data").glob("*.tsv*"))
    assert len(data_files) == 27
    assert not (ROOT / "results" / "release_data.zip").exists()
    redundant = {
        "figure_02_source_runs.tsv",
        "figure_04_base_robustness_increments.tsv",
        "figure_04_orientation_sensitivity_closure_results.tsv",
    }
    assert redundant.isdisjoint({path.name for path in data_files})


def test_figure_source_tables_are_rectangular() -> None:
    """Every release-safe TSV must retain the complete declared column schema."""
    for path in sorted((ROOT / "results" / "figure_data").glob("*.tsv*")):
        opener = gzip.open if path.suffix == ".gz" else open

        with opener(
            path,
            "rt",
            encoding="utf-8",
            newline="",
        ) as handle:
            rows = list(csv.reader(handle, delimiter="\t"))

        assert rows, f"Empty figure source: {path.name}"
        expected_width = len(rows[0])
        assert expected_width > 0
        assert all(len(row) == expected_width for row in rows), (
            f"Ragged figure source: {path.name}"
        )


def test_table_source_contains_complete_release() -> None:
    """The table source must contain five main and 15 supplementary tables."""
    payload = yaml.safe_load((ROOT / "results" / "table_data.yml").read_text(encoding="utf-8"))
    assert len(payload["tables"]) == 20
    assert (ROOT / "results" / "tables.md").is_file()


def test_release_scan_ignores_local_virtual_environment(tmp_path, monkeypatch) -> None:
    """Local virtual environments must not trigger machine-path privacy failures."""
    from gendep.commands import validate_release

    local_environment = tmp_path / ".venv" / "metadata"
    local_environment.mkdir(parents=True)
    (local_environment / "install.txt").write_text(
        "C:" + r"\Users\example\project", encoding="utf-8"
    )
    public_file = tmp_path / "README.md"
    public_file.write_text("release-safe text", encoding="utf-8")

    monkeypatch.setattr(validate_release, "ROOT", tmp_path)
    files = validate_release.repository_files()
    assert files == [public_file]
    assert validate_release.scan_text_files(files) == []


def test_table_markdown_is_human_readable_and_source_remains_exact() -> None:
    """Render rounded labels for readers while retaining unrounded YAML values."""
    markdown = (ROOT / "results" / "tables.md").read_text(encoding="utf-8")
    source = (ROOT / "results" / "table_data.yml").read_text(encoding="utf-8")
    assert "Adjusted OR per SD" in markdown
    assert "adjusted_odds_ratio_per_SD" not in markdown
    assert "0.934505548835992" in source
    assert "0.934505548835992" not in markdown
    assert "Displayed values are rounded for readability" in markdown


def test_citation_metadata_matches_release_and_custom_licence() -> None:
    """Keep GitHub citation metadata aligned with the package release."""
    citation = yaml.safe_load((ROOT / "CITATION.cff").read_text(encoding="utf-8"))
    assert citation["cff-version"] == "1.2.0"
    assert citation["version"] == "1.0.0-rc6"
    assert citation["repository-code"].endswith("gendep-genopred-remission")
    assert citation["license-url"].endswith("LICENSE.md")
    assert citation["preferred-citation"]["type"] == "thesis"
    assert citation["preferred-citation"]["thesis-type"] == "MSc dissertation"
    assert citation["preferred-citation"]["institution"]["name"] == "King's College London"
    assert "institution" not in {key for key in citation if key != "preferred-citation"}
    assert "license" not in citation
