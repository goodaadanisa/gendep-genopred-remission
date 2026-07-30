"""Validate the privacy, structure and release-safe outputs of the repository."""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
FORBIDDEN_SUFFIXES = {".rdata", ".rds", ".xlsx", ".docx", ".vcf", ".pgen", ".pvar", ".psam"}
PRIVATE_PATTERNS = [
    re.compile(r"/scratch/users/", re.I),
    re.compile(r"C:\\Users\\", re.I),
    re.compile(r"BEGIN (?:OPENSSH|RSA|EC|DSA) PRIVATE KEY"),
]
TEXT_SUFFIXES = {".py", ".R", ".r", ".sh", ".md", ".yml", ".yaml", ".tsv", ".toml", ".cff", ".txt"}


def parse_args() -> argparse.Namespace:
    """Parse release-validation options."""
    return argparse.ArgumentParser(description=__doc__).parse_args()


def repository_files() -> list[Path]:
    """Return repository files while excluding local Git and cache metadata."""
    excluded = {
        ".git",
        ".venv",
        "venv",
        "__pycache__",
        ".pytest_cache",
        ".mypy_cache",
        ".ruff_cache",
        ".tox",
        ".nox",
        "build",
        "dist",
    }
    return [path for path in ROOT.rglob("*") if path.is_file() and not excluded.intersection(path.parts)]


def scan_text_files(paths: list[Path]) -> list[str]:
    """Scan inspectable text files for machine-specific paths or private keys."""
    matches: list[str] = []
    for path in paths:
        if path.suffix not in TEXT_SUFFIXES and path.name not in {"README.md", "DATA_AVAILABILITY.md"}:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except (OSError, UnicodeDecodeError):
            continue
        for pattern in PRIVATE_PATTERNS:
            if pattern.search(text):
                matches.append(f"{path.relative_to(ROOT)}:{pattern.pattern}")
    return matches


def main() -> None:
    """Run release checks and stop at the first failed invariant."""
    parse_args()

    def check(name: str, condition: bool, detail: str) -> None:
        """Stop on a failed release invariant and report successful checks."""
        if not condition:
            raise SystemExit(f"{name}=FAIL: {detail}")
        print(f"{name}=PASS {detail}")

    scripts = sorted(path for path in (ROOT / "scripts").iterdir() if path.is_file())
    figure_inputs = sorted((ROOT / "results" / "figure_data").glob("*.tsv*"))
    pngs = sorted((ROOT / "results" / "figures").glob("*.png"))
    table_data = ROOT / "results" / "table_data.yml"

    check("public_script_count", len(scripts) == 8, f"observed={len(scripts)} expected=8")
    check("figure_input_count", len(figure_inputs) == 27, f"observed={len(figure_inputs)} expected=27")
    check("figure_count", len(pngs) == 13, f"observed={len(pngs)} expected=13")
    check("release_zip_removed", not (ROOT / "results" / "release_data.zip").exists(), "no opaque release archive")
    check("table_data", table_data.is_file(), str(table_data.relative_to(ROOT)))
    payload = yaml.safe_load(table_data.read_text(encoding="utf-8")) or {}
    check("table_count", len(payload.get("tables", [])) == 20, "expected=20")

    tracked = repository_files()
    scan_targets = [path for path in tracked if path.resolve() != Path(__file__).resolve()]
    forbidden = [str(path.relative_to(ROOT)) for path in tracked if path.suffix.lower() in FORBIDDEN_SUFFIXES]
    check("controlled_file_extensions", not forbidden, str(forbidden))
    check("private_path_or_key_scan", not scan_text_files(scan_targets), "inspectable text files scanned")
    print(f"RELEASE_VALIDATION=PASS files={len(tracked)}")


if __name__ == "__main__":
    main()
