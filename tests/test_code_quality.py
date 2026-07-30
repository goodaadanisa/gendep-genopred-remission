"""Maintainability and editorial checks for the public scientific codebase."""
from __future__ import annotations

import ast
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
CORE = ROOT / "src" / "gendep"
PLACEHOLDER_DOCSTRING = re.compile(r"Implement the .* workflow operation\.")
NONEXISTENT_SCRIPT = re.compile(r"scripts/(?!gendep\.py\b)[A-Za-z0-9_./-]+\.py")


def test_long_public_functions_have_docstrings() -> None:
    """Require documentation for public functions large enough to need explanation."""
    missing: list[str] = []
    for path in CORE.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            if node.name.startswith("_"):
                continue
            length = (node.end_lineno or node.lineno) - node.lineno + 1
            if length >= 20 and not ast.get_docstring(node):
                missing.append(f"{path.relative_to(ROOT)}:{node.lineno}:{node.name}")
    assert not missing, "Undocumented long public functions:\n" + "\n".join(missing)


def test_python_modules_have_module_docstrings() -> None:
    """Require every installable Python module to state its responsibility."""
    missing = []
    for path in CORE.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        if not ast.get_docstring(tree):
            missing.append(str(path.relative_to(ROOT)))
    assert not missing, "Missing module docstrings: " + ", ".join(missing)


def test_placeholder_docstrings_are_not_committed() -> None:
    """Reject mechanically generated docstrings that do not explain a contract."""
    matches: list[str] = []
    for path in CORE.rglob("*.py"):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in ast.walk(tree):
            if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                continue
            docstring = ast.get_docstring(node)
            if docstring and PLACEHOLDER_DOCSTRING.fullmatch(docstring):
                matches.append(f"{path.relative_to(ROOT)}:{node.lineno}:{node.name}")
    assert not matches, "Placeholder docstrings:\n" + "\n".join(matches)


def test_non_init_modules_have_no_obvious_unused_imports() -> None:
    """Catch straightforward unused imports without requiring a third-party linter."""
    unused: list[str] = []
    for path in CORE.rglob("*.py"):
        if path.name == "__init__.py":
            continue
        tree = ast.parse(path.read_text(encoding="utf-8"))
        used_names = {node.id for node in ast.walk(tree) if isinstance(node, ast.Name)}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported = [(alias.asname or alias.name.split(".")[0], node.lineno) for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                imported = [
                    (alias.asname or alias.name, node.lineno)
                    for alias in node.names
                    if alias.name != "*"
                ]
            else:
                continue
            for name, line in imported:
                if name != "annotations" and name not in used_names:
                    unused.append(f"{path.relative_to(ROOT)}:{line}:{name}")
    assert not unused, "Obvious unused imports:\n" + "\n".join(unused)


def test_documentation_does_not_reference_nonexistent_python_scripts() -> None:
    """Require documentation to use the installed CLI rather than invented wrappers."""
    matches: list[str] = []
    paths = [ROOT / "README.md", ROOT / "DATA_AVAILABILITY.md"]
    paths.extend((ROOT / "docs").rglob("*.md"))
    paths.extend(CORE.rglob("*.py"))
    for pattern in ("*.py", "*.R", "*.sh"):
        paths.extend((ROOT / "scripts").rglob(pattern))
    for path in paths:
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
            if NONEXISTENT_SCRIPT.search(line):
                matches.append(f"{path.relative_to(ROOT)}:{line_number}:{line.strip()}")
    assert not matches, "Nonexistent script examples:\n" + "\n".join(matches)


def test_github_workflow_rejects_duplicate_yaml_keys() -> None:
    """Parse the CI workflow with a loader that rejects duplicate mapping keys."""

    class UniqueKeyLoader(yaml.SafeLoader):
        pass

    def construct_mapping(loader, node, deep=False):
        mapping = {}
        for key_node, value_node in node.value:
            key = loader.construct_object(key_node, deep=deep)
            if key in mapping:
                raise AssertionError(f"Duplicate YAML key: {key}")
            mapping[key] = loader.construct_object(value_node, deep=deep)
        return mapping

    UniqueKeyLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping,
    )
    workflow = ROOT / ".github" / "workflows" / "tests.yml"
    payload = yaml.load(workflow.read_text(encoding="utf-8"), Loader=UniqueKeyLoader)
    assert payload["jobs"]["syntax"]["steps"]


def test_runtime_configuration_has_one_source_per_policy() -> None:
    """Keep runtime configuration distinct from generated policy snapshots."""
    removed_duplicates = {
        "clinical.yml",
        "clinical_schema.tsv",
        "predictors.yml",
        "secondary_analyses.yml",
    }
    present = {path.name for path in (ROOT / "config").iterdir() if path.is_file()}
    assert removed_duplicates.isdisjoint(present)
    resampling = yaml.safe_load((ROOT / "config" / "analyses.yml").read_text(encoding="utf-8"))
    assert set(resampling) == {
        "master_seed",
        "outer_folds",
        "outer_repeats",
        "inner_folds",
        "alpha_grid",
    }
