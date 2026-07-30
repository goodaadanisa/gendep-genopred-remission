"""Command-line interface for the Python stages of the GENDEP workflow."""
from __future__ import annotations

import importlib
import sys

COMMANDS = {
    "audit-clinical": "gendep.commands.audit_clinical",
    "define-predictors": "gendep.commands.define_predictors",
    "build-analysis": "gendep.commands.build_analysis_dataset",
    "reconstruct-target": "gendep.commands.reconstruct_target",
    "validate-target": "gendep.commands.validate_target",
    "prepare-gwas": "gendep.commands.prepare_gwas",
    "resampling": "gendep.commands.generate_resampling",
    "aggregate-models": "gendep.commands.aggregate_models",
    "model-diagnostics": "gendep.commands.model_diagnostics",
    "orientation-sensitivity": "gendep.commands.build_orientation_sensitivity",
    "tables": "gendep.commands.make_tables",
    "figures": "gendep.commands.make_figures",
    "validate": "gendep.commands.validate_release",
}


def usage() -> str:
    """Return concise help listing the available workflow subcommands."""
    names = "\n".join(f"  {name}" for name in COMMANDS)
    return (
        "Usage: gendep <command> [arguments]\n\n"
        "Commands:\n" + names + "\n\n"
        "Run a command with --help for its detailed inputs and outputs."
    )


def main() -> None:
    """Dispatch one named subcommand to its implementation module."""
    if len(sys.argv) < 2 or sys.argv[1] in {"-h", "--help"}:
        print(usage())
        return
    command = sys.argv[1]
    module_name = COMMANDS.get(command)
    if module_name is None:
        raise SystemExit(f"Unknown command: {command}\n\n{usage()}")
    module = importlib.import_module(module_name)
    sys.argv = [f"gendep {command}", *sys.argv[2:]]
    module.main()


if __name__ == "__main__":
    main()
