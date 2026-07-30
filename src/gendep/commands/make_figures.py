"""Regenerate 12 dissertation figures from release-safe aggregate inputs.

Supplementary Figure S4 remains as a committed PNG because its participant-level
plotting input is controlled. The command never deletes input or provenance
directories.
"""
from __future__ import annotations

import argparse
import importlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
DATA_DIR = ROOT / "results" / "figure_data"
OUTPUT_DIR = ROOT / "results" / "figures"
FIGURES = {
    "workflow": "gendep.figures.workflow",
    "target-reconstruction": "gendep.figures.target_reconstruction",
    "primary-performance": "gendep.figures.primary_performance",
    "robustness": "gendep.figures.robustness",
    "individual-prs": "gendep.figures.individual_prs",
    "orientation-support": "gendep.figures.orientation_support",
    "chromosome-retention": "gendep.figures.chromosome_retention",
    "prs-distributions": "gendep.figures.prs_distributions",
    "ancestry": "gendep.figures.ancestry",
    "repeat-increments": "gendep.figures.repeat_increments",
    "alpha-selection": "gendep.figures.alpha_selection",
    "individual-prs-increments": "gendep.figures.individual_prs_increments",
    "coefficient-stability": "gendep.figures.coefficient_stability",
}

PUBLIC_DEFAULTS = [key for key in FIGURES if key != "ancestry"]


def parse_args() -> argparse.Namespace:
    """Parse figure-selection and output-format options."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="+", choices=tuple(FIGURES))
    parser.add_argument(
        "--all-formats",
        action="store_true",
        help="Retain SVG and PDF outputs as well as PNG.",
    )
    return parser.parse_args()


def main() -> None:
    """Generate selected figures and verify the complete PNG release set."""
    args = parse_args()
    if not DATA_DIR.is_dir():
        raise SystemExit(f"Missing figure input directory: {DATA_DIR}")

    selected = args.only or PUBLIC_DEFAULTS
    for key in selected:
        module = importlib.import_module(FIGURES[key])
        original_argv = sys.argv[:]
        try:
            # Each plotting module exposes a small command-line main().
            sys.argv = [key]
            module.main()
        except Exception as exc:
            raise RuntimeError(f"Figure generation failed for '{key}'") from exc
        finally:
            sys.argv = original_argv

    if not args.all_formats:
        for path in OUTPUT_DIR.iterdir():
            if path.suffix.lower() in {".svg", ".pdf"}:
                path.unlink()

    pngs = list(OUTPUT_DIR.glob("*.png"))
    if args.only is None:
        committed_ancestry = OUTPUT_DIR / "Supplementary_Figure_S4_ancestry_probability.png"
        if not committed_ancestry.is_file():
            raise SystemExit("The controlled-data ancestry figure is missing from the committed release")
        if len(pngs) != 13:
            raise SystemExit(f"Expected 13 committed PNG figures, found {len(pngs)}")
        print("FIGURE_REPORT=PASS generated=12 controlled_input_required=1 png_total=13")
    else:
        print(f"FIGURE_REPORT=PASS generated={len(selected)} png_total={len(pngs)}")


if __name__ == "__main__":
    main()
