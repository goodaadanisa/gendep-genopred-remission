"""Render polished dissertation tables from the machine-readable release source."""
from __future__ import annotations

import argparse
import math
import re
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[3]
DEFAULT_DATA = ROOT / "results" / "table_data.yml"
DEFAULT_OUTPUT = ROOT / "results" / "tables.md"

COLUMN_LABELS = {
    "overall_n430": "Overall (n=430)",
    "treatment_group_1_n210": "Treatment group 1 (n=210)",
    "treatment_group_2_n220": "Treatment group 2 (n=220)",
    "stage_or_validation": "Stage or validation",
    "reference_or_denominator": "Reference or denominator",
    "analysis_label": "Analysis",
    "participants": "Participants",
    "comparator": "Comparator",
    "prs_model": "Comparator + PRS",
    "increment": "Increment",
    "lower_95": "95% CI lower",
    "upper_95": "95% CI upper",
    "adjusted_odds_ratio_per_SD": "Adjusted OR per SD",
    "or_lower_95": "OR 95% CI lower",
    "or_upper_95": "OR 95% CI upper",
    "wald_p_value": "Wald p-value",
    "bh_fdr_p_value": "BH-FDR p-value",
    "auc_increment": "AUC increment",
    "auc_lower_95": "AUC 95% CI lower",
    "auc_upper_95": "AUC 95% CI upper",
    "auc_positive_repetitions": "Positive AUC repeats (of 10)",
    "ratio_of_odds_ratios": "Ratio of odds ratios",
    "ratio_lower_95": "Ratio 95% CI lower",
    "ratio_upper_95": "Ratio 95% CI upper",
    "group_1_or": "Group 1 OR",
    "group_1_lower_95": "Group 1 OR 95% CI lower",
    "group_1_upper_95": "Group 1 OR 95% CI upper",
    "group_2_or": "Group 2 OR",
    "group_2_lower_95": "Group 2 OR 95% CI lower",
    "group_2_upper_95": "Group 2 OR 95% CI upper",
    "brier_improvement": "Brier improvement",
    "log_loss_improvement": "Log-loss improvement",
    "included_in_46_variable_summary_set": "In 46-variable summary set",
    "expected_storage_class": "Expected storage class",
    "expected_final_disposition": "Final disposition",
    "primary_model_role": "Primary model role",
    "penalty_factor": "Penalty factor",
    "roundtrip_pass": "Round-trip pass",
    "retention_percentage": "Retention (%)",
    "exact_source_filename": "Source file",
    "sample_size_strategy": "Sample-size strategy",
    "outcome_used_in_standardisation": "Outcome used in standardisation",
    "range_or_integrity": "Range or integrity",
    "fixed_design_or_seed": "Fixed design or seed",
    "pearson_correlation": "Pearson correlation",
    "spearman_correlation": "Spearman correlation",
    "mean_absolute_difference": "Mean absolute difference",
    "maximum_absolute_difference": "Maximum absolute difference",
    "participants_with_changed_score": "Participants with changed score",
    "accepted_nonzero_weights": "Primary non-zero weights",
    "nonzero_posterior_weights": "Non-zero posterior weights",
    "model_or_contrast": "Model or contrast",
    "calibration_measure": "Calibration measure",
    "ideal_value": "Ideal value",
    "nonzero_all_fits": "Non-zero (all fits)",
    "nonzero_alpha_gt_0_fits": "Non-zero (alpha > 0 fits)",
    "predominant_sign": "Predominant sign",
    "sign_consistency_percent": "Sign consistency (%)",
    "median_nonzero_coefficient": "Median non-zero coefficient",
    "q1_nonzero_coefficient": "Q1 non-zero coefficient",
    "q3_nonzero_coefficient": "Q3 non-zero coefficient",
}

HIDDEN_COLUMNS = {
    "table_01_baseline_characteristics": {"summary_type"},
    "table_03_model_performance": {"analysis_key"},
    "supplementary_table_s08_complete_model_performance": {"analysis_key"},
}

INTEGER_COLUMNS = {
    "position",
    "participants",
    "matched_variants",
    "retained_variants",
    "excluded_eur_frequency_ties",
    "genotype_comparisons",
    "genotype_mismatches",
    "input_rows",
    "autosomal_rows",
    "n",
    "auc_positive_repetitions",
    "accepted_nonzero_weights",
    "excluded_nonzero_weights",
    "participants_with_changed_score",
    "source_rows",
    "cleaned_harmonised_gwas_rows",
    "target_score_markers",
    "zero_weight_markers",
    "nonzero_posterior_weights",
}

CATEGORICAL_COLUMNS = {
    "expected_storage_class",
    "semantic_type",
    "primary_model_role",
    "sample_size_strategy",
    "analysis_role",
    "expected_final_disposition",
}


def parse_args() -> argparse.Namespace:
    """Parse table-source and Markdown-output paths."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def escape(value: object) -> str:
    """Escape a value for safe use in a GitHub Markdown table."""
    return str(value).replace("|", "\\|").replace("\n", " ")


def humanise_identifier(value: str) -> str:
    """Convert a machine-readable label to restrained human-readable text."""
    text = value.replace("_", " ")
    replacements = {
        "prs": "PRS",
        "gwas": "GWAS",
        "auc": "AUC",
        "fdr": "FDR",
        "or": "OR",
        "sd": "SD",
        "eur": "EUR",
        "pc": "PC",
        "q1": "Q1",
        "q3": "Q3",
        "bh": "BH",
    }
    words = [replacements.get(word.lower(), word) for word in text.split()]
    rendered = " ".join(words)
    return rendered[:1].upper() + rendered[1:]


def column_label(column: str) -> str:
    """Return a publication-style label for one source column."""
    return COLUMN_LABELS.get(column, humanise_identifier(column))


def parse_number(value: object) -> float | None:
    """Parse a finite scalar number without altering composite text values."""
    if isinstance(value, bool) or value is None:
        return None
    text = str(value).strip()
    if not re.fullmatch(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?", text):
        return None
    number = float(text)
    return number if math.isfinite(number) else None


def format_p_value(number: float) -> str:
    """Format a p-value to three significant figures."""
    if 0 < number < 0.001:
        return "<0.001"
    return f"{number:.3f}"


def format_cell(column: str, value: object) -> str:
    """Format one table cell while preserving exact source data in YAML."""
    text = str(value).strip()
    if text.upper() == "TRUE":
        return "Yes"
    if text.upper() == "FALSE":
        return "No"
    if text.lower() == "not_applicable":
        return "Not applicable"
    if column in CATEGORICAL_COLUMNS and "_" in text:
        return humanise_identifier(text).lower()

    number = parse_number(value)
    if number is None:
        return text
    if column.endswith("p_value"):
        return format_p_value(number)
    if column in INTEGER_COLUMNS or re.fullmatch(r"chromosome", column):
        return f"{int(round(number)):,}"
    if "percent" in column or column == "retention_percentage":
        return f"{number:.1f}"
    if column in {"comparator", "prs_model", "increment", "lower_95", "upper_95"} or any(
        token in column for token in ("auc", "brier", "log_loss", "calibration")
    ):
        rounded = 0.0 if abs(number) < 0.00005 else number
        return f"{rounded:.4f}"
    if any(
        token in column
        for token in (
            "odds_ratio",
            "_or",
            "correlation",
            "coefficient",
            "difference",
            "minimum",
            "maximum",
            "median",
            "mean",
            "q1",
            "q3",
            "sd",
            "estimate",
            "lower_95",
            "upper_95",
        )
    ):
        return f"{number:.3f}"
    if number.is_integer():
        return f"{int(number):,}"
    return f"{number:.3f}"


def render_table(table: dict[str, object]) -> list[str]:
    """Render one structured table record as polished Markdown lines."""
    table_id = str(table["id"])
    title = str(table["title"])
    hidden = HIDDEN_COLUMNS.get(table_id, set())
    source_columns = [str(value) for value in table.get("columns", [])]
    visible_indices = [index for index, column in enumerate(source_columns) if column not in hidden]
    columns = [source_columns[index] for index in visible_indices]
    rows = table.get("rows", [])
    lines = [f"### {title}", ""]
    if columns:
        lines.append("| " + " | ".join(escape(column_label(value)) for value in columns) + " |")
        lines.append("| " + " | ".join("---" for _ in columns) + " |")
        for row in rows:
            cells = [format_cell(source_columns[index], row[index]) for index in visible_indices]
            lines.append("| " + " | ".join(escape(value) for value in cells) + " |")
    lines.append("")
    return lines


def main() -> None:
    """Write all configured main and supplementary tables to one Markdown file."""
    args = parse_args()
    if not args.data.is_file():
        raise SystemExit(f"Missing table source: {args.data}")
    payload = yaml.safe_load(args.data.read_text(encoding="utf-8")) or {}
    tables = payload.get("tables", [])
    if len(tables) != 20:
        raise SystemExit(f"Expected 20 configured tables, found {len(tables)}")

    lines = [
        "# Dissertation tables",
        "",
        "Five main and fifteen supplementary tables from the final dissertation. "
        "Displayed values are rounded for readability; exact unrounded values remain in "
        "`table_data.yml`.",
        "",
        "## Main tables",
        "",
    ]
    for index, table in enumerate(tables):
        if index == 5:
            lines.extend(["## Supplementary tables", ""])
        lines.extend(render_table(table))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines), encoding="utf-8")
    print(f"TABLE_REPORT=PASS tables={len(tables)} output={args.output}")


if __name__ == "__main__":
    main()
