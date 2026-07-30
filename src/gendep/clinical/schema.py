"""Fixed, outcome-blind clinical schema and predictor policy.

The constants in this module are derived from the supplied workbook
inventory and the fixed methods specification. They do not contain participant-level
values. The policy retains source order, removes only seven later exact
duplicates, and keeps outcome and identifiers outside every predictor set.
"""

from __future__ import annotations

import re

IDENTIFIER_COLUMNS = ("subjectid", "Row.names", "bloodsampleid.x")
TREATMENT_COLUMN = "drug"
OUTCOME_COLUMN = "hdremit.all"
PRS_COLUMNS = (
    "PRS_MDD", "PRS_ANX", "PRS_BIP", "PRS_SCZ",
    "PRS_NEUR", "PRS_INSOM", "PRS_SWB", "PRS_EA",
)
PC_COLUMNS = tuple(f"ancestry_PC{index}" for index in range(1, 7))
EUR_INDICATOR = "model_based_EUR_keep"
ANCESTRY_POPULATIONS = ("AFR", "AMR", "CSA", "EAS", "EUR", "MID")

CLINICAL_COLUMNS = (
    "subjectid", "Row.names", "bloodsampleid.x", "drug", "educ",
    "marital.rec", "occup.rec", "children.rec", "smoker", "bmi_wk0",
    "age", "ageonset", "sex",
    *(f"madrs{index}wk0" for index in range(1, 11)),
    *(f"hamd{index}wk0" for index in range(1, 18)),
    *(f"bdi{index}wk0" for index in range(1, 19)),
    "bdi19awk0", "bdi19wk0", "bdi20wk0", "bdi21wk0",
    "f1score0", "f2score0", "f3score0",
    "f61score0", "f62score0", "f63score0", "f64score0", "f65score0", "f66score0",
    "mscore", "melanch", "atscore", "atyp", "anxdep", "anxsomdep",
    *(f"k{index}" for index in range(1, 33)),
    "newbleq1wk0", "newbleq2wk0", "newbleq4wk0", "newbleq5wk0",
    "newbleq6wk0", "newbleq7wk0", "newbleq8wk0", "newbleq9wk0",
    "nevents", "eventsbin", "mhssri", "everssri", "mhtca", "evertca",
    "mhdual", "everdual", "mhoth", "everoth", "mhantidep", "everantidep",
    "concssri", "conctca", "concantidep", "mcbenzo", "mczhypn", "mhmirta",
    "madrs.total", "hrsd.total", "bdi.total", "hdremit.all",
)

if len(CLINICAL_COLUMNS) != 139 or len(set(CLINICAL_COLUMNS)) != 139:
    raise RuntimeError("The fixed clinical schema must contain 139 unique columns")

EXACT_DUPLICATE_GROUPS = (
    ("hamd12wk0", "bdi15wk0"),
    ("hamd13wk0", "bdi16wk0"),
    ("hamd14wk0", "bdi17wk0"),
    ("k11", "k28"),
    ("k13", "k29", "k31"),
    ("k15", "k32"),
)
REMOVED_EXACT_DUPLICATES = (
    "bdi15wk0", "bdi16wk0", "bdi17wk0", "k28", "k29", "k31", "k32",
)

BASELINE_CLINICAL_CANDIDATES = tuple(CLINICAL_COLUMNS[3:-1])
COMPREHENSIVE_CLINICAL = tuple(
    column for column in BASELINE_CLINICAL_CANDIDATES
    if column not in REMOVED_EXACT_DUPLICATES
)

# The 46-variable sensitivity set is explicit rather than derived at runtime,
# so a workbook reordering cannot silently alter the scientific specification.
SUMMARY_CLINICAL = (
    "drug", "educ", "marital.rec", "occup.rec", "children.rec", "smoker",
    "bmi_wk0", "age", "ageonset", "sex",
    "f1score0", "f2score0", "f3score0", "f61score0", "f62score0",
    "f63score0", "f64score0", "f65score0", "f66score0", "mscore",
    "melanch", "atscore", "atyp", "anxdep", "anxsomdep",
    "nevents", "eventsbin",
    "mhssri", "everssri", "mhtca", "evertca", "mhdual", "everdual",
    "mhoth", "everoth", "mhantidep", "everantidep", "concssri", "conctca",
    "concantidep", "mcbenzo", "mczhypn", "mhmirta",
    "madrs.total", "hrsd.total", "bdi.total",
)

MODEL_PREDICTOR_SETS = {
    "primary_comprehensive_clinical": (*COMPREHENSIVE_CLINICAL, *PC_COLUMNS),
    "primary_comprehensive_combined": (*COMPREHENSIVE_CLINICAL, *PC_COLUMNS, *PRS_COLUMNS),
    "secondary_prs_focused": (TREATMENT_COLUMN, *PC_COLUMNS, *PRS_COLUMNS),
    "sensitivity_summary_clinical": (*SUMMARY_CLINICAL, *PC_COLUMNS),
    "sensitivity_summary_combined": (*SUMMARY_CLINICAL, *PC_COLUMNS, *PRS_COLUMNS),
}
EXPECTED_MODEL_COUNTS = {
    "primary_comprehensive_clinical": 134,
    "primary_comprehensive_combined": 142,
    "secondary_prs_focused": 15,
    "sensitivity_summary_clinical": 52,
    "sensitivity_summary_combined": 60,
}

BINARY_CATEGORICAL = (
    "drug", "marital.rec", "occup.rec", "smoker", "sex", "melanch", "atyp",
    "anxdep", "anxsomdep", "newbleq1wk0", "newbleq2wk0", "newbleq4wk0",
    "newbleq5wk0", "newbleq6wk0", "newbleq7wk0", "newbleq8wk0",
    "newbleq9wk0", "eventsbin", "everssri", "evertca", "everdual", "everoth",
    "everantidep", "concssri", "conctca", "concantidep", "mcbenzo",
    "mczhypn", "mhmirta",
)
COUNT_NUMERIC = (
    "educ", "children.rec", "mscore", "atscore", "nevents", "mhssri",
    "mhtca", "mhdual", "mhoth", "mhantidep",
)
CONTINUOUS_CLINICAL = (
    "bmi_wk0", "age", "ageonset", "f1score0", "f2score0", "f3score0",
    "f61score0", "f62score0", "f63score0", "f64score0", "f65score0",
    "f66score0", "madrs.total", "hrsd.total", "bdi.total",
)
EXPECTED_SEMANTIC_COUNTS = {
    "binary_categorical": 29,
    "count_numeric": 10,
    "ordinal_MADRS_item": 10,
    "ordinal_HAMD_item": 16,
    "ordinal_HAMD_item_source_fractional": 1,
    "ordinal_BDI_item": 19,
    "ordinal_SCAN_item": 28,
    "continuous_numeric": 15,
    "continuous_ancestry_PC": 6,
    "continuous_polygenic_score": 8,
}


def source_section(variable: str) -> str:
    """Return the source-workbook section for a clinical variable."""
    position = CLINICAL_COLUMNS.index(variable) + 1
    if position <= 3:
        return "Administrative identifiers"
    if position == 4:
        return "Treatment allocation"
    if position <= 13:
        return "Demographics and baseline characteristics"
    if position <= 23:
        return "Baseline MADRS items"
    if position <= 40:
        return "Baseline HAMD_HRSD items"
    if position <= 62:
        return "Baseline BDI items"
    if position <= 109:
        return "Depression subtypes symptom factors dimensions and SCAN variables"
    if position <= 119:
        return "Stressful life events"
    if position <= 135:
        return "Antidepressant medication history"
    if position <= 138:
        return "Derived baseline clinical totals"
    return "Outcome"


def semantic_type(variable: str) -> str:
    """Return the fixed modelling type assigned to a clinical variable."""
    if variable in BINARY_CATEGORICAL:
        return "binary_categorical"
    if variable in COUNT_NUMERIC:
        return "count_numeric"
    if variable in CONTINUOUS_CLINICAL:
        return "continuous_numeric"
    if variable in PC_COLUMNS:
        return "continuous_ancestry_PC"
    if variable in PRS_COLUMNS:
        return "continuous_polygenic_score"
    if re.fullmatch(r"madrs\d+wk0", variable):
        return "ordinal_MADRS_item"
    if variable == "hamd16wk0":
        return "ordinal_HAMD_item_source_fractional"
    if re.fullmatch(r"hamd\d+wk0", variable):
        return "ordinal_HAMD_item"
    if re.fullmatch(r"bdi(?:\d+|\d+[a-z])wk0", variable):
        return "ordinal_BDI_item"
    if re.fullmatch(r"k\d+", variable):
        return "ordinal_SCAN_item"
    return "unclassified"


def expected_storage_class(variable: str) -> str:
    """Return the expected source storage class for a clinical variable."""
    if variable in IDENTIFIER_COLUMNS:
        return "text_or_categorical"
    if variable == OUTCOME_COLUMN or variable in BINARY_CATEGORICAL:
        return "binary_integer"
    if variable == "hamd16wk0" or variable in {"bmi_wk0", "hrsd.total"}:
        return "continuous_numeric"
    return "integer"


def final_disposition(variable: str) -> str:
    """Describe whether a source variable is retained, excluded or duplicated."""
    if variable in IDENTIFIER_COLUMNS:
        return "Identifier; excluded"
    if variable == OUTCOME_COLUMN:
        return "Outcome; excluded"
    if variable in REMOVED_EXACT_DUPLICATES:
        return "Later exact duplicate; excluded"
    return "Retained predictor"


def model_role(variable: str, model: str) -> str:
    """Return the fixed model role assigned to a clinical variable."""
    if variable not in MODEL_PREDICTOR_SETS[model]:
        return "excluded"
    if variable in PC_COLUMNS or variable == TREATMENT_COLUMN:
        return "forced_unpenalised"
    if variable in PRS_COLUMNS:
        return "PRS_predictor"
    return "regularised_clinical_predictor"


def validate_static_policy() -> None:
    """Validate static policy against the fixed project invariants."""
    if len(BASELINE_CLINICAL_CANDIDATES) != 135:
        raise RuntimeError("Expected 135 baseline clinical candidates")
    if len(REMOVED_EXACT_DUPLICATES) != 7:
        raise RuntimeError("Expected seven removed exact duplicates")
    if len(COMPREHENSIVE_CLINICAL) != 128:
        raise RuntimeError("Expected 128 comprehensive clinical predictors")
    if len(SUMMARY_CLINICAL) != 46:
        raise RuntimeError("Expected 46 summary clinical predictors")
    if len(set(SUMMARY_CLINICAL)) != len(SUMMARY_CLINICAL):
        raise RuntimeError("Summary predictor set contains duplicates")
    for name, predictors in MODEL_PREDICTOR_SETS.items():
        if len(predictors) != EXPECTED_MODEL_COUNTS[name]:
            raise RuntimeError(f"Unexpected predictor count for {name}")
        if len(predictors) != len(set(predictors)):
            raise RuntimeError(f"Duplicate predictor in {name}")
        forbidden = set(IDENTIFIER_COLUMNS) | {OUTCOME_COLUMN, EUR_INDICATOR}
        if forbidden & set(predictors):
            raise RuntimeError(f"Forbidden variable in {name}")
    observed: dict[str, int] = {}
    for variable in MODEL_PREDICTOR_SETS["primary_comprehensive_combined"]:
        kind = semantic_type(variable)
        observed[kind] = observed.get(kind, 0) + 1
    if observed != EXPECTED_SEMANTIC_COUNTS:
        raise RuntimeError(f"Semantic count mismatch: {observed}")


validate_static_policy()
