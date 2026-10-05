from __future__ import annotations

from pathlib import Path

import pandas as pd

CASE_TERMS = {"case", "celiac", "coeliac", "celiac disease", "coeliac disease", "diagnosed celiac"}
CONTROL_TERMS = {"control", "no", "none", "negative", "no celiac", "no coeliac"}
EXCLUDE_TERMS = {"gluten sensitivity", "wheat allergy", "ibs", "unknown", "maybe", "unsure", ""}


def normalize_celiac_label(value: object) -> str | None:
    text = str(value).strip().lower()
    if text in CASE_TERMS:
        return "case"
    if text in CONTROL_TERMS:
        return "control"
    if text in EXCLUDE_TERMS:
        return None
    if "celiac" in text or "coeliac" in text:
        return "case"
    return None


def load_clean_phenotypes(path: str | Path, label_column: str = "celiac_status") -> pd.DataFrame:
    df = pd.read_csv(path)
    if "sample_id" not in df.columns:
        raise ValueError("Phenotype CSV must contain 'sample_id'.")
    if label_column not in df.columns:
        raise ValueError(f"Phenotype CSV must contain '{label_column}'.")

    cleaned = df[["sample_id", label_column]].copy()
    cleaned[label_column] = cleaned[label_column].map(normalize_celiac_label)
    cleaned = cleaned.dropna(subset=[label_column])
    return cleaned.drop_duplicates(subset=["sample_id"])

