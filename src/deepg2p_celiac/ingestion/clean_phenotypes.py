from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from deepg2p_celiac.ingestion.phenotype_loader import normalize_celiac_label
from deepg2p_celiac.ingestion.validate_data import find_suspicious_sample_ids


def clean_phenotype_table(
    input_path: str | Path,
    label_column: str,
    output_path: str | Path,
    excluded_path: str | Path | None = None,
) -> dict:
    raw = pd.read_csv(input_path)
    required = {"sample_id", label_column}
    missing = required.difference(raw.columns)
    if missing:
        raise ValueError(f"Phenotype file is missing required columns: {sorted(missing)}")

    suspicious = find_suspicious_sample_ids(raw["sample_id"])
    if suspicious:
        raise ValueError(
            "Suspicious sample IDs found. Refusing clean because synthetic/demo data is disabled: "
            + ", ".join(suspicious[:10])
        )

    working = raw[["sample_id", label_column]].copy()
    working["celiac_status"] = working[label_column].map(normalize_celiac_label)

    cleaned = working.dropna(subset=["celiac_status"])[["sample_id", "celiac_status"]]
    cleaned = cleaned.drop_duplicates(subset=["sample_id"])
    excluded = working[working["celiac_status"].isna()][["sample_id", label_column]]

    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    cleaned.to_csv(output, index=False)

    if excluded_path:
        excluded_output = Path(excluded_path)
        excluded_output.parent.mkdir(parents=True, exist_ok=True)
        excluded.to_csv(excluded_output, index=False)

    return {
        "input_rows": int(len(raw)),
        "cleaned_rows": int(len(cleaned)),
        "excluded_rows": int(len(excluded)),
        "case_count": int((cleaned["celiac_status"] == "case").sum()),
        "control_count": int((cleaned["celiac_status"] == "control").sum()),
        "output": str(output),
        "excluded_output": str(excluded_path) if excluded_path else None,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Clean real celiac phenotype labels.")
    parser.add_argument("--input", required=True, help="Raw phenotype CSV.")
    parser.add_argument("--label-column", default="celiac_status", help="Raw phenotype label column.")
    parser.add_argument("--output", required=True, help="Cleaned case/control CSV output.")
    parser.add_argument("--excluded-output", help="Optional CSV for excluded ambiguous rows.")
    parser.add_argument("--report", help="Optional JSON report output.")
    args = parser.parse_args()

    report = clean_phenotype_table(
        args.input,
        args.label_column,
        args.output,
        args.excluded_output,
    )
    text = json.dumps(report, indent=2)
    if args.report:
        report_path = Path(args.report)
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()

