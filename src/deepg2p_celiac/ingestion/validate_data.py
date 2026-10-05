from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

SUSPICIOUS_SAMPLE_TERMS = {
    "synthetic",
    "simulated",
    "fake",
    "dummy",
    "mock",
    "toy",
    "random",
    "example",
    "test",
}


def directory_size_bytes(path: str | Path) -> int:
    root = Path(path)
    if not root.exists():
        return 0
    return sum(file.stat().st_size for file in root.rglob("*") if file.is_file())


def find_suspicious_sample_ids(sample_ids: pd.Series) -> list[str]:
    suspicious: list[str] = []
    for sample_id in sample_ids.astype(str):
        lowered = sample_id.lower()
        if any(term in lowered for term in SUSPICIOUS_SAMPLE_TERMS):
            suspicious.append(sample_id)
    return suspicious


def validate_phenotype_file(path: str | Path, label_column: str = "celiac_status") -> dict:
    phenotype_path = Path(path)
    if not phenotype_path.exists():
        raise FileNotFoundError(f"Phenotype file not found: {phenotype_path}")

    df = pd.read_csv(phenotype_path)
    required = {"sample_id", label_column}
    missing = required.difference(df.columns)
    if missing:
        raise ValueError(f"Phenotype file is missing required columns: {sorted(missing)}")

    suspicious = find_suspicious_sample_ids(df["sample_id"])
    labels = set(df[label_column].dropna().astype(str).str.lower())
    invalid_labels = labels.difference({"case", "control"})
    if invalid_labels:
        raise ValueError(f"Phenotype file contains invalid labels: {sorted(invalid_labels)}")
    if suspicious:
        raise ValueError(
            "Suspicious sample IDs found. Refusing run because synthetic/demo data is disabled: "
            + ", ".join(suspicious[:10])
        )

    return {
        "phenotype_file": str(phenotype_path),
        "samples": int(len(df)),
        "case_count": int((df[label_column].str.lower() == "case").sum()),
        "control_count": int((df[label_column].str.lower() == "control").sum()),
    }


def validate_genotype_dir(path: str | Path) -> dict:
    genotype_dir = Path(path)
    if not genotype_dir.exists():
        raise FileNotFoundError(f"Genotype directory not found: {genotype_dir}")
    files = [file for file in genotype_dir.iterdir() if file.is_file() and not file.name.startswith(".")]
    if not files:
        raise ValueError(f"No genotype files found in: {genotype_dir}")

    suspicious_files = [
        file.name
        for file in files
        if any(term in file.name.lower() for term in SUSPICIOUS_SAMPLE_TERMS)
    ]
    if suspicious_files:
        raise ValueError(
            "Suspicious genotype filenames found. Refusing run because synthetic/demo data is disabled: "
            + ", ".join(suspicious_files[:10])
        )

    return {
        "genotype_dir": str(genotype_dir),
        "genotype_files": len(files),
        "size_mb": round(directory_size_bytes(genotype_dir) / 1024 / 1024, 3),
    }


def validate_storage(paths: list[str | Path], max_gb: float = 6.0) -> dict:
    total = sum(directory_size_bytes(path) for path in paths)
    total_gb = total / 1024 / 1024 / 1024
    if total_gb > max_gb:
        raise ValueError(f"Storage budget exceeded: {total_gb:.2f} GB used, max is {max_gb:.2f} GB.")
    return {"storage_used_gb": round(total_gb, 4), "storage_limit_gb": max_gb}


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate real-data inputs before a run.")
    parser.add_argument("--genotype-dir", default="data/raw/genotypes")
    parser.add_argument("--phenotypes", default="data/raw/phenotypes.csv")
    parser.add_argument("--max-gb", type=float, default=6.0)
    parser.add_argument("--output", help="Optional JSON validation report path.")
    args = parser.parse_args()

    report = {
        "synthetic_data_allowed": False,
        "genotypes": validate_genotype_dir(args.genotype_dir),
        "phenotypes": validate_phenotype_file(args.phenotypes),
        "storage": validate_storage(["data/raw", "data/interim", "data/processed"], args.max_gb),
    }

    text = json.dumps(report, indent=2)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()

