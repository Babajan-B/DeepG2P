from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import pandas as pd

from deepg2p_celiac.ingestion.validate_data import validate_phenotype_file, validate_storage
from deepg2p_celiac.ingestion.phenotype_loader import load_clean_phenotypes
from deepg2p_celiac.models.baseline import train_baselines
from deepg2p_celiac.qc.matrix_qc import filter_missingness, impute_panel_mean


def main() -> None:
    parser = argparse.ArgumentParser(description="Train celiac baseline classifiers.")
    parser.add_argument("--genotypes", required=True, help="Panel genotype matrix CSV.")
    parser.add_argument("--phenotypes", required=True, help="Phenotype CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for model outputs.")
    parser.add_argument("--max-gb", type=float, default=6.0, help="Maximum local project data size.")
    parser.add_argument("--min-samples-per-class", type=int, default=5)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    validate_phenotype_file(args.phenotypes)
    validate_storage(["data/raw", "data/interim", "data/processed"], args.max_gb)

    genotypes = pd.read_csv(args.genotypes)
    genotypes = impute_panel_mean(filter_missingness(genotypes))
    phenotypes = load_clean_phenotypes(args.phenotypes)
    results = train_baselines(
        genotypes,
        phenotypes,
        min_samples_per_class=args.min_samples_per_class,
    )

    summary = []
    best = max(results, key=lambda result: result.balanced_accuracy)
    for result in results:
        summary.append(
            {
                "model_name": result.model_name,
                "balanced_accuracy": result.balanced_accuracy,
                "auroc": result.auroc,
                "feature_columns": result.feature_columns,
            }
        )
    joblib.dump(
        {"model": best.model, "feature_columns": best.feature_columns, "model_name": best.model_name},
        output_dir / "best_model.joblib",
    )
    (output_dir / "metrics.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
