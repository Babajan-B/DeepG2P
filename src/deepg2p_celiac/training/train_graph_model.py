from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import pandas as pd

from deepg2p_celiac.ingestion.validate_data import validate_phenotype_file, validate_storage
from deepg2p_celiac.models.graph_smoothing import train_graph_smoothed_models


def main() -> None:
    parser = argparse.ArgumentParser(description="Train graph-smoothed celiac classifiers.")
    parser.add_argument("--gene-features", required=True, help="Gene feature matrix CSV.")
    parser.add_argument("--adjacency", required=True, help="Gene adjacency CSV.")
    parser.add_argument("--phenotypes", required=True, help="Clean phenotype CSV.")
    parser.add_argument("--output-dir", required=True, help="Directory for graph model outputs.")
    parser.add_argument("--alpha", type=float, default=0.5, help="Neighbor smoothing weight.")
    parser.add_argument("--max-gb", type=float, default=6.0)
    parser.add_argument("--min-samples-per-class", type=int, default=5)
    args = parser.parse_args()

    validate_phenotype_file(args.phenotypes)
    validate_storage(["data/raw", "data/interim", "data/processed"], args.max_gb)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_features = pd.read_csv(args.gene_features)
    adjacency = pd.read_csv(args.adjacency, index_col=0)
    phenotypes = pd.read_csv(args.phenotypes)

    results = train_graph_smoothed_models(
        gene_features,
        adjacency,
        phenotypes,
        alpha=args.alpha,
        min_samples_per_class=args.min_samples_per_class,
    )
    best = max(results, key=lambda result: result.balanced_accuracy)
    joblib.dump(
        {"model": best.model, "feature_columns": best.feature_columns, "model_name": best.model_name},
        output_dir / "best_graph_model.joblib",
    )
    metrics = [
        {
            "model_name": result.model_name,
            "balanced_accuracy": result.balanced_accuracy,
            "auroc": result.auroc,
            "feature_columns": result.feature_columns,
        }
        for result in results
    ]
    (output_dir / "graph_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

