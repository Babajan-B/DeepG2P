from __future__ import annotations

import argparse
import json
from pathlib import Path

import joblib
import pandas as pd

from deepg2p_celiac.analysis.validation import run_scientific_validation
from deepg2p_celiac.graph.export_graph_dataset import export_graph_dataset
from deepg2p_celiac.ingestion.clean_phenotypes import clean_phenotype_table
from deepg2p_celiac.ingestion.extract_panel import extract_panel_matrix
from deepg2p_celiac.ingestion.validate_data import (
    validate_genotype_dir,
    validate_phenotype_file,
    validate_storage,
)
from deepg2p_celiac.models.baseline import train_baselines
from deepg2p_celiac.models.graph_smoothing import train_graph_smoothed_models
from deepg2p_celiac.interpretation.report import render_html_report
from deepg2p_celiac.interpretation.top_snps import write_top_features
from deepg2p_celiac.pipeline.reporting import (
    build_run_metadata,
    summarize_matrix,
    summarize_storage,
    write_json_report,
)
from deepg2p_celiac.qc.matrix_qc import filter_missingness, impute_panel_mean


def run_pipeline(
    genotype_dir: str | Path,
    phenotype_file: str | Path,
    panel_file: str | Path,
    output_dir: str | Path,
    label_column: str = "celiac_status",
    max_gb: float = 6.0,
    min_samples_per_class: int = 5,
    validation_folds: int = 5,
    validation_permutations: int = 100,
) -> dict:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    validation = {
        "genotypes": validate_genotype_dir(genotype_dir),
        "storage_before": validate_storage(["data/raw", "data/interim", "data/processed"], max_gb),
    }

    clean_phenotypes_path = output / "phenotypes_clean.csv"
    excluded_phenotypes_path = output / "phenotypes_excluded.csv"
    phenotype_report = clean_phenotype_table(
        phenotype_file,
        label_column,
        clean_phenotypes_path,
        excluded_phenotypes_path,
    )
    validation["phenotypes"] = validate_phenotype_file(clean_phenotypes_path)

    raw_matrix_path = output / "celiac_panel_matrix_raw.csv"
    raw_matrix = extract_panel_matrix(genotype_dir, panel_file)
    raw_matrix.to_csv(raw_matrix_path, index=False)

    filtered_matrix = filter_missingness(raw_matrix)
    filtered_matrix = impute_panel_mean(filtered_matrix)
    filtered_matrix_path = output / "celiac_panel_matrix_qc.csv"
    filtered_matrix.to_csv(filtered_matrix_path, index=False)

    phenotypes = pd.read_csv(clean_phenotypes_path)
    results = train_baselines(
        filtered_matrix,
        phenotypes,
        min_samples_per_class=min_samples_per_class,
    )
    model_dir = output / "model"
    model_dir.mkdir(exist_ok=True)
    best = max(results, key=lambda result: result.balanced_accuracy)
    joblib.dump(
        {"model": best.model, "feature_columns": best.feature_columns, "model_name": best.model_name},
        model_dir / "best_model.joblib",
    )
    top_snps_path = model_dir / "top_snp_features.csv"
    write_top_features(model_dir / "best_model.joblib", top_snps_path, panel_file)
    metrics = [
        {
            "model_name": result.model_name,
            "balanced_accuracy": result.balanced_accuracy,
            "auroc": result.auroc,
            "feature_columns": result.feature_columns,
        }
        for result in results
    ]
    (model_dir / "metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")

    graph_dir = output / "graph"
    graph_summary = export_graph_dataset(filtered_matrix_path, panel_file, graph_dir)
    graph_results = train_graph_smoothed_models(
        pd.read_csv(graph_summary["gene_features"]),
        pd.read_csv(graph_summary["gene_adjacency"], index_col=0),
        phenotypes,
        min_samples_per_class=min_samples_per_class,
    )
    graph_model_dir = output / "graph_model"
    graph_model_dir.mkdir(exist_ok=True)
    best_graph = max(graph_results, key=lambda result: result.balanced_accuracy)
    joblib.dump(
        {
            "model": best_graph.model,
            "feature_columns": best_graph.feature_columns,
            "model_name": best_graph.model_name,
        },
        graph_model_dir / "best_graph_model.joblib",
    )
    top_genes_path = graph_model_dir / "top_gene_features.csv"
    write_top_features(graph_model_dir / "best_graph_model.joblib", top_genes_path)
    graph_metrics = [
        {
            "model_name": result.model_name,
            "balanced_accuracy": result.balanced_accuracy,
            "auroc": result.auroc,
            "feature_columns": result.feature_columns,
        }
        for result in graph_results
    ]
    (graph_model_dir / "graph_metrics.json").write_text(
        json.dumps(graph_metrics, indent=2),
        encoding="utf-8",
    )
    validation_report = run_scientific_validation(
        filtered_matrix_path,
        clean_phenotypes_path,
        output / "scientific_validation",
        folds=validation_folds,
        permutations=validation_permutations,
        max_gb=max_gb,
    )

    report = {
        **build_run_metadata(genotype_dir, phenotype_file, panel_file, output),
        "validation": validation,
        "phenotype_cleaning": phenotype_report,
        "raw_matrix": summarize_matrix(raw_matrix_path),
        "qc_matrix": summarize_matrix(filtered_matrix_path),
        "models": metrics,
        "best_model": best.model_name,
        "graph": graph_summary,
        "graph_models": graph_metrics,
        "best_graph_model": best_graph.model_name,
        "scientific_validation": validation_report,
        "storage_mb": summarize_storage(["data/raw", "data/interim", "data/processed"]),
    }
    run_report_path = output / "run_report.json"
    write_json_report(report, run_report_path)
    html_report_path = output / "research_report.html"
    render_html_report(run_report_path, html_report_path, top_snps_path, top_genes_path)
    report["html_report"] = str(html_report_path)
    write_json_report(report, run_report_path)
    validate_storage(["data/raw", "data/interim", "data/processed"], max_gb)
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the real-data DeepG2P-Celiac MVP pipeline.")
    parser.add_argument("--genotype-dir", default="data/raw/genotypes")
    parser.add_argument("--phenotypes", default="data/raw/phenotypes.csv")
    parser.add_argument("--panel", default="configs/celiac_snp_panel.csv")
    parser.add_argument("--output-dir", default="data/processed/run")
    parser.add_argument("--label-column", default="celiac_status")
    parser.add_argument("--max-gb", type=float, default=6.0)
    parser.add_argument("--min-samples-per-class", type=int, default=5)
    parser.add_argument("--validation-folds", type=int, default=5)
    parser.add_argument("--validation-permutations", type=int, default=100)
    args = parser.parse_args()

    report = run_pipeline(
        args.genotype_dir,
        args.phenotypes,
        args.panel,
        args.output_dir,
        args.label_column,
        args.max_gb,
        args.min_samples_per_class,
        args.validation_folds,
        args.validation_permutations,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
