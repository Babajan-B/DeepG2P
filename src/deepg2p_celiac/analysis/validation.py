from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    brier_score_loss,
    confusion_matrix,
    f1_score,
    roc_auc_score,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from deepg2p_celiac.ingestion.validate_data import validate_phenotype_file, validate_storage
from deepg2p_celiac.models.baseline import prepare_training_frame


@dataclass(frozen=True)
class FoldMetrics:
    fold: int
    accuracy: float
    balanced_accuracy: float
    f1: float
    auroc: float | None
    brier: float | None
    tn: int
    fp: int
    fn: int
    tp: int


def candidate_models(random_state: int = 42) -> dict[str, object]:
    return {
        "logistic_regression": Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=1000, class_weight="balanced")),
            ]
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=300,
            random_state=random_state,
            class_weight="balanced",
            min_samples_leaf=2,
        ),
    }


def class_balance(y: pd.Series) -> dict:
    counts = y.value_counts().sort_index()
    return {
        "control": int(counts.get(0, 0)),
        "case": int(counts.get(1, 0)),
        "total": int(len(y)),
        "case_fraction": float(counts.get(1, 0) / len(y)) if len(y) else 0.0,
    }


def calibration_bins(y_true: np.ndarray, probabilities: np.ndarray, bins: int = 10) -> list[dict]:
    edges = np.linspace(0, 1, bins + 1)
    rows = []
    for idx in range(bins):
        lower = edges[idx]
        upper = edges[idx + 1]
        if idx == bins - 1:
            mask = (probabilities >= lower) & (probabilities <= upper)
        else:
            mask = (probabilities >= lower) & (probabilities < upper)
        if not mask.any():
            continue
        rows.append(
            {
                "bin_lower": float(lower),
                "bin_upper": float(upper),
                "n": int(mask.sum()),
                "mean_predicted": float(probabilities[mask].mean()),
                "observed_case_fraction": float(y_true[mask].mean()),
            }
        )
    return rows


def evaluate_model_cv(
    x: pd.DataFrame,
    y: pd.Series,
    model: object,
    folds: int = 5,
    random_state: int = 42,
) -> tuple[list[FoldMetrics], list[dict]]:
    min_class = int(y.value_counts().min())
    if min_class < 2:
        raise ValueError("Cross-validation requires at least 2 samples per class.")
    n_splits = min(folds, min_class)
    splitter = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    metrics: list[FoldMetrics] = []
    prediction_rows: list[dict] = []

    for fold, (train_idx, test_idx) in enumerate(splitter.split(x, y), start=1):
        estimator = clone(model)
        x_train, x_test = x.iloc[train_idx], x.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        estimator.fit(x_train, y_train)
        predicted = estimator.predict(x_test)
        probabilities = (
            estimator.predict_proba(x_test)[:, 1] if hasattr(estimator, "predict_proba") else None
        )
        tn, fp, fn, tp = confusion_matrix(y_test, predicted, labels=[0, 1]).ravel()
        auroc = None
        brier = None
        if probabilities is not None and y_test.nunique() == 2:
            auroc = float(roc_auc_score(y_test, probabilities))
            brier = float(brier_score_loss(y_test, probabilities))
        metrics.append(
            FoldMetrics(
                fold=fold,
                accuracy=float(accuracy_score(y_test, predicted)),
                balanced_accuracy=float(balanced_accuracy_score(y_test, predicted)),
                f1=float(f1_score(y_test, predicted, zero_division=0)),
                auroc=auroc,
                brier=brier,
                tn=int(tn),
                fp=int(fp),
                fn=int(fn),
                tp=int(tp),
            )
        )
        if probabilities is not None:
            for sample_id, truth, probability, label in zip(
                x_test.index,
                y_test,
                probabilities,
                predicted,
                strict=False,
            ):
                prediction_rows.append(
                    {
                        "sample_index": int(sample_id),
                        "fold": fold,
                        "truth": int(truth),
                        "predicted": int(label),
                        "probability_case": float(probability),
                    }
                )
    return metrics, prediction_rows


def permutation_test(
    x: pd.DataFrame,
    y: pd.Series,
    model: object,
    observed_balanced_accuracy: float,
    permutations: int = 100,
    folds: int = 5,
    random_state: int = 42,
) -> dict:
    rng = np.random.default_rng(random_state)
    scores = []
    y_array = y.to_numpy()
    for _ in range(permutations):
        shuffled = pd.Series(rng.permutation(y_array), index=y.index)
        metrics, _ = evaluate_model_cv(x, shuffled, model, folds=folds, random_state=random_state)
        scores.append(float(np.mean([metric.balanced_accuracy for metric in metrics])))
    scores_array = np.array(scores)
    p_value = float((np.sum(scores_array >= observed_balanced_accuracy) + 1) / (permutations + 1))
    return {
        "permutations": permutations,
        "observed_balanced_accuracy": observed_balanced_accuracy,
        "null_mean_balanced_accuracy": float(scores_array.mean()) if len(scores_array) else None,
        "null_std_balanced_accuracy": float(scores_array.std()) if len(scores_array) else None,
        "empirical_p_value": p_value,
    }


def summarize_metrics(metrics: list[FoldMetrics]) -> dict:
    rows = pd.DataFrame([asdict(metric) for metric in metrics])
    summary = {}
    for column in ["accuracy", "balanced_accuracy", "f1", "auroc", "brier"]:
        values = rows[column].dropna()
        summary[column] = {
            "mean": float(values.mean()) if not values.empty else None,
            "std": float(values.std(ddof=0)) if len(values) > 1 else 0.0,
        }
    summary["confusion_matrix_sum"] = {
        "tn": int(rows["tn"].sum()),
        "fp": int(rows["fp"].sum()),
        "fn": int(rows["fn"].sum()),
        "tp": int(rows["tp"].sum()),
    }
    return summary


def run_scientific_validation(
    genotypes_path: str | Path,
    phenotypes_path: str | Path,
    output_dir: str | Path,
    folds: int = 5,
    permutations: int = 100,
    random_state: int = 42,
    max_gb: float = 6.0,
) -> dict:
    validate_phenotype_file(phenotypes_path)
    validate_storage(["data/raw", "data/interim", "data/processed"], max_gb)

    genotypes = pd.read_csv(genotypes_path)
    phenotypes = pd.read_csv(phenotypes_path)
    x, y = prepare_training_frame(genotypes, phenotypes)
    x = x.reset_index(drop=True)
    y = y.reset_index(drop=True)

    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "genotypes": str(genotypes_path),
        "phenotypes": str(phenotypes_path),
        "class_balance": class_balance(y),
        "folds_requested": folds,
        "permutations_requested": permutations,
        "models": {},
        "warning": "Research validation only. Not a medical diagnostic evaluation.",
    }

    for model_name, model in candidate_models(random_state).items():
        metrics, predictions = evaluate_model_cv(x, y, model, folds=folds, random_state=random_state)
        metrics_rows = [asdict(metric) for metric in metrics]
        pd.DataFrame(metrics_rows).to_csv(output / f"{model_name}_fold_metrics.csv", index=False)
        pd.DataFrame(predictions).to_csv(output / f"{model_name}_cv_predictions.csv", index=False)
        probabilities = np.array([row["probability_case"] for row in predictions])
        truths = np.array([row["truth"] for row in predictions])
        calibration = calibration_bins(truths, probabilities) if len(predictions) else []
        pd.DataFrame(calibration).to_csv(output / f"{model_name}_calibration.csv", index=False)
        summary = summarize_metrics(metrics)
        observed = summary["balanced_accuracy"]["mean"]
        permutation = permutation_test(
            x,
            y,
            model,
            observed_balanced_accuracy=observed or 0.0,
            permutations=permutations,
            folds=folds,
            random_state=random_state,
        )
        report["models"][model_name] = {
            "summary": summary,
            "permutation_test": permutation,
            "files": {
                "fold_metrics": str(output / f"{model_name}_fold_metrics.csv"),
                "cv_predictions": str(output / f"{model_name}_cv_predictions.csv"),
                "calibration": str(output / f"{model_name}_calibration.csv"),
            },
        }

    report_path = output / "scientific_validation_report.json"
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Run scientific validation on real celiac data.")
    parser.add_argument("--genotypes", required=True, help="QC genotype matrix CSV.")
    parser.add_argument("--phenotypes", required=True, help="Clean phenotype CSV.")
    parser.add_argument("--output-dir", required=True, help="Validation output directory.")
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--permutations", type=int, default=100)
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--max-gb", type=float, default=6.0)
    args = parser.parse_args()

    report = run_scientific_validation(
        args.genotypes,
        args.phenotypes,
        args.output_dir,
        args.folds,
        args.permutations,
        args.random_state,
        args.max_gb,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

