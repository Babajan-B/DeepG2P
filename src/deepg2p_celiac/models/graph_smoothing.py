from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class GraphTrainingResult:
    model_name: str
    balanced_accuracy: float
    auroc: float | None
    feature_columns: list[str]
    model: object


def smooth_gene_features(
    gene_features: pd.DataFrame,
    adjacency: pd.DataFrame,
    alpha: float = 0.5,
) -> pd.DataFrame:
    if "sample_id" not in gene_features.columns:
        raise ValueError("Gene feature matrix must contain 'sample_id'.")
    if not 0 <= alpha <= 1:
        raise ValueError("alpha must be between 0 and 1.")

    genes = [column for column in gene_features.columns if column != "sample_id"]
    adjacency = adjacency.reindex(index=genes, columns=genes).fillna(0)
    matrix = gene_features[genes].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy()
    adj = adjacency.to_numpy(dtype=float)
    degrees = adj.sum(axis=1)
    normalized = np.divide(adj, degrees[:, None], out=np.zeros_like(adj), where=degrees[:, None] != 0)
    smoothed = (1 - alpha) * matrix + alpha * matrix @ normalized.T

    output = pd.DataFrame(smoothed, columns=genes)
    output.insert(0, "sample_id", gene_features["sample_id"].to_numpy())
    return output


def train_graph_smoothed_models(
    gene_features: pd.DataFrame,
    adjacency: pd.DataFrame,
    phenotypes: pd.DataFrame,
    label_column: str = "celiac_status",
    test_size: float = 0.25,
    random_state: int = 42,
    min_samples_per_class: int = 5,
    alpha: float = 0.5,
) -> list[GraphTrainingResult]:
    smoothed = smooth_gene_features(gene_features, adjacency, alpha)
    merged = smoothed.merge(phenotypes[["sample_id", label_column]], on="sample_id", how="inner")
    if merged.empty:
        raise ValueError("No overlapping sample IDs between gene features and phenotypes.")

    feature_columns = [column for column in smoothed.columns if column != "sample_id"]
    x = merged[feature_columns].apply(pd.to_numeric, errors="coerce").fillna(0)
    y = merged[label_column].map({"control": 0, "case": 1})
    if y.isna().any() or y.nunique() < 2:
        raise ValueError("Graph training requires valid case and control labels.")
    class_counts = y.value_counts()
    too_small = class_counts[class_counts < min_samples_per_class]
    if not too_small.empty:
        raise ValueError(
            "Graph training set is too small for a meaningful split. "
            f"Minimum samples per class is {min_samples_per_class}; observed "
            f"{too_small.to_dict()}."
        )

    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y.astype(int),
        test_size=test_size,
        random_state=random_state,
        stratify=y,
    )

    models = {
        "graph_smoothed_logistic_regression": Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=1000, class_weight="balanced")),
            ]
        ),
        "graph_smoothed_random_forest": RandomForestClassifier(
            n_estimators=200,
            random_state=random_state,
            class_weight="balanced",
            min_samples_leaf=2,
        ),
    }

    results: list[GraphTrainingResult] = []
    for name, model in models.items():
        model.fit(x_train, y_train)
        predictions = model.predict(x_test)
        probabilities = (
            model.predict_proba(x_test)[:, 1] if hasattr(model, "predict_proba") else None
        )
        auroc = None
        if probabilities is not None and y_test.nunique() == 2:
            auroc = float(roc_auc_score(y_test, probabilities))
        results.append(
            GraphTrainingResult(
                model_name=name,
                balanced_accuracy=float(balanced_accuracy_score(y_test, predictions)),
                auroc=auroc,
                feature_columns=feature_columns,
                model=model,
            )
        )
    return results

