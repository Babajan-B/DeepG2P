from __future__ import annotations

from dataclasses import dataclass

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


@dataclass(frozen=True)
class TrainingResult:
    model_name: str
    balanced_accuracy: float
    auroc: float | None
    feature_columns: list[str]
    model: object


def prepare_training_frame(
    genotypes: pd.DataFrame,
    phenotypes: pd.DataFrame,
    label_column: str = "celiac_status",
) -> tuple[pd.DataFrame, pd.Series]:
    merged = genotypes.merge(phenotypes[["sample_id", label_column]], on="sample_id", how="inner")
    if merged.empty:
        raise ValueError("No overlapping sample IDs between genotypes and phenotypes.")

    feature_columns = [column for column in genotypes.columns if column != "sample_id"]
    x = merged[feature_columns].apply(pd.to_numeric, errors="coerce")
    x = x.fillna(x.mean()).fillna(0)
    y = merged[label_column].map({"control": 0, "case": 1})
    if y.isna().any():
        raise ValueError("Labels must normalize to 'case' or 'control'.")
    return x, y.astype(int)


def train_baselines(
    genotypes: pd.DataFrame,
    phenotypes: pd.DataFrame,
    label_column: str = "celiac_status",
    test_size: float = 0.25,
    random_state: int = 42,
    min_samples_per_class: int = 5,
) -> list[TrainingResult]:
    x, y = prepare_training_frame(genotypes, phenotypes, label_column)
    if y.nunique() < 2:
        raise ValueError("Training requires both case and control samples.")
    class_counts = y.value_counts()
    too_small = class_counts[class_counts < min_samples_per_class]
    if not too_small.empty:
        raise ValueError(
            "Training set is too small for a meaningful split. "
            f"Minimum samples per class is {min_samples_per_class}; observed "
            f"{too_small.to_dict()}."
        )

    stratify = y
    x_train, x_test, y_train, y_test = train_test_split(
        x,
        y,
        test_size=test_size,
        random_state=random_state,
        stratify=stratify,
    )

    models = {
        "logistic_regression": Pipeline(
            [
                ("scale", StandardScaler()),
                ("model", LogisticRegression(max_iter=1000, class_weight="balanced")),
            ]
        ),
        "random_forest": RandomForestClassifier(
            n_estimators=200,
            random_state=random_state,
            class_weight="balanced",
            min_samples_leaf=2,
        ),
    }

    results: list[TrainingResult] = []
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
            TrainingResult(
                model_name=name,
                balanced_accuracy=float(balanced_accuracy_score(y_test, predictions)),
                auroc=auroc,
                feature_columns=list(x.columns),
                model=model,
            )
        )
    return results
