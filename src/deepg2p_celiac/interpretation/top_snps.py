from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd


def _extract_estimator(model: object) -> object:
    return model.named_steps["model"] if hasattr(model, "named_steps") else model


def feature_importance_table(model_bundle: dict, limit: int = 10) -> pd.DataFrame:
    model = model_bundle["model"]
    feature_columns = model_bundle["feature_columns"]
    classifier = _extract_estimator(model)
    if hasattr(classifier, "coef_"):
        importances = classifier.coef_[0]
    elif hasattr(classifier, "feature_importances_"):
        importances = classifier.feature_importances_
    else:
        return pd.DataFrame(columns=["feature", "importance"])

    rows = [
        {"feature": feature, "importance": float(importance)}
        for feature, importance in zip(feature_columns, importances, strict=False)
    ]
    return (
        pd.DataFrame(rows)
        .assign(abs_importance=lambda frame: frame["importance"].abs())
        .sort_values("abs_importance", ascending=False)
        .drop(columns=["abs_importance"])
        .head(limit)
    )


def top_logistic_features(model_bundle: dict, limit: int = 5) -> pd.DataFrame:
    table = feature_importance_table(model_bundle, limit)
    return table.rename(columns={"feature": "rsid"})


def annotate_feature_importance(
    importance: pd.DataFrame,
    panel_path: str | Path,
    feature_column: str = "feature",
) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    if feature_column not in importance.columns:
        raise ValueError(f"Importance table must contain '{feature_column}'.")
    if "rsid" in panel.columns:
        annotated = importance.merge(
            panel,
            left_on=feature_column,
            right_on="rsid",
            how="left",
        )
    else:
        annotated = importance.copy()
    return annotated


def write_top_features(
    model_path: str | Path,
    output_path: str | Path,
    panel_path: str | Path | None = None,
    limit: int = 10,
) -> pd.DataFrame:
    bundle = joblib.load(model_path)
    importance = feature_importance_table(bundle, limit=limit)
    if panel_path:
        importance = annotate_feature_importance(importance, panel_path)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    importance.to_csv(output, index=False)
    return importance
