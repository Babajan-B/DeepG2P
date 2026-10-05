from __future__ import annotations

from pathlib import Path

import joblib
import pandas as pd

from deepg2p_celiac.annotation.snp_to_gene import aggregate_gene_features, load_snp_gene_map
from deepg2p_celiac.ingestion.extract_panel import extract_sample_panel, load_panel
from deepg2p_celiac.models.graph_smoothing import smooth_gene_features


def _predict_bundle(bundle: dict, features: pd.DataFrame) -> tuple[float, str]:
    feature_columns = bundle["feature_columns"]
    x = features.reindex(columns=feature_columns).apply(pd.to_numeric, errors="coerce").fillna(0)
    model = bundle["model"]
    score = float(model.predict_proba(x)[:, 1][0]) if hasattr(model, "predict_proba") else 0.0
    label = "case" if score >= 0.5 else "control"
    return score, label


def load_attribution(path: str | Path | None, limit: int = 5) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    table = pd.read_csv(path).head(limit)
    return table.fillna("").to_dict("records")


def predict_snp_panel(
    genotype_path: str | Path,
    panel_path: str | Path,
    model_path: str | Path,
) -> dict:
    panel = load_panel(panel_path)
    row = extract_sample_panel(genotype_path, panel)
    frame = pd.DataFrame([row])
    bundle = joblib.load(model_path)
    score, label = _predict_bundle(bundle, frame)
    return {
        "sample_id": str(row["sample_id"]),
        "model_name": str(bundle["model_name"]),
        "research_risk_score": score,
        "predicted_label": label,
    }


def predict_graph_smoothed(
    genotype_path: str | Path,
    panel_path: str | Path,
    graph_model_path: str | Path,
    adjacency_path: str | Path,
    alpha: float = 0.5,
) -> dict:
    panel = load_panel(panel_path)
    row = extract_sample_panel(genotype_path, panel)
    snp_frame = pd.DataFrame([row])
    snp_gene_map = load_snp_gene_map(panel_path)
    gene_features = aggregate_gene_features(snp_frame, snp_gene_map)
    adjacency = pd.read_csv(adjacency_path, index_col=0)
    smoothed = smooth_gene_features(gene_features, adjacency, alpha=alpha)
    bundle = joblib.load(graph_model_path)
    score, label = _predict_bundle(bundle, smoothed)
    return {
        "sample_id": str(row["sample_id"]),
        "model_name": str(bundle["model_name"]),
        "research_risk_score": score,
        "predicted_label": label,
    }

