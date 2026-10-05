import json
from pathlib import Path

import joblib
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from deepg2p_celiac.interpretation.report import render_html_report
from deepg2p_celiac.interpretation.top_snps import write_top_features


def test_write_top_features_with_panel_annotation(tmp_path: Path) -> None:
    x = pd.DataFrame({"rs1": [0, 0, 1, 2, 2, 1], "rs2": [0, 1, 1, 2, 2, 0]})
    y = [0, 0, 0, 1, 1, 1]
    model = Pipeline(
        [
            ("scale", StandardScaler()),
            ("model", LogisticRegression(max_iter=1000)),
        ]
    )
    model.fit(x, y)
    model_path = tmp_path / "model.joblib"
    output_path = tmp_path / "top.csv"
    panel_path = tmp_path / "panel.csv"
    joblib.dump({"model": model, "feature_columns": ["rs1", "rs2"]}, model_path)
    pd.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "gene": ["GENE1", "GENE2"],
            "category": ["immune", "immune"],
        }
    ).to_csv(panel_path, index=False)

    result = write_top_features(model_path, output_path, panel_path, limit=2)

    assert output_path.exists()
    assert {"feature", "importance", "rsid", "gene"}.issubset(result.columns)


def test_render_html_report(tmp_path: Path) -> None:
    run_report = tmp_path / "run_report.json"
    html_report = tmp_path / "report.html"
    top_snps = tmp_path / "top_snps.csv"
    top_genes = tmp_path / "top_genes.csv"
    run_report.write_text(
        json.dumps(
            {
                "created_at_utc": "2026-01-01T00:00:00+00:00",
                "project": "DeepG2P-Celiac",
                "synthetic_data_allowed": False,
                "raw_matrix": {"samples": 10, "snps": 2},
                "qc_matrix": {"samples": 10, "snps": 2},
                "graph": {"genes": 2, "edges": 1},
                "models": [{"model_name": "baseline", "balanced_accuracy": 0.5, "auroc": 0.5}],
                "graph_models": [{"model_name": "graph", "balanced_accuracy": 0.5, "auroc": 0.5}],
            }
        ),
        encoding="utf-8",
    )
    pd.DataFrame([{"feature": "rs1", "importance": 1.0, "gene": "GENE1"}]).to_csv(
        top_snps,
        index=False,
    )
    pd.DataFrame([{"feature": "GENE1", "importance": 1.0}]).to_csv(top_genes, index=False)

    html_text = render_html_report(run_report, html_report, top_snps, top_genes)

    assert html_report.exists()
    assert "Research-only output" in html_text
    assert "GENE1" in html_text

