from pathlib import Path

import joblib
import pandas as pd
from sklearn.linear_model import LogisticRegression

from deepg2p_celiac.api.prediction import predict_graph_smoothed, predict_snp_panel


def _write_genotype(path: Path) -> None:
    path.write_text(
        "rs1\t6\t1\tAA\n"
        "rs2\t6\t2\tGG\n",
        encoding="utf-8",
    )


def test_predict_snp_panel(tmp_path: Path) -> None:
    genotype = tmp_path / "real_upload.txt"
    panel = tmp_path / "panel.csv"
    model_path = tmp_path / "model.joblib"
    _write_genotype(genotype)
    pd.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "gene": ["GENE1", "GENE2"],
            "category": ["immune", "immune"],
        }
    ).to_csv(panel, index=False)

    x = pd.DataFrame({"rs1": [0, 0, 2, 2], "rs2": [0, 1, 2, 2]})
    y = [0, 0, 1, 1]
    model = LogisticRegression().fit(x, y)
    joblib.dump(
        {"model": model, "feature_columns": ["rs1", "rs2"], "model_name": "logistic"},
        model_path,
    )

    result = predict_snp_panel(genotype, panel, model_path)

    assert result["sample_id"] == "real_upload"
    assert 0 <= result["research_risk_score"] <= 1


def test_predict_graph_smoothed(tmp_path: Path) -> None:
    genotype = tmp_path / "real_upload.txt"
    panel = tmp_path / "panel.csv"
    model_path = tmp_path / "graph_model.joblib"
    adjacency = tmp_path / "adjacency.csv"
    _write_genotype(genotype)
    pd.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "gene": ["GENE1", "GENE2"],
            "category": ["immune", "immune"],
        }
    ).to_csv(panel, index=False)
    pd.DataFrame([[0, 1], [1, 0]], index=["GENE1", "GENE2"], columns=["GENE1", "GENE2"]).to_csv(
        adjacency
    )

    x = pd.DataFrame({"GENE1": [0, 0, 2, 2], "GENE2": [0, 1, 2, 2]})
    y = [0, 0, 1, 1]
    model = LogisticRegression().fit(x, y)
    joblib.dump(
        {"model": model, "feature_columns": ["GENE1", "GENE2"], "model_name": "graph_logistic"},
        model_path,
    )

    result = predict_graph_smoothed(genotype, panel, model_path, adjacency)

    assert result["sample_id"] == "real_upload"
    assert result["model_name"] == "graph_logistic"

