import pandas as pd

from deepg2p_celiac.models.graph_smoothing import smooth_gene_features, train_graph_smoothed_models


def test_smooth_gene_features() -> None:
    features = pd.DataFrame({"sample_id": ["real_a"], "GENE1": [1.0], "GENE2": [0.0]})
    adjacency = pd.DataFrame([[0, 1], [1, 0]], index=["GENE1", "GENE2"], columns=["GENE1", "GENE2"])

    smoothed = smooth_gene_features(features, adjacency, alpha=0.5)

    assert smoothed.loc[0, "GENE1"] == 0.5
    assert smoothed.loc[0, "GENE2"] == 0.5


def test_train_graph_smoothed_models() -> None:
    features = pd.DataFrame(
        {
            "sample_id": [f"real_{i}" for i in range(8)],
            "GENE1": [0, 0, 0, 1, 2, 2, 2, 1],
            "GENE2": [0, 1, 0, 1, 2, 1, 2, 2],
        }
    )
    adjacency = pd.DataFrame([[0, 1], [1, 0]], index=["GENE1", "GENE2"], columns=["GENE1", "GENE2"])
    phenotypes = pd.DataFrame(
        {
            "sample_id": [f"real_{i}" for i in range(8)],
            "celiac_status": ["control", "control", "control", "control", "case", "case", "case", "case"],
        }
    )

    results = train_graph_smoothed_models(
        features,
        adjacency,
        phenotypes,
        min_samples_per_class=4,
    )

    assert {result.model_name for result in results} == {
        "graph_smoothed_logistic_regression",
        "graph_smoothed_random_forest",
    }

