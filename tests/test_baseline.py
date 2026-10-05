import pandas as pd
import pytest

from deepg2p_celiac.models.baseline import train_baselines


def test_train_baselines() -> None:
    genotypes = pd.DataFrame(
        {
            "sample_id": [f"s{i}" for i in range(8)],
            "rs1": [0, 0, 0, 1, 2, 2, 2, 1],
            "rs2": [0, 1, 0, 1, 2, 1, 2, 2],
        }
    )
    phenotypes = pd.DataFrame(
        {
            "sample_id": [f"s{i}" for i in range(8)],
            "celiac_status": ["control", "control", "control", "control", "case", "case", "case", "case"],
        }
    )

    results = train_baselines(genotypes, phenotypes, test_size=0.25, min_samples_per_class=4)

    assert {result.model_name for result in results} == {"logistic_regression", "random_forest"}
    assert all(result.feature_columns == ["rs1", "rs2"] for result in results)


def test_train_baselines_rejects_tiny_class_counts() -> None:
    genotypes = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c", "d"],
            "rs1": [0, 1, 2, 2],
        }
    )
    phenotypes = pd.DataFrame(
        {
            "sample_id": ["a", "b", "c", "d"],
            "celiac_status": ["control", "control", "case", "case"],
        }
    )

    with pytest.raises(ValueError, match="too small"):
        train_baselines(genotypes, phenotypes)
