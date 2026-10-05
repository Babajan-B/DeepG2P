from pathlib import Path

import pandas as pd

from deepg2p_celiac.analysis.validation import (
    calibration_bins,
    class_balance,
    run_scientific_validation,
)


def test_class_balance() -> None:
    y = pd.Series([0, 0, 1, 1, 1])

    assert class_balance(y) == {"control": 2, "case": 3, "total": 5, "case_fraction": 0.6}


def test_calibration_bins() -> None:
    rows = calibration_bins(
        y_true=pd.Series([0, 1, 1]).to_numpy(),
        probabilities=pd.Series([0.1, 0.6, 0.9]).to_numpy(),
        bins=2,
    )

    assert len(rows) == 2
    assert rows[0]["n"] == 1


def test_run_scientific_validation(tmp_path: Path) -> None:
    genotypes = tmp_path / "genotypes.csv"
    phenotypes = tmp_path / "phenotypes.csv"
    output = tmp_path / "validation"
    pd.DataFrame(
        {
            "sample_id": [f"real_{i}" for i in range(12)],
            "rs1": [0, 0, 0, 1, 0, 1, 2, 2, 2, 1, 2, 1],
            "rs2": [0, 1, 0, 1, 0, 1, 2, 1, 2, 2, 1, 2],
        }
    ).to_csv(genotypes, index=False)
    pd.DataFrame(
        {
            "sample_id": [f"real_{i}" for i in range(12)],
            "celiac_status": [
                "control",
                "control",
                "control",
                "control",
                "control",
                "control",
                "case",
                "case",
                "case",
                "case",
                "case",
                "case",
            ],
        }
    ).to_csv(phenotypes, index=False)

    report = run_scientific_validation(
        genotypes,
        phenotypes,
        output,
        folds=3,
        permutations=2,
    )

    assert report["class_balance"]["case"] == 6
    assert (output / "scientific_validation_report.json").exists()
    assert (output / "logistic_regression_fold_metrics.csv").exists()
    assert "permutation_test" in report["models"]["logistic_regression"]

