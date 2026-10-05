from pathlib import Path

import pandas as pd

from deepg2p_celiac.analysis.evidence import (
    build_evidence_report,
    normalize_association_table,
    score_gene_evidence,
    split_gene_field,
)


def test_split_gene_field() -> None:
    assert split_gene_field("HLA-DQA1, HLA-DQB1; IL2") == ["HLA-DQA1", "HLA-DQB1", "IL2"]


def test_normalize_association_table_filters_trait(tmp_path: Path) -> None:
    associations = tmp_path / "associations.csv"
    pd.DataFrame(
        {
            "SNPS": ["rs1", "rs2"],
            "MAPPED_GENE": ["HLA-DQA1", "NOTCEL"],
            "DISEASE/TRAIT": ["celiac disease", "other disease"],
            "P-VALUE": [1e-8, 1e-4],
        }
    ).to_csv(associations, index=False)

    normalized = normalize_association_table(associations)

    assert list(normalized["gene"]) == ["HLA-DQA1"]
    assert normalized.loc[0, "p_value"] == 1e-8


def test_score_gene_evidence() -> None:
    table = pd.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "gene": ["HLA-DQA1", "HLA-DQA1"],
            "trait": ["celiac", "celiac"],
            "p_value": [1e-8, 1e-6],
        }
    )

    scores = score_gene_evidence(table)

    assert scores.loc[0, "gene"] == "HLA-DQA1"
    assert scores.loc[0, "association_count"] == 2


def test_build_evidence_report(tmp_path: Path) -> None:
    associations = tmp_path / "associations.csv"
    output = tmp_path / "evidence"
    pd.DataFrame(
        {
            "SNPS": ["rs1", "rs2"],
            "MAPPED_GENE": ["HLA-DQA1", "HLA-DQB1"],
            "DISEASE/TRAIT": ["celiac disease", "celiac disease"],
            "P-VALUE": [1e-8, 1e-7],
        }
    ).to_csv(associations, index=False)

    report = build_evidence_report([associations], output)

    assert report["association_rows"] == 2
    assert report["genes"] == 2
    assert (output / "open_source_evidence_report.json").exists()

