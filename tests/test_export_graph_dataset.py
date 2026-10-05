from pathlib import Path

import pandas as pd

from deepg2p_celiac.graph.export_graph_dataset import export_graph_dataset


def test_export_graph_dataset_writes_artifacts(tmp_path: Path) -> None:
    genotypes = tmp_path / "genotypes.csv"
    panel = tmp_path / "panel.csv"
    output = tmp_path / "graph"
    pd.DataFrame(
        {
            "sample_id": ["real_a", "real_b"],
            "rs1": [0, 2],
            "rs2": [1, 1],
        }
    ).to_csv(genotypes, index=False)
    pd.DataFrame(
        {
            "rsid": ["rs1", "rs2"],
            "gene": ["HLA-DQA1", "HLA-DQB1"],
            "category": ["HLA_tag", "HLA_tag"],
        }
    ).to_csv(panel, index=False)

    report = export_graph_dataset(genotypes, panel, output)

    assert report["samples"] == 2
    assert report["genes"] == 2
    assert Path(report["gene_features"]).exists()
    assert Path(report["gene_adjacency"]).exists()

