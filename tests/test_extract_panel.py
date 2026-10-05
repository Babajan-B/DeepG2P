from pathlib import Path

from deepg2p_celiac.ingestion.extract_panel import extract_panel_matrix


def test_extract_panel_matrix(tmp_path: Path) -> None:
    panel = tmp_path / "panel.csv"
    panel.write_text("rsid,gene\nrs2187668,HLA-DQA1\nrs7454108,HLA-DQB1\n", encoding="utf-8")
    input_dir = tmp_path / "genotypes"
    input_dir.mkdir()
    (input_dir / "sample_1.txt").write_text(
        "rs2187668\t6\t32605884\tAG\nrs7454108\t6\t32632665\tTT\n",
        encoding="utf-8",
    )

    matrix = extract_panel_matrix(input_dir, panel)

    assert matrix.shape == (1, 3)
    assert matrix.loc[0, "sample_id"] == "sample_1"
    assert set(matrix.columns) == {"sample_id", "rs2187668", "rs7454108"}

