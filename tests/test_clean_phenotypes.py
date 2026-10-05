from pathlib import Path

from deepg2p_celiac.ingestion.clean_phenotypes import clean_phenotype_table


def test_clean_phenotype_table_writes_clean_and_excluded_outputs(tmp_path: Path) -> None:
    raw = tmp_path / "phenotypes.csv"
    clean = tmp_path / "clean.csv"
    excluded = tmp_path / "excluded.csv"
    raw.write_text(
        "sample_id,celiac_status\n"
        "real_case,celiac disease\n"
        "real_control,no celiac\n"
        "real_excluded,gluten sensitivity\n",
        encoding="utf-8",
    )

    report = clean_phenotype_table(raw, "celiac_status", clean, excluded)

    assert report["cleaned_rows"] == 2
    assert report["excluded_rows"] == 1
    assert "real_case,case" in clean.read_text(encoding="utf-8")
    assert "real_excluded" in excluded.read_text(encoding="utf-8")

