from pathlib import Path

import pytest

from deepg2p_celiac.ingestion.validate_data import (
    validate_genotype_dir,
    validate_phenotype_file,
    validate_storage,
)


def test_validate_phenotype_rejects_fixture_like_sample_id(tmp_path: Path) -> None:
    phenotypes = tmp_path / "phenotypes.csv"
    phenotypes.write_text("sample_id,celiac_status\nsynthetic_001,case\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Suspicious sample IDs"):
        validate_phenotype_file(phenotypes)


def test_validate_genotype_dir_rejects_fixture_like_filename(tmp_path: Path) -> None:
    genotype_dir = tmp_path / "genotypes"
    genotype_dir.mkdir()
    (genotype_dir / "mock_sample.txt").write_text("rs1\t1\t1\tAA\n", encoding="utf-8")

    with pytest.raises(ValueError, match="Suspicious genotype filenames"):
        validate_genotype_dir(genotype_dir)


def test_validate_storage_reports_budget(tmp_path: Path) -> None:
    (tmp_path / "real_sample.txt").write_text("rs1\t1\t1\tAA\n", encoding="utf-8")

    report = validate_storage([tmp_path], max_gb=6)

    assert report["storage_limit_gb"] == 6
    assert report["storage_used_gb"] >= 0

