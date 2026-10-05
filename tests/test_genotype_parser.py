from pathlib import Path

from deepg2p_celiac.ingestion.genotype_parser import encode_genotype, read_consumer_genotype


def test_read_consumer_genotype_ignores_comments(tmp_path: Path) -> None:
    genotype = tmp_path / "sample.txt"
    genotype.write_text(
        "# comment\n"
        "rs2187668\t6\t32605884\tAG\n"
        "rs7454108 6 32632665 TT\n",
        encoding="utf-8",
    )

    df = read_consumer_genotype(genotype)

    assert list(df["rsid"]) == ["rs2187668", "rs7454108"]
    assert list(df["genotype"]) == ["AG", "TT"]


def test_encode_genotype_counts_effect_allele() -> None:
    assert encode_genotype("AG", "A") == 1.0
    assert encode_genotype("AA", "A") == 2.0
    assert encode_genotype("--", "A") is None

