from deepg2p_celiac.ingestion.phenotype_loader import normalize_celiac_label


def test_normalize_celiac_label() -> None:
    assert normalize_celiac_label("celiac disease") == "case"
    assert normalize_celiac_label("No celiac") == "control"
    assert normalize_celiac_label("gluten sensitivity") is None

