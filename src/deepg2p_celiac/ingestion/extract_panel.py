from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

from deepg2p_celiac.ingestion.genotype_parser import (
    encode_genotype,
    infer_sample_id,
    read_consumer_genotype,
)
from deepg2p_celiac.ingestion.validate_data import validate_genotype_dir


def load_panel(panel_path: str | Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    if "rsid" not in panel.columns:
        raise ValueError("Panel CSV must contain an 'rsid' column.")
    panel["rsid"] = panel["rsid"].astype(str)
    return panel


def extract_sample_panel(genotype_path: str | Path, panel: pd.DataFrame) -> dict[str, float | str | None]:
    genotype_file = Path(genotype_path)
    sample_id = infer_sample_id(genotype_file)
    genotypes = read_consumer_genotype(genotype_file)
    target = genotypes[genotypes["rsid"].isin(set(panel["rsid"]))]

    row: dict[str, float | str | None] = {"sample_id": sample_id}
    for rsid in panel["rsid"]:
        calls = target.loc[target["rsid"] == rsid, "genotype"]
        row[rsid] = encode_genotype(calls.iloc[0]) if not calls.empty else None
    return row


def extract_panel_matrix(input_dir: str | Path, panel_path: str | Path) -> pd.DataFrame:
    panel = load_panel(panel_path)
    paths = sorted(Path(input_dir).glob("*"))
    files = [path for path in paths if path.is_file() and not path.name.startswith(".")]
    rows = [extract_sample_panel(path, panel) for path in files]
    if not rows:
        columns = ["sample_id", *panel["rsid"].tolist()]
        return pd.DataFrame(columns=columns)
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description="Extract a celiac SNP panel from genotype files.")
    parser.add_argument("--input-dir", required=True, help="Directory of raw genotype text files.")
    parser.add_argument("--panel", required=True, help="CSV file containing target rsids.")
    parser.add_argument("--output", required=True, help="Output genotype matrix CSV.")
    parser.add_argument(
        "--allow-test-fixtures",
        action="store_true",
        help="Allow fixture-like filenames. Intended only for software tests, not real runs.",
    )
    args = parser.parse_args()

    if not args.allow_test_fixtures:
        validate_genotype_dir(args.input_dir)

    matrix = extract_panel_matrix(args.input_dir, args.panel)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    matrix.to_csv(output, index=False)


if __name__ == "__main__":
    main()
