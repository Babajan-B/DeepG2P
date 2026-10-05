"""Extract celiac disease associations from the full NHGRI-EBI GWAS Catalog association file.

Download the file first (about 60 MB zipped):
    curl -L -o data/raw/open_source/gwas-catalog-associations-full.zip \
        https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations-full.zip

Then run from the repository root:
    python scripts/extract_gwas_catalog.py

Rows are kept when DISEASE/TRAIT contains "celiac", "coeliac" or "gluten" (case-insensitive).
MAPPED_TRAIT is also searched when the file has that column (the ontology-annotated release). The manuscript used the release indexed 18 May 2026 (259 associations);
later releases may return more rows.
"""
from __future__ import annotations

import argparse
import io
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
PATTERN = r"celiac|coeliac|gluten"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--zip", type=Path, default=ROOT / "data/raw/open_source/gwas-catalog-associations-full.zip")
    ap.add_argument("--output", type=Path,
                    default=ROOT / "data/processed/open_source_evidence/gwas_catalog_celiac_associations.csv")
    args = ap.parse_args()

    kept = []
    with zipfile.ZipFile(args.zip) as zf:
        name = next(n for n in zf.namelist() if n.endswith(".tsv"))
        with zf.open(name) as handle:
            reader = pd.read_csv(io.TextIOWrapper(handle, encoding="utf-8"), sep="\t", dtype=str,
                                 chunksize=100_000, low_memory=False)
            for chunk in reader:
                hit = chunk["DISEASE/TRAIT"].fillna("").str.contains(PATTERN, case=False)
                if "MAPPED_TRAIT" in chunk:
                    hit |= chunk["MAPPED_TRAIT"].fillna("").str.contains(PATTERN, case=False)
                kept.append(chunk[hit])
    out = pd.concat(kept, ignore_index=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.output, index=False)
    print(f"{len(out)} associations from {out['PUBMEDID'].nunique()} studies -> {args.output}")


if __name__ == "__main__":
    main()
