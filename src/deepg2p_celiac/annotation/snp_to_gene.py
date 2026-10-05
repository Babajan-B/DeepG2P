from __future__ import annotations

from pathlib import Path

import pandas as pd


def load_snp_gene_map(panel_path: str | Path) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    required = {"rsid", "gene"}
    missing = required.difference(panel.columns)
    if missing:
        raise ValueError(f"Panel is missing required columns: {sorted(missing)}")
    return panel[["rsid", "gene", "category"]].dropna(subset=["rsid", "gene"]).drop_duplicates()


def aggregate_gene_features(genotypes: pd.DataFrame, snp_gene_map: pd.DataFrame) -> pd.DataFrame:
    if "sample_id" not in genotypes.columns:
        raise ValueError("Genotype matrix must contain 'sample_id'.")

    output = pd.DataFrame({"sample_id": genotypes["sample_id"]})
    for gene, group in snp_gene_map.groupby("gene"):
        snps = [rsid for rsid in group["rsid"] if rsid in genotypes.columns]
        if snps:
            output[gene] = genotypes[snps].mean(axis=1)
    return output

