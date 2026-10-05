from __future__ import annotations

import pandas as pd


def filter_missingness(
    matrix: pd.DataFrame,
    variant_missingness_max: float = 0.05,
    sample_missingness_max: float = 0.10,
) -> pd.DataFrame:
    if "sample_id" not in matrix.columns:
        raise ValueError("Genotype matrix must contain 'sample_id'.")

    snp_columns = [column for column in matrix.columns if column != "sample_id"]
    variant_missingness = matrix[snp_columns].isna().mean(axis=0)
    kept_snps = variant_missingness[variant_missingness <= variant_missingness_max].index.tolist()
    filtered = matrix[["sample_id", *kept_snps]].copy()

    sample_missingness = filtered[kept_snps].isna().mean(axis=1) if kept_snps else pd.Series(1, index=filtered.index)
    filtered = filtered.loc[sample_missingness <= sample_missingness_max].copy()
    return filtered


def impute_panel_mean(matrix: pd.DataFrame) -> pd.DataFrame:
    result = matrix.copy()
    snp_columns = [column for column in result.columns if column != "sample_id"]
    for column in snp_columns:
        result[column] = pd.to_numeric(result[column], errors="coerce")
        result[column] = result[column].fillna(result[column].mean())
    return result

