from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

from deepg2p_celiac.graph.build_graph import build_gene_graph

PVALUE_COLUMNS = ["p_value", "pvalue", "p-value", "P-VALUE", "P-VALUE_MLOG"]
RSID_COLUMNS = ["rsid", "SNPS", "SNP_ID_CURRENT", "variant_id", "variant"]
GENE_COLUMNS = ["gene", "MAPPED_GENE", "REPORTED GENE(S)", "mapped_gene"]
TRAIT_COLUMNS = ["trait", "DISEASE/TRAIT", "efo_trait", "disease_trait"]


def _first_existing(columns: list[str], candidates: list[str]) -> str | None:
    lower_map = {column.lower(): column for column in columns}
    for candidate in candidates:
        if candidate in columns:
            return candidate
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]
    return None


def _normalize_pvalue(value: object, column: str | None = None) -> float | None:
    parsed = pd.to_numeric(value, errors="coerce")
    if pd.isna(parsed):
        return None
    number = float(parsed)
    if column == "P-VALUE_MLOG":
        return 10 ** (-number)
    return number


def split_gene_field(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value)
    for separator in [";", ",", " - ", "/", "|"]:
        text = text.replace(separator, ",")
    genes = []
    for item in text.split(","):
        gene = item.strip()
        if gene and gene.upper() not in {"NR", "INTERGENIC", "PENDING"}:
            genes.append(gene)
    return sorted(set(genes))


def normalize_association_table(path: str | Path, trait_filter: str = "celiac") -> pd.DataFrame:
    table = pd.read_csv(path)
    rsid_col = _first_existing(list(table.columns), RSID_COLUMNS)
    gene_col = _first_existing(list(table.columns), GENE_COLUMNS)
    pvalue_col = _first_existing(list(table.columns), PVALUE_COLUMNS)
    trait_col = _first_existing(list(table.columns), TRAIT_COLUMNS)
    if not gene_col:
        raise ValueError("Association table must contain a gene column.")

    rows = []
    for record in table.to_dict("records"):
        trait = str(record.get(trait_col, "")) if trait_col else ""
        if trait_filter and trait and trait_filter.lower() not in trait.lower():
            continue
        p_value = _normalize_pvalue(record.get(pvalue_col), pvalue_col) if pvalue_col else None
        for gene in split_gene_field(record.get(gene_col)):
            rows.append(
                {
                    "rsid": str(record.get(rsid_col, "")).strip() if rsid_col else "",
                    "gene": gene,
                    "trait": trait,
                    "p_value": p_value,
                    "source_file": str(path),
                }
            )
    return pd.DataFrame(rows).drop_duplicates()


def score_gene_evidence(associations: pd.DataFrame) -> pd.DataFrame:
    if associations.empty:
        return pd.DataFrame(columns=["gene", "association_count", "best_p_value", "evidence_score"])

    working = associations.copy()
    working["p_value"] = pd.to_numeric(working["p_value"], errors="coerce")
    grouped = working.groupby("gene", dropna=True)
    rows = []
    for gene, group in grouped:
        best_p = group["p_value"].dropna().min()
        if pd.isna(best_p) or best_p <= 0:
            p_score = 0.0
        else:
            p_score = min(-np.log10(best_p), 300)
        count = int(len(group))
        rows.append(
            {
                "gene": gene,
                "association_count": count,
                "best_p_value": float(best_p) if not pd.isna(best_p) else None,
                "evidence_score": float(p_score + count),
            }
        )
    return pd.DataFrame(rows).sort_values("evidence_score", ascending=False)


def add_graph_context(gene_scores: pd.DataFrame) -> pd.DataFrame:
    genes = gene_scores["gene"].dropna().astype(str).tolist()
    graph = build_gene_graph(genes)
    centrality = nx.degree_centrality(graph) if graph.number_of_nodes() else {}
    result = gene_scores.copy()
    result["graph_degree"] = result["gene"].map(lambda gene: int(graph.degree(gene)) if gene in graph else 0)
    result["degree_centrality"] = result["gene"].map(lambda gene: float(centrality.get(gene, 0.0)))
    result["combined_evidence_graph_score"] = (
        result["evidence_score"] * (1 + result["degree_centrality"])
    )
    return result.sort_values("combined_evidence_graph_score", ascending=False)


def build_evidence_report(
    association_files: list[str | Path],
    output_dir: str | Path,
    trait_filter: str = "celiac",
) -> dict:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    normalized_tables = [
        normalize_association_table(path, trait_filter=trait_filter) for path in association_files
    ]
    associations = (
        pd.concat(normalized_tables, ignore_index=True).drop_duplicates()
        if normalized_tables
        else pd.DataFrame()
    )
    associations_path = output / "normalized_celiac_associations.csv"
    associations.to_csv(associations_path, index=False)

    gene_scores = score_gene_evidence(associations)
    graph_scores = add_graph_context(gene_scores) if not gene_scores.empty else gene_scores
    gene_scores_path = output / "gene_evidence_scores.csv"
    graph_scores.to_csv(gene_scores_path, index=False)

    report = {
        "input_files": [str(path) for path in association_files],
        "trait_filter": trait_filter,
        "association_rows": int(len(associations)),
        "genes": int(graph_scores["gene"].nunique()) if not graph_scores.empty else 0,
        "top_genes": graph_scores.head(20).to_dict("records") if not graph_scores.empty else [],
        "outputs": {
            "associations": str(associations_path),
            "gene_evidence_scores": str(gene_scores_path),
        },
        "interpretation_warning": (
            "Open-source association evidence ranks known loci; it does not validate predictive "
            "performance without individual-level genotype and phenotype data."
        ),
    }
    (output / "open_source_evidence_report.json").write_text(
        json.dumps(report, indent=2),
        encoding="utf-8",
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Build celiac evidence report from open association files.")
    parser.add_argument("--associations", nargs="+", required=True, help="GWAS/open association CSV files.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--trait-filter", default="celiac")
    args = parser.parse_args()

    report = build_evidence_report(args.associations, args.output_dir, args.trait_filter)
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
