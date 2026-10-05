from __future__ import annotations

import argparse
from pathlib import Path

import networkx as nx
import pandas as pd

from deepg2p_celiac.annotation.snp_to_gene import aggregate_gene_features, load_snp_gene_map
from deepg2p_celiac.graph.build_graph import build_gene_graph


def export_graph_dataset(
    genotype_matrix_path: str | Path,
    panel_path: str | Path,
    output_dir: str | Path,
    edge_path: str | Path | None = None,
) -> dict:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)

    genotypes = pd.read_csv(genotype_matrix_path)
    snp_gene_map = load_snp_gene_map(panel_path)
    gene_features = aggregate_gene_features(genotypes, snp_gene_map)
    gene_features_path = output / "gene_features.csv"
    gene_features.to_csv(gene_features_path, index=False)

    edge_table = pd.read_csv(edge_path) if edge_path else None
    genes = [column for column in gene_features.columns if column != "sample_id"]
    graph = build_gene_graph(genes, edge_table)

    edges_path = output / "gene_edges.csv"
    pd.DataFrame(
        [
            {
                "gene_a": gene_a,
                "gene_b": gene_b,
                "source": data.get("source", "unknown"),
            }
            for gene_a, gene_b, data in graph.edges(data=True)
        ]
    ).to_csv(edges_path, index=False)

    adjacency = nx.to_pandas_adjacency(graph, nodelist=genes, dtype=int)
    adjacency_path = output / "gene_adjacency.csv"
    adjacency.to_csv(adjacency_path)

    graphml_path = output / "gene_graph.graphml"
    nx.write_graphml(graph, graphml_path)

    return {
        "gene_features": str(gene_features_path),
        "gene_edges": str(edges_path),
        "gene_adjacency": str(adjacency_path),
        "gene_graphml": str(graphml_path),
        "samples": int(len(gene_features)),
        "genes": int(len(genes)),
        "edges": int(graph.number_of_edges()),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Export celiac gene-level graph dataset files.")
    parser.add_argument("--genotypes", required=True, help="QC genotype matrix CSV.")
    parser.add_argument("--panel", required=True, help="Celiac SNP panel CSV.")
    parser.add_argument("--output-dir", required=True, help="Output graph dataset directory.")
    parser.add_argument("--edges", help="Optional CSV with gene_a,gene_b columns.")
    args = parser.parse_args()

    report = export_graph_dataset(args.genotypes, args.panel, args.output_dir, args.edges)
    for key, value in report.items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()

