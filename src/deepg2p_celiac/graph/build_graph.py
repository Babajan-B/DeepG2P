from __future__ import annotations

import argparse
from itertools import combinations
from pathlib import Path

import networkx as nx
import pandas as pd


DEFAULT_IMMUNE_EDGES = [
    ("HLA-DQA1", "HLA-DQB1"),
    ("HLA-DQA1", "HLA-DRB1"),
    ("HLA-DQB1", "HLA-DRB1"),
    ("IL2", "IL18RAP"),
    ("IL2", "SH2B3"),
    ("CCR3", "IL18RAP"),
    ("ATXN2", "SH2B3"),
]


def build_gene_graph(genes: list[str], edge_table: pd.DataFrame | None = None) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(sorted(set(genes)))

    if edge_table is not None and {"gene_a", "gene_b"}.issubset(edge_table.columns):
        for record in edge_table.to_dict("records"):
            if record["gene_a"] in graph and record["gene_b"] in graph:
                graph.add_edge(record["gene_a"], record["gene_b"], source="external")
    else:
        for gene_a, gene_b in DEFAULT_IMMUNE_EDGES:
            if gene_a in graph and gene_b in graph:
                graph.add_edge(gene_a, gene_b, source="curated_seed")

    if graph.number_of_edges() == 0:
        for gene_a, gene_b in combinations(graph.nodes, 2):
            graph.add_edge(gene_a, gene_b, source="fallback_complete")
    return graph


def main() -> None:
    parser = argparse.ArgumentParser(description="Build a celiac immune gene graph.")
    parser.add_argument("--panel", required=True, help="Celiac SNP panel CSV.")
    parser.add_argument("--edges", help="Optional CSV with gene_a,gene_b columns.")
    parser.add_argument("--output", required=True, help="GraphML output path.")
    args = parser.parse_args()

    panel = pd.read_csv(args.panel)
    genes = sorted(panel["gene"].dropna().unique())
    edges = pd.read_csv(args.edges) if args.edges else None
    graph = build_gene_graph(genes, edges)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    nx.write_graphml(graph, output)


if __name__ == "__main__":
    main()

