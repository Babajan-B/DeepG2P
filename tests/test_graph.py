from deepg2p_celiac.graph.build_graph import build_gene_graph


def test_build_gene_graph_uses_seed_edges() -> None:
    graph = build_gene_graph(["HLA-DQA1", "HLA-DQB1", "IL2"])

    assert "HLA-DQA1" in graph
    assert graph.has_edge("HLA-DQA1", "HLA-DQB1")

