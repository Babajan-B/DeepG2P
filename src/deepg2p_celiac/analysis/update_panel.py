from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def update_panel_from_evidence(
    panel_path: str | Path,
    evidence_scores_path: str | Path,
    output_path: str | Path,
    max_new_genes: int = 25,
) -> pd.DataFrame:
    panel = pd.read_csv(panel_path)
    evidence = pd.read_csv(evidence_scores_path)
    if "gene" not in panel.columns or "gene" not in evidence.columns:
        raise ValueError("Both panel and evidence files must contain a 'gene' column.")

    panel_genes = set(panel["gene"].dropna().astype(str))
    evidence = evidence.sort_values(
        "combined_evidence_graph_score"
        if "combined_evidence_graph_score" in evidence.columns
        else "evidence_score",
        ascending=False,
    )
    new_rows = []
    for record in evidence.to_dict("records"):
        gene = str(record["gene"])
        if gene in panel_genes:
            continue
        new_rows.append(
            {
                "rsid": "",
                "gene": gene,
                "chromosome": "",
                "position_hg19": "",
                "category": "open_evidence_gene",
                "notes": (
                    "Gene prioritized from open-source association evidence; "
                    "requires SNP-level curation before genotype extraction."
                ),
            }
        )
        panel_genes.add(gene)
        if len(new_rows) >= max_new_genes:
            break

    updated = pd.concat([panel, pd.DataFrame(new_rows)], ignore_index=True)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    updated.to_csv(output, index=False)
    return updated


def main() -> None:
    parser = argparse.ArgumentParser(description="Create a panel candidate file from evidence scores.")
    parser.add_argument("--panel", default="configs/celiac_snp_panel.csv")
    parser.add_argument("--evidence-scores", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--max-new-genes", type=int, default=25)
    args = parser.parse_args()

    updated = update_panel_from_evidence(
        args.panel,
        args.evidence_scores,
        args.output,
        args.max_new_genes,
    )
    print(f"Wrote {len(updated)} rows to {args.output}")


if __name__ == "__main__":
    main()

