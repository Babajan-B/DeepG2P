from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd


def _table_html(rows: list[dict], columns: list[str]) -> str:
    header = "".join(f"<th>{html.escape(column)}</th>" for column in columns)
    body_rows = []
    for row in rows:
        body_rows.append(
            "<tr>"
            + "".join(f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns)
            + "</tr>"
        )
    if not body_rows:
        body_rows.append(f"<tr><td colspan='{len(columns)}'>No rows available</td></tr>")
    return f"<table><thead><tr>{header}</tr></thead><tbody>{''.join(body_rows)}</tbody></table>"


def render_html_report(
    run_report_path: str | Path,
    output_path: str | Path,
    top_snps_path: str | Path | None = None,
    top_genes_path: str | Path | None = None,
) -> str:
    report = json.loads(Path(run_report_path).read_text(encoding="utf-8"))
    top_snps = pd.read_csv(top_snps_path).head(10).to_dict("records") if top_snps_path else []
    top_genes = pd.read_csv(top_genes_path).head(10).to_dict("records") if top_genes_path else []
    model_rows = report.get("models", [])
    graph_model_rows = report.get("graph_models", [])

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>DeepG2P-Celiac Research Report</title>
  <style>
    body {{
      background: #10130f;
      color: #f4ead2;
      font-family: Menlo, Consolas, monospace;
      margin: 0;
      padding: 32px;
    }}
    main {{
      max-width: 1120px;
      margin: 0 auto;
    }}
    h1, h2 {{
      color: #f0b84f;
      letter-spacing: 0;
    }}
    section {{
      border-top: 1px solid #5b4b2e;
      padding: 20px 0;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      margin: 12px 0;
      background: #171b15;
    }}
    th, td {{
      border: 1px solid #3a3f33;
      padding: 8px;
      text-align: left;
      vertical-align: top;
      font-size: 13px;
    }}
    th {{
      color: #f0b84f;
      background: #202419;
    }}
    .warning {{
      color: #ffcf7a;
      border: 1px solid #6f5523;
      padding: 12px;
      background: #1e1a10;
    }}
    code {{
      color: #9be58f;
    }}
  </style>
</head>
<body>
<main>
  <h1>DeepG2P-Celiac Research Report</h1>
  <p class="warning">Research-only output. This report is not a medical diagnosis.</p>
  <section>
    <h2>Run Metadata</h2>
    <p><strong>Created:</strong> {html.escape(str(report.get("created_at_utc", "")))}</p>
    <p><strong>Project:</strong> {html.escape(str(report.get("project", "")))}</p>
    <p><strong>Synthetic data allowed:</strong> {html.escape(str(report.get("synthetic_data_allowed", "")))}</p>
  </section>
  <section>
    <h2>Data Summary</h2>
    <p><strong>Raw matrix:</strong> {html.escape(str(report.get("raw_matrix", {})))}</p>
    <p><strong>QC matrix:</strong> {html.escape(str(report.get("qc_matrix", {})))}</p>
    <p><strong>Graph:</strong> {html.escape(str(report.get("graph", {})))}</p>
  </section>
  <section>
    <h2>Baseline Models</h2>
    {_table_html(model_rows, ["model_name", "balanced_accuracy", "auroc"])}
  </section>
  <section>
    <h2>Graph Models</h2>
    {_table_html(graph_model_rows, ["model_name", "balanced_accuracy", "auroc"])}
  </section>
  <section>
    <h2>Top SNP Features</h2>
    {_table_html(top_snps, ["feature", "importance", "rsid", "gene", "category", "notes"])}
  </section>
  <section>
    <h2>Top Gene Features</h2>
    {_table_html(top_genes, ["feature", "importance"])}
  </section>
</main>
</body>
</html>
"""
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(html_text, encoding="utf-8")
    return html_text


def main() -> None:
    parser = argparse.ArgumentParser(description="Render an HTML DeepG2P-Celiac research report.")
    parser.add_argument("--run-report", required=True, help="run_report.json path.")
    parser.add_argument("--output", required=True, help="HTML output path.")
    parser.add_argument("--top-snps", help="Optional top SNP feature CSV.")
    parser.add_argument("--top-genes", help="Optional top gene feature CSV.")
    args = parser.parse_args()

    render_html_report(args.run_report, args.output, args.top_snps, args.top_genes)


if __name__ == "__main__":
    main()
