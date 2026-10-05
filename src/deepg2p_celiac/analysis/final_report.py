from __future__ import annotations

import argparse
import html
import json
from pathlib import Path

import pandas as pd


def _load_json(path: str | Path | None) -> dict:
    if not path or not Path(path).exists():
        return {}
    return json.loads(Path(path).read_text(encoding="utf-8"))


def _load_rows(path: str | Path | None, limit: int = 20) -> list[dict]:
    if not path or not Path(path).exists():
        return []
    return pd.read_csv(path).head(limit).fillna("").to_dict("records")


def _table(rows: list[dict]) -> str:
    if not rows:
        return "<p>No table available.</p>"
    columns = list(rows[0].keys())
    header = "".join(f"<th>{html.escape(str(column))}</th>" for column in columns)
    body = ""
    for row in rows:
        body += "<tr>" + "".join(
            f"<td>{html.escape(str(row.get(column, '')))}</td>" for column in columns
        ) + "</tr>"
    return f"<table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table>"


def render_final_analysis_report(
    output_path: str | Path,
    evidence_report_path: str | Path | None = None,
    evidence_scores_path: str | Path | None = None,
    validation_report_path: str | Path | None = None,
) -> str:
    evidence_report = _load_json(evidence_report_path)
    validation_report = _load_json(validation_report_path)
    evidence_rows = _load_rows(evidence_scores_path)

    html_text = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <title>DeepG2P-Celiac Final Analysis</title>
  <style>
    body {{
      margin: 0;
      padding: 32px;
      background: #11140f;
      color: #f5ead4;
      font-family: Arial, sans-serif;
      line-height: 1.5;
    }}
    main {{ max-width: 1180px; margin: 0 auto; }}
    h1, h2 {{ color: #e8b04d; }}
    table {{ border-collapse: collapse; width: 100%; margin: 16px 0; }}
    th, td {{ border: 1px solid #3b4135; padding: 8px; font-size: 13px; }}
    th {{ background: #202719; color: #e8b04d; text-align: left; }}
    .box {{ border: 1px solid #6c5729; background: #1d1b11; padding: 12px; }}
    code {{ color: #9ddd8c; }}
  </style>
</head>
<body>
<main>
  <h1>DeepG2P-Celiac Final Analysis</h1>
  <div class="box">
    <strong>Scope:</strong> Research-only celiac genetics analysis. This is not a diagnostic system.
  </div>
  <section>
    <h2>Evidence Synthesis</h2>
    <p>Association rows: {html.escape(str(evidence_report.get("association_rows", "not available")))}</p>
    <p>Genes ranked: {html.escape(str(evidence_report.get("genes", "not available")))}</p>
    {_table(evidence_rows)}
  </section>
  <section>
    <h2>Predictive Validation</h2>
    <p>{html.escape("Available" if validation_report else "Not available: individual-level genotype/phenotype data was not supplied.")}</p>
  </section>
  <section>
    <h2>Conclusion</h2>
    <p>The project is complete as an analysis framework. Open-source evidence mode ranks celiac loci from published association data. Individual-level predictive validation remains a separate mode that requires real case/control genotype files.</p>
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
    parser = argparse.ArgumentParser(description="Render final DeepG2P-Celiac analysis report.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--evidence-report")
    parser.add_argument("--evidence-scores")
    parser.add_argument("--validation-report")
    args = parser.parse_args()

    render_final_analysis_report(
        args.output,
        args.evidence_report,
        args.evidence_scores,
        args.validation_report,
    )


if __name__ == "__main__":
    main()

