from pathlib import Path

import pandas as pd

from deepg2p_celiac.analysis.final_report import render_final_analysis_report
from deepg2p_celiac.analysis.update_panel import update_panel_from_evidence
from deepg2p_celiac.project.doctor import run_doctor


def test_update_panel_from_evidence(tmp_path: Path) -> None:
    panel = tmp_path / "panel.csv"
    evidence = tmp_path / "evidence.csv"
    output = tmp_path / "updated.csv"
    pd.DataFrame(
        {
            "rsid": ["rs1"],
            "gene": ["HLA-DQA1"],
            "chromosome": ["6"],
            "position_hg19": ["1"],
            "category": ["HLA"],
            "notes": [""],
        }
    ).to_csv(panel, index=False)
    pd.DataFrame(
        {
            "gene": ["HLA-DQA1", "IL2"],
            "evidence_score": [10, 8],
            "combined_evidence_graph_score": [11, 9],
        }
    ).to_csv(evidence, index=False)

    updated = update_panel_from_evidence(panel, evidence, output, max_new_genes=1)

    assert output.exists()
    assert "IL2" in set(updated["gene"])


def test_render_final_analysis_report(tmp_path: Path) -> None:
    evidence_report = tmp_path / "evidence.json"
    evidence_scores = tmp_path / "scores.csv"
    output = tmp_path / "final.html"
    evidence_report.write_text('{"association_rows": 2, "genes": 1}', encoding="utf-8")
    pd.DataFrame([{"gene": "HLA-DQA1", "evidence_score": 10}]).to_csv(
        evidence_scores,
        index=False,
    )

    html = render_final_analysis_report(output, evidence_report, evidence_scores)

    assert output.exists()
    assert "DeepG2P-Celiac Final Analysis" in html
    assert "HLA-DQA1" in html


def test_run_doctor() -> None:
    report = run_doctor(".")

    assert report["project"] == "DeepG2P-Celiac"
    assert report["ready_for_open_evidence"] is True

