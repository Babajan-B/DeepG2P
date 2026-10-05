from __future__ import annotations

import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from deepg2p_celiac.ingestion.validate_data import directory_size_bytes


def write_json_report(report: dict, output_path: str | Path) -> None:
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2), encoding="utf-8")


def build_run_metadata(
    genotype_dir: str | Path,
    phenotype_file: str | Path,
    panel_file: str | Path,
    output_dir: str | Path,
) -> dict:
    return {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "project": "DeepG2P-Celiac",
        "synthetic_data_allowed": False,
        "clinical_use_allowed": False,
        "python_version": platform.python_version(),
        "inputs": {
            "genotype_dir": str(genotype_dir),
            "phenotype_file": str(phenotype_file),
            "panel_file": str(panel_file),
        },
        "outputs": {
            "output_dir": str(output_dir),
        },
    }


def summarize_matrix(path: str | Path) -> dict:
    df = pd.read_csv(path)
    snp_columns = [column for column in df.columns if column != "sample_id"]
    return {
        "path": str(path),
        "samples": int(len(df)),
        "snps": int(len(snp_columns)),
        "missing_rate": float(df[snp_columns].isna().mean().mean()) if snp_columns else None,
    }


def summarize_storage(paths: list[str | Path]) -> dict:
    return {
        str(path): round(directory_size_bytes(path) / 1024 / 1024, 3)
        for path in paths
    }

