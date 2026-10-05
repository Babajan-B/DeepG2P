from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


REQUIRED_PATHS = [
    "configs/celiac_snp_panel.csv",
    "configs/data_manifest.yaml",
    "configs/model.yaml",
    "configs/qc.yaml",
    "docs/data_policy.md",
    "docs/final_analysis_scope.md",
    "src/deepg2p_celiac",
]


def run_doctor(root: str | Path = ".") -> dict:
    project_root = Path(root)
    required = {
        path: (project_root / path).exists()
        for path in REQUIRED_PATHS
    }
    commands = {
        "plink": shutil.which("plink") is not None,
        "plink2": shutil.which("plink2") is not None,
    }
    data_dirs = {
        "raw": (project_root / "data/raw").exists(),
        "interim": (project_root / "data/interim").exists(),
        "processed": (project_root / "data/processed").exists(),
    }
    return {
        "project": "DeepG2P-Celiac",
        "required_paths": required,
        "data_dirs": data_dirs,
        "optional_commands": commands,
        "ready_for_open_evidence": all(required.values()),
        "ready_for_individual_validation": (project_root / "data/raw/genotypes").exists()
        and (project_root / "data/raw/phenotypes.csv").exists(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Check DeepG2P-Celiac project readiness.")
    parser.add_argument("--root", default=".")
    parser.add_argument("--output", help="Optional JSON output path.")
    args = parser.parse_args()

    report = run_doctor(args.root)
    text = json.dumps(report, indent=2)
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()

