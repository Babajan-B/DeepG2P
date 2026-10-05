from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import pandas as pd


def _open_text(path: str | Path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def parse_ieu_summary_vcf(
    vcf_path: str | Path,
    output_path: str | Path,
    top_n: int = 1000,
    panel_path: str | Path | None = None,
    panel_overlap_path: str | Path | None = None,
) -> dict:
    rows = []
    total = 0
    with _open_text(vcf_path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            total += 1
            chrom, pos, variant_id, ref, alt, _qual, _filter, _info, fmt, sample = line.rstrip().split(
                "\t"
            )
            keys = fmt.split(":")
            values = sample.split(":")
            parsed = dict(zip(keys, values, strict=False))
            lp = pd.to_numeric(parsed.get("LP"), errors="coerce")
            p_value = 10 ** (-float(lp)) if not pd.isna(lp) else None
            rows.append(
                {
                    "chromosome": chrom,
                    "position": int(pos),
                    "rsid": variant_id,
                    "ref": ref,
                    "alt": alt,
                    "effect_size": pd.to_numeric(parsed.get("ES"), errors="coerce"),
                    "standard_error": pd.to_numeric(parsed.get("SE"), errors="coerce"),
                    "p_value": p_value,
                    "minus_log10_p": float(lp) if not pd.isna(lp) else None,
                    "sample_size": pd.to_numeric(parsed.get("SS"), errors="coerce"),
                    "study": "ieu-a-276",
                    "trait": "Celiac disease",
                }
            )

    table = pd.DataFrame(rows)
    table = table.sort_values("p_value", na_position="last").head(top_n)
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(output, index=False)

    overlap_count = None
    if panel_path and panel_overlap_path:
        panel = pd.read_csv(panel_path)
        overlap = table.merge(panel, on="rsid", how="inner", suffixes=("", "_panel"))
        panel_output = Path(panel_overlap_path)
        panel_output.parent.mkdir(parents=True, exist_ok=True)
        overlap.to_csv(panel_output, index=False)
        overlap_count = int(len(overlap))

    report = {
        "vcf": str(vcf_path),
        "total_variants_scanned": total,
        "top_variants_written": int(len(table)),
        "output": str(output),
        "panel_overlap_output": str(panel_overlap_path) if panel_overlap_path else None,
        "panel_overlap_count": overlap_count,
        "warning": "Summary-statistics VCF; not individual-level genotype data.",
    }
    report_path = output.with_suffix(".report.json")
    report_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert IEU OpenGWAS summary-stat VCF to CSV.")
    parser.add_argument("--vcf", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--top-n", type=int, default=1000)
    parser.add_argument("--panel")
    parser.add_argument("--panel-overlap-output")
    args = parser.parse_args()

    report = parse_ieu_summary_vcf(
        args.vcf,
        args.output,
        args.top_n,
        args.panel,
        args.panel_overlap_output,
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()

