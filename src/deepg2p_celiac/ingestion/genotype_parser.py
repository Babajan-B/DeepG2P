from __future__ import annotations

from pathlib import Path

import pandas as pd

REQUIRED_COLUMNS = ["rsid", "chromosome", "position", "genotype"]


def infer_sample_id(path: Path) -> str:
    return path.name.split(".")[0]


def read_consumer_genotype(path: str | Path) -> pd.DataFrame:
    """Read a 23andMe-like raw genotype file.

    The parser supports comment-prefixed files and either whitespace or tab separators.
    """

    genotype_path = Path(path)
    rows: list[list[str]] = []
    with genotype_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split()
            if len(parts) < 4:
                continue
            rows.append(parts[:4])

    df = pd.DataFrame(rows, columns=REQUIRED_COLUMNS)
    if df.empty:
        return pd.DataFrame(columns=REQUIRED_COLUMNS)

    df["rsid"] = df["rsid"].astype(str)
    df["chromosome"] = df["chromosome"].astype(str)
    df["position"] = pd.to_numeric(df["position"], errors="coerce").astype("Int64")
    df["genotype"] = df["genotype"].astype(str).str.upper()
    return df.dropna(subset=["position"])


def encode_genotype(genotype: str, effect_allele: str | None = None) -> float | None:
    """Encode genotype as allele count.

    If an effect allele is supplied, count that allele. Otherwise, encode homozygous
    reference-like calls as 0, heterozygous calls as 1, and alternate homozygous calls as 2
    using alphabetical allele ordering as a deterministic fallback.
    """

    call = genotype.strip().upper()
    if not call or call in {"--", "00", "NN"}:
        return None
    if len(call) != 2:
        return None

    if effect_allele:
        allele = effect_allele.upper()
        return float(call.count(allele))

    alleles = sorted(call)
    if alleles[0] == alleles[1]:
        return 0.0 if alleles[0] in {"A", "C"} else 2.0
    return 1.0

