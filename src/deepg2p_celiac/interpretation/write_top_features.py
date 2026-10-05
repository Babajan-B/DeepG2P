from __future__ import annotations

import argparse

from deepg2p_celiac.interpretation.top_snps import write_top_features


def main() -> None:
    parser = argparse.ArgumentParser(description="Write top model feature attributions.")
    parser.add_argument("--model", required=True, help="Model bundle joblib path.")
    parser.add_argument("--output", required=True, help="Output CSV path.")
    parser.add_argument("--panel", help="Optional SNP panel CSV for annotation.")
    parser.add_argument("--limit", type=int, default=10)
    args = parser.parse_args()

    write_top_features(args.model, args.output, args.panel, args.limit)


if __name__ == "__main__":
    main()

