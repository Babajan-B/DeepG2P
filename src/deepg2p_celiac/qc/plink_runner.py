from __future__ import annotations

import argparse
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PlinkCommand:
    executable: str
    args: list[str]

    def as_list(self) -> list[str]:
        return [self.executable, *self.args]

    def as_shell_string(self) -> str:
        return " ".join(self.as_list())


def build_plink_qc_command(
    input_prefix: str | Path,
    output_prefix: str | Path,
    geno: float = 0.05,
    mind: float = 0.10,
    maf: float = 0.01,
    hwe: float = 1e-6,
    executable: str = "plink",
) -> PlinkCommand:
    return PlinkCommand(
        executable=executable,
        args=[
            "--bfile",
            str(input_prefix),
            "--geno",
            str(geno),
            "--mind",
            str(mind),
            "--maf",
            str(maf),
            "--hwe",
            str(hwe),
            "--make-bed",
            "--out",
            str(output_prefix),
        ],
    )


def run_plink_qc(command: PlinkCommand, dry_run: bool = False) -> subprocess.CompletedProcess | None:
    if dry_run:
        return None
    if shutil.which(command.executable) is None:
        raise FileNotFoundError(
            f"{command.executable!r} was not found on PATH. Install PLINK or use --dry-run."
        )
    return subprocess.run(command.as_list(), check=True, capture_output=True, text=True)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build or run a PLINK QC command.")
    parser.add_argument("--input-prefix", required=True, help="Input PLINK bed/bim/fam prefix.")
    parser.add_argument("--output-prefix", required=True, help="Output PLINK prefix.")
    parser.add_argument("--geno", type=float, default=0.05)
    parser.add_argument("--mind", type=float, default=0.10)
    parser.add_argument("--maf", type=float, default=0.01)
    parser.add_argument("--hwe", type=float, default=1e-6)
    parser.add_argument("--executable", default="plink")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    command = build_plink_qc_command(
        args.input_prefix,
        args.output_prefix,
        args.geno,
        args.mind,
        args.maf,
        args.hwe,
        args.executable,
    )
    print(command.as_shell_string())
    run_plink_qc(command, args.dry_run)


if __name__ == "__main__":
    main()

