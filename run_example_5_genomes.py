#!/usr/bin/env python3
"""
run_example_5_genomes.py

Example runner for the Layer-1 synthetic MAG simulator.

First download genomes:
    python download_reference_genomes.py --accessions examples/accessions_5_genomes.tsv --outdir .

Then run:
    python run_example_5_genomes.py
"""

import subprocess
import sys
from pathlib import Path


def main() -> None:
    root = Path(__file__).resolve().parent

    cmd = [
        sys.executable,
        str(root / "synthetic_mag_layer1.py"),
        "--metadata", str(root / "metadata.tsv"),
        "--outdir", str(root / "synthetic_MAGs_5genomes"),
        "--completeness", "0.95,0.80,0.60",
        "--contamination", "0,0.05,0.15",
        "--heterogeneity", "0,0.10,0.30",
        "--fragmentation", "medium,high",
        "--seeds", "1,2,3",
    ]

    print("Running:")
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)

    print("")
    print("Done. Results are in:")
    print(root / "synthetic_MAGs_5genomes")


if __name__ == "__main__":
    main()
