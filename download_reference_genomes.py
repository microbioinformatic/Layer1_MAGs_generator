#!/usr/bin/env python3
"""
download_reference_genomes.py

Pure-Python NCBI reference genome downloader for the Layer-1 synthetic MAG tool.

No Bash and no NCBI Datasets CLI are required. The script reads NCBI assembly
summary files, finds requested assembly accessions, downloads genomic FASTA
files, decompresses them, and writes metadata.tsv for synthetic_mag_layer1.py.

Required input TSV columns:
    genome_id
    assembly_accession
    species
    strain
    group
    output_fasta

Optional column:
    organism_label

Usage:
    python download_reference_genomes.py --accessions examples/accessions_5_genomes.tsv --outdir .
"""

from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
import urllib.request
from pathlib import Path
from typing import Dict, List


ASSEMBLY_SUMMARY_URLS = {
    "refseq": "https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt",
    "genbank": "https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt",
}


def download_url(url: str, output_path: Path, timeout: int = 180) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading: {url}")
    with urllib.request.urlopen(url, timeout=timeout) as response:
        with open(output_path, "wb") as out:
            shutil.copyfileobj(response, out)
    print(f"Wrote: {output_path}")


def read_assembly_summary(summary_path: Path) -> Dict[str, str]:
    """Return assembly accession -> FTP path."""
    ftp_by_accession: Dict[str, str] = {}
    header = None

    with open(summary_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue

            if line.startswith("#") and "assembly_accession" in line and "ftp_path" in line:
                header = line.lstrip("#").strip().split("\t")
                continue

            if line.startswith("#") or header is None:
                continue

            parts = line.split("\t")
            if len(parts) < len(header):
                continue

            row = dict(zip(header, parts))
            accession = row.get("assembly_accession")
            ftp_path = row.get("ftp_path")

            if accession and ftp_path and ftp_path != "na":
                ftp_by_accession[accession] = ftp_path

    print(f"Parsed {len(ftp_by_accession):,} assemblies from {summary_path.name}")
    return ftp_by_accession


def read_accessions_table(path: Path) -> List[dict]:
    required = {"genome_id", "assembly_accession", "species", "strain", "group", "output_fasta"}

    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Accessions table is missing columns: {sorted(missing)}")
        rows = list(reader)

    if not rows:
        raise ValueError("Accessions table is empty")

    return rows


def write_metadata(rows: List[dict], outdir: Path, metadata_name: str) -> Path:
    metadata_path = outdir / metadata_name

    with open(metadata_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=["genome_id", "fasta", "species", "strain", "group"],
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "genome_id": row["genome_id"],
                "fasta": row["output_fasta"],
                "species": row["species"],
                "strain": row["strain"],
                "group": row["group"],
            })

    print(f"Wrote metadata: {metadata_path}")
    return metadata_path


def download_genome(row: dict, ftp_path: str, outdir: Path) -> None:
    accession = row["assembly_accession"]
    genome_id = row["genome_id"]
    output_fasta = outdir / row["output_fasta"]

    folder_name = ftp_path.rstrip("/").split("/")[-1]
    gz_name = folder_name + "_genomic.fna.gz"
    url = ftp_path.rstrip("/") + "/" + gz_name
    gz_out = outdir / "ncbi_downloads" / gz_name

    print("")
    print(f"Genome: {genome_id}")
    print(f"Accession: {accession}")
    print(f"FTP path: {ftp_path}")

    download_url(url, gz_out)

    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(gz_out, "rb") as fin:
        with open(output_fasta, "wb") as fout:
            shutil.copyfileobj(fin, fout)

    print(f"Decompressed FASTA: {output_fasta}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Download reference genome FASTAs from NCBI.")
    parser.add_argument("--accessions", required=True, type=Path)
    parser.add_argument("--outdir", default=Path("."), type=Path)
    parser.add_argument("--metadata-name", default="metadata.tsv")
    parser.add_argument("--source", choices=["refseq", "genbank", "both"], default="both")
    args = parser.parse_args()

    outdir = args.outdir.resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    rows = read_accessions_table(args.accessions)
    write_metadata(rows, outdir, args.metadata_name)

    sources = ["refseq", "genbank"] if args.source == "both" else [args.source]

    ftp_by_accession: Dict[str, str] = {}
    for source in sources:
        url = ASSEMBLY_SUMMARY_URLS[source]
        summary_path = outdir / "ncbi_downloads" / url.split("/")[-1]
        download_url(url, summary_path)
        ftp_by_accession.update(read_assembly_summary(summary_path))

    missing = [row["assembly_accession"] for row in rows if row["assembly_accession"] not in ftp_by_accession]
    if missing:
        raise RuntimeError(
            "These accessions were not found in the selected NCBI assembly summaries:\n"
            + "\n".join(missing)
        )

    for row in rows:
        download_genome(row, ftp_by_accession[row["assembly_accession"]], outdir)

    print("")
    print("Done.")
    print(f"Project directory: {outdir}")
    print(f"Metadata for MAG simulation: {outdir / args.metadata_name}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print("")
        print("ERROR:")
        print(exc)
        sys.exit(1)
