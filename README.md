# Layer-1 Synthetic MAG Generator

A lightweight Python tool for generating **synthetic draft MAG-like genomes** directly from complete reference genomes.

This is a **Layer-1 simulator**: it does not simulate sequencing reads, metagenomic assembly, or binning. Instead, it directly fragments known complete genomes into synthetic MAGs with controlled levels of:

- genome fragmentation
- true completeness
- contaminant DNA
- same-group / strain-level heterogeneity
- coordinate-level truth labels

The goal is to create controlled benchmark datasets where every contig has a known origin.

## Why this tool exists

Many metagenomic simulators work through read simulation, assembly, and binning. That is useful for realism, but it makes it harder to isolate individual effects.

This tool is designed for controlled experiments such as:

- How does marker-gene completeness compare with true genome completeness?
- How does contamination affect MAG quality estimates?
- Can QC tools distinguish distant contamination from same-species strain heterogeneity?
- How much does fragmentation affect functional annotation?
- When does a mixed-strain MAG become biologically misleading?

Because every synthetic contig is derived from a known reference coordinate, the output includes exact ground truth.

## Repository contents

```text
.
├── synthetic_mag_layer1.py
├── download_reference_genomes.py
├── run_example_5_genomes.py
├── examples/
│   └── accessions_5_genomes.tsv
├── README.md
└── SUMMARY.md
```

## Requirements

Python 3.8 or newer.

No non-standard Python packages are required.

## Quick start

### 1. Download example reference genomes

```bash
python download_reference_genomes.py --accessions examples/accessions_5_genomes.tsv --outdir .
```

This creates:

```text
genomes/
ncbi_downloads/
metadata.tsv
```

### 2. Generate synthetic MAGs

```bash
python synthetic_mag_layer1.py ^
  --metadata metadata.tsv ^
  --outdir synthetic_MAGs_5genomes ^
  --completeness 0.95,0.80,0.60 ^
  --contamination 0,0.05,0.15 ^
  --heterogeneity 0,0.10,0.30 ^
  --fragmentation medium,high ^
  --seeds 1,2,3
```

On Linux/macOS, replace `^` with `\`.

Or run:

```bash
python run_example_5_genomes.py
```

## Input metadata format

The simulator requires a tab-separated metadata table with these columns:

```text
genome_id    fasta    species    strain    group
```

The `group` column defines which genomes may be treated as related strains for heterogeneity simulations.

## Input accessions format for downloader

The downloader requires a tab-separated table with these columns:

```text
genome_id    assembly_accession    species    strain    group    output_fasta
```

Optional:

```text
organism_label
```

## Output structure

```text
synthetic_MAGs_5genomes/
├── mags/
│   └── MAG_*.fna
├── truth/
│   └── MAG_*.truth.tsv
├── summary/
│   └── MAG_*.summary.json
└── manifest.tsv
```

## Definitions

```text
target_bp = bp in the MAG derived from the intended reference genome
heterogeneity_bp = bp in the MAG derived from related same-group genomes
contaminant_bp = bp in the MAG derived from unrelated genomes
total_MAG_bp = target_bp + heterogeneity_bp + contaminant_bp
```

The simulator reports:

```text
true_completeness = target_bp / full_target_genome_bp
true_contamination = contaminant_bp / total_MAG_bp
true_heterogeneity = heterogeneity_bp / (target_bp + heterogeneity_bp)
```

## Heterogeneity modes

```bash
--heterogeneity-mode additive
```

Adds related-strain contigs to the MAG.

```bash
--heterogeneity-mode replacement
```

Removes approximately equivalent target contigs and replaces them with related-strain contigs. This creates a more genome-sized mixed-strain/mosaic MAG.

## Notes

The current implementation uses unrestricted random fragmentation, so contig breaks can occur inside coding regions or intergenic regions.

Future extensions could add GFF-aware fragmentation modes:

- unrestricted
- intergenic-preferred
- intergenic-only
- gene-disruptive

## Suggested validation tools

The synthetic MAGs can be analyzed with common MAG QC or annotation tools, for example:

- CheckM / CheckM2
- GUNC
- GTDB-Tk
- QUAST / MetaQUAST
- functional annotation workflows

The value of this simulator is that estimated quality can be compared against exact coordinate-level truth.
