# Conceptual Summary: Layer-1 Synthetic MAG Generator

## Overview

The Layer-1 Synthetic MAG Generator creates controlled draft MAG-like genome bins directly from complete reference genomes. It is intended for benchmarking and method development in metagenomics, microbial genomics, and MAG quality assessment.

Unlike read-based metagenome simulators, this tool does not simulate sequencing reads, assembly, or binning. It starts from complete reference FASTA files and directly generates fragmented, incomplete, contaminated, and/or strain-heterogeneous synthetic MAGs.

This makes it possible to isolate specific factors that affect MAG interpretation.

## Main idea

Given a set of complete genomes, the tool creates artificial MAGs by:

1. Fragmenting complete genomes into contigs.
2. Retaining only a fraction of target-genome contigs to simulate incompleteness.
3. Adding contigs from unrelated genomes to simulate contamination.
4. Adding or replacing contigs with related same-group genomes to simulate strain heterogeneity.
5. Writing exact truth labels for every contig.

The result is a collection of MAG FASTA files with known coordinate-level origin.

## Why Layer 1?

MAG simulation can be separated into two conceptual layers.

### Layer 1: direct synthetic MAG generation

Input:

```text
complete genomes
```

Output:

```text
synthetic MAG FASTA files + truth labels
```

Layer 1 is useful when the scientific question requires precise control over genome completeness, contamination, fragmentation, and heterogeneity.

### Layer 2: read-based metagenomic simulation

Input:

```text
genomes + abundance profiles
```

Intermediate steps:

```text
read simulation -> assembly -> binning
```

Output:

```text
recovered MAGs
```

Layer 2 is more realistic but less controlled, because assembly and binning introduce additional sources of variation.

This repository focuses on Layer 1.

## Core parameters

### Fragmentation

The complete genome is cut into contigs. The current implementation uses unrestricted random fragmentation, meaning cuts can occur anywhere in the genome.

### Completeness

Completeness is simulated by retaining only a subset of target-genome contigs.

```text
true_completeness = target_bp / full_target_genome_bp
```

This is a coordinate-level quantity, not a marker-gene estimate.

### Contamination

Contamination is simulated by adding contigs from genomes outside the target group.

```text
true_contamination = contaminant_bp / total_MAG_bp
```

### Strain heterogeneity

Strain heterogeneity is simulated using genomes in the same user-defined group.

```text
true_heterogeneity = heterogeneity_bp / (target_bp + heterogeneity_bp)
```

This is intentionally separated from contamination. Same-species or same-group heterogeneity is biologically different from contamination by unrelated organisms.

## Heterogeneity modes

### Additive mode

Related-strain contigs are added on top of the retained target contigs.

### Replacement mode

Some target contigs are removed and replaced by related-strain contigs. This creates a more genome-sized synthetic MAG and can mimic a mosaic/mixed-strain reconstruction.

## Output truth

Every synthetic MAG has three layers of output.

### FASTA

The MAG sequence itself.

### Truth table

A tab-separated file where every contig is assigned:

- MAG ID
- contig ID
- source genome
- source sequence
- start coordinate
- end coordinate
- length
- role: target, contaminant, or heterogeneity

### Summary JSON

A machine-readable file containing the requested and observed simulation parameters.

The observed values may differ slightly from requested values because the simulator selects whole contigs rather than fractional contigs.

## Scientific use cases

This tool can support studies asking:

1. How well do estimated MAG completeness values match true coordinate-level completeness?
2. How do contamination estimates behave when contaminant contigs come from unrelated genomes?
3. How do tools behave when the signal is same-species strain heterogeneity rather than distant contamination?
4. How much fragmentation is needed before annotation quality degrades?
5. How robust are genome classification and functional annotation tools under controlled MAG degradation?
6. Can current MAG QC tools distinguish incomplete genomes, contaminated genomes, and mixed-strain genomes?

## Novelty and positioning

Synthetic metagenome and MAG simulation tools already exist, especially tools that simulate reads and then perform assembly and binning. The contribution of this tool is narrower and more controlled.

The distinctive feature is:

```text
direct coordinate-level generation of MAG-like FASTA bins from complete references,
with independent control of fragmentation, true completeness, contamination,
and same-group strain heterogeneity.
```

This makes the tool useful as a controlled benchmark layer before moving to more realistic read-based simulations.
