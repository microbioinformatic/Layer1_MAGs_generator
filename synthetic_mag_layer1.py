#!/usr/bin/env python3
"""
synthetic_mag_layer1.py

Layer-1 synthetic MAG generator:
- Input: complete reference genomes + metadata
- Output: fragmented synthetic MAG FASTA files with truth tables and summaries
- Controls: completeness, contamination, fragmentation, and same-group strain heterogeneity

This does NOT simulate reads, assembly, or binning. It directly creates controlled draft MAGs
from known references, so every contig has exact coordinate-level truth.

Metadata TSV required columns:
    genome_id    fasta    species    strain    group

Example:
    G01    genomes/SL1344.fna    Salmonella_enterica    SL1344    Salmonella_enterica
    G02    genomes/LT2.fna       Salmonella_enterica    LT2       Salmonella_enterica
    G03    genomes/EcoliK12.fna  Escherichia_coli       K12       Escherichia_coli

Usage example:
    python synthetic_mag_layer1.py \
      --metadata metadata.tsv \
      --outdir synthetic_MAGs \
      --completeness 0.95,0.80,0.60 \
      --contamination 0,0.05,0.15 \
      --heterogeneity 0,0.10,0.30 \
      --fragmentation medium,high \
      --seeds 1,2,3

Author: generated with ChatGPT for Leyden Fernandez
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple


@dataclass
class GenomeMeta:
    genome_id: str
    fasta: Path
    species: str
    strain: str
    group: str


@dataclass
class SeqRecord:
    genome_id: str
    seq_id: str
    seq: str


@dataclass
class Segment:
    genome_id: str
    seq_id: str
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    seq: str
    role: str   # target | contaminant | heterogeneity

    @property
    def length(self) -> int:
        return self.end - self.start + 1


FRAGMENTATION_PRESETS = {
    "low":    {"target_contigs": 30,  "min_contig": 5000, "mean_contig": 150000},
    "medium": {"target_contigs": 120, "min_contig": 2000, "mean_contig": 50000},
    "high":   {"target_contigs": 400, "min_contig": 1000, "mean_contig": 15000},
}


def parse_float_list(text: str) -> List[float]:
    return [float(x.strip()) for x in text.split(",") if x.strip()]


def parse_int_list(text: str) -> List[int]:
    return [int(x.strip()) for x in text.split(",") if x.strip()]


def parse_str_list(text: str) -> List[str]:
    return [x.strip() for x in text.split(",") if x.strip()]


def read_fasta(path: Path, genome_id: str) -> List[SeqRecord]:
    records: List[SeqRecord] = []
    name = None
    chunks: List[str] = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(SeqRecord(genome_id, name, "".join(chunks).upper()))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
    if name is not None:
        records.append(SeqRecord(genome_id, name, "".join(chunks).upper()))
    if not records:
        raise ValueError(f"No FASTA records found in {path}")
    return records


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))


def safe_name(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", text)


def deterministic_seed(*items: object) -> int:
    raw = "|".join(str(x) for x in items).encode()
    return int(hashlib.md5(raw).hexdigest()[:8], 16)


def load_metadata(path: Path) -> Dict[str, GenomeMeta]:
    out: Dict[str, GenomeMeta] = {}
    with open(path, "r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"genome_id", "fasta", "species", "strain", "group"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Metadata missing columns: {sorted(missing)}")
        base = path.parent
        for row in reader:
            gid = row["genome_id"].strip()
            fasta = Path(row["fasta"].strip())
            if not fasta.is_absolute():
                fasta = base / fasta
            out[gid] = GenomeMeta(
                genome_id=gid,
                fasta=fasta,
                species=row["species"].strip(),
                strain=row["strain"].strip(),
                group=row["group"].strip(),
            )
    if not out:
        raise ValueError("Metadata table is empty")
    return out


def genome_lengths(records_by_genome: Dict[str, List[SeqRecord]]) -> Dict[str, int]:
    return {gid: sum(len(r.seq) for r in recs) for gid, recs in records_by_genome.items()}


def random_contig_lengths(total_len: int, target_contigs: int, min_contig: int,
                          mean_contig: int, rng: random.Random,
                          distribution: str = "lognormal") -> List[int]:
    """Generate approximate contig lengths summing to total_len."""
    if total_len <= min_contig:
        return [total_len]

    target_contigs = max(1, min(target_contigs, max(1, total_len // min_contig)))

    if distribution == "uniform":
        weights = [rng.random() + 0.01 for _ in range(target_contigs)]
    elif distribution == "lognormal":
        # many small contigs, few larger contigs
        weights = [rng.lognormvariate(math.log(max(mean_contig, 1)), 1.2) for _ in range(target_contigs)]
    else:
        raise ValueError(f"Unknown distribution: {distribution}")

    s = sum(weights)
    lengths = [max(min_contig, int(total_len * w / s)) for w in weights]

    # Adjust sum to total_len
    diff = total_len - sum(lengths)
    i = 0
    while diff != 0 and lengths:
        j = i % len(lengths)
        if diff > 0:
            lengths[j] += 1
            diff -= 1
        else:
            if lengths[j] > min_contig:
                lengths[j] -= 1
                diff += 1
        i += 1
        if i > 10_000_000:
            break

    # If min contig constraints made too many bp, merge from end
    while sum(lengths) > total_len and len(lengths) > 1:
        a = lengths.pop()
        lengths[-1] += a

    if sum(lengths) != total_len:
        lengths[-1] += total_len - sum(lengths)

    return [x for x in lengths if x > 0]


def fragment_records(records: List[SeqRecord], target_contigs: int, min_contig: int,
                     mean_contig: int, rng: random.Random,
                     distribution: str = "lognormal") -> List[Segment]:
    """Fragment each FASTA record independently, with unrestricted breakpoints."""
    total_len = sum(len(r.seq) for r in records)
    segments: List[Segment] = []

    for rec in records:
        # allocate contigs proportional to replicon length
        n = max(1, round(target_contigs * len(rec.seq) / total_len))
        lengths = random_contig_lengths(
            total_len=len(rec.seq),
            target_contigs=n,
            min_contig=min(min_contig, max(1, len(rec.seq))),
            mean_contig=mean_contig,
            rng=rng,
            distribution=distribution,
        )

        pos = 0
        for L in lengths:
            start0 = pos
            end0 = min(pos + L, len(rec.seq))
            if end0 <= start0:
                continue
            segments.append(
                Segment(
                    genome_id=rec.genome_id,
                    seq_id=rec.seq_id,
                    start=start0 + 1,
                    end=end0,
                    seq=rec.seq[start0:end0],
                    role="target",
                )
            )
            pos = end0

    rng.shuffle(segments)
    return segments


def select_bp(segments: List[Segment], target_bp: int, rng: random.Random,
              role: str) -> List[Segment]:
    """Randomly select whole contigs until at least target_bp is reached."""
    shuffled = list(segments)
    rng.shuffle(shuffled)
    selected: List[Segment] = []
    total = 0
    for s in shuffled:
        if total >= target_bp:
            break
        selected.append(Segment(s.genome_id, s.seq_id, s.start, s.end, s.seq, role))
        total += s.length
    return selected


def n50(lengths: Sequence[int]) -> Tuple[int, int]:
    """Return N50 and L50."""
    if not lengths:
        return 0, 0
    total = sum(lengths)
    half = total / 2
    acc = 0
    for i, L in enumerate(sorted(lengths, reverse=True), start=1):
        acc += L
        if acc >= half:
            return L, i
    return 0, 0


def write_mag_fasta(path: Path, mag_id: str, segments: List[Segment],
                    meta: Dict[str, GenomeMeta]) -> List[Dict[str, object]]:
    truth_rows: List[Dict[str, object]] = []
    with open(path, "w", encoding="utf-8") as out:
        for idx, seg in enumerate(segments, start=1):
            contig_id = f"{mag_id}_contig_{idx:05d}"
            gm = meta[seg.genome_id]
            header = (
                f">{contig_id} source_genome={seg.genome_id} "
                f"source_seq={seg.seq_id} start={seg.start} end={seg.end} "
                f"role={seg.role} species={safe_name(gm.species)} strain={safe_name(gm.strain)}"
            )
            out.write(header + "\n")
            out.write(wrap_fasta(seg.seq) + "\n")

            truth_rows.append({
                "mag_id": mag_id,
                "contig_id": contig_id,
                "source_genome": seg.genome_id,
                "source_seq": seg.seq_id,
                "source_species": gm.species,
                "source_strain": gm.strain,
                "source_group": gm.group,
                "start": seg.start,
                "end": seg.end,
                "length": seg.length,
                "role": seg.role,
            })
    return truth_rows


def write_truth(path: Path, rows: List[Dict[str, object]]) -> None:
    if not rows:
        return
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def choose_contaminant_genomes(target_gid: str, meta: Dict[str, GenomeMeta],
                               same_group: bool, rng: random.Random) -> List[str]:
    target_group = meta[target_gid].group
    candidates = []
    for gid, gm in meta.items():
        if gid == target_gid:
            continue
        if same_group:
            if gm.group == target_group:
                candidates.append(gid)
        else:
            if gm.group != target_group:
                candidates.append(gid)
    rng.shuffle(candidates)
    return candidates


def make_mag_for_one_condition(
    target_gid: str,
    meta: Dict[str, GenomeMeta],
    records_by_genome: Dict[str, List[SeqRecord]],
    lengths_by_genome: Dict[str, int],
    fragmented_by_genome: Dict[Tuple[str, str, int], List[Segment]],
    completeness: float,
    contamination: float,
    heterogeneity: float,
    fragmentation: str,
    target_contigs: int,
    min_contig: int,
    mean_contig: int,
    distribution: str,
    seed: int,
    heterogeneity_mode: str,
    contamination_same_group: bool,
) -> Tuple[str, List[Segment], Dict[str, object]]:

    rng = random.Random(deterministic_seed(target_gid, completeness, contamination,
                                           heterogeneity, fragmentation, seed))

    # Fragment target genome
    key = (target_gid, fragmentation, seed)
    if key not in fragmented_by_genome:
        fragmented_by_genome[key] = fragment_records(
            records_by_genome[target_gid],
            target_contigs=target_contigs,
            min_contig=min_contig,
            mean_contig=mean_contig,
            rng=random.Random(deterministic_seed("frag", target_gid, fragmentation, seed)),
            distribution=distribution,
        )
    target_segments_all = fragmented_by_genome[key]

    target_bp_requested = int(lengths_by_genome[target_gid] * completeness)
    target_segments = select_bp(target_segments_all, target_bp_requested, rng, "target")
    target_bp = sum(s.length for s in target_segments)

    # Heterogeneity: add or replace using related strains in same group
    hetero_segments: List[Segment] = []
    related = choose_contaminant_genomes(target_gid, meta, same_group=True, rng=rng)

    if heterogeneity > 0 and related:
        # bp target: heterogeneity / (target + heterogeneity)
        hetero_bp_requested = int((heterogeneity / max(1e-9, 1 - heterogeneity)) * target_bp)

        pool: List[Segment] = []
        for gid in related:
            k = (gid, fragmentation, seed)
            if k not in fragmented_by_genome:
                fragmented_by_genome[k] = fragment_records(
                    records_by_genome[gid],
                    target_contigs=target_contigs,
                    min_contig=min_contig,
                    mean_contig=mean_contig,
                    rng=random.Random(deterministic_seed("frag", gid, fragmentation, seed)),
                    distribution=distribution,
                )
            pool.extend(fragmented_by_genome[k])

        hetero_segments = select_bp(pool, hetero_bp_requested, rng, "heterogeneity")

        if heterogeneity_mode == "replacement":
            # Remove target contigs approximately equal to heterogeneity bp,
            # then add related-strain contigs. This keeps MAG size more genome-like.
            remove_bp = sum(s.length for s in hetero_segments)
            target_segments_sorted = list(target_segments)
            rng.shuffle(target_segments_sorted)
            kept = []
            removed = 0
            for s in target_segments_sorted:
                if removed < remove_bp:
                    removed += s.length
                else:
                    kept.append(s)
            target_segments = kept
            target_bp = sum(s.length for s in target_segments)
        elif heterogeneity_mode == "additive":
            pass
        else:
            raise ValueError("--heterogeneity-mode must be additive or replacement")

    # Contamination: add unrelated or same-group contigs
    contaminant_segments: List[Segment] = []
    clean_bp = target_bp + sum(s.length for s in hetero_segments)
    if contamination > 0:
        # contaminant / total = C; contaminant = C * clean / (1-C)
        contaminant_bp_requested = int((contamination / max(1e-9, 1 - contamination)) * clean_bp)

        contaminant_gids = choose_contaminant_genomes(
            target_gid, meta, same_group=contamination_same_group, rng=rng
        )
        # If no non-same-group genomes exist, fall back to any non-target genome.
        if not contaminant_gids:
            contaminant_gids = [gid for gid in meta if gid != target_gid]
            rng.shuffle(contaminant_gids)

        pool: List[Segment] = []
        for gid in contaminant_gids:
            k = (gid, fragmentation, seed)
            if k not in fragmented_by_genome:
                fragmented_by_genome[k] = fragment_records(
                    records_by_genome[gid],
                    target_contigs=target_contigs,
                    min_contig=min_contig,
                    mean_contig=mean_contig,
                    rng=random.Random(deterministic_seed("frag", gid, fragmentation, seed)),
                    distribution=distribution,
                )
            pool.extend(fragmented_by_genome[k])

        contaminant_segments = select_bp(pool, contaminant_bp_requested, rng, "contaminant")

    segments = target_segments + hetero_segments + contaminant_segments
    rng.shuffle(segments)

    total_bp = sum(s.length for s in segments)
    target_bp_final = sum(s.length for s in segments if s.role == "target")
    hetero_bp_final = sum(s.length for s in segments if s.role == "heterogeneity")
    contaminant_bp_final = sum(s.length for s in segments if s.role == "contaminant")
    N50, L50 = n50([s.length for s in segments])

    mag_id = (
        f"MAG_{safe_name(target_gid)}"
        f"_C{int(round(completeness*100)):03d}"
        f"_X{int(round(contamination*100)):03d}"
        f"_H{int(round(heterogeneity*100)):03d}"
        f"_F{fragmentation}"
        f"_S{seed}"
    )

    summary = {
        "mag_id": mag_id,
        "target_genome": target_gid,
        "target_species": meta[target_gid].species,
        "target_strain": meta[target_gid].strain,
        "target_group": meta[target_gid].group,
        "seed": seed,
        "fragmentation": fragmentation,
        "heterogeneity_mode": heterogeneity_mode,
        "completeness_requested": completeness,
        "contamination_requested": contamination,
        "heterogeneity_requested": heterogeneity,
        "target_genome_length_bp": lengths_by_genome[target_gid],
        "total_mag_bp": total_bp,
        "target_bp": target_bp_final,
        "heterogeneity_bp": hetero_bp_final,
        "contaminant_bp": contaminant_bp_final,
        "true_completeness": target_bp_final / lengths_by_genome[target_gid] if lengths_by_genome[target_gid] else 0,
        "true_contamination": contaminant_bp_final / total_bp if total_bp else 0,
        "true_heterogeneity": hetero_bp_final / (target_bp_final + hetero_bp_final) if (target_bp_final + hetero_bp_final) else 0,
        "contig_count": len(segments),
        "N50": N50,
        "L50": L50,
        "largest_contig": max([s.length for s in segments], default=0),
        "smallest_contig": min([s.length for s in segments], default=0),
    }
    return mag_id, segments, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Layer-1 synthetic MAG generator from complete genomes.")
    parser.add_argument("--metadata", required=True, type=Path, help="Metadata TSV with genome_id,fasta,species,strain,group")
    parser.add_argument("--outdir", required=True, type=Path, help="Output directory")
    parser.add_argument("--completeness", default="0.95,0.80,0.60", help="Comma-separated completeness fractions")
    parser.add_argument("--contamination", default="0,0.05,0.15", help="Comma-separated contamination fractions")
    parser.add_argument("--heterogeneity", default="0,0.10,0.30", help="Comma-separated heterogeneity fractions")
    parser.add_argument("--fragmentation", default="medium,high", help="Comma-separated: low,medium,high")
    parser.add_argument("--seeds", default="1,2,3", help="Comma-separated integer seeds")
    parser.add_argument("--distribution", default="lognormal", choices=["lognormal", "uniform"], help="Contig size distribution")
    parser.add_argument("--heterogeneity-mode", default="additive", choices=["additive", "replacement"],
                        help="additive adds related-strain contigs; replacement swaps target contigs for related-strain contigs")
    parser.add_argument("--contamination-same-group", action="store_true",
                        help="Use same-group genomes as contaminants. Default uses different groups when possible.")
    parser.add_argument("--only-genomes", default=None,
                        help="Optional comma-separated genome_ids to simulate instead of all genomes")
    parser.add_argument("--target-contigs", type=int, default=None,
                        help="Override target contig count for all fragmentation presets")
    parser.add_argument("--min-contig", type=int, default=None,
                        help="Override minimum contig size for all fragmentation presets")
    parser.add_argument("--mean-contig", type=int, default=None,
                        help="Override lognormal mean contig size for all fragmentation presets")
    args = parser.parse_args()

    meta = load_metadata(args.metadata)
    selected_gids = parse_str_list(args.only_genomes) if args.only_genomes else list(meta.keys())

    missing = [gid for gid in selected_gids if gid not in meta]
    if missing:
        raise ValueError(f"Genome IDs requested with --only-genomes not in metadata: {missing}")

    records_by_genome: Dict[str, List[SeqRecord]] = {}
    for gid, gm in meta.items():
        records_by_genome[gid] = read_fasta(gm.fasta, gid)

    lengths_by_genome = genome_lengths(records_by_genome)

    outdir = args.outdir
    mag_dir = outdir / "mags"
    truth_dir = outdir / "truth"
    summary_dir = outdir / "summary"
    for d in (mag_dir, truth_dir, summary_dir):
        d.mkdir(parents=True, exist_ok=True)

    completeness_values = parse_float_list(args.completeness)
    contamination_values = parse_float_list(args.contamination)
    heterogeneity_values = parse_float_list(args.heterogeneity)
    fragmentation_values = parse_str_list(args.fragmentation)
    seeds = parse_int_list(args.seeds)

    for f in fragmentation_values:
        if f not in FRAGMENTATION_PRESETS:
            raise ValueError(f"Unknown fragmentation preset '{f}'. Use: {sorted(FRAGMENTATION_PRESETS)}")

    fragmented_by_genome: Dict[Tuple[str, str, int], List[Segment]] = {}
    manifest_rows: List[Dict[str, object]] = []

    for gid in selected_gids:
        for frag in fragmentation_values:
            preset = dict(FRAGMENTATION_PRESETS[frag])
            target_contigs = args.target_contigs or preset["target_contigs"]
            min_contig = args.min_contig or preset["min_contig"]
            mean_contig = args.mean_contig or preset["mean_contig"]

            for seed in seeds:
                for comp in completeness_values:
                    for cont in contamination_values:
                        for het in heterogeneity_values:
                            mag_id, segments, summary = make_mag_for_one_condition(
                                target_gid=gid,
                                meta=meta,
                                records_by_genome=records_by_genome,
                                lengths_by_genome=lengths_by_genome,
                                fragmented_by_genome=fragmented_by_genome,
                                completeness=comp,
                                contamination=cont,
                                heterogeneity=het,
                                fragmentation=frag,
                                target_contigs=target_contigs,
                                min_contig=min_contig,
                                mean_contig=mean_contig,
                                distribution=args.distribution,
                                seed=seed,
                                heterogeneity_mode=args.heterogeneity_mode,
                                contamination_same_group=args.contamination_same_group,
                            )

                            fasta_path = mag_dir / f"{mag_id}.fna"
                            truth_path = truth_dir / f"{mag_id}.truth.tsv"
                            summary_path = summary_dir / f"{mag_id}.summary.json"

                            truth_rows = write_mag_fasta(fasta_path, mag_id, segments, meta)
                            write_truth(truth_path, truth_rows)
                            with open(summary_path, "w", encoding="utf-8") as handle:
                                json.dump(summary, handle, indent=2)

                            manifest = dict(summary)
                            manifest["fasta"] = str(fasta_path)
                            manifest["truth"] = str(truth_path)
                            manifest["summary"] = str(summary_path)
                            manifest_rows.append(manifest)

    manifest_path = outdir / "manifest.tsv"
    if manifest_rows:
        fields = list(manifest_rows[0].keys())
        with open(manifest_path, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
            writer.writeheader()
            writer.writerows(manifest_rows)

    print(f"Done. Generated {len(manifest_rows)} synthetic MAGs.")
    print(f"Manifest: {manifest_path}")
    print(f"MAG FASTA files: {mag_dir}")
    print(f"Truth tables: {truth_dir}")
    print(f"Summaries: {summary_dir}")


if __name__ == "__main__":
    main()
