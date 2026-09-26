#!/usr/bin/env python3
"""Draw a PCAM-aligned, multi-species JCVI macro-synteny ribbon figure.

The figure is derived exclusively from existing JCVI BED and anchors.simple
files. PZI, PSP and PPA are displayed as A/B subgenome tracks (8 + 8
chromosomes); diploid taxa retain their major chromosomes. Adjacent-track
ribbons are matched through their shared PZI/PCAM-projected intervals so that
no new sequence alignment is required.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

from PIL import Image, ImageDraw, ImageFont


SPECIES = [
    "pzi", "psp", "ppa", "pmal", "pav", "pcam", "pspe",
    "ppe", "pdu", "psal", "pmu_refseq", "par", "gtr", "rch",
]

SPECIES_NAMES = {
    "pzi": "Prunus zippeliana",
    "psp": "Prunus spinulosa",
    "ppa": "Prunus padus",
    "pmal": "Prunus mahaleb",
    "pav": "Prunus avium",
    "pcam": "Prunus campanulata",
    "pspe": "Prunus speciosa",
    "ppe": "Prunus persica",
    "pdu": "Prunus dulcis",
    "psal": "Prunus salicina",
    "pmu_refseq": "Prunus mume",
    "par": "Prunus armeniaca",
    "gtr": "Gillenia trifoliata",
    "rch": "Rosa chinensis",
}

POLYPLOID_COMPONENTS = {
    "pzi": ("pzi_hapA1", "pzi_hapB1"),
    "psp": ("psp_hapA1", "psp_hapB1"),
    "ppa": ("ppa_hapA", "ppa_hapB"),
}

GENERIC_LABEL = {
    "pmal": "pmal",
    "pav": "pav",
    "pcam": "pcam",
    "pspe": "pspe",
    "ppe": "ppe",
    "pdu": "pdu",
    "psal": "psal",
    "pmu_refseq": "pmu_refseq",
    "par": "par",
    "gtr": "gtr",
    "rch": "rch",
}


@dataclass
class Gene:
    seqid: str
    start: int
    end: int

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass
class BedData:
    genes: dict[str, Gene]
    lengths: dict[str, int]
    gene_counts: dict[str, int]


@dataclass
class Projection:
    species: str
    canonical: str
    canonical_start: float
    canonical_end: float
    target_seqid: str
    target_start: float
    target_end: float
    score: int
    source: str


@dataclass
class TrackChromosome:
    species: str
    seqid: str
    display: str
    canonical: str
    pcam_bin: int
    length: int
    reverse: bool
    group: str
    x1: float = 0.0
    x2: float = 0.0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir", type=Path,
        default=Path("downstream_analysis/local_A1B1_ribbon_inputs"),
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=Path(
            "analysis/karyotype/jcvi/figures/custom_svg/ribbons/"
            "multispecies_A1B1_pcam_aligned"
        ),
    )
    parser.add_argument("--min-chromosome-mb", type=float, default=5.0)
    parser.add_argument("--min-block-score", type=int, default=10)
    parser.add_argument("--max-ribbons-per-homolog", type=int, default=3)
    parser.add_argument("--png-scale", type=float, default=0.75)
    return parser.parse_args()


def natural_key(value: str) -> list[object]:
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", value)]


def read_bed(path: Path) -> BedData:
    genes: dict[str, Gene] = {}
    lengths: dict[str, int] = defaultdict(int)
    counts: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED has fewer than four columns: {path}:{line_number}")
            seqid, start_text, end_text, gene_id = fields[:4]
            start, end = int(start_text), int(end_text)
            if gene_id in genes:
                continue
            genes[gene_id] = Gene(seqid, start, end)
            lengths[seqid] = max(lengths[seqid], end)
            counts[seqid] += 1
    if not genes:
        raise ValueError(f"No BED records: {path}")
    return BedData(genes, dict(lengths), dict(counts))


def read_copy1_map(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("copy1") and row.get("display_seqid"):
                mapping[row["copy1"]] = row["display_seqid"]
    if len(mapping) != 16:
        raise ValueError(f"Expected 16 copy1 chromosomes in {path}, found {len(mapping)}")
    return mapping


def canonical_label(seqid: str) -> str | None:
    match = re.search(r"(C\d{2})_([AB])", seqid)
    if not match:
        return None
    return f"{match.group(1)}_{match.group(2)}"


def read_simple(path: Path) -> list[tuple[str, str, str, str, int, str]]:
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("##") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            elif line.startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 6:
                continue
            a, b, c, d, score_text, orientation = fields
            try:
                score = int(float(score_text))
            except ValueError:
                raise ValueError(f"Invalid score: {path}:{line_number}") from None
            rows.append((a, b, c, d, score, orientation))
    return rows


def fraction(gene: Gene, bed: BedData) -> float:
    return max(0.0, min(1.0, gene.midpoint / max(1, bed.lengths[gene.seqid])))


def pcam_reference_model(
    input_dir: Path,
    beds: dict[str, BedData],
    copy1_map: dict[str, str],
) -> tuple[dict[str, bool], dict[str, int], dict[str, str]]:
    """Infer PZI chromosome direction and homolog order from direct PCAM anchors."""
    pcam_seqids = sorted(
        [s for s, length in beds["pcam"].lengths.items() if length >= 5_000_000],
        key=natural_key,
    )
    pcam_index = {seqid: i for i, seqid in enumerate(pcam_seqids)}
    orientation_votes: dict[str, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
    target_votes: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))

    sources = [
        (
            "pzi_copy1", beds["pzi_copy1"], copy1_map,
            input_dir / "generic_anchors" / "pzi_copy1.pcam.anchors.simple",
        ),
        (
            "pzi_hapA1", beds["pzi_hapA1"], None,
            input_dir / "subgenome_anchors" / "pzi_hapA1.pcam.anchors.simple",
        ),
        (
            "pzi_hapB1", beds["pzi_hapB1"], None,
            input_dir / "subgenome_anchors" / "pzi_hapB1.pcam.anchors.simple",
        ),
    ]
    for _, qbed, seq_map, path in sources:
        for a, b, c, d, score, orientation in read_simple(path):
            if a not in qbed.genes or b not in qbed.genes or c not in beds["pcam"].genes:
                continue
            qseq = qbed.genes[a].seqid
            rseq = beds["pcam"].genes[c].seqid
            mapped = seq_map.get(qseq) if seq_map else qseq
            canonical = canonical_label(mapped or "")
            if not canonical or rseq not in pcam_index:
                continue
            orientation_votes[canonical][orientation == "-"] += max(1, score)
            target_votes[canonical][rseq] += max(1, score)

    reverse = {
        key: votes.get(True, 0) > votes.get(False, 0)
        for key, votes in orientation_votes.items()
    }
    best_target = {
        key: max(votes, key=votes.get) for key, votes in target_votes.items()
    }
    bins = {key: pcam_index[target] for key, target in best_target.items()}
    return reverse, bins, best_target


def make_projections(
    species: str,
    qbed: BedData,
    rbed: BedData,
    simple: Path,
    qseq_to_display: dict[str, str] | None,
    query_reverse: dict[str, bool],
    target_reverse_votes: dict[str, dict[bool, int]],
) -> list[Projection]:
    projections: list[Projection] = []
    for a, b, c, d, score, orientation in read_simple(simple):
        if not all(x in qbed.genes for x in (a, b)) or not all(x in rbed.genes for x in (c, d)):
            continue
        qa, qb, rc, rd = qbed.genes[a], qbed.genes[b], rbed.genes[c], rbed.genes[d]
        if qa.seqid != qb.seqid or rc.seqid != rd.seqid:
            continue
        qdisplay = qseq_to_display.get(qa.seqid) if qseq_to_display else qa.seqid
        canonical = canonical_label(qdisplay or "")
        if not canonical:
            continue
        q0, q1 = sorted((fraction(qa, qbed), fraction(qb, qbed)))
        query_is_reverse = query_reverse.get(canonical, False)
        if query_is_reverse:
            q0, q1 = 1.0 - q1, 1.0 - q0
        r0, r1 = sorted((fraction(rc, rbed), fraction(rd, rbed)))
        target_is_reverse = query_is_reverse ^ (orientation == "-")
        target_reverse_votes[rc.seqid][target_is_reverse] += max(1, score)
        projections.append(
            Projection(
                species, canonical, q0, q1, rc.seqid, r0, r1,
                max(1, score), simple.name,
            )
        )
    return projections


def major_seqids(bed: BedData, threshold: int, maximum: int = 8) -> list[str]:
    retained = [seqid for seqid, length in bed.lengths.items() if length >= threshold]
    if len(retained) > maximum:
        retained = sorted(retained, key=lambda x: bed.lengths[x], reverse=True)[:maximum]
    return retained


def native_labels(seqids: Iterable[str]) -> dict[str, str]:
    ordered = sorted(seqids, key=natural_key)
    return {seqid: f"Chr{i:02d}" for i, seqid in enumerate(ordered, 1)}


def build_tracks_and_projections(
    input_dir: Path, min_chr_bp: int
) -> tuple[
    dict[str, list[TrackChromosome]], dict[str, list[Projection]],
    dict[str, BedData], dict[str, int], dict[str, str],
]:
    bed_paths: dict[str, Path] = {}
    for path in (input_dir / "generic_bed").glob("*.bed"):
        bed_paths[path.stem] = path
    for path in (input_dir / "subgenome_bed").glob("*.bed"):
        bed_paths[path.stem] = path
    beds = {label: read_bed(path) for label, path in bed_paths.items()}
    required = {
        "pzi_copy1", "pzi_hapA1", "pzi_hapB1", "psp_hapA1", "psp_hapB1",
        "ppa_hapA", "ppa_hapB", "pcam", "pmal", "pav", "pspe", "ppe",
        "pdu", "psal", "pmu_refseq", "par", "gtr", "rch",
    }
    missing = sorted(required - beds.keys())
    if missing:
        raise ValueError(f"Missing BED inputs: {', '.join(missing)}")

    copy1_map = read_copy1_map(input_dir / "maps" / "chromosome_rename_map.tsv")
    query_reverse, pcam_bins, pcam_targets = pcam_reference_model(input_dir, beds, copy1_map)

    projections: dict[str, list[Projection]] = defaultdict(list)
    reverse_votes: dict[str, dict[str, dict[bool, int]]] = defaultdict(
        lambda: defaultdict(lambda: defaultdict(int))
    )

    # PZI projections are assembled from all direct A/B comparisons. This
    # preserves block-scale intervals while keeping PZI itself as the target.
    pzi_projection_sources = [
        ("pzi_hapA1", "pzi_hapA1.pcam.anchors.simple"),
        ("pzi_hapB1", "pzi_hapB1.pcam.anchors.simple"),
    ]
    for qlabel, filename in pzi_projection_sources:
        qbed = beds[qlabel]
        path = input_dir / "subgenome_anchors" / filename
        for a, b, _, _, score, _ in read_simple(path):
            if a not in qbed.genes or b not in qbed.genes:
                continue
            ga, gb = qbed.genes[a], qbed.genes[b]
            if ga.seqid != gb.seqid:
                continue
            canonical = canonical_label(ga.seqid)
            if not canonical:
                continue
            q0, q1 = sorted((fraction(ga, qbed), fraction(gb, qbed)))
            if query_reverse.get(canonical, False):
                q0, q1 = 1.0 - q1, 1.0 - q0
            projections["pzi"].append(
                Projection("pzi", canonical, q0, q1, ga.seqid, q0, q1, score, filename)
            )

    for species, components in (("psp", POLYPLOID_COMPONENTS["psp"]), ("ppa", POLYPLOID_COMPONENTS["ppa"])):
        for component in components:
            if component.endswith("A1"):
                query, filename = "pzi_hapA1", f"pzi_hapA1.{component}.anchors.simple"
            elif component.endswith("B1"):
                query, filename = "pzi_hapB1", f"pzi_hapB1.{component}.anchors.simple"
            elif component.endswith("hapA"):
                query, filename = "pzi_hapA1", "pzi_hapA1.ppa_hapA.anchors.simple"
            else:
                query, filename = "pzi_hapB1", "pzi_hapB1.ppa_hapB.anchors.simple"
            projections[species].extend(
                make_projections(
                    species, beds[query], beds[component],
                    input_dir / "subgenome_anchors" / filename,
                    None, query_reverse, reverse_votes[species],
                )
            )

    for species, label in GENERIC_LABEL.items():
        # GTR is represented by a single nine-chromosome BED, but its
        # anchors were generated separately against the PZI A1/B1 haplotypes.
        # Do not attempt to read a nonexistent pzi_copy1.gtr file here.
        if species == "gtr":
            continue
        filename = f"pzi_copy1.{label}.anchors.simple"
        projections[species].extend(
            make_projections(
                species, beds["pzi_copy1"], beds[label],
                input_dir / "generic_anchors" / filename,
                copy1_map, query_reverse, reverse_votes[species],
            )
        )

    if "gtr" in SPECIES:
        for query, filename in (
            ("pzi_hapA1", "pzi_hapA1.gtr.anchors.simple"),
            ("pzi_hapB1", "pzi_hapB1.gtr.anchors.simple"),
        ):
            projections["gtr"].extend(
                make_projections(
                    "gtr", beds[query], beds["gtr"],
                    input_dir / "subgenome_anchors" / filename,
                    None, query_reverse, reverse_votes["gtr"],
                )
            )

    tracks: dict[str, list[TrackChromosome]] = {}
    for species in SPECIES:
        if species in POLYPLOID_COMPONENTS:
            chroms: list[TrackChromosome] = []
            for group, component in zip(("A", "B"), POLYPLOID_COMPONENTS[species]):
                for seqid, length in beds[component].lengths.items():
                    canonical = canonical_label(seqid)
                    if not canonical:
                        continue
                    votes = reverse_votes[species].get(seqid, {})
                    reverse = query_reverse.get(canonical, False) if species == "pzi" else (
                        votes.get(True, 0) > votes.get(False, 0)
                    )
                    chroms.append(
                        TrackChromosome(
                            species, seqid, seqid, canonical,
                            pcam_bins.get(canonical, 99), length, reverse, group,
                        )
                    )
            tracks[species] = sorted(
                chroms, key=lambda c: (0 if c.group == "A" else 1, c.pcam_bin, natural_key(c.seqid))
            )
        else:
            label = GENERIC_LABEL[species]
            selected = major_seqids(
                beds[label], min_chr_bp, 9 if species == "gtr" else 8
            )
            native = native_labels(selected)
            canonical_scores: dict[str, dict[str, int]] = defaultdict(lambda: defaultdict(int))
            for block in projections[species]:
                canonical_scores[block.target_seqid][block.canonical] += block.score
            chroms = []
            for seqid in selected:
                scores = canonical_scores.get(seqid, {})
                canonical = max(scores, key=scores.get) if scores else "unassigned"
                pcam_bin = pcam_bins.get(canonical, 99)
                votes = reverse_votes[species].get(seqid, {})
                reverse = votes.get(True, 0) > votes.get(False, 0)
                chroms.append(
                    TrackChromosome(
                        species, seqid, native[seqid], canonical, pcam_bin,
                        beds[label].lengths[seqid], reverse, "diploid",
                    )
                )
            tracks[species] = sorted(chroms, key=lambda c: (c.pcam_bin, natural_key(c.seqid)))

    return tracks, projections, beds, pcam_bins, pcam_targets


def assign_geometry(tracks: dict[str, list[TrackChromosome]]) -> None:
    left, right, center_gap = 270.0, 2180.0, 88.0
    half_width = (right - left - center_gap) / 2.0
    cell_poly = half_width / 8.0
    diploid_count = max(
        8,
        max(
            (len(chromosomes) for species, chromosomes in tracks.items()
             if species not in POLYPLOID_COMPONENTS),
            default=8,
        ),
    )
    cell_dip = (right - left) / diploid_count
    global_max = max(c.length for values in tracks.values() for c in values)

    for species, chromosomes in tracks.items():
        if species in POLYPLOID_COMPONENTS:
            groups = {
                "A": [c for c in chromosomes if c.group == "A"],
                "B": [c for c in chromosomes if c.group == "B"],
            }
            for group, values in groups.items():
                group_left = left if group == "A" else left + half_width + center_gap
                for index, chrom in enumerate(values):
                    center = group_left + (index + 0.5) * cell_poly
                    width = cell_poly * (0.58 + 0.34 * chrom.length / global_max)
                    chrom.x1, chrom.x2 = center - width / 2.0, center + width / 2.0
        else:
            used_bins: set[int] = set()
            for fallback, chrom in enumerate(chromosomes):
                index = chrom.pcam_bin if 0 <= chrom.pcam_bin < 8 and chrom.pcam_bin not in used_bins else fallback
                used_bins.add(index)
                center = left + (index + 0.5) * cell_dip
                width = cell_dip * (0.56 + 0.36 * chrom.length / global_max)
                chrom.x1, chrom.x2 = center - width / 2.0, center + width / 2.0


def track_lookup(tracks: dict[str, list[TrackChromosome]]) -> dict[str, dict[str, TrackChromosome]]:
    return {species: {chrom.seqid: chrom for chrom in values} for species, values in tracks.items()}


def map_fraction(chrom: TrackChromosome, value: float) -> float:
    value = max(0.0, min(1.0, value))
    if chrom.reverse:
        value = 1.0 - value
    return chrom.x1 + value * (chrom.x2 - chrom.x1)


def match_ribbons(
    upper: str,
    lower: str,
    projections: dict[str, list[Projection]],
    lookup: dict[str, dict[str, TrackChromosome]],
    min_score: int,
    max_per_homolog: int,
) -> list[dict[str, object]]:
    candidates: dict[str, list[tuple[float, Projection, Projection, float, float]]] = defaultdict(list)
    upper_blocks = [b for b in projections[upper] if b.score >= min_score and b.target_seqid in lookup[upper]]
    lower_blocks = [b for b in projections[lower] if b.score >= min_score and b.target_seqid in lookup[lower]]
    lower_by_key: dict[str, list[Projection]] = defaultdict(list)
    for block in lower_blocks:
        lower_by_key[block.canonical].append(block)
    for a in upper_blocks:
        for b in lower_by_key.get(a.canonical, []):
            overlap0 = max(a.canonical_start, b.canonical_start)
            overlap1 = min(a.canonical_end, b.canonical_end)
            overlap = overlap1 - overlap0
            if overlap <= 0:
                continue
            shorter = min(a.canonical_end - a.canonical_start, b.canonical_end - b.canonical_start)
            if overlap < 0.008 and overlap / max(shorter, 1e-9) < 0.18:
                continue
            quality = overlap * math.sqrt(max(1, a.score) * max(1, b.score))
            candidates[a.canonical].append((quality, a, b, overlap0, overlap1))

    selected: list[dict[str, object]] = []
    for canonical, values in candidates.items():
        values.sort(key=lambda row: row[0], reverse=True)
        used_pairs: set[tuple[str, str, int, int]] = set()
        for quality, a, b, c0, c1 in values:
            pair_key = (
                a.target_seqid, b.target_seqid,
                round(c0 * 1000), round(c1 * 1000),
            )
            if pair_key in used_pairs:
                continue
            used_pairs.add(pair_key)

            def project(block: Projection, c: float) -> float:
                span = max(1e-9, block.canonical_end - block.canonical_start)
                t = (c - block.canonical_start) / span
                return block.target_start + t * (block.target_end - block.target_start)

            upper_chrom = lookup[upper][a.target_seqid]
            lower_chrom = lookup[lower][b.target_seqid]
            ua = sorted((map_fraction(upper_chrom, project(a, c0)), map_fraction(upper_chrom, project(a, c1))))
            lb = sorted((map_fraction(lower_chrom, project(b, c0)), map_fraction(lower_chrom, project(b, c1))))
            selected.append({
                "upper": upper, "lower": lower, "canonical": canonical,
                "upper_seqid": a.target_seqid, "lower_seqid": b.target_seqid,
                "upper_x1": ua[0], "upper_x2": ua[1],
                "lower_x1": lb[0], "lower_x2": lb[1],
                "canonical_start": c0, "canonical_end": c1,
                "quality": quality, "upper_source": a.source, "lower_source": b.source,
            })
            if len(used_pairs) >= max_per_homolog:
                break
    return selected


def ribbon_path(x1: float, x2: float, y1: float, z1: float, z2: float, y2: float) -> str:
    control = (y2 - y1) * 0.43
    return (
        f"M {x1:.2f},{y1:.2f} "
        f"C {x1:.2f},{y1 + control:.2f} {z1:.2f},{y2 - control:.2f} {z1:.2f},{y2:.2f} "
        f"L {z2:.2f},{y2:.2f} "
        f"C {z2:.2f},{y2 - control:.2f} {x2:.2f},{y1 + control:.2f} {x2:.2f},{y1:.2f} Z"
    )


def bezier(p0: tuple[float, float], p1: tuple[float, float], p2: tuple[float, float], p3: tuple[float, float], n: int = 18):
    for i in range(n + 1):
        t = i / n
        u = 1 - t
        yield (
            u**3 * p0[0] + 3 * u*u*t * p1[0] + 3*u*t*t*p2[0] + t**3*p3[0],
            u**3 * p0[1] + 3 * u*u*t * p1[1] + 3*u*t*t*p2[1] + t**3*p3[1],
        )


def draw_outputs(
    outdir: Path,
    tracks: dict[str, list[TrackChromosome]],
    projections: dict[str, list[Projection]],
    pcam_targets: dict[str, str],
    min_score: int,
    max_per_homolog: int,
    png_scale: float,
) -> tuple[Path, Path, list[dict[str, object]]]:
    outdir.mkdir(parents=True, exist_ok=True)
    width = 2320
    top = 150
    row_gap = 118
    bar_h = 18
    bottom = 95
    height = top + row_gap * (len(SPECIES) - 1) + bottom
    y_positions = {species: top + i * row_gap for i, species in enumerate(SPECIES)}
    lookup = track_lookup(tracks)
    all_ribbons: list[dict[str, object]] = []
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        all_ribbons.extend(
            match_ribbons(upper, lower, projections, lookup, min_score, max_per_homolog)
        )

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif}.title{font-size:27px;font-weight:700;fill:#18212A}.subtitle{font-size:15px;fill:#5B6670}.species{font-size:18px;font-style:italic;fill:#17212B}.chrom{font-size:11px;font-weight:600;fill:#FFFFFF;letter-spacing:0}.group{font-size:13px;font-weight:700;fill:#53606B}.note{font-size:12px;fill:#66717B}</style>',
        '<text class="title" x="270" y="43">PCAM-aligned macro-synteny across Prunus genomes</text>',
        '<text class="subtitle" x="270" y="69">PZI, PSP and PPA are partitioned into A (left) and B (right) subgenomes; ribbons are derived from screened JCVI blocks</text>',
        '<line x1="1181" y1="92" x2="1181" y2="125" stroke="#AAB3BB" stroke-width="1"/>',
        '<text class="group" x="725" y="106" text-anchor="middle">A subgenome</text>',
        '<text class="group" x="1640" y="106" text-anchor="middle">B subgenome</text>',
    ]

    # Ribbons are drawn first so chromosome bars and labels remain crisp.
    for ribbon in sorted(all_ribbons, key=lambda r: float(r["quality"])):
        y1 = y_positions[str(ribbon["upper"])] + bar_h / 2
        y2 = y_positions[str(ribbon["lower"])] - bar_h / 2
        path = ribbon_path(
            float(ribbon["upper_x1"]), float(ribbon["upper_x2"]), y1,
            float(ribbon["lower_x1"]), float(ribbon["lower_x2"]), y2,
        )
        title = html.escape(
            f'{ribbon["upper"]}-{ribbon["lower"]} {ribbon["canonical"]}: '
            f'{ribbon["upper_seqid"]} to {ribbon["lower_seqid"]}'
        )
        svg.append(
            f'<path d="{path}" fill="#8798A8" fill-opacity="0.24" '
            f'stroke="#718291" stroke-opacity="0.18" stroke-width="0.45"><title>{title}</title></path>'
        )

    for species in SPECIES:
        y = y_positions[species]
        reference = " (reference)" if species == "pcam" else ""
        svg.append(
            f'<text class="species" x="245" y="{y + 6:.2f}" text-anchor="end">'
            f'{html.escape(SPECIES_NAMES[species])}{reference}</text>'
        )
        if species in POLYPLOID_COMPONENTS:
            svg.append(f'<text class="group" x="278" y="{y - 17:.2f}">A1</text>')
            ppa_suffix = "B" if species == "ppa" else "B1"
            a_suffix = "A" if species == "ppa" else "A1"
            svg[-1] = f'<text class="group" x="278" y="{y - 17:.2f}">{a_suffix}</text>'
            svg.append(f'<text class="group" x="1205" y="{y - 17:.2f}">{ppa_suffix}</text>')
        for chrom in tracks[species]:
            fill = "#26323D" if chrom.group != "B" else "#34434F"
            svg.append(
                f'<rect x="{chrom.x1:.2f}" y="{y - bar_h/2:.2f}" '
                f'width="{chrom.x2 - chrom.x1:.2f}" height="{bar_h}" rx="8" '
                f'fill="{fill}" stroke="#17212A" stroke-width="0.7"><title>'
                f'{html.escape(chrom.seqid)}; {chrom.length:,} bp; '
                f'PCAM homolog {chrom.pcam_bin + 1 if chrom.pcam_bin < 8 else "NA"}; '
                f'orientation {"reverse" if chrom.reverse else "forward"}</title></rect>'
            )
            svg.append(
                f'<text class="chrom" x="{(chrom.x1 + chrom.x2)/2:.2f}" '
                f'y="{y + 4:.2f}" text-anchor="middle">{html.escape(chrom.display)}</text>'
            )

    svg.extend([
        f'<text class="note" x="270" y="{height - 34}">Only chromosomes at least 5 Mb are shown. '
        f'Block score threshold: {min_score}; maximum {max_per_homolog} ribbons per projected homolog and adjacent species pair.</text>',
        f'<text class="note" x="2180" y="{height - 34}" text-anchor="end">Chromosome direction is oriented to P. campanulata without altering internal block structure.</text>',
        '</svg>',
    ])
    svg_path = outdir / "multispecies_A1B1_pcam_aligned.ribbon.svg"
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")

    # Python/Pillow preview generated from the same geometry.
    scale = png_scale
    image = Image.new("RGB", (round(width * scale), round(height * scale)), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    try:
        font = ImageFont.truetype("arial.ttf", max(8, round(15 * scale)))
        small = ImageFont.truetype("arial.ttf", max(7, round(11 * scale)))
        title_font = ImageFont.truetype("arialbd.ttf", max(12, round(26 * scale)))
    except OSError:
        font = small = title_font = ImageFont.load_default()

    def S(point: tuple[float, float]) -> tuple[int, int]:
        return round(point[0] * scale), round(point[1] * scale)

    draw.text(S((270, 30)), "PCAM-aligned macro-synteny across Prunus genomes", fill=(24, 33, 42, 255), font=title_font)
    for ribbon in sorted(all_ribbons, key=lambda r: float(r["quality"])):
        y1 = y_positions[str(ribbon["upper"])] + bar_h / 2
        y2 = y_positions[str(ribbon["lower"])] - bar_h / 2
        control = (y2 - y1) * 0.43
        x1, x2 = float(ribbon["upper_x1"]), float(ribbon["upper_x2"])
        z1, z2 = float(ribbon["lower_x1"]), float(ribbon["lower_x2"])
        left_curve = list(bezier((x1, y1), (x1, y1 + control), (z1, y2 - control), (z1, y2)))
        right_curve = list(bezier((z2, y2), (z2, y2 - control), (x2, y1 + control), (x2, y1)))
        draw.polygon([S(p) for p in left_curve + right_curve], fill=(135, 152, 168, 62))
    for species in SPECIES:
        y = y_positions[species]
        label = SPECIES_NAMES[species] + (" (reference)" if species == "pcam" else "")
        bbox = draw.textbbox((0, 0), label, font=font)
        draw.text(S((245 - (bbox[2] - bbox[0]) / scale, y - 7)), label, fill=(23, 33, 43, 255), font=font)
        for chrom in tracks[species]:
            fill = (38, 50, 61, 255) if chrom.group != "B" else (52, 67, 79, 255)
            draw.rounded_rectangle(
                [S((chrom.x1, y - bar_h/2)), S((chrom.x2, y + bar_h/2))],
                radius=max(2, round(7 * scale)), fill=fill, outline=(23, 33, 42, 255), width=1,
            )
            bbox = draw.textbbox((0, 0), chrom.display, font=small)
            tx = (chrom.x1 + chrom.x2) / 2 - (bbox[2] - bbox[0]) / (2 * scale)
            draw.text(S((tx, y - 6)), chrom.display, fill=(255, 255, 255, 255), font=small)
    png_path = outdir / "multispecies_A1B1_pcam_aligned.ribbon.png"
    image.save(png_path, dpi=(300, 300))
    return svg_path, png_path, all_ribbons


def write_audits(
    outdir: Path,
    tracks: dict[str, list[TrackChromosome]],
    ribbons: list[dict[str, object]],
    pcam_targets: dict[str, str],
) -> None:
    with (outdir / "chromosome_order_orientation.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "species", "row_order", "subgenome", "display_order", "display_label",
            "source_seqid", "length_bp", "pzi_canonical_homolog", "pcam_homolog",
            "pcam_source_seqid", "orientation_to_pcam",
        ])
        for row_order, species in enumerate(SPECIES, 1):
            for display_order, chrom in enumerate(tracks[species], 1):
                writer.writerow([
                    species, row_order, chrom.group, display_order, chrom.display,
                    chrom.seqid, chrom.length, chrom.canonical,
                    chrom.pcam_bin + 1 if chrom.pcam_bin < 8 else "NA",
                    pcam_targets.get(chrom.canonical, "NA"),
                    "reverse" if chrom.reverse else "forward",
                ])

    with (outdir / "ribbon_projection_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        columns = [
            "upper_species", "lower_species", "pzi_canonical_homolog",
            "upper_source_seqid", "lower_source_seqid", "canonical_start_fraction",
            "canonical_end_fraction", "match_quality", "upper_jcvi_source", "lower_jcvi_source",
        ]
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in ribbons:
            writer.writerow({
                "upper_species": row["upper"],
                "lower_species": row["lower"],
                "pzi_canonical_homolog": row["canonical"],
                "upper_source_seqid": row["upper_seqid"],
                "lower_source_seqid": row["lower_seqid"],
                "canonical_start_fraction": f'{float(row["canonical_start"]):.6f}',
                "canonical_end_fraction": f'{float(row["canonical_end"]):.6f}',
                "match_quality": f'{float(row["quality"]):.4f}',
                "upper_jcvi_source": row["upper_source"],
                "lower_jcvi_source": row["lower_source"],
            })

    (outdir / "README.txt").write_text(
        "\n".join([
            "PCAM-aligned multi-species A/B-subgenome JCVI ribbon plot",
            "",
            "Rows: pzi, psp, ppa, pmal, pav, pcam, pspe, ppe, pdu, psal, pmu, par, rch.",
            "PZI and PSP show hapA1 (left) and hapB1 (right). PPA shows hapA (left) and hapB (right).",
            "PZI chromosome names follow chromosome_rename_map.tsv (C01-C08_A1/B1).",
            "All chromosome orders and whole-chromosome directions are aligned to PCAM.",
            "No chromosome-internal rearrangement or block direction was altered.",
            "",
            "Important interpretation:",
            "Ribbons between adjacent rows are reference-projected macro-synteny intervals.",
            "Each interval is supported by screened JCVI blocks sharing a PZI/PCAM-projected locus.",
            "They are not newly calculated direct pairwise JCVI blocks between every adjacent taxon.",
            "See ribbon_projection_audit.tsv for the source anchors used by every ribbon.",
        ]) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    args = parse_args()
    tracks, projections, _, _, pcam_targets = build_tracks_and_projections(
        args.input_dir, round(args.min_chromosome_mb * 1_000_000)
    )
    assign_geometry(tracks)
    svg, png, ribbons = draw_outputs(
        args.outdir, tracks, projections, pcam_targets,
        args.min_block_score, args.max_ribbons_per_homolog, args.png_scale,
    )
    write_audits(args.outdir, tracks, ribbons, pcam_targets)
    print(f"[OK] SVG: {svg}")
    print(f"[OK] PNG: {png}")
    print(f"[OK] chromosomes: {sum(len(x) for x in tracks.values())}")
    print(f"[OK] ribbons: {len(ribbons)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
