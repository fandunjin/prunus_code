#!/usr/bin/env python3
"""Prepare and render a PCAM-ordered seven-species A1/B1 synteny ribbon figure.

The script uses the 0903 JCVI BED/anchors for PZI, PSP, PPA, PCAM, PPE and
GTR. RCH is the only taxon absent from the 0903 pair plan; its existing JCVI
RCH anchors are used only to obtain the final GTR--RCH projected edge, while
the RCH track itself is retained as the seven chromosome reference track.

The output includes JCVI-compatible BED, seqids, layout and simple-anchor
files. The ``--prepare-only`` mode stages these inputs without drawing; the
companion native JCVI job then renders the final figure with
``jcvi.graphics.karyotype``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import html
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor
from reportlab.pdfgen import canvas


SPECIES = ["pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"]
POLYPLOID = {"pzi", "psp", "ppa"}
PALETTE = {
    "C01": "#0072B2", "C02": "#E69F00", "C03": "#009E73", "C04": "#D55E00",
    "C05": "#CC79A7", "C06": "#56B4E9", "C07": "#BEA800", "C08": "#505050",
}

# PCAM order used by the supplied native ribbon figure. C01_A and C03_A both
# map to PCAM chromosome 1 because they are the two retained PZI homologs.
GROUP_ORDER = [
    ("C01_A", 1), ("C03_A", 1), ("C07_A", 2), ("C08_A", 3),
    ("C04_A", 4), ("C05_A", 5), ("C06_A", 6), ("C02_A", 8),
    ("C03_B", 1), ("C07_B", 2), ("C08_B", 3), ("C04_B", 4),
    ("C05_B", 5), ("C06_B", 6), ("C01_B", 7), ("C02_B", 8),
]
GROUP_INDEX = {group: i for i, (group, _) in enumerate(GROUP_ORDER)}
GROUP_BIN = {group: pcam - 1 for group, pcam in GROUP_ORDER}
DEFAULT_REVERSE = {
    "C01_A": False, "C03_A": False, "C07_A": True, "C08_A": True,
    "C04_A": False, "C05_A": False, "C06_A": True, "C02_A": True,
    "C03_B": False, "C07_B": True, "C08_B": False, "C04_B": False,
    "C05_B": False, "C06_B": True, "C01_B": True, "C02_B": True,
}
RCH_ORDER = [
    ("NC_037091.1", 1, True), ("NC_037088.1", 2, False),
    ("NC_037089.1", 3, True), ("NC_037092.1", 4, False),
    ("NC_037094.1", 5, True), ("NC_037090.1", 8, True),
    ("NC_037093.1", 8, True),
]


@dataclass(frozen=True)
class Gene:
    seqid: str
    start: int
    end: int

    @property
    def midpoint(self) -> float:
        return (self.start + self.end) / 2.0


@dataclass
class Bed:
    genes: dict[str, Gene]
    seq_genes: dict[str, list[str]]
    lengths: dict[str, int]


@dataclass(frozen=True)
class Projection:
    canonical: str
    seqid: str
    canonical_start: float
    canonical_end: float
    target_start: float
    target_end: float
    score: int
    gene1: str
    gene2: str
    source: str


@dataclass
class Chromosome:
    species: str
    source_seqid: str
    display: str
    canonical: str
    pcam_bin: int
    reverse: bool
    length: int
    group: str
    x1: float = 0.0
    x2: float = 0.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hap-input", type=Path, default=Path("0903/haploid_repartition/inputs"))
    p.add_argument(
        "--jcvi-input", type=Path,
        default=Path("0903/downstream_analysis/figure2_repartition_ks/inputs"),
    )
    p.add_argument(
        "--pair-work", type=Path,
        default=Path("0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work"),
    )
    p.add_argument(
        "--rch-input", type=Path,
        default=Path("downstream_analysis/local_A1B1_ribbon_inputs"),
    )
    p.add_argument(
        "--outdir", type=Path,
        default=Path("analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_0903_7species_A1B1_pcam_pairwise"),
    )
    p.add_argument("--min-score", type=int, default=4)
    p.add_argument("--max-ribbons-per-homolog", type=int, default=4)
    p.add_argument("--basepair", action="store_true", help="Render only the base-pair-scaled outputs")
    p.add_argument("--prepare-only", action="store_true", help="Write JCVI inputs without rendering a figure")
    return p.parse_args()


def natural_key(value: str) -> list[object]:
    return [int(x) if x.isdigit() else x.lower() for x in re.split(r"(\d+)", value)]


def read_bed(path: Path) -> Bed:
    genes: dict[str, Gene] = {}
    lengths: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_no, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 4:
                raise ValueError(f"BED has fewer than four columns: {path}:{line_no}")
            seqid, start, end, gene = fields[:4]
            if gene in genes:
                continue
            genes[gene] = Gene(seqid, int(start), int(end))
            lengths[seqid] = max(lengths[seqid], int(end))
    seq_genes: dict[str, list[str]] = defaultdict(list)
    for gene, value in genes.items():
        seq_genes[value.seqid].append(gene)
    for seqid in seq_genes:
        seq_genes[seqid].sort(key=lambda g: (genes[g].midpoint, g))
    if not genes:
        raise ValueError(f"No BED records: {path}")
    return Bed(genes, dict(seq_genes), dict(lengths))


def parse_simple(path: Path) -> list[tuple[str, str, str, str, int, str]]:
    rows = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6:
                continue
            try:
                rows.append((*fields[:4], int(float(fields[4])), fields[5]))
            except ValueError:
                continue
    return rows


def source_canonical(seqid: str) -> str | None:
    m = re.search(r"C0*(\d+)_?([AB])(?:1)?$", seqid)
    if not m:
        m = re.search(r"chr0*(\d+)_([AB])(?:1)?$", seqid)
    if not m:
        return None
    return f"C{int(m.group(1)):02d}_{m.group(2)}"


def fraction(gene: Gene, bed: Bed) -> float:
    return max(0.0, min(1.0, gene.midpoint / max(1, bed.lengths[gene.seqid])))


def anchor_path(pair_work: Path, query: str, target: str) -> Path:
    return pair_work / f"{query}_vs_{target}" / f"{query}.{target}.anchors.simple"


def require(path: Path) -> None:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty input: {path}")


def fasta_lengths(path: Path) -> dict[str, int]:
    """Return sequence lengths from a FASTA without loading sequences in memory."""
    lengths: dict[str, int] = {}
    if not path.is_file() or path.stat().st_size == 0:
        return lengths
    opener = gzip.open if path.suffix == ".gz" else open
    current: str | None = None
    length = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if current is not None:
                    lengths[current] = length
                current = raw[1:].strip().split(None, 1)[0]
                length = 0
            elif current is not None:
                length += len(raw.strip())
        if current is not None:
            lengths[current] = length
    return lengths


def genome_lengths(args: argparse.Namespace) -> dict[str, dict[str, int]]:
    """Load chromosome lengths where the 0903 bundle contains a genome FASTA.

    JCVI's native ``--basepair`` mode infers chromosome size from the largest
    BED end coordinate.  The staged BEDs therefore receive an end sentinel
    based on these FASTA lengths, preserving the true assembly coordinate
    scale without changing any anchor gene IDs.
    """
    root = args.hap_input.resolve().parents[2]
    fasta_paths = {
        "pzi_hapA1": args.hap_input / "pzi" / "genome" / "pzi_hapA1.fa",
        "pzi_hapB1": args.hap_input / "pzi" / "genome" / "pzi_hapB1.fa",
        "psp_hapA1": args.hap_input / "psp" / "genome" / "psp_hapA1.fa",
        "psp_hapB1": args.hap_input / "psp" / "genome" / "psp_hapB1.fa",
        "ppa_hapA": args.hap_input / "ppa" / "genome" / "ppa_hapA.fa",
        "ppa_hapB": args.hap_input / "ppa" / "genome" / "ppa_hapB.fa",
        "pcam": root / "0903" / "TE" / "pcam" / "pcam.fa",
        "ppe": root / "0903" / "TE" / "ppe" / "ppe.fa",
        "rch": root / "data" / "rch" / "GCF_002994745.2_RchiOBHm-V2_genomic.fna.gz",
    }
    result: dict[str, dict[str, int]] = {}
    for component, path in fasta_paths.items():
        lengths = fasta_lengths(path)
        if lengths:
            result[component] = lengths
    return result


def source_orientation(paths: list[tuple[Path, Bed, Bed]], fallback: dict[str, bool]) -> dict[str, bool]:
    votes: dict[str, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
    for path, qbed, rbed in paths:
        require(path)
        for a, b, c, d, score, orient in parse_simple(path):
            if a not in qbed.genes or b not in qbed.genes or c not in rbed.genes:
                continue
            if qbed.genes[a].seqid != qbed.genes[b].seqid:
                continue
            key = source_canonical(qbed.genes[a].seqid)
            if key:
                votes[qbed.genes[a].seqid][orient == "-"] += max(1, score)
    out = dict(fallback)
    for seqid, values in votes.items():
        out[seqid] = values.get(True, 0) > values.get(False, 0)
    return out


def target_orientation(
    paths: list[tuple[Path, Bed, Bed]], source_orient: dict[str, bool],
    target_bed: Bed, fallback: dict[str, bool] | None = None,
) -> dict[str, bool]:
    votes: dict[str, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
    for path, qbed, rbed in paths:
        require(path)
        for a, b, c, d, score, orient in parse_simple(path):
            if a not in qbed.genes or b not in qbed.genes or c not in rbed.genes:
                continue
            qseq = qbed.genes[a].seqid
            tseq = rbed.genes[c].seqid
            key = source_canonical(qseq)
            if not key or tseq not in target_bed.lengths:
                continue
            qrev = source_orient.get(qseq, DEFAULT_REVERSE.get(key, False))
            votes[tseq][qrev ^ (orient == "-")] += max(1, score)
    out = dict(fallback or {})
    for seqid, values in votes.items():
        out[seqid] = values.get(True, 0) > values.get(False, 0)
    return out


def make_projection(
    path: Path, qbed: Bed, rbed: Bed, source_orient: dict[str, bool],
    use_query: bool,
) -> list[Projection]:
    out: list[Projection] = []
    for a, b, c, d, score, _orient in parse_simple(path):
        if score < 1 or not all(x in qbed.genes for x in (a, b)) or not all(x in rbed.genes for x in (c, d)):
            continue
        qa, qb, rc, rd = qbed.genes[a], qbed.genes[b], rbed.genes[c], rbed.genes[d]
        if qa.seqid != qb.seqid or rc.seqid != rd.seqid:
            continue
        key = source_canonical(qa.seqid)
        if not key:
            continue
        q0, q1 = sorted((fraction(qa, qbed), fraction(qb, qbed)))
        if source_orient.get(qa.seqid, DEFAULT_REVERSE.get(key, False)):
            q0, q1 = 1.0 - q1, 1.0 - q0
        if use_query:
            s0, s1 = sorted((fraction(qa, qbed), fraction(qb, qbed)))
            seqid, g1, g2 = qa.seqid, a, b
        else:
            s0, s1 = sorted((fraction(rc, rbed), fraction(rd, rbed)))
            seqid, g1, g2 = rc.seqid, c, d
        out.append(Projection(key, seqid, q0, q1, s0, s1, max(1, score), g1, g2, path.name))
    return out


def pcam_track_projections(path: Path, qbed: Bed, rbed: Bed, source_orient: dict[str, bool]) -> list[Projection]:
    """Use the target PCAM interval but the source PZI canonical coordinate."""
    return make_projection(path, qbed, rbed, source_orient, use_query=False)


def build_tracks_and_projections(args: argparse.Namespace):
    jcvi = args.jcvi_input
    work = args.pair_work
    beds = {
        "pzi_hapA1": read_bed(jcvi / "pzi_hapA1" / "pzi_hapA1.bed"),
        "pzi_hapB1": read_bed(jcvi / "pzi_hapB1" / "pzi_hapB1.bed"),
        "psp_hapA1": read_bed(jcvi / "psp_hapA1" / "psp_hapA1.bed"),
        "psp_hapB1": read_bed(jcvi / "psp_hapB1" / "psp_hapB1.bed"),
        "ppa_hapA": read_bed(jcvi / "ppa_hapA" / "ppa_hapA.bed"),
        "ppa_hapB": read_bed(jcvi / "ppa_hapB" / "ppa_hapB.bed"),
        "pcam": read_bed(jcvi / "pcam" / "pcam.bed"),
        "ppe": read_bed(jcvi / "ppe" / "ppe.bed"),
        "gtr": read_bed(jcvi / "gtr" / "gtr.bed"),
        "rch": read_bed(args.rch_input / "generic_bed" / "rch.bed"),
    }
    genome_size = genome_lengths(args)
    specs = {
        "pzi": [("pzi_hapA1", "pcam"), ("pzi_hapB1", "pcam")],
        "psp": [("psp_hapA1", "pcam"), ("psp_hapB1", "pcam")],
        "ppa": [("ppa_hapA", "pcam"), ("ppa_hapB", "pcam")],
    }
    paths: dict[str, list[tuple[Path, Bed, Bed]]] = defaultdict(list)
    for species, pairs in specs.items():
        for q, r in pairs:
            path = anchor_path(work, q, r)
            paths[species].append((path, beds[q], beds[r]))

    pzi_fallback = {seq: DEFAULT_REVERSE.get(source_canonical(seq) or "", False)
                    for q in ("pzi_hapA1", "pzi_hapB1") for seq in beds[q].lengths}
    pzi_orient = source_orientation(paths["pzi"], pzi_fallback)
    psp_orient = source_orientation(paths["psp"], {})
    ppa_orient = source_orientation(paths["ppa"], {})

    pzi_to_pcam = paths["pzi"]
    pzi_to_ppe = [(anchor_path(work, q, "ppe"), beds[q], beds["ppe"])
                  for q in ("pzi_hapA1", "pzi_hapB1")]
    pzi_to_gtr = [(anchor_path(work, q, "gtr"), beds[q], beds["gtr"])
                  for q in ("pzi_hapA1", "pzi_hapB1")]
    rch_base = args.rch_input / ".." / ".." / "analysis" / "karyotype" / "jcvi" / "figures" / "jcvi_native" / "reference_external"
    # Resolve the repository-level reference directory directly from the known
    # local input bundle instead of assuming the current working directory.
    rch_ref = Path("analysis/karyotype/jcvi/figures/jcvi_native/reference_external")
    rch_pairs = []
    for q in ("pzi_hapA1", "pzi_hapB1"):
        qbed = read_bed(rch_ref / f"{q}_vs_rch" / f"{q}.bed")
        path = rch_ref / f"{q}_vs_rch" / f"{q}.rch.anchors.simple"
        rch_pairs.append((path, qbed, beds["rch"]))

    ppe_orient = target_orientation(pzi_to_ppe, pzi_orient, beds["ppe"])
    gtr_orient = target_orientation(pzi_to_gtr, pzi_orient, beds["gtr"])
    rch_orient = target_orientation(rch_pairs, pzi_orient, beds["rch"], {s: rev for s, _, rev in RCH_ORDER})

    projections: dict[str, list[Projection]] = defaultdict(list)
    projections["pzi"] = [p for path, qbed, rbed in pzi_to_pcam for p in make_projection(path, qbed, rbed, pzi_orient, True)]
    projections["psp"] = [p for path, qbed, rbed in paths["psp"] for p in make_projection(path, qbed, rbed, psp_orient, True)]
    projections["ppa"] = [p for path, qbed, rbed in paths["ppa"] for p in make_projection(path, qbed, rbed, ppa_orient, True)]
    projections["pcam"] = [p for path, qbed, rbed in pzi_to_pcam for p in pcam_track_projections(path, qbed, rbed, pzi_orient)]
    projections["ppe"] = [p for path, qbed, rbed in pzi_to_ppe for p in make_projection(path, qbed, rbed, pzi_orient, False)]
    projections["gtr"] = [p for path, qbed, rbed in pzi_to_gtr for p in make_projection(path, qbed, rbed, pzi_orient, False)]
    projections["rch"] = [p for path, qbed, rbed in rch_pairs for p in make_projection(path, qbed, rbed, pzi_orient, False)]

    # Derive GTR order from the new PZI-GTR anchors and PCAM homolog bins.
    gtr_scores: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    for path, qbed, rbed in pzi_to_gtr:
        for a, b, c, d, score, _orient in parse_simple(path):
            if a not in qbed.genes or b not in qbed.genes or c not in rbed.genes:
                continue
            key = source_canonical(qbed.genes[a].seqid)
            if key and key in GROUP_BIN:
                gtr_scores[rbed.genes[c].seqid][GROUP_BIN[key]] += max(1, score)
    gtr_order = sorted(beds["gtr"].lengths, key=lambda s: (min(gtr_scores[s]) if gtr_scores[s] else 99, natural_key(s)))

    tracks: dict[str, list[Chromosome]] = {}
    for species, components, orient in (
        ("pzi", ("pzi_hapA1", "pzi_hapB1"), pzi_orient),
        ("psp", ("psp_hapA1", "psp_hapB1"), {**psp_orient}),
        ("ppa", ("ppa_hapA", "ppa_hapB"), {**ppa_orient}),
    ):
        chroms = []
        for group, component in zip(("A", "B"), components):
            for seqid in beds[component].lengths:
                key = source_canonical(seqid)
                if not key:
                    continue
                # 0903 haplotype annotations use the compact C1A/C1B style.
                # Keep the internal C01_A/C01_B key for mapping, but expose
                # the compact chromosome label in the rendered figure.
                compact = f"C{int(key[1:3])}{group}"
                length = max(
                    beds[component].lengths[seqid],
                    genome_size.get(component, {}).get(seqid, 0),
                )
                chroms.append(Chromosome(species, seqid, compact, key, GROUP_BIN.get(key, 99),
                                         orient.get(seqid, DEFAULT_REVERSE.get(key, False)), length, group))
        tracks[species] = sorted(chroms, key=lambda c: (0 if c.group == "A" else 1, GROUP_INDEX.get(c.canonical, 99)))

    for species, bedname, order, orient in (
        ("pcam", "pcam", sorted(beds["pcam"].lengths, key=natural_key), {}),
        ("ppe", "ppe", sorted(beds["ppe"].lengths, key=natural_key), ppe_orient),
        ("gtr", "gtr", gtr_order, gtr_orient),
    ):
        chroms = []
        for i, seqid in enumerate(order, 1):
            if seqid not in beds[bedname].lengths:
                continue
            length = max(
                beds[bedname].lengths[seqid],
                genome_size.get(bedname, {}).get(seqid, 0),
            )
            chroms.append(Chromosome(species, seqid, f"Chr{i:02d}", "", i - 1,
                                     orient.get(seqid, False), length, "diploid"))
        tracks[species] = chroms
    tracks["rch"] = [Chromosome("rch", seqid, f"Chr{i:02d}", "", pcam, rch_orient.get(seqid, rev),
                                 max(beds["rch"].lengths[seqid], genome_size.get("rch", {}).get(seqid, 0)), "diploid")
                     for i, (seqid, pcam, rev) in enumerate(RCH_ORDER, 1) if seqid in beds["rch"].lengths]
    return beds, tracks, projections, paths, rch_pairs, genome_size


def assign_geometry(tracks: dict[str, list[Chromosome]], basepair: bool) -> tuple[int, int]:
    width, height = 1296, 1728
    left, right, gap = 120.0, 1260.0, 80.0
    for species, chroms in tracks.items():
        if species in POLYPLOID:
            groups = {g: [c for c in chroms if c.group == g] for g in ("A", "B")}
            for group, vals in groups.items():
                lo = left if group == "A" else left + (right - left + gap) / 2
                hi = lo + (right - left - gap) / 2
                if basepair:
                    total = sum(c.length for c in vals) or 1
                    cursor = lo
                    for c in vals:
                        span = (hi - lo) * c.length / total
                        c.x1, c.x2 = cursor + span * 0.05, cursor + span * 0.95
                        cursor += span
                else:
                    cell = (hi - lo) / max(1, len(vals))
                    for i, c in enumerate(vals):
                        c.x1, c.x2 = lo + i * cell + cell * 0.08, lo + (i + 1) * cell - cell * 0.08
        else:
            if basepair:
                total = sum(c.length for c in chroms) or 1
                cursor = left
                for c in chroms:
                    span = (right - left) * c.length / total
                    c.x1, c.x2 = cursor + span * 0.025, cursor + span * 0.975
                    cursor += span
            else:
                cell = (right - left) / max(1, len(chroms))
                for i, c in enumerate(chroms):
                    c.x1, c.x2 = left + i * cell + cell * 0.05, left + (i + 1) * cell - cell * 0.05
    return width, height


def map_x(chrom: Chromosome, value: float) -> float:
    value = max(0.0, min(1.0, value))
    if chrom.reverse:
        value = 1.0 - value
    return chrom.x1 + value * (chrom.x2 - chrom.x1)


def match_pair(upper: str, lower: str, projections: dict[str, list[Projection]], tracks: dict[str, list[Chromosome]], min_score: int, cap: int):
    up = {c.source_seqid: c for c in tracks[upper]}
    lo = {c.source_seqid: c for c in tracks[lower]}
    by_key: dict[str, list[Projection]] = defaultdict(list)
    upper_by_key: dict[str, list[Projection]] = defaultdict(list)
    for b in projections[lower]:
        if b.score >= min_score and b.seqid in lo:
            by_key[b.canonical].append(b)
    for a in projections[upper]:
        if a.score >= min_score and a.seqid in up:
            upper_by_key[a.canonical].append(a)
    out = []
    for canonical, upper_values in upper_by_key.items():
        candidates = []
        for a in upper_values:
            for b in by_key.get(canonical, []):
                c0, c1 = max(a.canonical_start, b.canonical_start), min(a.canonical_end, b.canonical_end)
                if c1 <= c0:
                    continue
                quality = (c1 - c0) * math.sqrt(a.score * b.score)
                candidates.append((quality, a, b, c0, c1))
        candidates.sort(key=lambda x: x[0], reverse=True)
        used = 0
        seen = set()
        for quality, aa, bb, c0, c1 in candidates:
            pair = (aa.seqid, bb.seqid, round(c0, 5), round(c1, 5))
            if pair in seen:
                continue
            seen.add(pair)
            ta0 = aa.target_start + (c0 - aa.canonical_start) / max(1e-9, aa.canonical_end - aa.canonical_start) * (aa.target_end - aa.target_start)
            ta1 = aa.target_start + (c1 - aa.canonical_start) / max(1e-9, aa.canonical_end - aa.canonical_start) * (aa.target_end - aa.target_start)
            tb0 = bb.target_start + (c0 - bb.canonical_start) / max(1e-9, bb.canonical_end - bb.canonical_start) * (bb.target_end - bb.target_start)
            tb1 = bb.target_start + (c1 - bb.canonical_start) / max(1e-9, bb.canonical_end - bb.canonical_start) * (bb.target_end - bb.target_start)
            out.append({"canonical": aa.canonical, "upper_seqid": aa.seqid, "lower_seqid": bb.seqid,
                        "upper_x1": map_x(up[aa.seqid], ta0), "upper_x2": map_x(up[aa.seqid], ta1),
                        "lower_x1": map_x(lo[bb.seqid], tb0), "lower_x2": map_x(lo[bb.seqid], tb1),
                        "quality": quality})
            used += 1
            if used >= cap:
                break
    return out


def ribbon_points(x1: float, x2: float, y1: float, z1: float, z2: float, y2: float):
    control = (y2 - y1) * 0.43
    def cubic(a, b, c, d, n=18):
        for i in range(n + 1):
            t = i / n; u = 1 - t
            yield (u**3*a[0] + 3*u*u*t*b[0] + 3*u*t*t*c[0] + t**3*d[0],
                   u**3*a[1] + 3*u*u*t*b[1] + 3*u*t*t*c[1] + t**3*d[1])
    return list(cubic((x1, y1), (x1, y1 + control), (z1, y2 - control), (z1, y2))) + list(cubic((z2, y2), (z2, y2 - control), (x2, y1 + control), (x2, y1)))[::-1]


def write_staged_inputs(out: Path, beds: dict[str, Bed], tracks: dict[str, list[Chromosome]], projections: dict[str, list[Projection]], paths, rch_pairs, min_score: int, cap: int, genome_size: dict[str, dict[str, int]]):
    (out / "bed").mkdir(parents=True, exist_ok=True)
    (out / "anchors").mkdir(parents=True, exist_ok=True)
    lookup = {s: {c.source_seqid: c for c in cs} for s, cs in tracks.items()}
    bed_sources = {
        "pzi": ("pzi_hapA1", "pzi_hapB1"), "psp": ("psp_hapA1", "psp_hapB1"),
        "ppa": ("ppa_hapA", "ppa_hapB"), "pcam": ("pcam",), "ppe": ("ppe",),
        "gtr": ("gtr",), "rch": ("rch",),
    }
    for species, components in bed_sources.items():
        lines = []
        for comp in components:
            bed = beds[comp]
            for gene, g in bed.genes.items():
                if g.seqid not in lookup[species]:
                    continue
                c = lookup[species][g.seqid]
                lines.append(f"{c.display}\t{g.start}\t{g.end}\t{gene}\t0\t+")
            # Make JCVI's inferred base-pair size equal to the source genome
            # length.  The sentinel is never referenced by an anchor file.
            for seqid in bed.lengths:
                target_length = max(
                    bed.lengths[seqid], genome_size.get(comp, {}).get(seqid, 0)
                )
                if target_length <= 0 or seqid not in lookup[species]:
                    continue
                display = lookup[species][seqid].display
                sentinel = f"__jcvi_length__{comp}__{seqid}"
                lines.append(f"{display}\t{target_length - 1}\t{target_length}\t{sentinel}\t0\t+")
        (out / "bed" / f"{species}.bed").write_text("\n".join(lines) + "\n", encoding="utf-8")

    y = {s: 0.95 - i * (0.84 / (len(SPECIES) - 1)) for i, s in enumerate(SPECIES)}
    seqids = []
    layout = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
    basepair_layout = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
    # JCVI's --basepair mode uses the layout track span together with BED
    # coordinates.  Give each species a track width proportional to the sum
    # of its displayed chromosome spans; otherwise every row is forced into
    # the same .08-.92 interval and the figure appears equal-length.
    displayed_bp = {
        species: sum(chrom.length for chrom in tracks[species])
        for species in SPECIES
    }
    max_displayed_bp = max(displayed_bp.values()) or 1
    track_center = 0.50
    max_track_width = 0.84
    basepair_scale_rows = []
    for species in SPECIES:
        labels = []
        for c in tracks[species]:
            labels.append(c.display + ("-" if c.reverse else ""))
        seqids.append(",".join(labels))
        label = species
        layout.append(f" {y[species]:.5f}, .08, .92, 0, #4C5660, {label}, center, bed/{species}.bed, center")
        relative_width = displayed_bp[species] / max_displayed_bp
        track_width = max_track_width * relative_width
        xstart = track_center - track_width / 2.0
        xend = track_center + track_width / 2.0
        basepair_layout.append(
            f" {y[species]:.5f}, {xstart:.5f}, {xend:.5f}, 0, #4C5660, {label}, center, bed/{species}.bed, center"
        )
        basepair_scale_rows.append((species, displayed_bp[species], relative_width, xstart, xend))
    layout.append("# edges")
    basepair_layout.append("# edges")
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        layout.append(f" e, {SPECIES.index(upper)}, {SPECIES.index(lower)}, anchors/{upper}.{lower}.blocks.simple")
        basepair_layout.append(f" e, {SPECIES.index(upper)}, {SPECIES.index(lower)}, anchors/{upper}.{lower}.blocks.simple")
    (out / "multispecies_0903_A1B1_pcam.seqids").write_text("\n".join(seqids) + "\n", encoding="utf-8")
    (out / "multispecies_0903_A1B1_pcam.layout").write_text("\n".join(layout) + "\n", encoding="utf-8")
    (out / "multispecies_0903_A1B1_pcam_basepair.layout").write_text("\n".join(basepair_layout) + "\n", encoding="utf-8")
    with (out / "basepair_track_scale.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\tdisplayed_chromosome_span_bp\tdisplayed_chromosome_span_mb\trelative_to_max\txstart\txend\n")
        for species, span, relative, xstart, xend in basepair_scale_rows:
            handle.write(
                f"{species}\t{span}\t{span / 1_000_000:.6f}\t{relative:.6f}\t"
                f"{xstart:.5f}\t{xend:.5f}\n"
            )
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        upper_by_key: dict[str, list[Projection]] = defaultdict(list)
        lower_by_key: dict[str, list[Projection]] = defaultdict(list)
        for p in projections[upper]:
            if p.score >= min_score:
                upper_by_key[p.canonical].append(p)
        for p in projections[lower]:
            if p.score >= min_score:
                lower_by_key[p.canonical].append(p)
        rows = []
        for key in sorted(set(upper_by_key) & set(lower_by_key)):
            candidates = []
            for a in upper_by_key[key]:
                for b in lower_by_key[key]:
                    overlap = min(a.canonical_end, b.canonical_end) - max(a.canonical_start, b.canonical_start)
                    if overlap > 0:
                        candidates.append((overlap * math.sqrt(a.score * b.score), a, b))
            for _, a, b in sorted(candidates, key=lambda x: x[0], reverse=True)[:cap]:
                # JCVI's native SimpleFile parser treats the text before '*'
                # as the Shade highlight colour.  Use the same eight
                # chromosome colours as the 0903 Circos links so the native
                # JCVI ribbons and the Circos figure share one colour key.
                color = PALETTE.get(key.split("_", 1)[0], "#8798A8")
                rows.append(
                    f"{color}*{a.gene1}\t{a.gene2}\t{b.gene1}\t{b.gene2}\t"
                    f"{min(a.score, b.score)}\t+"
                )
        (out / "anchors" / f"{upper}.{lower}.blocks.simple").write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def write_figure(out: Path, tracks: dict[str, list[Chromosome]], projections: dict[str, list[Projection]], min_score: int, cap: int, basepair: bool):
    width, height = assign_geometry(tracks, basepair)
    rows = {s: 180 + i * 235 for i, s in enumerate(SPECIES)}
    bar_h = 26
    ribbons = []
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        ribbons.extend(match_pair(upper, lower, projections, tracks, min_score, cap))
    prefix = "multispecies_0903_A1B1_pcam.native.basepair" if basepair else "multispecies_0903_A1B1_pcam.native"
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" width="18in" height="24in" viewBox="0 0 {width} {height}">', '<rect width="100%" height="100%" fill="white"/>', '<style>text{font-family:Arial,Helvetica,sans-serif}.species{font-size:22px;font-style:italic;fill:#26323D}.chrom{font-size:12px;font-weight:600;fill:#FFFFFF}.tick{font-size:11px;fill:#53606B}</style>']
    ribbons_by_pair = []
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        for r in match_pair(upper, lower, projections, tracks, min_score, cap):
            r["upper"] = upper; r["lower"] = lower; ribbons_by_pair.append(r)
    for r in sorted(ribbons_by_pair, key=lambda x: x["quality"]):
        y1 = rows[r["upper"]] + bar_h / 2; y2 = rows[r["lower"]] - bar_h / 2
        pts = ribbon_points(r["upper_x1"], r["upper_x2"], y1, r["lower_x1"], r["lower_x2"], y2)
        d = "M " + " ".join(f"{x:.2f},{y:.2f}" for x, y in pts) + " Z"
        color = PALETTE.get(r["canonical"].split("_", 1)[0], "#8798A8")
        svg.append(f'<path d="{d}" fill="{color}" fill-opacity="0.43" stroke="{color}" stroke-opacity="0.25" stroke-width="0.6"/>')
    for species in SPECIES:
        y = rows[species]
        svg.append(f'<text class="species" x="103" y="{y + 7}" text-anchor="end">{species}</text>')
        for c in tracks[species]:
            svg.append(f'<rect x="{c.x1:.2f}" y="{y-bar_h/2:.2f}" width="{c.x2-c.x1:.2f}" height="{bar_h}" rx="10" fill="#4C5660" stroke="#1E2730" stroke-width="0.9"/>')
            svg.append(f'<text class="chrom" x="{(c.x1+c.x2)/2:.2f}" y="{y+4}" text-anchor="middle">{html.escape(c.display)}</text>')
    svg.append('</svg>')
    svg_path = out / f"{prefix}.svg"; svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")

    pdf_path = out / f"{prefix}.pdf"
    pdf = canvas.Canvas(str(pdf_path), pagesize=(width, height))
    pdf.setFillColor(HexColor("#FFFFFF")); pdf.rect(0, 0, width, height, fill=1, stroke=0)
    for r in sorted(ribbons_by_pair, key=lambda x: x["quality"]):
        y1 = rows[r["upper"]] + bar_h / 2; y2 = rows[r["lower"]] - bar_h / 2
        pts = ribbon_points(r["upper_x1"], r["upper_x2"], y1, r["lower_x1"], r["lower_x2"], y2)
        path = pdf.beginPath(); path.moveTo(*pts[0])
        for p in pts[1:]: path.lineTo(*p)
        path.close()
        color = PALETTE.get(r["canonical"].split("_", 1)[0], "#8798A8")
        pdf.setFillColor(HexColor(color)); pdf.setFillAlpha(0.43); pdf.drawPath(path, fill=1, stroke=0)
    pdf.setFillAlpha(1.0)
    pdf.setFont("Helvetica-Oblique", 16)
    for species in SPECIES:
        y = rows[species]; pdf.setFillColor(HexColor("#26323D")); pdf.drawRightString(103, y - 5, species)
        for c in tracks[species]:
            pdf.setFillColor(HexColor("#4C5660")); pdf.roundRect(c.x1, y-bar_h/2, c.x2-c.x1, bar_h, 10, fill=1, stroke=1)
            pdf.setFillColor(HexColor("#FFFFFF")); pdf.setFont("Helvetica-Bold", 8); pdf.drawCentredString((c.x1+c.x2)/2, y-3, c.display)
    pdf.save()

    scale = 2
    image = Image.new("RGB", (width*scale, height*scale), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    try:
        species_font = ImageFont.truetype("ariali.ttf", 30)
        chrom_font = ImageFont.truetype("arialbd.ttf", 16)
    except OSError:
        species_font = chrom_font = ImageFont.load_default()
    for r in sorted(ribbons_by_pair, key=lambda x: x["quality"]):
        y1 = rows[r["upper"]] + bar_h / 2; y2 = rows[r["lower"]] - bar_h / 2
        color = PALETTE.get(r["canonical"].split("_", 1)[0], "#8798A8")
        rgb = tuple(int(color[i:i+2], 16) for i in (1, 3, 5))
        draw.polygon([(round(x*scale), round(y*scale)) for x, y in ribbon_points(r["upper_x1"], r["upper_x2"], y1, r["lower_x1"], r["lower_x2"], y2)], fill=(*rgb, 110))
    for species in SPECIES:
        y = rows[species]
        draw.text((100*scale, (y-12)*scale), species, fill=(38,50,61,255), font=species_font, anchor="ra")
        for c in tracks[species]:
            draw.rounded_rectangle((round(c.x1*scale), round((y-bar_h/2)*scale), round(c.x2*scale), round((y+bar_h/2)*scale)), radius=18, fill=(76,86,96,255), outline=(30,39,48,255), width=2)
            draw.text((round((c.x1+c.x2)*scale/2), round((y-8)*scale)), c.display, fill=(255,255,255,255), font=chrom_font, anchor="ma")
    png_path = out / f"{prefix}.png"; image.save(png_path, dpi=(300, 300))
    return svg_path, pdf_path, png_path, ribbons_by_pair


def main() -> int:
    args = parse_args(); args.outdir.mkdir(parents=True, exist_ok=True)
    beds, tracks, projections, paths, rch_pairs, genome_size = build_tracks_and_projections(args)
    write_staged_inputs(args.outdir, beds, tracks, projections, paths, rch_pairs, args.min_score, args.max_ribbons_per_homolog, genome_size)
    if args.prepare_only:
        print(f"[OK] JCVI inputs prepared: {args.outdir}")
        return 0
    modes = [True] if args.basepair else [False, True]
    audits = []
    for basepair in modes:
        svg, pdf, png, ribbons = write_figure(args.outdir, tracks, projections, args.min_score, args.max_ribbons_per_homolog, basepair)
        audits.append(("basepair" if basepair else "gene_normalized", len(ribbons), svg.name, pdf.name, png.name))
    with (args.outdir / "ribbon_render_audit.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["scale", "ribbons", "svg", "pdf", "png"]); writer.writerows(audits)
    (args.outdir / "README.txt").write_text(
        "PCAM-ordered seven-species A1/B1 ribbon figure\n"
        "Species order: pzi, psp, ppa, pcam, ppe, gtr, rch.\n"
        "PZI, PSP and PPA each display eight A1 chromosomes on the left and eight B1 chromosomes on the right.\n"
        "JCVI pair parameters represented by the supplied 0903 inputs: cscore=0.70, min_size=4, minspan=30.\n"
        "The *.native.basepair.* outputs are the JCVI-native base-pair-scaled figures; the ordinary *.native.* outputs are gene-normalized.\n"
        "The basepair layout uses species-specific xstart/xend spans proportional to the displayed BED chromosome lengths.\n"
        "Ribbons are projected through the common PCAM/PZI homolog coordinate; no new alignment is inferred by the renderer.\n"
        "RCH is absent from the 0903 pair plan and uses the existing RCH JCVI evidence for the terminal projected edge only.\n",
        encoding="utf-8",
    )
    print(f"[OK] output: {args.outdir}")
    for row in audits: print(f"[OK] {row[0]} ribbons={row[1]} svg={row[2]} pdf={row[3]} png={row[4]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
