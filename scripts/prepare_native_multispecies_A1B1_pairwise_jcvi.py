#!/usr/bin/env python3
"""Build PCAM-ordered adjacent-track JCVI inputs from existing JCVI blocks.

The supplied local package has direct PZI-to-taxon blocks, rather than every
adjacent-taxon comparison. Adjacent native JCVI blocks are therefore projected
only where the two source blocks overlap on the same normalized PZI homolog.
No sequence alignment is run and all displayed endpoints remain real genes.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


SPECIES = [
    "pzi", "psp", "ppa", "pmal", "pav", "pcam", "pspe", "ppe", "pdu",
    "psal", "pmu_refseq", "par", "gtr", "rch",
]
COMPONENTS = {
    "pzi": ("pzi_hapA1", "pzi_hapB1"),
    "psp": ("psp_hapA1", "psp_hapB1"),
    "ppa": ("ppa_hapA", "ppa_hapB"),
}
# Exact C01-C08 A1-B1 link palette from the four-haplotype Circos figure.
# A and B members of one homologous C group intentionally share a color.
PZI_CIRCOS_COLORS = {
    "C01": "#0072B2",  # blue
    "C02": "#E69F00",  # orange
    "C03": "#009E73",  # bluish green
    "C04": "#D55E00",  # vermilion
    "C05": "#CC79A7",  # reddish purple
    "C06": "#56B4E9",  # sky blue
    "C07": "#BEA800",  # mustard yellow
    "C08": "#505050",  # charcoal
}


@dataclass(frozen=True)
class Proj:
    canonical: str
    start: float
    end: float
    seqid: str
    lo_gene: str
    hi_gene: str
    reverse: bool
    score: int


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--order-table", required=True, type=Path)
    p.add_argument("--outdir", required=True, type=Path)
    p.add_argument("--min-score", type=int, default=0)
    p.add_argument(
        "--max-blocks-per-homolog",
        type=int,
        default=0,
        help="Retain at most this many projected blocks per PZI homolog; 0 keeps all.",
    )
    return p.parse_args()


def canon(seqid: str) -> str | None:
    m = re.search(r"(C\d{2})_([AB])", seqid)
    return f"{m.group(1)}_{m.group(2)}" if m else None


def pzi_circos_color(seqid: str) -> str:
    """Return the Circos color of a PZI C##_A1/B1 chromosome."""
    match = re.search(r"(C\d{2})_[AB]1$", seqid)
    if not match or match.group(1) not in PZI_CIRCOS_COLORS:
        raise ValueError(f"No Circos color is defined for PZI chromosome {seqid!r}")
    return PZI_CIRCOS_COLORS[match.group(1)]


def pzi_group_color(canonical: str) -> str:
    """Return a C01-C08 Circos color for a normalized PZI homolog."""
    return PZI_CIRCOS_COLORS.get(canonical.split("_", 1)[0], "#B8B8B8")


def read_bed(path: Path) -> tuple[list[list[str]], dict[str, str], dict[str, int], dict[str, list[str]]]:
    rows: list[list[str]] = []
    gene_seq: dict[str, str] = {}
    gene_rank: dict[str, int] = {}
    seq_genes: dict[str, list[str]] = defaultdict(list)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[3] in gene_seq:
                continue
            gene = fields[3]
            gene_seq[gene] = fields[0]
            gene_rank[gene] = len(seq_genes[fields[0]])
            seq_genes[fields[0]].append(gene)
            rows.append(fields)
    return rows, gene_seq, gene_rank, dict(seq_genes)


def parse_simple(path: Path):
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
                yield fields[:4], int(float(fields[4])), fields[5]
            except ValueError:
                continue


def read_order(path: Path) -> dict[str, list[dict[str, str]]]:
    rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows[row["species"]].append(row)
    for species in rows:
        rows[species].sort(key=lambda x: int(x["display_order"]))
    return rows


def copy1_map(path: Path) -> dict[str, str]:
    out = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("copy1") and row.get("display_seqid"):
                out[row["copy1"]] = row["display_seqid"]
    return out


def derive_gtr_order(order: dict[str, list[dict[str, str]]], beds: dict, inp: Path) -> None:
    """Place G. triflorum chromosomes against the existing PCAM/PZI model."""
    pzi_orientation = {
        row["display_label"]: row["orientation_to_pcam"] == "reverse"
        for row in order["pzi"]
    }
    pzi_pcam = {row["display_label"]: int(row["pcam_homolog"]) for row in order["pzi"]}
    rank = {f"Chr{i:02d}": i for i in range(1, 10)}
    scores: dict[str, dict[int, int]] = defaultdict(lambda: defaultdict(int))
    reverse_votes: dict[str, dict[bool, int]] = defaultdict(lambda: defaultdict(int))
    for qlabel in ("pzi_hapA1", "pzi_hapB1"):
        _, qseq, _, _ = beds[qlabel]
        _, rseq, _, _ = beds["gtr"]
        path = inp / "subgenome_anchors" / f"{qlabel}.gtr.anchors.simple"
        for genes, score, orientation in parse_simple(path):
            if not all(g in qseq for g in genes[:2]) or not all(g in rseq for g in genes[2:]):
                continue
            if qseq[genes[0]] != qseq[genes[1]] or rseq[genes[2]] != rseq[genes[3]]:
                continue
            target = rseq[genes[2]]
            if target not in rank:
                continue
            source = qseq[genes[0]]
            key = canon(source)
            if key not in pzi_pcam:
                continue
            scores[target][pzi_pcam[key]] += score
            reverse_votes[target][pzi_orientation[key] ^ (orientation == "-")] += score
    rows = []
    for seqid in sorted(rank, key=rank.get):
        if seqid not in beds["gtr"][3]:
            continue
        best = max(scores[seqid], key=scores[seqid].get) if scores[seqid] else 99
        reverse = reverse_votes[seqid].get(True, 0) > reverse_votes[seqid].get(False, 0)
        rows.append({
            "species": "gtr", "source_seqid": seqid, "display_label": seqid,
            "pcam_homolog": str(best), "orientation_to_pcam": "reverse" if reverse else "forward",
        })
    rows.sort(key=lambda row: (int(row["pcam_homolog"]), rank[row["source_seqid"]]))
    for i, row in enumerate(rows, 1):
        row["display_order"] = str(i)
    order["gtr"] = rows


def projection_blocks(
    source_label: str,
    target_species: str,
    target_label: str,
    simple: Path,
    beds: dict,
    copy_map: dict[str, str],
    min_score: int,
) -> list[Proj]:
    _, qseq, qrank, qgenes = beds[source_label]
    _, rseq, rrank, rgenes = beds[target_label]
    out: list[Proj] = []
    for genes, score, _ in parse_simple(simple):
        if score < min_score or not all(g in qseq for g in genes[:2]) or not all(g in rseq for g in genes[2:]):
            continue
        qa, qb, rc, rd = genes
        if qseq[qa] != qseq[qb] or rseq[rc] != rseq[rd]:
            continue
        qid = copy_map.get(qseq[qa], qseq[qa]) if source_label == "pzi_copy1" else qseq[qa]
        key = canon(qid)
        if not key:
            continue
        qn = max(1, len(qgenes[qseq[qa]]) - 1)
        p0, p1 = sorted((qrank[qa] / qn, qrank[qb] / qn))
        rlo, rhi = sorted((rrank[rc], rrank[rd]))
        rlist = rgenes[rseq[rc]]
        slope_reverse = (qrank[qb] - qrank[qa]) * (rrank[rd] - rrank[rc]) < 0
        out.append(Proj(key, p0, p1, rseq[rc], rlist[rlo], rlist[rhi], slope_reverse, score))
    return out


def write_pair_projection(upper: list[Proj], lower: list[Proj], out: Path, cap: int) -> int:
    by_upper: dict[str, list[Proj]] = defaultdict(list)
    by_lower: dict[str, list[Proj]] = defaultdict(list)
    for block in upper:
        by_upper[block.canonical].append(block)
    for block in lower:
        by_lower[block.canonical].append(block)
    lines: list[str] = []
    seen: set[tuple[str, str, str, str]] = set()
    for key in sorted(set(by_upper) & set(by_lower)):
        candidates = []
        for a in by_upper[key]:
            for b in by_lower[key]:
                overlap = min(a.end, b.end) - max(a.start, b.start)
                if overlap <= 0:
                    continue
                quality = overlap * math.sqrt(a.score * b.score)
                candidates.append((quality, a, b))
        used = 0
        for _, a, b in sorted(candidates, key=lambda x: x[0], reverse=True):
            reverse = a.reverse ^ b.reverse
            c, d = (b.hi_gene, b.lo_gene) if reverse else (b.lo_gene, b.hi_gene)
            row = (a.lo_gene, a.hi_gene, c, d)
            if row in seen:
                continue
            seen.add(row)
            color = pzi_group_color(key)
            lines.append(
                f"{color}*" + "\t".join((*row, str(min(a.score, b.score)), "-" if reverse else "+"))
            )
            used += 1
            if cap > 0 and used >= cap:
                break
    out.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")
    return len(lines)


def main() -> int:
    a = parse_args()
    inp, out = a.input_dir, a.outdir
    out.mkdir(parents=True, exist_ok=True)
    (out / "bed").mkdir(exist_ok=True)
    (out / "anchors").mkdir(exist_ok=True)
    order = read_order(a.order_table)
    copy_map = copy1_map(inp / "maps" / "chromosome_rename_map.tsv")

    source_paths = {p.stem: p for p in (inp / "generic_bed").glob("*.bed")}
    source_paths.update({p.stem: p for p in (inp / "subgenome_bed").glob("*.bed")})
    beds = {label: read_bed(path) for label, path in source_paths.items()}
    derive_gtr_order(order, beds, inp)
    valid: dict[str, set[str]] = {}

    for species in SPECIES:
        labels = COMPONENTS.get(species, (species,))
        mapping = {r["source_seqid"]: r["display_label"] for r in order[species]}
        seen: set[str] = set()
        staged: list[str] = []
        for label in labels:
            for fields in beds[label][0]:
                if fields[3] in seen:
                    continue
                seen.add(fields[3])
                fields = list(fields)
                fields[0] = mapping.get(fields[0], fields[0])
                staged.append("\t".join(fields))
        (out / "bed" / f"{species}.bed").write_text("\n".join(staged) + "\n", encoding="utf-8")
        valid[species] = seen

    # Blocks projected to every displayed target track through its PZI homolog.
    projections: dict[str, list[Proj]] = {}
    source_specs = {
        "psp": [("pzi_hapA1", "psp_hapA1", "pzi_hapA1.psp_hapA1.anchors.simple"), ("pzi_hapB1", "psp_hapB1", "pzi_hapB1.psp_hapB1.anchors.simple")],
        "ppa": [("pzi_hapA1", "ppa_hapA", "pzi_hapA1.ppa_hapA.anchors.simple"), ("pzi_hapB1", "ppa_hapB", "pzi_hapB1.ppa_hapB.anchors.simple")],
    }
    for species in ("pmal", "pav", "pcam", "pspe", "ppe", "pdu", "psal", "pmu_refseq", "par", "rch"):
        source_specs[species] = [("pzi_copy1", species, f"pzi_copy1.{species}.anchors.simple")]
    source_specs["gtr"] = [
        ("pzi_hapA1", "gtr", "pzi_hapA1.gtr.anchors.simple"),
        ("pzi_hapB1", "gtr", "pzi_hapB1.gtr.anchors.simple"),
    ]
    for species, specs in source_specs.items():
        blocks: list[Proj] = []
        for q, r, filename in specs:
            folder = "subgenome_anchors" if q != "pzi_copy1" else "generic_anchors"
            blocks.extend(projection_blocks(q, species, r, inp / folder / filename, beds, copy_map, a.min_score))
        projections[species] = blocks

    counts: list[tuple[str, int, str]] = []
    # The first pair is genuine direct A1/B1 JCVI synteny, not a projection.
    # Keep every valid simple block.  In particular, do not collapse blocks by
    # chromosome pair: the full PZI-PSP input contains 60 A1 blocks plus 55 B1
    # blocks, and retaining them is needed to show local rearrangements and
    # discontinuous homology along the same chromosome pair.
    first = out / "anchors" / "pzi.psp.blocks.simple"
    direct: list[str] = []
    for filename in ("pzi_hapA1.psp_hapA1.anchors.simple", "pzi_hapB1.psp_hapB1.anchors.simple"):
        for genes, score, orient in parse_simple(inp / "subgenome_anchors" / filename):
            if score >= a.min_score and all(g in valid["pzi"] for g in genes[:2]) and all(g in valid["psp"] for g in genes[2:]):
                qseq = beds["pzi_hapA1"][1].get(genes[0]) or beds["pzi_hapB1"][1].get(genes[0])
                rseq = beds["psp_hapA1"][1].get(genes[2]) or beds["psp_hapB1"][1].get(genes[2])
                if qseq and rseq:
                    # JCVI ``color*`` prefix colors this block only.  This
                    # retains the exact Circos C01-C08 palette on all 115
                    # direct PZI-PSP blocks while later species-pair layers
                    # stay neutral gray.
                    direct.append(
                        f"{pzi_circos_color(qseq)}*"
                        + "\t".join((*genes, str(score), orient))
                    )
    first.write_text("\n".join(direct) + "\n", encoding="utf-8")
    counts.append(("pzi-psp", len(direct), first.name))
    for upper, lower in zip(SPECIES[1:], SPECIES[2:]):
        destination = out / "anchors" / f"{upper}.{lower}.blocks.simple"
        n = write_pair_projection(projections[upper], projections[lower], destination, a.max_blocks_per_homolog)
        counts.append((f"{upper}-{lower}", n, destination.name))

    # Match the established tree-ribbon base-pair view: all tracks are centred,
    # but the width of each track is proportional to the total span of the
    # displayed chromosomes. The supplied local bundle contains gene BED files
    # rather than FASTA indexes, so chromosome span is measured from the first
    # to last annotated gene on each displayed chromosome.
    displayed_bp: dict[str, int] = {}
    for species in SPECIES:
        starts: dict[str, int] = {}
        ends: dict[str, int] = {}
        for component in COMPONENTS.get(species, (species,)):
            for fields in beds[component][0]:
                seqid = fields[0]
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                starts[seqid] = min(starts.get(seqid, start), start)
                ends[seqid] = max(ends.get(seqid, end), end)
        displayed_bp[species] = sum(ends[seqid] - starts[seqid] + 1 for seqid in starts)

    seqids = []
    layout = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
    basepair_layout = ["# y, xstart, xend, rotation, color, label, va, bed, label_va"]
    ys = [0.96 - i * (0.88 / (len(SPECIES) - 1)) for i in range(len(SPECIES))]
    max_displayed_bp = max(displayed_bp.values())
    track_center = 0.50
    max_track_width = 0.84
    for i, species in enumerate(SPECIES):
        names = [r["display_label"] + ("-" if r["orientation_to_pcam"] == "reverse" else "") for r in order[species]]
        seqids.append(",".join(names))
        label = "pmu" if species == "pmu_refseq" else species
        layout.append(f" {ys[i]:.5f}, .08, .92, 0, #4C5660, {label}, center, bed/{species}.bed, center")
        relative_width = displayed_bp[species] / max_displayed_bp
        track_width = max_track_width * relative_width
        xstart = track_center - track_width / 2.0
        xend = track_center + track_width / 2.0
        basepair_layout.append(
            f" {ys[i]:.5f}, {xstart:.5f}, {xend:.5f}, 0, #4C5660, {label}, center, bed/{species}.bed, center"
        )
    layout.append("# edges")
    basepair_layout.append("# edges")
    for i, (_, n, filename) in enumerate(counts):
        if n:
            layout.append(f" e, {i}, {i + 1}, anchors/{filename}")
            basepair_layout.append(f" e, {i}, {i + 1}, anchors/{filename}")
    (out / "multispecies_A1B1_pairwise.seqids").write_text("\n".join(seqids) + "\n", encoding="utf-8")
    (out / "multispecies_A1B1_pairwise.layout").write_text("\n".join(layout) + "\n", encoding="utf-8")
    (out / "multispecies_A1B1_pairwise_basepair.layout").write_text(
        "\n".join(basepair_layout) + "\n", encoding="utf-8"
    )
    with (out / "basepair_track_scale.tsv").open("w", encoding="utf-8") as handle:
        handle.write("species\tdisplayed_gene_span_bp\tdisplayed_gene_span_mb\trelative_to_max\txstart\txend\n")
        for species in SPECIES:
            relative_width = displayed_bp[species] / max_displayed_bp
            track_width = max_track_width * relative_width
            xstart = track_center - track_width / 2.0
            xend = track_center + track_width / 2.0
            handle.write(
                f"{species}\t{displayed_bp[species]}\t{displayed_bp[species] / 1_000_000:.6f}\t"
                f"{relative_width:.6f}\t{xstart:.5f}\t{xend:.5f}\n"
            )
    with (out / "pairwise_block_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("adjacent_pair\tnative_jcvi_blocks\tblock_file\n")
        for pair, n, filename in counts:
            handle.write(f"{pair}\t{n}\tanchors/{filename}\n")
    with (out / "pzi_psp_circos_color_legend.tsv").open("w", encoding="ascii") as handle:
        handle.write("pzi_chromosome\tcircos_link_color\tdescription\n")
        for group, color in PZI_CIRCOS_COLORS.items():
            for hap in ("A1", "B1"):
                handle.write(f"{group}_{hap}\t{color}\t{group} PZI-PSP JCVI blocks\n")
    print(f"[OK] Wrote native JCVI inputs: {out}")
    for pair, n, _ in counts:
        print(f"[OK] {pair}: {n} blocks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
