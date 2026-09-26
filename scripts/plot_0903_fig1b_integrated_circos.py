#!/usr/bin/env python3
"""Create an integrated Figure-1b-style radial synteny plot for the 0903 data.

The figure is deliberately assembled from the frozen 0903 JCVI beds/blocks and
the PZI Circos tracks.  It does not infer new alignments.  The inner circle is
PZI (A1 on the left, B1 on the right) with its existing density tracks and
within-genome A-B links.  The outer bands show PZI, PSP, PPA, PCAM, PPE, GTR
and RCH; ribbons are restricted to the six adjacent block files available in
the native pairwise result.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MplPath


SPECIES = ("pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch")
POLYPLOID = {"pzi", "psp", "ppa"}
PAIR_COLORS = {
    "C1": "#0072B2", "C2": "#E69F00", "C3": "#009E73", "C4": "#D55E00",
    "C5": "#CC79A7", "C6": "#56B4E9", "C7": "#BEA800", "C8": "#505050",
}
SPECIES_TINT = {
    "pzi": "#586875", "psp": "#5F8076", "ppa": "#9B7B43",
    "pcam": "#6B7184", "ppe": "#81718A", "gtr": "#66836B", "rch": "#6C747A",
}
DIPLOID_GROUPS = ("C1", "C3", "C7", "C8", "C4", "C5", "C6", "C2")

# A is the left half (90--270 degrees); B is the right half (270--450 degrees).
GROUP_START = {"A": math.radians(90.0), "B": math.radians(270.0)}
GROUP_END = {"A": math.radians(270.0), "B": math.radians(450.0)}
GROUP_GAP = math.radians(4.0)
CHROM_GAP = math.radians(0.8)


@dataclass(frozen=True)
class Chrom:
    species: str
    seqid: str
    label: str
    group: str
    index: int
    length: int
    start: float = 0.0
    end: float = 0.0


@dataclass(frozen=True)
class Gene:
    seqid: str
    start: int
    end: int


@dataclass
class Bed:
    genes: dict[str, Gene]
    lengths: dict[str, int]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("0903"))
    p.add_argument(
        "--native", type=Path,
        default=Path("0903/native_multispecies_A1B1_pcam_pairwise"),
    )
    p.add_argument(
        "--pzi-circos", type=Path,
        default=Path("0903/haploid_repartition/circos/pzi"),
    )
    p.add_argument(
        "--out", type=Path,
        default=Path("0903/haploid_repartition/circos/fig1b_integrated"),
    )
    p.add_argument("--dpi", type=int, default=600)
    p.add_argument("--size-mm", type=float, default=180.0)
    p.add_argument("--ribbon-alpha", type=float, default=0.16)
    return p.parse_args()


def canonical(value: str) -> str | None:
    m = re.search(r"C0*(\d+)_?([AB])", value)
    if not m:
        return None
    return f"C{int(m.group(1))}{m.group(2)}"


def pair_key(value: str) -> str | None:
    c = canonical(value)
    if c:
        return c[:2]
    m = re.search(r"(?:Chr|chr)0*(\d+)$", value)
    if m:
        n = int(m.group(1))
        return DIPLOID_GROUPS[n - 1] if 1 <= n <= len(DIPLOID_GROUPS) else "C8"
    return None


def fasta_lengths(path: Path) -> dict[str, int]:
    if not path.is_file():
        return {}
    opener = gzip.open if path.suffix == ".gz" else open
    lengths: dict[str, int] = {}
    current = None
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if current is not None:
                    lengths[current] = n
                current = raw[1:].strip().split(None, 1)[0]
                n = 0
            elif current is not None:
                n += len(raw.strip())
        if current is not None:
            lengths[current] = n
    return lengths


def read_bed(path: Path) -> Bed:
    genes: dict[str, Gene] = {}
    lengths: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            f = raw.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            try:
                start, end = int(f[1]), int(f[2])
            except ValueError:
                continue
            if f[3] not in genes:
                genes[f[3]] = Gene(f[0], start, end)
            lengths[f[0]] = max(lengths[f[0]], end)
    if not genes:
        raise ValueError(f"Empty BED: {path}")
    return Bed(genes, dict(lengths))


def karyotype_lengths(path: Path, phase: str) -> dict[str, int]:
    out: dict[str, int] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            f = raw.split()
            if len(f) < 6 or f[0] != "chr":
                continue
            m = re.search(r"(C\d+[AB])", f[2])
            # The single-species karyotypes contain two records per canonical
            # chromosome (e.g. C1A1/C1A2).  Select the requested phase by the
            # trailing haplotype digit while retaining the canonical C1A key.
            if m and ((phase in {"A1", "B1"} and f[2].endswith("1")) or
                      (phase not in {"A1", "B1"} and f[2].startswith(phase))):
                out[m.group(1)] = int(f[5])
    return out


def load_all_beds(native: Path) -> dict[str, Bed]:
    return {s: read_bed(native / "bed" / f"{s}.bed") for s in SPECIES}


def species_chromosomes(native: Path, pzi_circos: Path) -> dict[str, list[Chrom]]:
    beds = load_all_beds(native)
    pzi_lengths = karyotype_lengths(pzi_circos / "data" / "karyotype.txt", "A1")
    pzi_b_lengths = karyotype_lengths(pzi_circos / "data" / "karyotype.txt", "B1")
    pzi_lengths.update(pzi_b_lengths)
    chroms: dict[str, list[Chrom]] = {}
    for species in ("pzi", "psp", "ppa"):
        values = []
        # Native beds use the canonical C1A--C8B names.  Use A1/B1 lengths
        # where available and BED maxima as a conservative fallback.
        for group in ("A", "B"):
            for i in range(1, 9):
                seq = f"C{i}{group}"
                bed_len = beds[species].lengths.get(seq, 0)
                if species == "pzi":
                    length = pzi_lengths.get(seq, bed_len)
                else:
                    # The single-species Circos karyotypes contain the exact
                    # A1/B1 lengths for PSP/PPA.
                    kary = pzi_circos.parent / species / "data" / "karyotype.txt"
                    phase = "A1" if group == "A" else "B1"
                    lookup = karyotype_lengths(kary, phase) if kary.is_file() else {}
                    length = lookup.get(seq, bed_len)
                values.append(Chrom(species, seq, seq, group, i, max(1, length)))
        chroms[species] = values

    # The remaining species retain their own chromosome number and length.
    fasta_candidates = {
        "pcam": pzi_circos.parent.parent.parent / "TE" / "pcam" / "pcam.fa",
        "ppe": pzi_circos.parent.parent.parent / "TE" / "ppe" / "ppe.fa",
        "rch": pzi_circos.parent.parent.parent.parent / "data" / "rch" / "GCF_002994745.2_RchiOBHm-V2_genomic.fna.gz",
    }
    for species in ("pcam", "ppe", "gtr", "rch"):
        names = sorted(beds[species].lengths, key=lambda x: int(re.search(r"\d+", x).group()))
        fasta = fasta_lengths(fasta_candidates[species]) if species in fasta_candidates else {}
        vals = []
        if species in ("pcam", "ppe"):
            for i, seq in enumerate(names, 1):
                vals.append(Chrom(species, seq, seq, "B", i, max(1, fasta.get(seq, beds[species].lengths[seq]))))
        elif species == "gtr":
            for i, seq in enumerate(names, 1):
                vals.append(Chrom(species, seq, seq, "A", i, beds[species].lengths[seq]))
        else:
            for i, seq in enumerate(names, 1):
                vals.append(Chrom(species, seq, seq, "A", i, max(1, fasta.get(seq, beds[species].lengths[seq]))))
        chroms[species] = vals
    return chroms


def assign_angles(chroms: dict[str, list[Chrom]]) -> dict[str, list[Chrom]]:
    result: dict[str, list[Chrom]] = {}
    for species, values in chroms.items():
        by_group = {g: [c for c in values if c.group == g] for g in ("A", "B")}
        out = []
        for group, items in by_group.items():
            if not items:
                continue
            start = GROUP_START[group] + GROUP_GAP / 2
            end = GROUP_END[group] - GROUP_GAP / 2
            usable = end - start - CHROM_GAP * len(items)
            total = sum(c.length for c in items) or 1
            cursor = start
            for c in items:
                span = usable * c.length / total
                out.append(Chrom(c.species, c.seqid, c.label, c.group, c.index, c.length,
                                 cursor + CHROM_GAP / 2, cursor + span + CHROM_GAP / 2))
                cursor += span + CHROM_GAP
        result[species] = sorted(out, key=lambda c: c.start)
    return result


def parse_blocks(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6:
                continue
            try:
                yield fields[0], fields[1], fields[2], fields[3], int(float(fields[4])), fields[5]
            except ValueError:
                continue


def map_chrom(chroms: dict[str, list[Chrom]]) -> dict[str, dict[str, Chrom]]:
    return {species: {c.seqid: c for c in values} for species, values in chroms.items()}


def block_ribbons(native: Path, chroms: dict[str, list[Chrom]], beds: dict[str, Bed], upper: str, lower: str, filename: str):
    lookup = map_chrom(chroms)
    for g1, g2, g3, g4, score, orient in parse_blocks(native / "anchors" / filename):
        if not all(g in beds[upper].genes for g in (g1, g2)) or not all(g in beds[lower].genes for g in (g3, g4)):
            continue
        a, b = beds[upper].genes[g1], beds[upper].genes[g2]
        c, d = beds[lower].genes[g3], beds[lower].genes[g4]
        if a.seqid != b.seqid or c.seqid != d.seqid:
            continue
        if a.seqid not in lookup[upper] or c.seqid not in lookup[lower]:
            continue
        ca, cb = lookup[upper][a.seqid], lookup[lower][c.seqid]
        alen, blen = beds[upper].lengths[a.seqid], beds[lower].lengths[c.seqid]
        u0 = min(a.start, b.start) / max(1, alen); u1 = max(a.end, b.end) / max(1, alen)
        l0 = min(c.start, d.start) / max(1, blen); l1 = max(c.end, d.end) / max(1, blen)
        theta = (ca.start + u0 * (ca.end - ca.start), ca.start + u1 * (ca.end - ca.start))
        phi = (cb.start + l0 * (cb.end - cb.start), cb.start + l1 * (cb.end - cb.start))
        key = pair_key(ca.seqid) or pair_key(cb.seqid) or "C8"
        yield theta, phi, key, score, orient, upper, lower


def pzi_internal_links(pzi_circos: Path, chroms: dict[str, list[Chrom]]):
    lookup = map_chrom(chroms)["pzi"]
    path = pzi_circos / "data" / "links_ab.txt"
    if not path.is_file():
        return
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        f = raw.split()
        if len(f) < 6 or not f[0].startswith("A1_") or not f[3].startswith("B1_"):
            continue
        q = canonical(f[0].split("_", 1)[1]); s = canonical(f[3].split("_", 1)[1])
        if not q or not s or q not in lookup or s not in lookup:
            continue
        try:
            q0, q1, s0, s1 = map(int, (f[1], f[2], f[4], f[5]))
        except ValueError:
            continue
        qlen, slen = lookup[q].length, lookup[s].length
        qa = lookup[q]; sb = lookup[s]
        yield ((qa.start + q0 / max(1, qlen) * (qa.end - qa.start), qa.start + q1 / max(1, qlen) * (qa.end - qa.start)),
               (sb.start + s0 / max(1, slen) * (sb.end - sb.start), sb.start + s1 / max(1, slen) * (sb.end - sb.start)),
               pair_key(q) or "C8")


def xy(theta: float, radius: float) -> tuple[float, float]:
    return radius * math.cos(theta), radius * math.sin(theta)


def ribbon_path(theta: tuple[float, float], phi: tuple[float, float], r0: float, r1: float) -> MplPath:
    # Arc boundaries are sampled; connectors are cubic Beziers directed through
    # the radial midpoint, matching the smooth ribbons in Figure 1b.
    t0, t1 = theta; p0, p1 = phi
    vertices = [xy(t0, r0)]
    codes = [MplPath.MOVETO]
    for t in [t0 + (t1 - t0) * i / 8 for i in range(1, 9)]:
        vertices.append(xy(t, r0)); codes.append(MplPath.LINETO)
    mid = (r0 + r1) / 2
    c1 = xy(t1, mid); c2 = xy(p1, mid); q = xy(p1, r1)
    vertices.extend([c1, c2, q]); codes.extend([MplPath.CURVE4] * 3)
    for t in [p1 - (p1 - p0) * i / 8 for i in range(1, 9)]:
        vertices.append(xy(t, r1)); codes.append(MplPath.LINETO)
    c1 = xy(p0, mid); c2 = xy(t0, mid); q = xy(t0, r0)
    vertices.extend([c1, c2, q]); codes.extend([MplPath.CURVE4] * 3)
    codes[-1] = MplPath.CLOSEPOLY
    return MplPath(vertices, codes)


def draw_annular_value(ax, chrom: Chrom, start: int, end: int, value: float, maximum: float, r0: float, r1: float, color, alpha=1.0):
    if end <= start:
        return
    frac = max(0.0, min(1.0, value / maximum if maximum else 0.0))
    a0 = math.degrees(chrom.start + start / chrom.length * (chrom.end - chrom.start))
    a1 = math.degrees(chrom.start + end / chrom.length * (chrom.end - chrom.start))
    if a1 <= a0:
        a1 += 360
    ax.add_patch(Wedge((0, 0), r0 + frac * (r1 - r0), a0, a1,
                       width=frac * (r1 - r0), facecolor=color, edgecolor="none", alpha=alpha))


def read_track(path: Path) -> dict[str, list[tuple[int, int, float]]]:
    out: dict[str, list[tuple[int, int, float]]] = defaultdict(list)
    if not path.is_file():
        return out
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            f = raw.split()
            if len(f) < 4:
                continue
            try:
                out[f[0]].append((int(f[1]), int(f[2]), float(f[3])))
            except ValueError:
                continue
    return out


def draw_pzi_tracks(ax, pzi_circos: Path, chroms: dict[str, list[Chrom]]):
    lookup = map_chrom(chroms)["pzi"]
    data = pzi_circos / "data"
    specs = [
        ("ncrna_density.txt", 0.305, 0.330, "#6F9E9E", "hist", 1.0),
        ("gc_fraction.txt", 0.340, 0.370, "#2A586C", "line", 1.0),
        ("helitron_coverage.txt", 0.380, 0.410, "#289691", "hist", 0.060),
        ("tir_coverage.txt", 0.420, 0.450, "#7E58A0", "hist", 0.41),
        ("copia_coverage.txt", 0.460, 0.490, "#4C78A8", "hist", 0.16),
        ("gypsy_coverage.txt", 0.500, 0.530, "#C45B52", "hist", 0.32),
        ("te_coverage.txt", 0.540, 0.570, "#708732", "hist", 1.0),
    ]
    for filename, r0, r1, color, kind, maximum in specs:
        tracks = read_track(data / filename)
        if not tracks:
            continue
        if kind == "line":
            for seq, rows in tracks.items():
                seq_key = canonical(seq)
                if seq_key not in lookup:
                    continue
                c = lookup[seq_key]
                xs, ys = [], []
                for start, end, value in rows:
                    theta = c.start + ((start + end) / 2) / c.length * (c.end - c.start)
                    xs.append(theta); ys.append(r0 + max(0, min(1, (value - 0.35) / 0.22)) * (r1 - r0))
                if xs:
                    ax.plot([y * math.cos(x) for x, y in zip(xs, ys)], [y * math.sin(x) for x, y in zip(xs, ys)], color=color, lw=0.45, alpha=0.9)
        else:
            for seq, rows in tracks.items():
                seq_key = canonical(seq)
                if seq_key not in lookup:
                    continue
                c = lookup[seq_key]
                for start, end, value in rows:
                    draw_annular_value(ax, c, start, end, value, maximum, r0, r1, color, 0.88)

    # Gene density is deliberately a heatmap, matching the existing PZI Circos.
    genes = read_track(data / "gene_density.txt")
    norm = Normalize(vmin=0, vmax=max([v for rows in genes.values() for _, _, v in rows] or [1]))
    cmap = mpl.colormaps.get_cmap("Spectral_r")
    for seq, rows in genes.items():
        seq_key = canonical(seq)
        if seq_key not in lookup:
            continue
        c = lookup[seq_key]
        for start, end, value in rows:
            a0 = math.degrees(c.start + start / c.length * (c.end - c.start))
            a1 = math.degrees(c.start + end / c.length * (c.end - c.start))
            ax.add_patch(Wedge((0, 0), 0.615, a0, a1, width=0.035, facecolor=cmap(norm(value)), edgecolor="none"))


def draw_chrom_ring(ax, chroms: dict[str, list[Chrom]], species: str, r0: float, r1: float, labels=True):
    for c in chroms[species]:
        if species in POLYPLOID:
            color = PAIR_COLORS.get(f"C{c.index}", "#68737A")
        else:
            color = SPECIES_TINT[species]
        a0, a1 = math.degrees(c.start), math.degrees(c.end)
        ax.add_patch(Wedge((0, 0), r1, a0, a1, width=r1 - r0, facecolor=color, edgecolor="white", linewidth=0.35))
        if labels and (c.end - c.start) > math.radians(3.0):
            theta = (c.start + c.end) / 2
            x, y = xy(theta, r1 + 0.012)
            deg = math.degrees(theta)
            rotation = deg - 90
            if 90 < deg < 270:
                rotation += 180
            ax.text(x, y, c.label, fontsize=5.2, rotation=rotation, rotation_mode="anchor",
                    ha="center", va="center", color="#26323D")


def draw_species_label(ax, species: str, chroms: dict[str, list[Chrom]], r: float):
    groups = sorted({c.group for c in chroms[species]})
    for group in groups:
        # Place labels in the open upper sector so they do not stack on the
        # horizontal A/B boundary where the ribbon density is highest.
        theta = math.radians(112 if group == "A" else 68)
        x, y = xy(theta, r)
        ax.text(x, y, species, fontsize=6.0, fontstyle="italic", color="#26323D", ha="center", va="center")


def draw_ribbons(ax, ribbons, r0: float, r1: float, alpha: float, zorder: int):
    for theta, phi, key, score, orient, upper, lower in sorted(ribbons, key=lambda x: x[3]):
        patch = PathPatch(ribbon_path(theta, phi, r0, r1), facecolor=PAIR_COLORS.get(key, "#8798A8"), edgecolor="none", alpha=alpha, zorder=zorder)
        ax.add_patch(patch)


def main() -> int:
    args = parse_args()
    native = args.native.resolve(); pzi_circos = args.pzi_circos.resolve(); out = args.out.resolve()
    (out / "figures").mkdir(parents=True, exist_ok=True)
    beds = load_all_beds(native)
    chroms = assign_angles(species_chromosomes(native, pzi_circos))

    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "svg.fonttype": "none", "pdf.fonttype": 42, "font.size": 7,
    })
    fig = plt.figure(figsize=(args.size_mm / 25.4, args.size_mm / 25.4), facecolor="white")
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_aspect("equal"); ax.axis("off")
    ax.set_xlim(-1.02, 1.02); ax.set_ylim(-1.02, 1.02)

    # PZI's existing inner tracks and chromosome ring.
    draw_pzi_tracks(ax, pzi_circos, chroms)
    draw_chrom_ring(ax, chroms, "pzi", 0.580, 0.605, labels=True)
    for theta, phi, key in pzi_internal_links(pzi_circos, chroms):
        ax.add_patch(PathPatch(ribbon_path(theta, phi, 0.585, 0.105), facecolor=PAIR_COLORS[key], edgecolor="none", alpha=0.14, zorder=2))

    # Outer species bands.  PCAM/PPE are placed on B; GTR/RCH on A.
    ring = {"pzi": (0.635, 0.655), "psp": (0.675, 0.695), "ppa": (0.715, 0.735),
            "pcam": (0.755, 0.775), "ppe": (0.795, 0.815), "gtr": (0.755, 0.775), "rch": (0.795, 0.815)}
    for species in SPECIES:
        draw_chrom_ring(ax, chroms, species, *ring[species], labels=True)
        draw_species_label(ax, species, chroms, ring[species][1] + 0.035)

    edge_specs = [
        ("pzi", "psp", "pzi.psp.blocks.simple", 0.655, 0.675),
        ("psp", "ppa", "psp.ppa.blocks.simple", 0.695, 0.715),
        ("ppa", "pcam", "ppa.pcam.blocks.simple", 0.735, 0.755),
        ("pcam", "ppe", "pcam.ppe.blocks.simple", 0.775, 0.795),
        ("ppe", "gtr", "ppe.gtr.blocks.simple", 0.815, 0.755),
        ("gtr", "rch", "gtr.rch.blocks.simple", 0.775, 0.795),
    ]
    audit = []
    for upper, lower, filename, r0, r1 in edge_specs:
        values = list(block_ribbons(native, chroms, beds, upper, lower, filename))
        # For the cross-side PPE-GTR edge, route through the gap outside PPE.
        draw_ribbons(ax, values, r0, r1, args.ribbon_alpha, 1)
        audit.append((upper, lower, filename, len(values)))

    # A/B labels and a compact track key; labels are direct and publication sized.
    ax.text(-0.97, 0.02, "A", fontsize=8, fontweight="bold", color="#26323D", ha="center", va="center")
    ax.text(0.97, 0.02, "B", fontsize=8, fontweight="bold", color="#26323D", ha="center", va="center")
    legend_x, legend_y = -0.98, -0.94
    for i, (name, color) in enumerate((("Gene density", "#D95F02"), ("GC", "#2A586C"), ("Gypsy", "#C45B52"), ("Copia", "#4C78A8"), ("TE", "#708732"))):
        x = legend_x + (i % 3) * 0.19; y = legend_y - (i // 3) * 0.035
        ax.plot([x, x + 0.025], [y, y], color=color, lw=2.0)
        ax.text(x + 0.032, y, name, fontsize=5.5, color="#53606B", va="center")

    base = out / "figures" / "fig1b_integrated_pzi_multispecies"
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white")
    fig.savefig(base.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    with (out / "figure_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("item", "value"))
        writer.writerow(("inner", "PZI A1/B1 with existing Circos tracks and A-B links"))
        writer.writerow(("outer_order", "PZI, PSP, PPA; B-side PCAM, PPE; A-side GTR, RCH"))
        writer.writerow(("scale", "base-pair chromosome lengths; FASTA lengths where available, BED maxima otherwise"))
        writer.writerow(("ribbon_alpha", args.ribbon_alpha))
        for row in audit:
            writer.writerow((f"{row[0]}_{row[1]}_blocks", row[3]))
    (out / "README_methods.txt").write_text(
        "Figure 1b-style radial plot assembled from the 0903 native JCVI beds and blocks.\n"
        "The inner PZI ring uses A1 chromosomes on the left and B1 chromosomes on the right,\n"
        "with the existing gene/TE/GC/ncRNA tracks and A1-B1 links. Outer ribbons are only\n"
        "PZI-PSP, PSP-PPA, PPA-PCAM, PCAM-PPE, PPE-GTR and GTR-RCH. Chromosome arcs are\n"
        "scaled by base-pair length; labels and source coordinates are retained in the\n"
        "figure_manifest.tsv and the staged native input bundle.\n",
        encoding="ascii",
    )
    print(f"[OK] {base}.svg")
    for row in audit:
        print(f"[OK] {row[0]}-{row[1]}: {row[3]} blocks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
