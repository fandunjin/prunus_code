#!/usr/bin/env python3
"""Plot PZI copy1 self-synteny highlights with a restricted taxon panel.

The PZI copy1 row is retained as the reference because the highlighted events
are internal PZI synteny blocks.  The external rows are restricted to PPA, PSP,
PPE, PCAM and GTR; their ribbons are the existing PCAM-aligned JCVI projections.
No new sequence alignment is performed by this script.
"""

from __future__ import annotations

import argparse
import csv
import html
import sys
from collections import defaultdict
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import plot_multispecies_A1B1_pcam_ribbon as base  # noqa: E402


SPECIES = ["pzi", "ppa", "psp", "ppe", "pcam", "gtr"]
HIGHLIGHT_COLORS = {
    "chr12_mid_chr11_tail": "#0072B2",   # blue
    "chr12_head_chr16_head": "#E69F00",  # orange
    "chr05_head_chr15_tail": "#009E73",  # bluish green
}
HIGHLIGHT_LABELS = {
    "chr12_mid_chr11_tail": "chr12 middle - chr11 distal block",
    "chr12_head_chr16_head": "chr12 head - chr16 head block",
    "chr05_head_chr15_tail": "chr05 head - chr15 distal fragments",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir", type=Path,
        default=ROOT / "downstream_analysis/local_A1B1_ribbon_inputs",
    )
    parser.add_argument(
        "--self-bed", type=Path,
        default=ROOT / (
            "downstream_analysis/karyotype/"
            "pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/pzi_copy1.bed"
        ),
    )
    parser.add_argument(
        "--self-anchors", type=Path,
        default=ROOT / (
            "downstream_analysis/karyotype/"
            "pzi_copy1_ordered_self_jcvi_local/pairs/pzi_copy1_self/"
            "pzi_copy1.pzi_copy1.anchors.simple"
        ),
    )
    parser.add_argument(
        "--outdir", type=Path,
        default=ROOT / (
            "analysis/karyotype/jcvi/figures/custom_svg/ribbons/"
            "pzi_self_highlight_ppa_psp_ppe_pcam_gtr"
        ),
    )
    parser.add_argument("--min-chromosome-mb", type=float, default=5.0)
    parser.add_argument("--min-block-score", type=int, default=10)
    parser.add_argument("--max-ribbons-per-homolog", type=int, default=5)
    parser.add_argument("--png-scale", type=float, default=1.0)
    return parser.parse_args()


def read_self_bed(path: Path) -> tuple[dict[str, tuple[str, int, int]], dict[str, int]]:
    genes: dict[str, tuple[str, int, int]] = {}
    lengths: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            seqid, start_text, end_text, gene = fields[:4]
            if gene in genes:
                continue
            start, end = int(start_text), int(end_text)
            genes[gene] = (seqid, start, end)
            lengths[seqid] = max(lengths[seqid], end)
    if not genes:
        raise ValueError(f"No self-comparison BED records: {path}")
    return genes, dict(lengths)


def self_fraction(record: tuple[str, int, int], lengths: dict[str, int]) -> float:
    seqid, start, end = record
    return ((start + end) / 2.0) / max(1, lengths[seqid])


def classify_block(
    qchr: str, q0: float, q1: float, tchr: str, t0: float, t1: float,
) -> str | None:
    """Select only the three user-specified copy1 regions.

    Coordinates are normalized to the named chromosome pair, so a reversed
    orientation in the JCVI anchor file does not change the biological label.
    """
    if {qchr, tchr} == {"chr11", "chr12"}:
        if qchr == "chr11":
            qchr, tchr, q0, q1, t0, t1 = tchr, qchr, t0, t1, q0, q1
        if qchr == "chr12" and tchr == "chr11":
            # One JCVI block falls in the middle/distal portion of chr12 and
            # the terminal portion of chr11; the other chr11-chr12 block does
            # not meet this distal criterion and remains grey.
            if q0 >= 0.50 and q1 <= 0.85 and t0 >= 0.70:
                return "chr12_mid_chr11_tail"
    if {qchr, tchr} == {"chr12", "chr16"}:
        if qchr == "chr16":
            qchr, tchr, q0, q1, t0, t1 = tchr, qchr, t0, t1, q0, q1
        if qchr == "chr12" and tchr == "chr16":
            if q1 <= 0.25 and t1 <= 0.25:
                return "chr12_head_chr16_head"
    if {qchr, tchr} == {"chr05", "chr15"}:
        if qchr == "chr15":
            qchr, tchr, q0, q1, t0, t1 = tchr, qchr, t0, t1, q0, q1
        if qchr == "chr05" and tchr == "chr15":
            if q1 <= 0.25 and t0 >= 0.70:
                return "chr05_head_chr15_tail"
    return None


def read_highlights(
    anchors: Path, self_bed: Path,
) -> list[dict[str, object]]:
    genes, lengths = read_self_bed(self_bed)
    selected: list[dict[str, object]] = []
    with anchors.open(encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            line = raw.strip()
            if not line or line.startswith("#") or line.startswith("StartGeneA"):
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6:
                continue
            qa, qb, ta, tb, score_text, orientation = fields
            if not all(gene in genes for gene in (qa, qb, ta, tb)):
                continue
            try:
                score = int(float(score_text))
            except ValueError:
                continue
            q = genes[qa]
            q_end = genes[qb]
            t = genes[ta]
            t_end = genes[tb]
            if q[0] != q_end[0] or t[0] != t_end[0]:
                continue
            q0, q1 = sorted((self_fraction(q, lengths), self_fraction(q_end, lengths)))
            t0, t1 = sorted((self_fraction(t, lengths), self_fraction(t_end, lengths)))
            label = classify_block(q[0], q0, q1, t[0], t0, t1)
            if not label:
                continue
            # Keep named chromosome direction and endpoints in the audit.
            if q[0] == "chr15" and t[0] == "chr05":
                q, q_end, t, t_end = t, t_end, q, q_end
                q0, q1, t0, t1 = t0, t1, q0, q1
            elif q[0] == "chr11" and t[0] == "chr12":
                q, q_end, t, t_end = t, t_end, q, q_end
                q0, q1, t0, t1 = t0, t1, q0, q1
            elif q[0] == "chr16" and t[0] == "chr12":
                q, q_end, t, t_end = t, t_end, q, q_end
                q0, q1, t0, t1 = t0, t1, q0, q1
            selected.append({
                "line": line_number,
                "class": label,
                "color": HIGHLIGHT_COLORS[label],
                "query_chr": q[0], "target_chr": t[0],
                "query_start": q0, "query_end": q1,
                "target_start": t0, "target_end": t1,
                "query_gene_a": q[0],  # replaced below with stable IDs
                "target_gene_a": t[0],
                "score": score, "orientation": orientation,
                "genes": (qa, qb, ta, tb),
            })
    # Correct the audit gene columns after the orientation normalization.
    for row in selected:
        qa, qb, ta, tb = row["genes"]
        if row["query_chr"] == "chr12" and genes[qa][0] == "chr11":
            qa, qb, ta, tb = ta, tb, qa, qb
        elif row["query_chr"] == "chr05" and genes[qa][0] == "chr15":
            qa, qb, ta, tb = ta, tb, qa, qb
        elif row["query_chr"] == "chr12" and genes[qa][0] == "chr16":
            qa, qb, ta, tb = ta, tb, qa, qb
        row["query_gene_a"], row["query_gene_b"] = qa, qb
        row["target_gene_a"], row["target_gene_b"] = ta, tb
    return selected


def self_ribbon_path(x1: float, x2: float, z1: float, z2: float, y: float, arch: float) -> str:
    x1, x2 = sorted((x1, x2))
    z1, z2 = sorted((z1, z2))
    # A narrow filled arch, with the two chromosome intervals as its bases.
    return (
        f"M {x1:.2f},{y:.2f} C {x1:.2f},{y-arch:.2f} {z1:.2f},{y-arch:.2f} {z1:.2f},{y:.2f} "
        f"L {z2:.2f},{y:.2f} C {z2:.2f},{y-arch:.2f} {x2:.2f},{y-arch:.2f} {x2:.2f},{y:.2f} Z"
    )


def draw_figure(
    outdir: Path,
    tracks: dict[str, list[base.TrackChromosome]],
    projections: dict[str, list[base.Projection]],
    pcam_targets: dict[str, str],
    highlights: list[dict[str, object]],
    copy1_map: dict[str, str],
    min_score: int,
    max_per_homolog: int,
    png_scale: float,
) -> tuple[Path, Path, Path]:
    outdir.mkdir(parents=True, exist_ok=True)
    width, top, row_gap, bar_h, bottom = 2320, 250, 154, 18, 118
    height = top + row_gap * (len(SPECIES) - 1) + bottom
    y_positions = {species: top + i * row_gap for i, species in enumerate(SPECIES)}
    lookup = base.track_lookup(tracks)
    ribbons: list[dict[str, object]] = []
    for upper, lower in zip(SPECIES, SPECIES[1:]):
        ribbons.extend(base.match_ribbons(
            upper, lower, projections, lookup, min_score, max_per_homolog
        ))

    # map copy1 chromosome names (chr01...) to displayed PZI C##_A1/B1 bars
    pzi_by_display = lookup["pzi"]
    pzi_map = {copy1: display for copy1, display in copy1_map.items()}
    pzi_highlights: list[dict[str, object]] = []
    for index, row in enumerate(highlights):
        qtrack = pzi_by_display.get(pzi_map.get(str(row["query_chr"]), ""))
        ttrack = pzi_by_display.get(pzi_map.get(str(row["target_chr"]), ""))
        if not qtrack or not ttrack:
            raise ValueError(f"Cannot map highlighted chromosome pair: {row}")
        qx = sorted((base.map_fraction(qtrack, float(row["query_start"])),
                     base.map_fraction(qtrack, float(row["query_end"]))))
        tx = sorted((base.map_fraction(ttrack, float(row["target_start"])),
                     base.map_fraction(ttrack, float(row["target_end"]))))
        pzi_highlights.append({**row, "query_track": qtrack.seqid, "target_track": ttrack.seqid,
                               "query_x": qx, "target_x": tx, "arch": 62 + 15 * index})

    svg: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif}.title{font-size:27px;font-weight:700;fill:#18212A}.subtitle{font-size:15px;fill:#5B6670}.species{font-size:18px;font-style:italic;fill:#17212B}.chrom{font-size:11px;font-weight:600;fill:#FFFFFF}.group{font-size:13px;font-weight:700;fill:#53606B}.legend{font-size:13px;fill:#2F3B45}.note{font-size:12px;fill:#66717B}</style>',
        '<text class="title" x="270" y="43">PZI copy1 internal synteny in a selected Prunus comparison panel</text>',
        '<text class="subtitle" x="270" y="69">PZI self-synteny highlights are coloured; interspecific JCVI ribbons are shown in grey</text>',
        '<line x1="1181" y1="92" x2="1181" y2="125" stroke="#AAB3BB" stroke-width="1"/>',
        '<text class="group" x="725" y="106" text-anchor="middle">PZI A1</text>',
        '<text class="group" x="1640" y="106" text-anchor="middle">PZI B1</text>',
    ]
    # Legend kept in the open upper-right area and intentionally concise.
    legend_x, legend_y = 1620, 34
    for i, label in enumerate(HIGHLIGHT_COLORS):
        y = legend_y + i * 20
        svg.append(f'<rect x="{legend_x}" y="{y-10}" width="12" height="12" rx="2" fill="{HIGHLIGHT_COLORS[label]}" fill-opacity="0.70"/>')
        svg.append(f'<text class="legend" x="{legend_x+18}" y="{y}">{html.escape(HIGHLIGHT_LABELS[label])}</text>')

    # Interspecific ribbons first, then colored PZI self-arches, then bars.
    for ribbon in sorted(ribbons, key=lambda r: float(r["quality"])):
        y1 = y_positions[str(ribbon["upper"])] + bar_h / 2
        y2 = y_positions[str(ribbon["lower"])] - bar_h / 2
        path = base.ribbon_path(float(ribbon["upper_x1"]), float(ribbon["upper_x2"]), y1,
                                float(ribbon["lower_x1"]), float(ribbon["lower_x2"]), y2)
        title = html.escape(f'{ribbon["upper"]}-{ribbon["lower"]} {ribbon["canonical"]}: {ribbon["upper_seqid"]} to {ribbon["lower_seqid"]}')
        svg.append(f'<path d="{path}" fill="#8798A8" fill-opacity="0.22" stroke="#718291" stroke-opacity="0.14" stroke-width="0.45"><title>{title}</title></path>')
    pzi_y = y_positions["pzi"]
    for row in pzi_highlights:
        path = self_ribbon_path(row["query_x"][0], row["query_x"][1], row["target_x"][0], row["target_x"][1], pzi_y - bar_h / 2, row["arch"])
        title = html.escape(f'{row["class"]}; {row["query_chr"]} to {row["target_chr"]}; score={row["score"]}; {row["orientation"]}')
        svg.append(f'<path d="{path}" fill="{row["color"]}" fill-opacity="0.72" stroke="{row["color"]}" stroke-opacity="0.94" stroke-width="3.2" stroke-linejoin="round"><title>{title}</title></path>')

    display_copy1 = {display: copy1 for copy1, display in copy1_map.items()}
    for species in SPECIES:
        y = y_positions[species]
        reference = " (reference)" if species == "pcam" else ""
        svg.append(f'<text class="species" x="245" y="{y+6:.2f}" text-anchor="end">{html.escape(base.SPECIES_NAMES[species])}{reference}</text>')
        if species in base.POLYPLOID_COMPONENTS:
            svg.append(f'<text class="group" x="278" y="{y-17:.2f}">{"A1" if species != "ppa" else "A"}</text>')
            svg.append(f'<text class="group" x="1205" y="{y-17:.2f}">{"B1" if species != "ppa" else "B"}</text>')
        for chrom in tracks[species]:
            fill = "#26323D" if chrom.group != "B" else "#34434F"
            svg.append(f'<rect x="{chrom.x1:.2f}" y="{y-bar_h/2:.2f}" width="{chrom.x2-chrom.x1:.2f}" height="{bar_h}" rx="8" fill="{fill}" stroke="#17212A" stroke-width="0.7"><title>{html.escape(chrom.seqid)}; {chrom.length:,} bp; orientation {"reverse" if chrom.reverse else "forward"}</title></rect>')
            label = display_copy1.get(chrom.seqid, chrom.display) if species == "pzi" else chrom.display
            svg.append(f'<text class="chrom" x="{(chrom.x1+chrom.x2)/2:.2f}" y="{y+4:.2f}" text-anchor="middle">{html.escape(label)}</text>')
    svg.extend([
        f'<text class="note" x="270" y="{height-39}">PZI internal blocks: JCVI self anchors; interspecific ribbons: PCAM-aligned projections; score threshold {min_score}.</text>',
        f'<text class="note" x="270" y="{height-20}">PZI chromosome labels are copy1 chr01-chr16; colored blocks are selected by the stated chromosome-end intervals.</text>',
        '</svg>',
    ])
    svg_path = outdir / "pzi_self_highlight_ppa_psp_ppe_pcam_gtr.svg"
    svg_path.write_text("\n".join(svg) + "\n", encoding="utf-8")

    # Raster preview is generated from the same coordinates in Python/Pillow.
    scale = png_scale
    image = Image.new("RGB", (round(width*scale), round(height*scale)), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    try:
        font_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"
        bold_path = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
        font = ImageFont.truetype(font_path, max(8, round(15*scale)))
        small = ImageFont.truetype(font_path, max(7, round(11*scale)))
        title_font = ImageFont.truetype(bold_path, max(12, round(26*scale)))
    except OSError:
        font = small = title_font = ImageFont.load_default()
    def S(point: tuple[float, float]) -> tuple[int, int]:
        return round(point[0]*scale), round(point[1]*scale)
    draw.text(S((270, 30)), "PZI copy1 internal synteny in a selected Prunus comparison panel", fill=(24,33,42,255), font=title_font)
    for ribbon in sorted(ribbons, key=lambda r: float(r["quality"])):
        y1 = y_positions[str(ribbon["upper"])] + bar_h / 2
        y2 = y_positions[str(ribbon["lower"])] - bar_h / 2
        control = (y2-y1)*0.43
        x1, x2 = float(ribbon["upper_x1"]), float(ribbon["upper_x2"])
        z1, z2 = float(ribbon["lower_x1"]), float(ribbon["lower_x2"])
        left_curve = list(base.bezier((x1,y1),(x1,y1+control),(z1,y2-control),(z1,y2)))
        right_curve = list(base.bezier((z2,y2),(z2,y2-control),(x2,y1+control),(x2,y1)))
        draw.polygon([S(p) for p in left_curve+right_curve], fill=(135,152,168,56))
    for row in pzi_highlights:
        col = row["color"]
        rgb = tuple(int(col[i:i+2],16) for i in (1,3,5))
        qx1, qx2 = row["query_x"]
        tx1, tx2 = row["target_x"]
        arch = float(row["arch"])
        left_curve = list(base.bezier((qx1, pzi_y-bar_h/2), (qx1, pzi_y-bar_h/2-arch), (tx1, pzi_y-bar_h/2-arch), (tx1, pzi_y-bar_h/2)))
        right_curve = list(base.bezier((tx2, pzi_y-bar_h/2), (tx2, pzi_y-bar_h/2-arch), (qx2, pzi_y-bar_h/2-arch), (qx2, pzi_y-bar_h/2)))
        draw.polygon([S(p) for p in left_curve + right_curve], fill=(*rgb,175))
    for species in SPECIES:
        y = y_positions[species]
        label = base.SPECIES_NAMES[species] + (" (reference)" if species == "pcam" else "")
        bbox = draw.textbbox((0,0), label, font=font)
        draw.text(S((245-(bbox[2]-bbox[0])/scale,y-7)), label, fill=(23,33,43,255), font=font)
        for chrom in tracks[species]:
            fill = (38,50,61,255) if chrom.group != "B" else (52,67,79,255)
            draw.rounded_rectangle([S((chrom.x1,y-bar_h/2)),S((chrom.x2,y+bar_h/2))], radius=max(2,round(7*scale)), fill=fill, outline=(23,33,42,255), width=1)
            label = display_copy1.get(chrom.seqid, chrom.display) if species == "pzi" else chrom.display
            bbox = draw.textbbox((0,0), label, font=small)
            draw.text(S(((chrom.x1+chrom.x2)/2-(bbox[2]-bbox[0])/(2*scale),y-6)), label, fill=(255,255,255,255), font=small)
    png_path = outdir / "pzi_self_highlight_ppa_psp_ppe_pcam_gtr.png"
    image.save(png_path, dpi=(300,300))

    # Keep a PDF deliverable in the selected Python backend.  The SVG remains
    # the editable/vector master; the PDF is a high-resolution raster export.
    import matplotlib.pyplot as plt
    fig = plt.figure(figsize=(width/100, height/100), dpi=100)
    fig.figimage(image, xo=0, yo=0)
    fig.savefig(outdir / "pzi_self_highlight_ppa_psp_ppe_pcam_gtr.pdf", dpi=300)
    plt.close(fig)
    return svg_path, png_path, outdir / "pzi_self_highlight_ppa_psp_ppe_pcam_gtr.pdf"


def write_audits(outdir: Path, highlights: list[dict[str, object]]) -> None:
    columns = [
        "highlight_class", "source_anchor_line", "query_chr", "target_chr",
        "query_start_fraction", "query_end_fraction", "target_start_fraction",
        "target_end_fraction", "query_gene_a", "query_gene_b", "target_gene_a",
        "target_gene_b", "score", "orientation",
    ]
    with (outdir / "highlight_blocks.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in highlights:
            writer.writerow({
                "highlight_class": row["class"], "source_anchor_line": row["line"],
                "query_chr": row["query_chr"], "target_chr": row["target_chr"],
                "query_start_fraction": f'{float(row["query_start"]):.6f}',
                "query_end_fraction": f'{float(row["query_end"]):.6f}',
                "target_start_fraction": f'{float(row["target_start"]):.6f}',
                "target_end_fraction": f'{float(row["target_end"]):.6f}',
                "query_gene_a": row["query_gene_a"], "query_gene_b": row["query_gene_b"],
                "target_gene_a": row["target_gene_a"], "target_gene_b": row["target_gene_b"],
                "score": row["score"], "orientation": row["orientation"],
            })


def main() -> int:
    args = parse_args()
    for path in (args.input_dir, args.self_bed, args.self_anchors):
        if not path.exists():
            raise FileNotFoundError(path)
    base.SPECIES = SPECIES
    tracks, projections, _, _, pcam_targets = base.build_tracks_and_projections(
        args.input_dir, round(args.min_chromosome_mb * 1_000_000)
    )
    base.assign_geometry(tracks)
    highlights = read_highlights(args.self_anchors, args.self_bed)
    expected = {
        "chr12_mid_chr11_tail", "chr12_head_chr16_head", "chr05_head_chr15_tail",
    }
    observed = {str(row["class"]) for row in highlights}
    if not expected <= observed:
        raise ValueError(f"Missing highlight class(es): {sorted(expected-observed)}")
    copy1_map = base.read_copy1_map(args.input_dir / "maps" / "chromosome_rename_map.tsv")
    svg, png, pdf = draw_figure(
        args.outdir, tracks, projections, pcam_targets, highlights, copy1_map,
        args.min_block_score, args.max_ribbons_per_homolog, args.png_scale,
    )
    write_audits(args.outdir, highlights)
    (args.outdir / "README.txt").write_text(
        "PZI copy1 self-highlight figure\n\n"
        "Rows: PZI reference plus PPA, PSP, PPE, PCAM and GTR.\n"
        "Grey ribbons are existing PCAM-aligned JCVI projections.\n"
        "Colored arches are selected PZI copy1 self-anchor blocks; see highlight_blocks.tsv.\n"
        "The PZI copy1 track is labeled chr01-chr16 and ordered by the existing PCAM model.\n",
        encoding="utf-8",
    )
    print(f"[OK] SVG: {svg}")
    print(f"[OK] PNG: {png}")
    print(f"[OK] PDF: {pdf}")
    print(f"[OK] highlighted blocks: {len(highlights)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
