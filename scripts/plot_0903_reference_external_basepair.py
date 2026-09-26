#!/usr/bin/env python3
"""Render the 0903 JCVI ribbons in the reference-external visual style.

This renderer reuses the BED files and screened JCVI projections prepared by
``plot_0903_multispecies_A1B1_ribbon.py``.  It changes only the presentation:
chromosome bars are laid out with widths proportional to their BED lengths,
and the ribbons use the compact external-reference style.  No alignments or
synteny blocks are recalculated here.
"""

from __future__ import annotations

import argparse
import csv
import html
import math
import sys
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.colors import HexColor
from reportlab.pdfgen import canvas

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

import plot_0903_multispecies_A1B1_ribbon as source  # noqa: E402


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hap-input", type=Path, default=Path("0903/haploid_repartition/inputs"))
    parser.add_argument("--jcvi-input", type=Path, default=Path("0903/downstream_analysis/figure2_repartition_ks/inputs"))
    parser.add_argument("--pair-work", type=Path, default=Path("0903/downstream_analysis/figure2_repartition_ks/01_jcvi_pairs/work"))
    parser.add_argument("--rch-input", type=Path, default=Path("downstream_analysis/local_A1B1_ribbon_inputs"))
    parser.add_argument(
        "--outdir", type=Path,
        default=Path("analysis/karyotype/jcvi/figures/custom_svg/ribbons/native_multispecies_A1B1_pcam_pairwise"),
    )
    parser.add_argument("--min-score", type=int, default=4)
    parser.add_argument("--max-ribbons-per-homolog", type=int, default=4)
    parser.add_argument("--width", type=int, default=1500)
    parser.add_argument("--height", type=int, default=1900)
    parser.add_argument("--left", type=float, default=205.0)
    parser.add_argument("--right", type=float, default=1445.0)
    parser.add_argument("--row-gap", type=float, default=270.0)
    parser.add_argument("--bar-height", type=float, default=9.0)
    parser.add_argument("--chrom-gap", type=float, default=4.0)
    parser.add_argument("--ribbon-alpha", type=float, default=0.48)
    parser.add_argument("--png-scale", type=int, default=3)
    return parser.parse_args()


def configure_tracks(tracks, width: int, left: float, right: float, chrom_gap: float) -> None:
    """Assign proportional x coordinates to every chromosome bar.

    Polyploid tracks retain an A/B split with a central gap.  Within each
    group, the available width is allocated by the actual BED chromosome
    lengths; diploid tracks use the complete horizontal span in the same way.
    """
    polyploid = source.POLYPLOID
    center_gap = 74.0
    for species, chroms in tracks.items():
        if species in polyploid:
            groups = [("A", [c for c in chroms if c.group == "A"]),
                      ("B", [c for c in chroms if c.group == "B"])]
            half = (right - left - center_gap) / 2.0
            for group, values in groups:
                x0 = left if group == "A" else left + half + center_gap
                available = max(1.0, half - chrom_gap * max(0, len(values) - 1))
                total = max(1, sum(c.length for c in values))
                cursor = x0
                for index, chrom in enumerate(values):
                    span = available * chrom.length / total
                    chrom.x1 = cursor
                    chrom.x2 = cursor + span
                    cursor += span + chrom_gap
        else:
            available = max(1.0, right - left - chrom_gap * max(0, len(chroms) - 1))
            total = max(1, sum(c.length for c in chroms))
            cursor = left
            for chrom in chroms:
                span = available * chrom.length / total
                chrom.x1 = cursor
                chrom.x2 = cursor + span
                cursor += span + chrom_gap


def ribbon_polygon(ribbon: dict[str, object], y1: float, y2: float) -> list[tuple[float, float]]:
    """Return a smooth quadrilateral polygon matching the reference style."""
    x1, x2 = float(ribbon["upper_x1"]), float(ribbon["upper_x2"])
    z1, z2 = float(ribbon["lower_x1"]), float(ribbon["lower_x2"])
    control = (y2 - y1) * 0.43

    def cubic(p0, p1, p2, p3, n=22):
        for i in range(n + 1):
            t = i / n
            u = 1.0 - t
            yield (
                u**3 * p0[0] + 3 * u * u * t * p1[0] + 3 * u * t * t * p2[0] + t**3 * p3[0],
                u**3 * p0[1] + 3 * u * u * t * p1[1] + 3 * u * t * t * p2[1] + t**3 * p3[1],
            )

    left_curve = list(cubic((x1, y1), (x1, y1 + control), (z1, y2 - control), (z1, y2)))
    right_curve = list(cubic((z2, y2), (z2, y2 - control), (x2, y1 + control), (x2, y1)))
    return left_curve + right_curve


def collect_ribbons(tracks, projections, min_score: int, cap: int):
    ribbons = []
    for upper, lower in zip(source.SPECIES, source.SPECIES[1:]):
        for ribbon in source.match_pair(upper, lower, projections, tracks, min_score, cap):
            ribbon = dict(ribbon)
            ribbon["upper"] = upper
            ribbon["lower"] = lower
            ribbons.append(ribbon)
    return ribbons


def write_audit(path: Path, tracks, ribbons, args: argparse.Namespace) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["record_type", "species", "group", "display", "source_seqid", "length_bp", "x1", "x2", "upper", "lower", "canonical", "quality"])
        for species in source.SPECIES:
            for chrom in tracks[species]:
                writer.writerow(["chromosome", species, chrom.group, chrom.display, chrom.source_seqid, chrom.length, f"{chrom.x1:.3f}", f"{chrom.x2:.3f}", "", "", "", ""])
        for ribbon in ribbons:
            writer.writerow(["ribbon", "", "", "", "", "", f"{float(ribbon['upper_x1']):.3f}", f"{float(ribbon['upper_x2']):.3f}", ribbon["upper"], ribbon["lower"], ribbon["canonical"], f"{float(ribbon['quality']):.5f}"])


def render_svg(path: Path, tracks, ribbons, args: argparse.Namespace) -> None:
    top = 130.0
    y = {species: top + i * args.row_gap for i, species in enumerate(source.SPECIES)}
    bar_y = {species: y[species] for species in source.SPECIES}
    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{args.width}" height="{args.height}" viewBox="0 0 {args.width} {args.height}">',
        '<rect width="100%" height="100%" fill="#FFFFFF"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif}.species{font-size:18px;font-style:italic;fill:#26323D}.chrom{font-size:10px;font-weight:600;fill:#FFFFFF}.group{font-size:12px;font-weight:700;fill:#5D6872}.note{font-size:11px;fill:#68747E}</style>',
    ]
    # Ribbons are behind the chromosome bars, as in the reference_external figures.
    for ribbon in sorted(ribbons, key=lambda r: float(r["quality"])):
        upper_y = bar_y[str(ribbon["upper"])] + args.bar_height / 2
        lower_y = bar_y[str(ribbon["lower"])] - args.bar_height / 2
        points = ribbon_polygon(ribbon, upper_y, lower_y)
        d = "M " + " ".join(f"{x:.2f},{yy:.2f}" for x, yy in points) + " Z"
        color = source.PALETTE.get(str(ribbon["canonical"]).split("_", 1)[0], "#8798A8")
        title = html.escape(f'{ribbon["upper"]}-{ribbon["lower"]} {ribbon["canonical"]}')
        lines.append(f'<path d="{d}" fill="{color}" fill-opacity="{args.ribbon_alpha:.3f}" stroke="{color}" stroke-opacity="0.22" stroke-width="0.55"><title>{title}</title></path>')

    for species in source.SPECIES:
        yy = bar_y[species]
        lines.append(f'<text class="species" x="{args.left - 22:.2f}" y="{yy + 6:.2f}" text-anchor="end">{html.escape(species)}</text>')
        if species in source.POLYPLOID:
            a = [c for c in tracks[species] if c.group == "A"]
            b = [c for c in tracks[species] if c.group == "B"]
            lines.append(f'<text class="group" x="{(a[0].x1 + a[-1].x2)/2:.2f}" y="{yy - 20:.2f}" text-anchor="middle">A1</text>')
            lines.append(f'<text class="group" x="{(b[0].x1 + b[-1].x2)/2:.2f}" y="{yy - 20:.2f}" text-anchor="middle">B1</text>')
        for chrom in tracks[species]:
            width = max(1.0, chrom.x2 - chrom.x1)
            lines.append(f'<rect x="{chrom.x1:.2f}" y="{yy - args.bar_height/2:.2f}" width="{width:.2f}" height="{args.bar_height:.2f}" rx="4" fill="#4C5660" stroke="#1E2730" stroke-width="0.8"><title>{html.escape(chrom.source_seqid)}; {chrom.length:,} bp</title></rect>')
            lines.append(f'<text class="chrom" x="{(chrom.x1 + chrom.x2)/2:.2f}" y="{yy + 3.5:.2f}" text-anchor="middle">{html.escape(chrom.display)}</text>')
    lines.append(f'<text class="note" x="{args.left:.2f}" y="{args.height - 28:.2f}">Chromosome widths are proportional to BED chromosome length; ribbons are screened JCVI projections.</text>')
    lines.append('</svg>')
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_pdf(path: Path, tracks, ribbons, args: argparse.Namespace) -> None:
    top = 130.0
    y = {species: top + i * args.row_gap for i, species in enumerate(source.SPECIES)}
    pdf = canvas.Canvas(str(path), pagesize=(args.width, args.height))
    pdf.setFillColor(HexColor("#FFFFFF"))
    pdf.rect(0, 0, args.width, args.height, fill=1, stroke=0)
    for ribbon in sorted(ribbons, key=lambda r: float(r["quality"])):
        upper_y = y[str(ribbon["upper"])] + args.bar_height / 2
        lower_y = y[str(ribbon["lower"])] - args.bar_height / 2
        points = ribbon_polygon(ribbon, upper_y, lower_y)
        p = pdf.beginPath()
        p.moveTo(points[0][0], args.height - points[0][1])
        for px, py in points[1:]:
            p.lineTo(px, args.height - py)
        p.close()
        color = source.PALETTE.get(str(ribbon["canonical"]).split("_", 1)[0], "#8798A8")
        pdf.setFillColor(HexColor(color))
        pdf.setFillAlpha(args.ribbon_alpha)
        pdf.drawPath(p, fill=1, stroke=0)
    pdf.setFillAlpha(1.0)
    for species in source.SPECIES:
        yy = y[species]
        pdf.setFillColor(HexColor("#26323D"))
        pdf.setFont("Helvetica-Oblique", 12)
        pdf.drawRightString(args.left - 22, args.height - yy + 4, species)
        for chrom in tracks[species]:
            pdf.setFillColor(HexColor("#4C5660"))
            pdf.setStrokeColor(HexColor("#1E2730"))
            pdf.roundRect(chrom.x1, args.height - yy - args.bar_height/2, chrom.x2 - chrom.x1, args.bar_height, 4, fill=1, stroke=1)
            pdf.setFillColor(HexColor("#FFFFFF"))
            pdf.setFont("Helvetica-Bold", 6.5)
            pdf.drawCentredString((chrom.x1 + chrom.x2) / 2, args.height - yy - 2.3, chrom.display)
    pdf.save()


def render_png(path: Path, tracks, ribbons, args: argparse.Namespace) -> None:
    scale = args.png_scale
    image = Image.new("RGBA", (args.width * scale, args.height * scale), (255, 255, 255, 255))
    draw = ImageDraw.Draw(image, "RGBA")
    try:
        species_font = ImageFont.truetype("ariali.ttf", 17 * scale)
        chrom_font = ImageFont.truetype("arialbd.ttf", 10 * scale)
    except OSError:
        species_font = chrom_font = ImageFont.load_default()
    top = 130.0
    y = {species: top + i * args.row_gap for i, species in enumerate(source.SPECIES)}
    for ribbon in sorted(ribbons, key=lambda r: float(r["quality"])):
        points = ribbon_polygon(ribbon, y[str(ribbon["upper"])] + args.bar_height/2, y[str(ribbon["lower"])] - args.bar_height/2)
        color = source.PALETTE.get(str(ribbon["canonical"]).split("_", 1)[0], "#8798A8")
        rgb = tuple(int(color[i:i+2], 16) for i in (1, 3, 5))
        draw.polygon([(round(px * scale), round(py * scale)) for px, py in points], fill=(*rgb, round(args.ribbon_alpha * 255)))
    for species in source.SPECIES:
        yy = y[species]
        draw.text(((args.left - 22) * scale, (yy - 8) * scale), species, fill=(38, 50, 61, 255), font=species_font, anchor="ra")
        for chrom in tracks[species]:
            box = (round(chrom.x1 * scale), round((yy - args.bar_height/2) * scale), round(chrom.x2 * scale), round((yy + args.bar_height/2) * scale))
            draw.rounded_rectangle(box, radius=4 * scale, fill=(76, 86, 96, 255), outline=(30, 39, 48, 255), width=max(1, scale))
            draw.text((round((chrom.x1 + chrom.x2) * scale / 2), round((yy - 5) * scale)), chrom.display, fill=(255, 255, 255, 255), font=chrom_font, anchor="mm")
    image.convert("RGB").save(path, dpi=(300, 300))


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    beds, tracks, projections, _paths, _rch_pairs = source.build_tracks_and_projections(args)
    configure_tracks(tracks, args.width, args.left, args.right, args.chrom_gap)
    ribbons = collect_ribbons(tracks, projections, args.min_score, args.max_ribbons_per_homolog)
    prefix = args.outdir / "multispecies_A1B1_pcam.pairwise.reference_external.basepair"
    # ``Path.with_suffix`` would replace the ``.basepair`` part of this
    # descriptive stem.  Append the extension explicitly so all three files
    # retain the scale/style identifier in their names.
    svg_path = Path(str(prefix) + ".svg")
    pdf_path = Path(str(prefix) + ".pdf")
    png_path = Path(str(prefix) + ".png")
    render_svg(svg_path, tracks, ribbons, args)
    render_pdf(pdf_path, tracks, ribbons, args)
    render_png(png_path, tracks, ribbons, args)
    write_audit(args.outdir / "reference_external_basepair_audit.tsv", tracks, ribbons, args)
    (args.outdir / "reference_external_basepair_README.txt").write_text(
        "Reference-external style, base-pair-scaled multispecies JCVI ribbon figure\n"
        "Species order: pzi, psp, ppa, pcam, ppe, gtr, rch.\n"
        "Chromosome bar widths are proportional to the chromosome lengths read from the JCVI BED files.\n"
        "Ribbons reuse the screened JCVI projections from the 0903 preparation; no new alignment was inferred.\n"
        f"min_score={args.min_score}; max_ribbons_per_homolog={args.max_ribbons_per_homolog}; chrom_gap={args.chrom_gap}px.\n",
        encoding="utf-8",
    )
    print(f"[OK] SVG: {svg_path}")
    print(f"[OK] PDF: {pdf_path}")
    print(f"[OK] PNG: {png_path}")
    print(f"[OK] ribbons: {len(ribbons)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
