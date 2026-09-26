#!/usr/bin/env python3
"""Draw a simple SVG heatmap from minimap2 chromosome-pair summary tables."""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def natural_key(text: str):
    parts = []
    buf = ""
    is_digit = False
    for ch in text:
        if ch.isdigit() != is_digit and buf:
            parts.append(int(buf) if is_digit else buf)
            buf = ch
            is_digit = ch.isdigit()
        else:
            buf += ch
            is_digit = ch.isdigit()
    if buf:
        parts.append(int(buf) if is_digit else buf)
    return parts


def get_first(row, *names):
    for name in names:
        if name in row and row[name] != "":
            return row[name]
    raise KeyError(names[0])


def read_summary(path: Path):
    matrix = defaultdict(float)
    ref_order = []
    query_order = []
    seen_ref = set()
    seen_query = set()
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            ref = get_first(
                row,
                "reference_chr",
                "target_chr",
                "pcam_chr",
                "ref_chr",
                "chromosome",
            )
            query = get_first(
                row,
                "query_chr",
                "pzi_chr",
                "query",
                "chromosome_query",
            )
            value = float(get_first(
                row,
                "ppa_union_covered_bp",
                "aligned_block_bp",
                "aligned_bp",
                "block_bp",
                "bp",
            ))
            matrix[(ref, query)] += value
            if ref not in seen_ref:
                ref_order.append(ref)
                seen_ref.add(ref)
            if query not in seen_query:
                query_order.append(query)
                seen_query.add(query)
    ref_order = sorted(ref_order, key=natural_key)
    query_order = sorted(query_order, key=natural_key)
    return matrix, ref_order, query_order


def color_for(value: float, vmax: float) -> str:
    if value <= 0 or vmax <= 0:
        return "#F7FAFC"
    frac = math.log10(value + 1) / math.log10(vmax + 1)
    frac = max(0.0, min(1.0, frac))
    # White -> yellow -> red.
    if frac < 0.5:
        t = frac / 0.5
        r = int(255)
        g = int(255 - (255 - 214) * t)
        b = int(255 - (255 - 102) * t)
    else:
        t = (frac - 0.5) / 0.5
        r = int(255 - (255 - 165) * t)
        g = int(214 - (214 - 0) * t)
        b = int(102 - (102 - 38) * t)
    return f"#{r:02X}{g:02X}{b:02X}"


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def draw_svg(
    matrix,
    ref_order,
    query_order,
    out_svg: Path,
    title: str,
    x_label: str,
    y_label: str,
    ref_colors=None,
    query_colors=None,
):
    cell = 30
    band = 12
    has_axis_tracks = bool(ref_colors) or bool(query_colors)
    left = 225 if has_axis_tracks else 150
    top = 185 if has_axis_tracks else 120
    right = 210
    bottom = 110
    width = left + len(query_order) * cell + right
    height = top + len(ref_order) * cell + bottom
    vmax = max(matrix.values()) if matrix else 0

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="32" font-family="Arial, sans-serif" font-size="20" font-weight="700">{esc(title)}</text>',
        f'<text x="{left}" y="58" font-family="Arial, sans-serif" font-size="12" fill="#4A5568">Cell color shows log-scaled PPA chromosome coverage by primary minimap2 alignment chains.</text>',
    ]

    query_band_y = top - band - 4
    query_label_y = top - band - 28
    ref_band_x = left - band - 4
    ref_label_x = left - band - 24

    for j, query in enumerate(query_order):
        x = left + j * cell + cell / 2
        lines.append(
            f'<text x="{x}" y="{query_label_y}" font-family="Arial, sans-serif" font-size="11" '
            f'text-anchor="end" transform="rotate(-45 {x} {query_label_y})">{esc(query)}</text>'
        )
        if query_colors:
            color = query_colors.get(query, "#CBD5E0")
            lines.append(
                f'<rect x="{left + j * cell}" y="{query_band_y}" width="{cell}" height="{band}" fill="{color}" stroke="#FFFFFF" stroke-width="0.2"/>'
            )

    for i, ref in enumerate(ref_order):
        y = top + i * cell + cell / 2 + 4
        lines.append(
            f'<text x="{ref_label_x}" y="{y}" font-family="Arial, sans-serif" font-size="11" '
            f'text-anchor="end">{esc(ref)}</text>'
        )
        if ref_colors:
            color = ref_colors.get(ref, "#CBD5E0")
            lines.append(
                f'<rect x="{ref_band_x}" y="{top + i * cell}" width="{band}" height="{cell}" fill="{color}" stroke="#FFFFFF" stroke-width="0.2"/>'
            )

    for i, ref in enumerate(ref_order):
        for j, query in enumerate(query_order):
            value = matrix.get((ref, query), 0.0)
            x = left + j * cell
            y = top + i * cell
            fill = color_for(value, vmax)
            label = f"{value / 1_000_000:.2f} Mb" if value else "0"
            lines.append(
                f'<rect x="{x}" y="{y}" width="{cell}" height="{cell}" fill="{fill}" '
                f'stroke="#CBD5E0" stroke-width="0.5"><title>{esc(ref)} vs {esc(query)}: {label}</title></rect>'
            )

    lines.extend(
        [
            f'<text x="{left + len(query_order) * cell / 2}" y="{height - 30}" font-family="Arial, sans-serif" font-size="13" text-anchor="middle">{esc(x_label)}</text>',
            f'<text x="24" y="{top + len(ref_order) * cell / 2}" font-family="Arial, sans-serif" font-size="13" text-anchor="middle" transform="rotate(-90 24 {top + len(ref_order) * cell / 2})">{esc(y_label)}</text>',
        ]
    )

    legend_x = left + len(query_order) * cell + 40
    legend_y = top
    lines.append(f'<text x="{legend_x}" y="{legend_y - 14}" font-family="Arial, sans-serif" font-size="12" font-weight="700">PPA covered bp</text>')
    for k in range(8):
        frac = k / 7
        val = 10 ** (frac * math.log10(vmax + 1)) - 1 if vmax > 0 else 0
        lines.append(
            f'<rect x="{legend_x}" y="{legend_y + k * 20}" width="22" height="18" fill="{color_for(val, vmax)}" stroke="#CBD5E0" stroke-width="0.5"/>'
        )
    lines.append(f'<text x="{legend_x + 30}" y="{legend_y + 12}" font-family="Arial, sans-serif" font-size="11">{vmax / 1_000_000:.1f} Mb</text>')
    lines.append(f'<text x="{legend_x + 30}" y="{legend_y + 7 * 20 + 12}" font-family="Arial, sans-serif" font-size="11">0</text>')

    lines.append("</svg>")
    out_svg.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--out-svg", required=True, type=Path)
    parser.add_argument("--title", required=True)
    parser.add_argument("--x-label", default="Prunus zippeliana chromosomes")
    parser.add_argument("--y-label", default="Peach reference sequences")
    parser.add_argument("--axis-colors", type=Path, help="Optional TSV with chr and color columns for both axes")
    args = parser.parse_args()

    matrix, ref_order, query_order = read_summary(args.summary)
    ref_colors = None
    query_colors = None
    if args.axis_colors and args.axis_colors.is_file():
        ref_colors = {}
        query_colors = {}
        with args.axis_colors.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                chr_name = row.get("chr", "")
                color = row.get("color", "")
                axis = row.get("axis", "both")
                if not chr_name or not color:
                    continue
                if axis in {"both", "ref", "target"}:
                    ref_colors[chr_name] = color
                if axis in {"both", "query"}:
                    query_colors[chr_name] = color
    draw_svg(
        matrix,
        ref_order,
        query_order,
        args.out_svg,
        args.title,
        args.x_label,
        args.y_label,
        ref_colors=ref_colors,
        query_colors=query_colors,
    )


if __name__ == "__main__":
    main()
