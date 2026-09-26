#!/usr/bin/env python3
"""Draw a minimap2 PAF genome synteny dotplot as SVG.

The x-axis is the target/reference genome and the y-axis is the query genome.
Each retained PAF alignment is drawn as a line segment colored by alignment
identity. This is intended for macro-synteny inspection, not publication-grade
gene-level collinearity.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path
from typing import Optional


def natural_key(text: str):
    parts = []
    buf = ""
    was_digit = None
    for ch in text:
        is_digit = ch.isdigit()
        if was_digit is not None and is_digit != was_digit:
            parts.append(int(buf) if was_digit else buf)
            buf = ch
        else:
            buf += ch
        was_digit = is_digit
    if buf:
        parts.append(int(buf) if was_digit else buf)
    return parts


def esc(text: str) -> str:
    return (
        text.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def parse_paf(path: Path):
    records = []
    target_lengths = {}
    query_lengths = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            qlen = int(fields[1])
            qstart = int(fields[2])
            qend = int(fields[3])
            strand = fields[4]
            tname = fields[5]
            tlen = int(fields[6])
            tstart = int(fields[7])
            tend = int(fields[8])
            matches = int(fields[9])
            block_len = int(fields[10])
            mapq = int(fields[11])
            ident = matches / block_len if block_len else 0.0
            records.append(
                {
                    "query": qname,
                    "query_len": qlen,
                    "query_start": qstart,
                    "query_end": qend,
                    "strand": strand,
                    "target": tname,
                    "target_len": tlen,
                    "target_start": tstart,
                    "target_end": tend,
                    "matches": matches,
                    "block_len": block_len,
                    "mapq": mapq,
                    "identity": ident,
                }
            )
            target_lengths[tname] = max(target_lengths.get(tname, 0), tlen)
            query_lengths[qname] = max(query_lengths.get(qname, 0), qlen)
    return records, target_lengths, query_lengths


def aggregate(records, key_name):
    agg = defaultdict(int)
    for rec in records:
        agg[rec[key_name]] += rec["block_len"]
    return agg


def read_axis_segments(path: Optional[Path]):
    if not path or not path.is_file():
        return {"target": defaultdict(list), "query": defaultdict(list)}
    segments = {"target": defaultdict(list), "query": defaultdict(list)}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            axis = row.get("axis", "")
            chr_name = row.get("chr", "")
            if axis not in segments or not chr_name:
                continue
            try:
                start = int(float(row.get("start", "0")))
                end = int(float(row.get("end", "0")))
            except ValueError:
                continue
            if end <= start:
                continue
            segments[axis][chr_name].append({
                "start": start,
                "end": end,
                "color": row.get("color", "#CBD5E0"),
                "label": row.get("label", ""),
            })
    for axis in segments:
        for chr_name in segments[axis]:
            segments[axis][chr_name].sort(key=lambda r: (r["start"], r["end"], r["label"]))
    return segments


def build_offsets(names, lengths, gap):
    offsets = {}
    pos = 0
    for name in names:
        offsets[name] = pos
        pos += lengths[name] + gap
    return offsets, pos - gap if names else 0


def color_for_identity(identity: float, vmin: float, vmax: float) -> str:
    if vmax <= vmin:
        frac = 1.0
    else:
        frac = max(0.0, min(1.0, (identity - vmin) / (vmax - vmin)))
    # Blue -> green -> yellow -> red.
    stops = [
        (0.0, (49, 130, 189)),
        (0.33, (102, 194, 165)),
        (0.66, (255, 237, 111)),
        (1.0, (215, 48, 39)),
    ]
    for (x0, c0), (x1, c1) in zip(stops, stops[1:]):
        if frac <= x1:
            t = 0 if x1 == x0 else (frac - x0) / (x1 - x0)
            rgb = tuple(round(c0[i] + (c1[i] - c0[i]) * t) for i in range(3))
            return f"#{rgb[0]:02X}{rgb[1]:02X}{rgb[2]:02X}"
    return "#D73027"


def tick_label_identity(value: float) -> str:
    return f"{value:.2f}"


def draw_dotplot(
    records,
    target_lengths,
    query_lengths,
    out_svg: Path,
    title: str,
    min_alignment_len: int,
    min_query_aggregate: int,
    min_target_aggregate: int,
    target_colors=None,
    query_colors=None,
    axis_segments=None,
    min_visible_px: float = 0.0,
):
    records = [rec for rec in records if rec["block_len"] >= min_alignment_len]
    query_agg = aggregate(records, "query")
    target_agg = aggregate(records, "target")
    records = [
        rec
        for rec in records
        if query_agg[rec["query"]] >= min_query_aggregate
        and target_agg[rec["target"]] >= min_target_aggregate
    ]

    query_names = sorted({rec["query"] for rec in records}, key=natural_key)
    target_names = sorted({rec["target"] for rec in records}, key=natural_key)
    query_lengths = {name: query_lengths[name] for name in query_names}
    target_lengths = {name: target_lengths[name] for name in target_names}

    gap = 3_000_000
    target_offsets, target_total = build_offsets(target_names, target_lengths, gap)
    query_offsets, query_total = build_offsets(query_names, query_lengths, gap)

    plot_w = 980
    plot_h = max(760, int(plot_w * (query_total / target_total))) if target_total else 760
    has_axis_tracks = bool(axis_segments) or bool(target_colors) or bool(query_colors)
    left = 240 if has_axis_tracks else 170
    top = 160 if has_axis_tracks else 100
    right = 180
    bottom = 130
    width = left + plot_w + right
    height = top + plot_h + bottom

    x_scale = plot_w / target_total if target_total else 1
    y_scale = plot_h / query_total if query_total else 1

    identities = [rec["identity"] for rec in records]
    vmin = min(identities) if identities else 0
    vmax = max(identities) if identities else 1

    def x_of(target, pos):
        return left + (target_offsets[target] + pos) * x_scale

    def y_of(query, pos):
        return top + plot_h - (query_offsets[query] + pos) * y_scale

    lines = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="{left}" y="24" font-family="Arial, sans-serif" font-size="13">Post-filtering number of alignments: {len(records)}    minimum alignment length (-m): {min_alignment_len:.0e}</text>',
        f'<text x="{left}" y="43" font-family="Arial, sans-serif" font-size="13">Post-filtering number of queries: {len(query_names)}    minimum query aggregate alignment length (-q): {min_query_aggregate:.0e}</text>',
        f'<text x="{left}" y="62" font-family="Arial, sans-serif" font-size="13">Post-filtering number of targets: {len(target_names)}    minimum target aggregate alignment length: {min_target_aggregate:.0e}</text>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="white" stroke="#9CA3AF" stroke-width="1"/>',
    ]
    axis_segments = axis_segments or {"target": {}, "query": {}}
    target_band_y = top - 30
    target_label_y = top - 50
    query_band_x = left - 26
    query_label_x = left - 38

    for name in target_names:
        x = x_of(name, 0)
        lines.append(f'<line x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}" stroke="#E5E7EB" stroke-width="1"/>')
        target_segment_rows = axis_segments.get("target", {}).get(name, [])
        if target_segment_rows:
            for seg in target_segment_rows:
                sx = x_of(name, max(0, min(seg["start"], target_lengths[name])))
                ex = x_of(name, max(0, min(seg["end"], target_lengths[name])))
                label = seg.get("label", "")
                lines.append(
                    f'<rect x="{sx:.2f}" y="{target_band_y}" width="{max(0.4, ex - sx):.2f}" height="12" '
                    f'fill="{seg.get("color", "#CBD5E0")}" stroke="#FFFFFF" stroke-width="0.15">'
                    f'<title>{esc(name)}:{seg["start"]}-{seg["end"]} {esc(label)}</title></rect>'
                )
        elif target_colors:
            color = target_colors.get(name, "#CBD5E0")
            lines.append(
                f'<rect x="{x:.2f}" y="{target_band_y}" width="{target_lengths[name] * x_scale:.2f}" height="12" fill="{color}" stroke="#FFFFFF" stroke-width="0.2"/>'
            )
        cx = x_of(name, target_lengths[name] / 2)
        lines.append(
            f'<text x="{cx:.2f}" y="{target_label_y}" font-family="Arial, sans-serif" font-size="10" '
            f'text-anchor="end" transform="rotate(-35 {cx:.2f} {target_label_y})">{esc(name)}</text>'
        )

    for name in query_names:
        y = y_of(name, query_lengths[name] / 2)
        query_segment_rows = axis_segments.get("query", {}).get(name, [])
        if query_segment_rows:
            for seg in query_segment_rows:
                sy1 = y_of(name, max(0, min(seg["start"], query_lengths[name])))
                sy2 = y_of(name, max(0, min(seg["end"], query_lengths[name])))
                y_top = min(sy1, sy2)
                y_height = max(0.4, abs(sy2 - sy1))
                label = seg.get("label", "")
                lines.append(
                    f'<rect x="{query_band_x}" y="{y_top:.2f}" width="13" height="{y_height:.2f}" '
                    f'fill="{seg.get("color", "#CBD5E0")}" stroke="#FFFFFF" stroke-width="0.15">'
                    f'<title>{esc(name)}:{seg["start"]}-{seg["end"]} {esc(label)}</title></rect>'
                )
        elif query_colors:
            color = query_colors.get(name, "#CBD5E0")
            lines.append(
                f'<rect x="{query_band_x}" y="{y - (query_lengths[name] * y_scale) / 2:.2f}" width="13" height="{query_lengths[name] * y_scale:.2f}" fill="{color}" stroke="#FFFFFF" stroke-width="0.2"/>'
            )
        lines.append(
            f'<text x="{query_label_x}" y="{y + 3:.2f}" font-family="Arial, sans-serif" font-size="9" '
            f'text-anchor="end">{esc(name)}</text>'
        )
        y0 = y_of(name, 0)
        lines.append(f'<line x1="{left}" y1="{y0:.2f}" x2="{left + plot_w}" y2="{y0:.2f}" stroke="#F3F4F6" stroke-width="0.6"/>')

    records_sorted = sorted(records, key=lambda r: r["block_len"])
    for rec in records_sorted:
        x1 = x_of(rec["target"], rec["target_start"])
        x2 = x_of(rec["target"], rec["target_end"])
        if rec["strand"] == "-":
            y1 = y_of(rec["query"], rec["query_end"])
            y2 = y_of(rec["query"], rec["query_start"])
        else:
            y1 = y_of(rec["query"], rec["query_start"])
            y2 = y_of(rec["query"], rec["query_end"])
        color = color_for_identity(rec["identity"], vmin, vmax)
        width_px = 0.6 if rec["block_len"] < 500_000 else 1.1
        label = (
            f'{rec["target"]}:{rec["target_start"]}-{rec["target_end"]} vs '
            f'{rec["query"]}:{rec["query_start"]}-{rec["query_end"]}; '
            f'identity={rec["identity"]:.3f}; block={rec["block_len"]}'
        )
        visible_len = math.hypot(x2 - x1, y2 - y1)
        if min_visible_px > 0 and visible_len < min_visible_px:
            cx = (x1 + x2) / 2
            cy = (y1 + y2) / 2
            radius = max(0.7, min_visible_px / 2)
            lines.append(
                f'<circle cx="{cx:.2f}" cy="{cy:.2f}" r="{radius:.2f}" '
                f'fill="{color}" stroke="{color}" stroke-width="0.35" opacity="0.92">'
                f'<title>{esc(label)}</title></circle>'
            )
        else:
            lines.append(
                f'<line x1="{x1:.2f}" y1="{y1:.2f}" x2="{x2:.2f}" y2="{y2:.2f}" '
                f'stroke="{color}" stroke-width="{width_px}" stroke-linecap="round" opacity="0.82">'
                f'<title>{esc(label)}</title></line>'
            )

    lines.extend(
        [
            f'<text x="{left + plot_w / 2}" y="{height - 28}" font-family="Arial, sans-serif" font-size="14" text-anchor="middle">Target</text>',
            f'<text x="34" y="{top + plot_h / 2}" font-family="Arial, sans-serif" font-size="14" text-anchor="middle" transform="rotate(-90 34 {top + plot_h / 2})">Query</text>',
        ]
    )

    legend_x = left + plot_w + 60
    legend_y = top + plot_h * 0.44
    lines.append(f'<text x="{legend_x}" y="{legend_y - 18}" font-family="Arial, sans-serif" font-size="11" font-weight="700">Alignment identity</text>')
    grad_h = 120
    steps = 30
    for i in range(steps):
        frac = i / (steps - 1)
        val = vmax - frac * (vmax - vmin)
        y = legend_y + frac * grad_h
        lines.append(
            f'<rect x="{legend_x}" y="{y:.2f}" width="18" height="{grad_h / steps + 1:.2f}" fill="{color_for_identity(val, vmin, vmax)}"/>'
        )
    lines.append(f'<rect x="{legend_x}" y="{legend_y}" width="18" height="{grad_h}" fill="none" stroke="#9CA3AF" stroke-width="0.6"/>')
    lines.append(f'<text x="{legend_x + 26}" y="{legend_y + 5}" font-family="Arial, sans-serif" font-size="10">{tick_label_identity(vmax)}</text>')
    lines.append(f'<text x="{legend_x + 26}" y="{legend_y + grad_h}" font-family="Arial, sans-serif" font-size="10">{tick_label_identity(vmin)}</text>')

    if title:
        lines.append(f'<text x="{left + plot_w / 2}" y="{height - 8}" font-family="Arial, sans-serif" font-size="11" fill="#4B5563" text-anchor="middle">{esc(title)}</text>')

    lines.append("</svg>")
    out_svg.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", required=True, type=Path)
    parser.add_argument("--out-svg", required=True, type=Path)
    parser.add_argument("--title", default="")
    parser.add_argument("--min-alignment-len", type=int, default=100_000)
    parser.add_argument("--min-query-aggregate", type=int, default=500_000)
    parser.add_argument("--min-target-aggregate", type=int, default=500_000)
    parser.add_argument("--min-visible-px", type=float, default=0.0, help="Render shorter projected alignments as visible points with this minimum diameter")
    parser.add_argument("--axis-colors", type=Path, help="Optional TSV with chr and color columns for both axes")
    parser.add_argument("--axis-segments", type=Path, help="Optional TSV with axis, chr, start, end, color, label columns")
    args = parser.parse_args()

    records, target_lengths, query_lengths = parse_paf(args.paf)
    axis_segments = read_axis_segments(args.axis_segments)
    target_colors = None
    query_colors = None
    if args.axis_colors and args.axis_colors.is_file():
        target_colors = {}
        query_colors = {}
        with args.axis_colors.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                chr_name = row.get("chr", "")
                color = row.get("color", "")
                axis = row.get("axis", "both")
                if not chr_name or not color:
                    continue
                if axis in {"both", "target"}:
                    target_colors[chr_name] = color
                if axis in {"both", "query"}:
                    query_colors[chr_name] = color
    draw_dotplot(
        records,
        target_lengths,
        query_lengths,
        args.out_svg,
        args.title,
        args.min_alignment_len,
        args.min_query_aggregate,
        args.min_target_aggregate,
        target_colors=target_colors,
        query_colors=query_colors,
        axis_segments=axis_segments,
        min_visible_px=args.min_visible_px,
    )


if __name__ == "__main__":
    main()
