#!/usr/bin/env python3
"""Draw JCVI gene-level synteny figures from downstream_09 outputs."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Pair:
    comparison: str
    comparison_type: str
    block_id: str
    query_label: str
    ref_label: str
    query_gene: str
    ref_gene: str
    query_seqid: str
    query_start: int
    query_end: int
    query_strand: str
    ref_seqid: str
    ref_start: int
    ref_end: int
    ref_strand: str


@dataclass(frozen=True)
class Block:
    comparison: str
    comparison_type: str
    block_id: str
    query_label: str
    ref_label: str
    n_pairs: int
    query_seqids: str
    query_start_min: int | None
    query_end_max: int | None
    ref_seqids: str
    ref_start_min: int | None
    ref_end_max: int | None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jcvi-dir", required=True, help="analysis/karyotype/jcvi directory")
    parser.add_argument("--outdir", required=True, help="Output figure directory")
    parser.add_argument("--min-block-pairs", type=int, default=10, help="Minimum gene pairs for ribbon blocks")
    parser.add_argument("--max-dotplot-points", type=int, default=120000, help="Downsample very dense dotplots")
    parser.add_argument("--query-order", nargs="*", default=["pzi_A2B2", "pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2"])
    parser.add_argument("--ref-order", nargs="*", default=["ppe", "pdu", "pmu_refseq", "pcam", "pspe", "rch", "fve"])
    return parser.parse_args()


def natural_key(text: str) -> list[object]:
    parts: list[object] = []
    buf = ""
    was_digit: bool | None = None
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


def esc(text: object) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def to_int(value: str | None) -> int | None:
    if value in {"", "NA", "None", None}:
        return None
    try:
        return int(float(value))
    except ValueError:
        return None


def read_pairs(path: Path) -> list[Pair]:
    pairs: list[Pair] = []
    if not path.is_file():
        return pairs
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            qstart = to_int(row.get("query_start0"))
            qend = to_int(row.get("query_end"))
            rstart = to_int(row.get("ref_start0"))
            rend = to_int(row.get("ref_end"))
            if qstart is None or qend is None or rstart is None or rend is None:
                continue
            pairs.append(
                Pair(
                    comparison=row["comparison"],
                    comparison_type=row.get("comparison_type", "external_reference"),
                    block_id=row["block_id"],
                    query_label=row["query_label"],
                    ref_label=row["ref_label"],
                    query_gene=row["query_gene"],
                    ref_gene=row["ref_gene"],
                    query_seqid=row["query_seqid"],
                    query_start=qstart,
                    query_end=qend,
                    query_strand=row.get("query_strand", "."),
                    ref_seqid=row["ref_seqid"],
                    ref_start=rstart,
                    ref_end=rend,
                    ref_strand=row.get("ref_strand", "."),
                )
            )
    return pairs


def read_blocks(path: Path) -> list[Block]:
    blocks: list[Block] = []
    if not path.is_file():
        return blocks
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            blocks.append(
                Block(
                    comparison=row["comparison"],
                    comparison_type=row.get("comparison_type", "external_reference"),
                    block_id=row["block_id"],
                    query_label=row["query_label"],
                    ref_label=row["ref_label"],
                    n_pairs=int(row.get("n_pairs", "0") or 0),
                    query_seqids=row.get("query_seqids", "NA"),
                    query_start_min=to_int(row.get("query_start0_min")),
                    query_end_max=to_int(row.get("query_end_max")),
                    ref_seqids=row.get("ref_seqids", "NA"),
                    ref_start_min=to_int(row.get("ref_start0_min")),
                    ref_end_max=to_int(row.get("ref_end_max")),
                )
            )
    return blocks


def read_status(path: Path) -> dict[str, dict[str, str]]:
    status: dict[str, dict[str, str]] = {}
    if not path.is_file():
        return status
    with path.open(newline="", encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            status[row["comparison"]] = row
    return status


def read_bed_lengths(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = defaultdict(int)
    if not path.is_file():
        return lengths
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            end = to_int(fields[2])
            if end is not None:
                lengths[fields[0]] = max(lengths[fields[0]], end)
    return dict(lengths)


def chromosome_offsets(lengths: dict[str, int]) -> tuple[dict[str, int], int, list[str]]:
    order = sorted(lengths, key=natural_key)
    offsets: dict[str, int] = {}
    cursor = 0
    gap = max(100000, int(sum(lengths.values()) * 0.002)) if lengths else 1
    for seqid in order:
        offsets[seqid] = cursor
        cursor += lengths[seqid] + gap
    return offsets, max(cursor - gap, 1), order


def block_color(block_id: str) -> str:
    palette = ["#1F77B4", "#D62728", "#2CA02C", "#9467BD", "#FF7F0E", "#17BECF", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"]
    return palette[sum(ord(ch) for ch in block_id) % len(palette)]


def svg_header(width: int, height: int) -> list[str]:
    return [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect x="0" y="0" width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif}.grid{stroke:#D6DEE8;stroke-width:.7}.tick{font-size:9px;fill:#333}.title{font-size:16px;font-weight:700;fill:#111}.subtitle{font-size:11px;fill:#555}</style>',
    ]


def write_svg(path: Path, lines: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def dotplot(comparison: str, pairs: list[Pair], q_lengths: dict[str, int], r_lengths: dict[str, int], out: Path, max_points: int) -> None:
    width, height = 1300, 1100
    left, right, top, bottom = 120, 40, 80, 140
    plot_w = width - left - right
    plot_h = height - top - bottom
    q_offsets, q_total, q_order = chromosome_offsets(q_lengths)
    r_offsets, r_total, r_order = chromosome_offsets(r_lengths)
    step = max(1, math.ceil(len(pairs) / max_points)) if max_points > 0 else 1
    sampled = pairs[::step]
    lines = svg_header(width, height)
    lines.append(f'<text class="title" x="{left}" y="36">{esc(comparison)} JCVI syntenic gene-pair dotplot</text>')
    lines.append(f'<text class="subtitle" x="{left}" y="56">x: reference genome, y: PZI query genome; points are screened JCVI syntenic gene pairs</text>')
    lines.append(f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#FBFCFE" stroke="#222"/>')
    for seqid in r_order:
        x = left + plot_w * r_offsets[seqid] / r_total
        lines.append(f'<line class="grid" x1="{x:.2f}" y1="{top}" x2="{x:.2f}" y2="{top + plot_h}"/>')
        mid = left + plot_w * (r_offsets[seqid] + r_lengths[seqid] / 2) / r_total
        lines.append(f'<text class="tick" x="{mid:.2f}" y="{top + plot_h + 18}" text-anchor="middle" transform="rotate(45 {mid:.2f},{top + plot_h + 18})">{esc(seqid)}</text>')
    for seqid in q_order:
        y = top + plot_h - plot_h * q_offsets[seqid] / q_total
        lines.append(f'<line class="grid" x1="{left}" y1="{y:.2f}" x2="{left + plot_w}" y2="{y:.2f}"/>')
        mid = top + plot_h - plot_h * (q_offsets[seqid] + q_lengths[seqid] / 2) / q_total
        lines.append(f'<text class="tick" x="{left - 8}" y="{mid + 3:.2f}" text-anchor="end">{esc(seqid)}</text>')
    for p in sampled:
        x = left + plot_w * (r_offsets[p.ref_seqid] + (p.ref_start + p.ref_end) / 2) / r_total
        y = top + plot_h - plot_h * (q_offsets[p.query_seqid] + (p.query_start + p.query_end) / 2) / q_total
        lines.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="1.15" fill="{block_color(p.block_id)}" opacity="0.62"/>')
    lines.append(f'<text class="subtitle" x="{left}" y="{height - 18}">pairs={len(pairs)}; plotted={len(sampled)}; downsample_step={step}</text>')
    lines.append("</svg>")
    write_svg(out, lines)


def ribbonplot(comparison: str, blocks: list[Block], q_lengths: dict[str, int], r_lengths: dict[str, int], min_block_pairs: int, out: Path) -> None:
    kept = [b for b in blocks if b.n_pairs >= min_block_pairs and b.query_start_min is not None and b.ref_start_min is not None]
    width, height = 1200, 1100
    left_x, right_x = 280, 920
    top, bottom = 90, 70
    panel_h = height - top - bottom
    q_offsets, q_total, q_order = chromosome_offsets(q_lengths)
    r_offsets, r_total, r_order = chromosome_offsets(r_lengths)

    def q_y(seqid: str, pos: int) -> float:
        return top + panel_h * (q_offsets[seqid] + pos) / q_total

    def r_y(seqid: str, pos: int) -> float:
        return top + panel_h * (r_offsets[seqid] + pos) / r_total

    lines = svg_header(width, height)
    lines.append(f'<text class="title" x="90" y="36">{esc(comparison)} JCVI macro-synteny ribbon plot</text>')
    lines.append(f'<text class="subtitle" x="90" y="56">blocks with at least {min_block_pairs} screened syntenic gene pairs</text>')
    for seqid in q_order:
        y1 = q_y(seqid, 0)
        y2 = q_y(seqid, q_lengths[seqid])
        lines.append(f'<line x1="{left_x}" y1="{y1:.2f}" x2="{left_x}" y2="{y2:.2f}" stroke="#1D2330" stroke-width="6"/>')
        lines.append(f'<text class="tick" x="{left_x - 12}" y="{(y1 + y2) / 2 + 3:.2f}" text-anchor="end">{esc(seqid)}</text>')
    for seqid in r_order:
        y1 = r_y(seqid, 0)
        y2 = r_y(seqid, r_lengths[seqid])
        lines.append(f'<line x1="{right_x}" y1="{y1:.2f}" x2="{right_x}" y2="{y2:.2f}" stroke="#1D2330" stroke-width="6"/>')
        lines.append(f'<text class="tick" x="{right_x + 12}" y="{(y1 + y2) / 2 + 3:.2f}" text-anchor="start">{esc(seqid)}</text>')
    for b in kept:
        qseq = b.query_seqids.split(",")[0]
        rseq = b.ref_seqids.split(",")[0]
        if qseq not in q_lengths or rseq not in r_lengths or b.query_end_max is None or b.ref_end_max is None:
            continue
        yq1 = q_y(qseq, b.query_start_min or 0)
        yq2 = q_y(qseq, b.query_end_max)
        yr1 = r_y(rseq, b.ref_start_min or 0)
        yr2 = r_y(rseq, b.ref_end_max)
        color = block_color(b.block_id)
        opacity = min(0.78, 0.18 + math.log10(b.n_pairs + 1) / 3.2)
        path = f'M {left_x + 5},{yq1:.2f} C {left_x + 260},{yq1:.2f} {right_x - 260},{yr1:.2f} {right_x - 5},{yr1:.2f} L {right_x - 5},{yr2:.2f} C {right_x - 260},{yr2:.2f} {left_x + 260},{yq2:.2f} {left_x + 5},{yq2:.2f} Z'
        lines.append(f'<path d="{path}" fill="{color}" opacity="{opacity:.2f}" stroke="{color}" stroke-width="0.5"><title>{esc(b.block_id)} pairs={b.n_pairs}</title></path>')
    lines.append(f'<text class="subtitle" x="90" y="{height - 22}">blocks_total={len(blocks)}; blocks_plotted={len(kept)}</text>')
    lines.append("</svg>")
    write_svg(out, lines)


def heat_color(value: float, vmax: float) -> str:
    if value <= 0 or vmax <= 0:
        return "#F7FAFC"
    frac = max(0.0, min(1.0, math.log10(value + 1) / math.log10(vmax + 1)))
    if frac < 0.5:
        t = frac / 0.5
        r, g, b = 255, int(255 - 45 * t), int(255 - 185 * t)
    else:
        t = (frac - 0.5) / 0.5
        r, g, b = int(255 - 95 * t), int(210 - 130 * t), int(70 - 55 * t)
    return f"#{r:02X}{g:02X}{b:02X}"


def heatmap(title: str, matrix: dict[tuple[str, str], float], rows: list[str], cols: list[str], out: Path) -> None:
    cell_w, cell_h = 120, 42
    left, top = 150, 90
    width = left + cell_w * len(cols) + 70
    height = top + cell_h * len(rows) + 90
    vmax = max(matrix.values()) if matrix else 0
    lines = svg_header(width, height)
    lines.append(f'<text class="title" x="{left}" y="36">{esc(title)}</text>')
    for i, col in enumerate(cols):
        x = left + i * cell_w + cell_w / 2
        lines.append(f'<text class="tick" x="{x:.2f}" y="{top - 14}" text-anchor="middle">{esc(col)}</text>')
    for r, row in enumerate(rows):
        y = top + r * cell_h
        lines.append(f'<text class="tick" x="{left - 10}" y="{y + cell_h / 2 + 4:.2f}" text-anchor="end">{esc(row)}</text>')
        for c, col in enumerate(cols):
            value = matrix.get((row, col), 0)
            x = left + c * cell_w
            lines.append(f'<rect x="{x}" y="{y}" width="{cell_w}" height="{cell_h}" fill="{heat_color(value, vmax)}" stroke="#D1D5DB"/>')
            label = str(int(value)) if abs(value - int(value)) < 1e-9 else f"{value:.1f}"
            lines.append(f'<text class="tick" x="{x + cell_w / 2:.2f}" y="{y + cell_h / 2 + 4:.2f}" text-anchor="middle">{esc(label)}</text>')
    lines.append("</svg>")
    write_svg(out, lines)


def main() -> int:
    args = parse_args()
    jcvi_dir = Path(args.jcvi_dir)
    outdir = Path(args.outdir)
    summary_dir = jcvi_dir / "summary"
    inputs_dir = jcvi_dir / "inputs"
    pairs = read_pairs(summary_dir / "jcvi_syntenic_pairs.tsv")
    blocks = read_blocks(summary_dir / "jcvi_block_summary.tsv")
    status = read_status(summary_dir / "jcvi_pair_status.tsv")
    if not status:
        raise SystemExit(f"No JCVI pair status table found under {summary_dir}")

    pair_by_comp: dict[str, list[Pair]] = defaultdict(list)
    for p in pairs:
        pair_by_comp[p.comparison].append(p)
    block_by_comp: dict[str, list[Block]] = defaultdict(list)
    for b in blocks:
        block_by_comp[b.comparison].append(b)

    for comparison, row in sorted(status.items()):
        query = row["query_label"]
        ref = row["ref_label"]
        comparison_type = row.get("comparison_type", "external_reference")
        q_lengths = read_bed_lengths(inputs_dir / query / f"{query}.bed")
        r_lengths = read_bed_lengths(inputs_dir / ref / f"{ref}.bed")
        comp_pairs = pair_by_comp.get(comparison, [])
        comp_blocks = block_by_comp.get(comparison, [])
        plot_group = "haplotype_internal" if comparison_type == "internal_haplotype" else "reference_external"
        if q_lengths and r_lengths and comp_pairs:
            dotplot(comparison, comp_pairs, q_lengths, r_lengths, outdir / "dotplots" / plot_group / f"{comparison}.jcvi_dotplot.svg", args.max_dotplot_points)
        if q_lengths and r_lengths and comp_blocks:
            ribbonplot(comparison, comp_blocks, q_lengths, r_lengths, args.min_block_pairs, outdir / "ribbons" / plot_group / f"{comparison}.jcvi_ribbon.svg")

    pair_count: dict[tuple[str, str], float] = {}
    block_count: dict[tuple[str, str], float] = {}
    mean_pairs: dict[tuple[str, str], float] = {}
    for row in status.values():
        key = (row["query_label"], row["ref_label"])
        n_pairs = float(row.get("n_syntenic_pairs", "0") or 0)
        n_blocks = float(row.get("n_blocks", "0") or 0)
        pair_count[key] = n_pairs
        block_count[key] = n_blocks
        mean_pairs[key] = n_pairs / n_blocks if n_blocks > 0 else 0

    external_status = {k: v for k, v in status.items() if v.get("comparison_type", "external_reference") != "internal_haplotype"}
    internal_status = {k: v for k, v in status.items() if v.get("comparison_type", "external_reference") == "internal_haplotype"}

    rows = [x for x in args.query_order if any(row["query_label"] == x for row in external_status.values())]
    cols = [x for x in args.ref_order if any(row["ref_label"] == x for row in external_status.values())]
    external_keys = {(row["query_label"], row["ref_label"]) for row in external_status.values()}
    heatmap("JCVI screened syntenic gene-pair count", {k: v for k, v in pair_count.items() if k in external_keys}, rows, cols, outdir / "summary_heatmaps" / "reference_external" / "jcvi_syntenic_pair_count_heatmap.svg")
    heatmap("JCVI screened syntenic block count", {k: v for k, v in block_count.items() if k in external_keys}, rows, cols, outdir / "summary_heatmaps" / "reference_external" / "jcvi_block_count_heatmap.svg")
    heatmap("Mean syntenic gene pairs per block", {k: v for k, v in mean_pairs.items() if k in external_keys}, rows, cols, outdir / "summary_heatmaps" / "reference_external" / "jcvi_mean_pairs_per_block_heatmap.svg")

    if internal_status:
        hap_order = ["pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2"]
        internal_keys = {(row["query_label"], row["ref_label"]) for row in internal_status.values()}
        int_rows = [x for x in hap_order if any(k[0] == x for k in internal_keys)]
        int_cols = [x for x in hap_order if any(k[1] == x for k in internal_keys)]
        heatmap("JCVI internal haplotype syntenic gene-pair count", {k: v for k, v in pair_count.items() if k in internal_keys}, int_rows, int_cols, outdir / "summary_heatmaps" / "haplotype_internal" / "jcvi_haplotype_pair_count_heatmap.svg")
        heatmap("JCVI internal haplotype syntenic block count", {k: v for k, v in block_count.items() if k in internal_keys}, int_rows, int_cols, outdir / "summary_heatmaps" / "haplotype_internal" / "jcvi_haplotype_block_count_heatmap.svg")

    manifest = outdir / "figure_manifest.tsv"
    manifest.parent.mkdir(parents=True, exist_ok=True)
    with manifest.open("w", encoding="utf-8") as handle:
        print("figure_type\tpath\tn_records", file=handle)
        for path in sorted((outdir / "dotplots").glob("**/*.svg")):
            comparison = path.name.replace(".jcvi_dotplot.svg", "")
            print(f"dotplot\t{path}\t{len(pair_by_comp.get(comparison, []))}", file=handle)
        for path in sorted((outdir / "ribbons").glob("**/*.svg")):
            comparison = path.name.replace(".jcvi_ribbon.svg", "")
            n = sum(1 for b in block_by_comp.get(comparison, []) if b.n_pairs >= args.min_block_pairs)
            print(f"ribbon\t{path}\t{n}", file=handle)
        for path in sorted((outdir / "summary_heatmaps").glob("**/*.svg")):
            print(f"summary_heatmap\t{path}\tNA", file=handle)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
