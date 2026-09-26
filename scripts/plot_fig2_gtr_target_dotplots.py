#!/usr/bin/env python3
"""Draw GTR-reference JCVI dotplots for all Figure 2 target genomes."""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path


MAIN = [
    "pzi_hapA1_vs_gtr", "pzi_hapB1_vs_gtr",
    "psp_hapA1_vs_gtr", "psp_hapB1_vs_gtr",
    "ppa_hapA_vs_gtr", "ppa_hapB_vs_gtr",
    "pcam_vs_gtr", "ppe_vs_gtr",
]
SENSITIVITY = [
    "pzi_hapA2_vs_gtr", "pzi_hapB2_vs_gtr",
    "psp_hapA2_vs_gtr", "psp_hapB2_vs_gtr",
]
DISPLAY = {
    "pzi_hapA1": "PZI-A1", "pzi_hapA2": "PZI-A2",
    "pzi_hapB1": "PZI-B1", "pzi_hapB2": "PZI-B2",
    "psp_hapA1": "PSP-A1", "psp_hapA2": "PSP-A2",
    "psp_hapB1": "PSP-B1", "psp_hapB2": "PSP-B2",
    "ppa_hapA": "PPA-A", "ppa_hapB": "PPA-B",
    "pcam": "PCAM", "ppe": "PPE", "gtr": "GTR",
}
PALETTE = ["#3B6FB6", "#D4554B", "#3F8C75", "#9165A8", "#D08A3C", "#5B8FA8", "#B36A8C", "#7C7C7C", "#8AA34B"]


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True)
    parser.add_argument("--inputs", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--max-grid-points", type=int, default=7000)
    return parser.parse_args()


def esc(value):
    return str(value).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def read_seqids(path):
    return [value for value in path.read_text(encoding="utf-8").strip().split(",") if value]


def primary_seqids(label, seqids):
    if label == "gtr":
        return [seqid for seqid in seqids if re.fullmatch(r"Chr0[1-9]", seqid)]
    if label == "ppe":
        return [seqid for seqid in seqids if re.fullmatch(r"NC_0340(09|10|11|12|13|14|15|16)\.1", seqid)]
    if label.startswith(("pzi_hap", "psp_hap", "ppa_hap")):
        return seqids[:8]
    if label == "pcam":
        return seqids[:8]
    return seqids


def read_lengths(path, allowed):
    lengths = defaultdict(int)
    allowed_set = set(allowed)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3 or fields[0] not in allowed_set:
                continue
            try:
                lengths[fields[0]] = max(lengths[fields[0]], int(fields[2]))
            except ValueError:
                continue
    return {seqid: lengths[seqid] for seqid in allowed if lengths[seqid] > 0}


def axis_layout(lengths, order):
    total_bp = sum(lengths.values())
    gap = max(1, int(total_bp * 0.006))
    offsets = {}
    cursor = 0
    retained_order = [seqid for seqid in order if seqid in lengths]
    for seqid in retained_order:
        offsets[seqid] = cursor
        cursor += lengths[seqid] + gap
    return offsets, max(1, cursor - gap), retained_order


def read_pairs(path):
    grouped = defaultdict(list)
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["comparison"]].append(row)
    return grouped


def svg_start(width, height):
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (width, height, width, height),
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Helvetica,sans-serif;letter-spacing:0}.title{font-size:17px;font-weight:700}.label{font-size:11px;fill:#222}.tick{font-size:8px;fill:#444}.meta{font-size:9px;fill:#555}.frame{fill:#FCFCFC;stroke:#333;stroke-width:.8}.grid{stroke:#CBD1D8;stroke-width:.55}</style>',
    ]


def panel_svg(lines, x0, y0, width, height, comparison, records, inputs, max_points, panel_label=None):
    query = comparison.rsplit("_vs_gtr", 1)[0]
    ref = "gtr"
    q_seqids = primary_seqids(query, read_seqids(inputs / query / (query + ".seqids")))
    r_seqids = primary_seqids(ref, read_seqids(inputs / ref / (ref + ".seqids")))
    q_lengths = read_lengths(inputs / query / (query + ".bed"), q_seqids)
    r_lengths = read_lengths(inputs / ref / (ref + ".bed"), r_seqids)
    q_offsets, q_total, q_order = axis_layout(q_lengths, q_seqids)
    r_offsets, r_total, r_order = axis_layout(r_lengths, r_seqids)
    primary = [row for row in records if row["query_seqid"] in q_offsets and row["ref_seqid"] in r_offsets]
    step = max(1, int(math.ceil(len(primary) / float(max_points)))) if max_points else 1
    sampled = primary[::step]

    left, right, top, bottom = 64, 18, 48, 65
    pw, ph = width - left - right, height - top - bottom
    if panel_label:
        lines.append('<text class="title" x="%.1f" y="%.1f">%s</text>' % (x0 + 8, y0 + 20, esc(panel_label)))
    lines.append('<text class="title" x="%.1f" y="%.1f">%s vs GTR</text>' % (x0 + left, y0 + 24, esc(DISPLAY.get(query, query))))
    lines.append('<rect class="frame" x="%.1f" y="%.1f" width="%.1f" height="%.1f"/>' % (x0 + left, y0 + top, pw, ph))
    for seqid in r_order:
        x = x0 + left + pw * r_offsets[seqid] / r_total
        lines.append('<line class="grid" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f"/>' % (x, y0 + top, x, y0 + top + ph))
        mid = x0 + left + pw * (r_offsets[seqid] + r_lengths[seqid] / 2.0) / r_total
        lines.append('<text class="tick" x="%.2f" y="%.2f" text-anchor="middle">%s</text>' % (mid, y0 + top + ph + 14, esc(seqid.replace("Chr", ""))))
    color_by_seqid = {seqid: PALETTE[index % len(PALETTE)] for index, seqid in enumerate(q_order)}
    for seqid in q_order:
        y = y0 + top + ph - ph * q_offsets[seqid] / q_total
        lines.append('<line class="grid" x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f"/>' % (x0 + left, y, x0 + left + pw, y))
        mid = y0 + top + ph - ph * (q_offsets[seqid] + q_lengths[seqid] / 2.0) / q_total
        short = seqid.replace("C0", "").replace("_A1", "A1").replace("_A2", "A2").replace("_B1", "B1").replace("_B2", "B2").replace("_A", "A").replace("_B", "B")
        lines.append('<text class="tick" x="%.2f" y="%.2f" text-anchor="end">%s</text>' % (x0 + left - 5, mid + 2.5, esc(short)))
    for row in sampled:
        qpos = (float(row["query_start0"]) + float(row["query_end"])) / 2.0
        rpos = (float(row["ref_start0"]) + float(row["ref_end"])) / 2.0
        x = x0 + left + pw * (r_offsets[row["ref_seqid"]] + rpos) / r_total
        y = y0 + top + ph - ph * (q_offsets[row["query_seqid"]] + qpos) / q_total
        lines.append('<circle cx="%.2f" cy="%.2f" r="0.75" fill="%s" opacity="0.52"/>' % (x, y, color_by_seqid[row["query_seqid"]]))
    lines.append('<text class="label" x="%.2f" y="%.2f" text-anchor="middle">GTR chromosome</text>' % (x0 + left + pw / 2.0, y0 + height - 20))
    lines.append('<text class="meta" x="%.2f" y="%.2f">primary pairs=%d/%d; plotted=%d</text>' % (x0 + left, y0 + height - 5, len(primary), len(records), len(sampled)))
    return {"comparison": comparison, "query_label": query, "total_pairs": len(records), "primary_pairs": len(primary), "plotted_pairs": len(sampled)}


def write_individual(path, comparison, records, inputs):
    width, height = 1100, 900
    lines = svg_start(width, height)
    result = panel_svg(lines, 0, 0, width, height, comparison, records, inputs, 60000)
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def write_grid(path, comparisons, grouped, inputs, columns, max_points):
    panel_w, panel_h = 520, 430
    rows = int(math.ceil(len(comparisons) / float(columns)))
    width, height = panel_w * columns, panel_h * rows
    lines = svg_start(width, height)
    results = []
    for index, comparison in enumerate(comparisons):
        x0 = (index % columns) * panel_w
        y0 = (index // columns) * panel_h
        results.append(panel_svg(lines, x0, y0, panel_w, panel_h, comparison, grouped[comparison], inputs, max_points))
    lines.append("</svg>")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return results


def main():
    args = parse_args()
    pairs_path = Path(args.pairs)
    inputs = Path(args.inputs)
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: " + str(outdir))
    grouped = read_pairs(pairs_path)
    missing = [comparison for comparison in MAIN + SENSITIVITY if comparison not in grouped]
    if missing:
        raise SystemExit("ERROR: missing comparisons: " + ",".join(missing))
    manifest_rows = []
    for group_name, comparisons in (("main", MAIN), ("sensitivity", SENSITIVITY)):
        for comparison in comparisons:
            result = write_individual(outdir / "individual" / group_name / (comparison + ".dotplot.svg"), comparison, grouped[comparison], inputs)
            result["analysis_role"] = group_name
            manifest_rows.append(result)
    write_grid(outdir / "gtr_target_main_grid.svg", MAIN, grouped, inputs, 4, args.max_grid_points)
    write_grid(outdir / "gtr_target_sensitivity_grid.svg", SENSITIVITY, grouped, inputs, 2, args.max_grid_points)
    with (outdir / "dotplot_manifest.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = ["comparison", "query_label", "analysis_role", "total_pairs", "primary_pairs", "plotted_pairs"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    (outdir / "README.txt").write_text(
        "All plots use screened JCVI syntenic gene pairs. X=GTR Chr01-Chr09; Y=target primary nuclear chromosomes.\n"
        "PPE chloroplast NC_014697.1 and unplaced scaffolds are excluded from plotting only.\n"
        "Individual plots retain all primary-chromosome pairs; summary grids may be evenly subsampled for file size.\n"
        "Colors identify target chromosomes and do not encode subgenome affinity or evolutionary distance.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("GTR-target dotplots completed.\n", encoding="ascii")
    print("[OK] Wrote 12 dotplots and two summary grids: " + str(outdir))


if __name__ == "__main__":
    main()
