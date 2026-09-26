#!/usr/bin/env python3
"""Plot a chromosome-pair heatmap from an asm10 PAF alignment."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


def fasta_lengths(path: Path) -> tuple[list[str], dict[str, int]]:
    order = []
    lengths = {}
    name = None
    length = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = raw[1:].split(None, 1)[0]
                order.append(name)
                length = 0
            else:
                length += len(raw.strip())
    if name is not None:
        lengths[name] = length
    if not order or any(lengths.get(x, 0) <= 0 for x in order):
        raise SystemExit(f"Invalid FASTA: {path}")
    return order, lengths


def union_length(intervals: list[tuple[int, int]]) -> int:
    if not intervals:
        return 0
    merged = []
    for start, end in sorted(intervals):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return sum(end - start for start, end in merged)


def parse_paf(path: Path, query: set[str], target: set[str], min_mapq: int, min_aln: int):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 12:
                continue
            qname = fields[0]
            tname = fields[5]
            if qname not in query or tname not in target:
                continue
            qstart, qend = int(fields[2]), int(fields[3])
            mapq = int(fields[11])
            if qend - qstart < min_aln or mapq < min_mapq:
                continue
            yield qname, tname, qstart, qend


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.linewidth": 0.7,
    })


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--query-fasta", type=Path, required=True)
    parser.add_argument("--target-fasta", type=Path, required=True)
    parser.add_argument("--out-stem", type=Path, required=True)
    parser.add_argument("--title", default="")
    parser.add_argument("--min-mapq", type=int, default=0)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    args = parser.parse_args()

    q_order, q_lengths = fasta_lengths(args.query_fasta)
    t_order, t_lengths = fasta_lengths(args.target_fasta)
    intervals: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    n = 0
    for qname, tname, qstart, qend in parse_paf(
        args.paf, set(q_order), set(t_order), args.min_mapq, args.min_aln_len
    ):
        intervals[(qname, tname)].append((qstart, qend))
        n += 1
    if not n:
        raise SystemExit(f"No PAF records passed filters: {args.paf}")

    matrix = []
    rows = []
    for tname in t_order:
        row = []
        for qname in q_order:
            value = union_length(intervals.get((qname, tname), []))
            row.append(value)
            rows.append({
                "target_chr": tname,
                "query_chr": qname,
                "aligned_union_bp": value,
                "query_length_bp": q_lengths[qname],
                "target_length_bp": t_lengths[tname],
                "query_coverage": f"{value / q_lengths[qname]:.6f}",
                "target_coverage": f"{value / t_lengths[tname]:.6f}",
            })
        matrix.append(row)

    args.out_stem.parent.mkdir(parents=True, exist_ok=True)
    with Path(f"{args.out_stem}.matrix.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    configure()
    values = [[math.log10(1 + x) for x in row] for row in matrix]
    vmax = max(1.0, max(max(row) for row in values))
    fig, ax = plt.subplots(figsize=(10.2, 8.6))
    image = ax.imshow(values, aspect="auto", interpolation="nearest", cmap="viridis", norm=Normalize(0, vmax))
    ax.set_xticks(range(len(q_order)), q_order, rotation=90, fontsize=7)
    ax.set_yticks(range(len(t_order)), t_order, fontsize=7)
    ax.set_xlabel("Query chromosome")
    ax.set_ylabel("Target chromosome")
    if args.title:
        ax.set_title(args.title)
    if len(q_order) >= 16:
        ax.axvline(7.5, color="white", lw=0.8)
    if len(t_order) >= 16:
        ax.axhline(7.5, color="white", lw=0.8)
    cbar = fig.colorbar(image, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("log10(1 + union-aligned bp)")
    fig.savefig(f"{args.out_stem}.svg", bbox_inches="tight")
    fig.savefig(f"{args.out_stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{args.out_stem}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"[OK] {args.out_stem}: {n} PAF records, {len(t_order)} x {len(q_order)} chromosome matrix")


if __name__ == "__main__":
    main()
