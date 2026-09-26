#!/usr/bin/env python3
"""Redraw existing JCVI dotplots with complete haplotype chromosome labels.

The JCVI anchors and BED coordinates are used unchanged.  This only replaces
JCVI's shortened chromosome labels (for example ``C1``) with the full IDs
already present in the BED files, such as ``C01_A1`` or ``C01_A``.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from xml.sax.saxutils import escape


def args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair-root", type=Path, required=True)
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--include", nargs="*", help="Pair directory names; default discovers PZI/PSP/PPA pairs")
    p.add_argument("--dpi", type=int, default=300)
    return p.parse_args()


def read_bed(path: Path):
    genes = {}
    order = []
    by_chr = {}
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4:
                continue
            chrom, gene = f[0], f[3]
            if chrom not in by_chr:
                by_chr[chrom] = []
                order.append(chrom)
            if gene not in genes:
                genes[gene] = (chrom, len(by_chr[chrom]))
                by_chr[chrom].append(gene)
    return genes, order, by_chr


def parse_anchors(path: Path, qgenes, rgenes):
    q = []
    r = []
    with path.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("###"):
                continue
            f = line.split()
            if len(f) < 2 or f[0] not in qgenes or f[1] not in rgenes:
                continue
            q.append(qgenes[f[0]])
            r.append(rgenes[f[1]])
    return q, r


def positions(order, by_chr):
    # Keep chromosomes separated while preserving their gene-order scale.
    gap = 0.018
    n = sum(max(1, len(by_chr[c])) for c in order)
    total = n + gap * max(0, len(order) - 1) * n
    scale = 1.0 / total
    start = 0.0
    out = {}
    centers = []
    boundaries = []
    for chrom in order:
        size = max(1, len(by_chr[chrom])) * scale
        out[chrom] = (start, size / max(1, len(by_chr[chrom])))
        centers.append(start + size / 2)
        boundaries.append(start)
        start += size + gap * n * scale
    boundaries.append(start - gap * n * scale)
    return out, centers, boundaries


def draw(pair_dir: Path, outdir: Path, dpi: int):
    bed_files = sorted(pair_dir.glob("*.bed"))
    anchor_files = [p for p in pair_dir.glob("*.anchors") if not p.name.endswith(".new")]
    if len(bed_files) < 2 or not anchor_files:
        return None
    # Names are query/ref encoded by the anchor filename; BEDs are symlinked
    # into the same pair directory by the existing JCVI workflow.
    anchor = anchor_files[0]
    stem = anchor.name.removesuffix(".anchors")
    m = re.match(r"(.+)\.(.+)$", stem)
    if not m:
        return None
    qlabel, rlabel = m.groups()
    qbed = pair_dir / f"{qlabel}.bed"
    rbed = pair_dir / f"{rlabel}.bed"
    if not qbed.exists() or not rbed.exists():
        if len(bed_files) != 2:
            return None
        qbed, rbed = bed_files
    qgenes, qorder, qby = read_bed(qbed)
    rgenes, rorder, rby = read_bed(rbed)
    qpoints, rpoints = parse_anchors(anchor, qgenes, rgenes)
    if not qpoints:
        return None
    qpos, qcenters, qbounds = positions(qorder, qby)
    rpos, rcenters, rbounds = positions(rorder, rby)
    x = [qpos[ch][0] + (rank + 0.5) * qpos[ch][1] for ch, rank in qpoints]
    y = [rpos[ch][0] + (rank + 0.5) * rpos[ch][1] for ch, rank in rpoints]

    # Direct SVG emission keeps the local workflow dependency-free. Coordinates
    # match the JCVI dotplot convention: query on x, reference on y, y inverted.
    width, height = 1500, 1500
    left, top, plot_w, plot_h = 190, 150, 1120, 1120
    sx = lambda value: left + value / qbounds[-1] * plot_w
    sy = lambda value: top + value / rbounds[-1] * plot_h
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<rect x="{left}" y="{top}" width="{plot_w}" height="{plot_h}" fill="#E9E9F2"/>',
        '<g stroke="#B8DDBF" stroke-width="1" opacity="0.9">',
    ]
    for b in qbounds[1:-1]:
        xx = sx(b)
        parts.append(f'<line x1="{xx:.2f}" y1="{top}" x2="{xx:.2f}" y2="{top + plot_h}"/>')
    for b in rbounds[1:-1]:
        yy = sy(b)
        parts.append(f'<line x1="{left}" y1="{yy:.2f}" x2="{left + plot_w}" y2="{yy:.2f}"/>')
    parts.append('</g>')
    parts.append('<g fill="#111111">')
    for xx, yy in zip(x, y):
        parts.append(f'<circle cx="{sx(xx):.2f}" cy="{sy(yy):.2f}" r="1.35"/>')
    parts.append('</g>')
    parts.append('<g font-family="Arial, Helvetica, sans-serif" fill="#666666" font-size="15">')
    for center, label in zip(qcenters, qorder):
        xx = sx(center)
        parts.append(f'<text x="{xx:.2f}" y="{top - 18}" text-anchor="start" transform="rotate(45 {xx:.2f} {top - 18})">{escape(label)}</text>')
    for center, label in zip(rcenters, rorder):
        yy = sy(center)
        parts.append(f'<text x="{left - 14}" y="{yy + 5:.2f}" text-anchor="end">{escape(label)}</text>')
    parts.append('</g>')
    parts.append(f'<text x="{left + plot_w / 2:.1f}" y="{height - 45}" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="22" fill="#222">{escape(qlabel)}</text>')
    parts.append(f'<text x="48" y="{top + plot_h / 2:.1f}" text-anchor="middle" transform="rotate(-90 48 {top + plot_h / 2:.1f})" font-family="Arial, Helvetica, sans-serif" font-size="22" fill="#222">{escape(rlabel)}</text>')
    parts.append(f'<text x="{width / 2:.1f}" y="62" text-anchor="middle" font-family="Arial, Helvetica, sans-serif" font-size="24" fill="#222">JCVI synteny: {escape(qlabel)} vs {escape(rlabel)} ({len(x):,} gene pairs)</text>')
    parts.append('</svg>')
    outdir.mkdir(parents=True, exist_ok=True)
    svg = outdir / f"{stem}.haplotype_labels.dotplot.svg"
    svg.write_text("\n".join(parts) + "\n", encoding="utf-8")
    return stem, len(x), svg


def main() -> int:
    a = args()
    if a.include:
        dirs = [a.pair_root / name for name in a.include]
    else:
        dirs = []
        for d in sorted(a.pair_root.iterdir()):
            if not d.is_dir():
                continue
            pair_names = d.name.split("_vs_")
            if len(pair_names) != 2:
                continue
            # Limit the default set to comparisons where both tracks belong
            # to PZI, PSP or PPA; reference taxa (PCAM/PPE/GTR) are excluded.
            if not all(x.startswith(("pzi_", "psp_", "ppa_")) for x in pair_names):
                continue
            if list(d.glob("*.anchors")):
                dirs.append(d)
    rows = []
    for d in dirs:
        result = draw(d, a.outdir, a.dpi)
        if result:
            stem, n, svg = result
            rows.append((d.name, stem, n, svg))
            print(f"[OK] {d.name}: {n:,} pairs")
        else:
            print(f"[SKIP] no usable JCVI pair in {d}")
    a.outdir.mkdir(parents=True, exist_ok=True)
    with (a.outdir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(("pair_directory", "anchor_stem", "gene_pairs", "svg"))
        w.writerows(rows)
    print(f"[OK] Wrote {len(rows)} labelled dotplots to {a.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
