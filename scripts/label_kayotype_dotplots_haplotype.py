#!/usr/bin/env python3
"""Replace chromosome labels on existing JCVI dotplot PNGs.

The raster plot, points, ordering and filtering are preserved.  Only the
query/reference chromosome-label margins are painted over and redrawn using
the PZI, PSP and PPA haplotype naming maps.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


PZI_COPY1 = {
    "chr01": "C06_A1", "chr02": "C02_A1", "chr03": "C03_A1",
    "chr04": "C04_A1", "chr05": "C05_A1", "chr06": "C01_A1",
    "chr07": "C04_B1", "chr08": "C07_A1", "chr09": "C08_A1",
    "chr10": "C07_B1", "chr11": "C02_B1", "chr12": "C03_B1",
    "chr13": "C06_B1", "chr14": "C08_B1", "chr15": "C01_B1",
    "chr16": "C05_B1",
}

# The generic PSP dotplots use the 16 chromosomes from
# ``psp_2n32.copy1.hap_provenance.chr16.fa``.  Despite the historical
# ``copy1`` filename, the sequence lengths match the even scaffold records in
# ``psp_A1A2B1B2.rename_map.tsv``; this is the A2+B2 set.  Keep this explicit
# mapping so that relabelling never changes the plotted coordinates.
PSP_GENERIC = {
    "chr01": "C03_A2", "chr02": "C02_A2", "chr03": "C06_A2",
    "chr04": "C01_A2", "chr05": "C06_B2", "chr06": "C08_A2",
    "chr07": "C08_B2", "chr08": "C07_B2", "chr09": "C07_A2",
    "chr10": "C04_B2", "chr11": "C01_B2", "chr12": "C03_B2",
    "chr13": "C04_A2", "chr14": "C05_A2", "chr15": "C05_B2",
    "chr16": "C02_B2",
}


def read_map(path: Path, key: str, value: str) -> dict[str, str]:
    out: dict[str, str] = {}
    with path.open(encoding="utf-8-sig", errors="replace", newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            if row.get(key) and row.get(value):
                out[row[key]] = row[value]
    return out


def read_chromosomes(bed: Path) -> list[str]:
    seen: list[str] = []
    with bed.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split()[0]
            if chrom not in seen:
                seen.append(chrom)
    return seen


def gene_counts(bed: Path) -> tuple[list[str], dict[str, int]]:
    order: list[str] = []
    counts: dict[str, int] = {}
    with bed.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            chrom = line.split()[0]
            if chrom not in counts:
                order.append(chrom)
                counts[chrom] = 0
            counts[chrom] += 1
    return order, counts


def font(size: int) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        r"C:\Windows\Fonts\arial.ttf",
        r"C:\Windows\Fonts\Arial.ttf",
        r"C:\Windows\Fonts\segoeui.ttf",
    ]
    for candidate in candidates:
        try:
            return ImageFont.truetype(candidate, size)
        except OSError:
            pass
    return ImageFont.load_default()


def source_bed(figdir: Path, stem: str, label: str) -> Path | None:
    direct = figdir / f"{label}.bed"
    if direct.is_file():
        return direct
    # Native JCVI inputs for haplotype plots are staged with the plot inputs.
    root = figdir.parents[2] / "analysis" / "karyotype" / "jcvi" / "figures" / "jcvi_native" / "reference_external"
    candidate = root / stem / f"{label}.bed"
    if candidate.is_file():
        return candidate
    # The local ks pair inputs contain the same BED coordinates.
    ks = figdir.parents[2] / "downstream_analysis" / "ks" / "jcvi" / "pairs"
    for pair in ks.glob(f"{label}_vs_* /{label}.bed"):
        if pair.is_file():
            return pair
    for pairdir in ks.glob(f"{label}_vs_*"):
        candidate = pairdir / f"{label}.bed"
        if candidate.is_file():
            return candidate
    return None


def hap_labels(label: str, chroms: list[str]) -> dict[str, str] | None:
    if label == "pzi_copy1":
        return {c: PZI_COPY1.get(c, c) for c in chroms}
    if label == "psp":
        return {c: PSP_GENERIC.get(c, c) for c in chroms}
    if label == "ppa":
        return None
    match = re.fullmatch(r"pzi_hap(A1|A2|B1|B2)", label)
    if match:
        suffix = match.group(1)
        return {c: f"C{i:02d}_{suffix}" for i, c in enumerate(chroms, 1)}
    match = re.fullmatch(r"psp_hap(A1|A2|B1|B2)", label)
    if match:
        suffix = match.group(1)
        return {c: f"C{i:02d}_{suffix}" for i, c in enumerate(chroms, 1)}
    return None


def ppa_labels(bed: Path, rename: dict[str, str]) -> dict[str, str]:
    # PPA BEDs also contain unplaced contigs.  They are intentionally left
    # unlabeled in the relabeled margin rather than displaying long accession
    # strings that overlap one another; only the 16 assigned chromosomes are
    # chromosome labels requested for the figure.
    return {c: rename[c] for c in read_chromosomes(bed) if c in rename}


def centers(order: list[str], counts: dict[str, int], start: float, end: float) -> dict[str, float]:
    total = sum(max(1, counts.get(c, 1)) for c in order)
    out: dict[str, float] = {}
    cursor = start
    span = end - start
    for chrom in order:
        width = span * max(1, counts.get(chrom, 1)) / total
        out[chrom] = cursor + width / 2
        cursor += width
    return out


def rotated_label(draw: ImageDraw.ImageDraw, xy: tuple[float, float], text: str, fnt, fill: str) -> None:
    bbox = draw.textbbox((0, 0), text, font=fnt)
    w, h = bbox[2] - bbox[0] + 12, bbox[3] - bbox[1] + 12
    layer = Image.new("RGBA", (w, h), (255, 255, 255, 0))
    ImageDraw.Draw(layer).text((6, 6), text, font=fnt, fill=fill)
    layer = layer.rotate(45, expand=True, resample=Image.Resampling.BICUBIC)
    draw._image.paste(layer, (int(xy[0] - layer.width / 2), int(xy[1] - layer.height / 2)), layer)


def relabel(path: Path, out: Path, query_labels: dict[str, str] | None,
            ref_labels: dict[str, str] | None, qbed: Path | None, rbed: Path | None) -> tuple[int, int]:
    image = Image.open(path).convert("RGB")
    width, height = image.size
    draw = ImageDraw.Draw(image)
    # These are the fixed JCVI plot margins for the 3600 x 3600 PNGs.
    sx0, sx1 = 360, min(width - 350, 3230)
    sy0, sy1 = 365, min(height - 350, 3230)
    top_changed = 0
    right_changed = 0
    if query_labels and qbed:
        qorder, qcounts = gene_counts(qbed)
        qcenters = centers(qorder, qcounts, sx0, sx1)
        # Paint only the chromosome-label strip, leaving title and plot intact.
        draw.rectangle((sx0 - 35, 170, sx1 + 35, sy0 - 2), fill="white")
        fnt = font(34 if len(qorder) <= 16 else 27)
        for chrom in qorder:
            text = query_labels.get(chrom, chrom)
            rotated_label(draw, (qcenters[chrom], 285), text, fnt, "#888888")
            top_changed += 1
        draw.line((sx0, sy0, sx1, sy0), fill="#B8DDBF", width=3)
    if ref_labels and rbed:
        rorder, rcounts = gene_counts(rbed)
        rcenters = centers(rorder, rcounts, sy0, sy1)
        draw.rectangle((sx1 + 2, sy0 - 5, width, sy1 + 8), fill="white")
        fnt = font(30 if len(rorder) <= 32 else 22)
        for chrom in rorder:
            if chrom not in ref_labels:
                continue
            text = ref_labels[chrom]
            bbox = draw.textbbox((0, 0), text, font=fnt)
            draw.text((sx1 + 38, rcenters[chrom] - (bbox[3] - bbox[1]) / 2), text, font=fnt, fill="#888888")
            right_changed += 1
    out.parent.mkdir(parents=True, exist_ok=True)
    image.save(out, format="PNG", optimize=True)
    return top_changed, right_changed


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figdir", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--ppa-map", type=Path, required=True)
    args = ap.parse_args()
    ppa_map = read_map(args.ppa_map, "ppa_source_seqid", "renamed_seqid")
    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in sorted(args.figdir.glob("*_vs_*.dotplot.png")):
        stem = path.name.removesuffix(".dotplot.png")
        query, ref = stem.split("_vs_", 1)
        qbed = source_bed(args.figdir, stem, query)
        rbed = source_bed(args.figdir, stem, ref)
        qchroms = read_chromosomes(qbed) if qbed else []
        rchroms = read_chromosomes(rbed) if rbed else []
        qlabels = hap_labels(query, qchroms)
        rlabels = hap_labels(ref, rchroms)
        if ref == "ppa" and rbed:
            rlabels = ppa_labels(rbed, ppa_map)
        if query == "ppa" and qbed:
            qlabels = ppa_labels(qbed, ppa_map)
        # Only copy/relabel images where a requested PZI/PSP/PPA label exists.
        if not qlabels and not rlabels:
            continue
        out = args.outdir / f"{stem}.haplotype_labels.png"
        ntop, nright = relabel(path, out, qlabels, rlabels, qbed, rbed)
        rows.append((path.name, out.name, query, ref, ntop, nright, str(qbed or ""), str(rbed or "")))
        print(f"[OK] {path.name}: query_labels={ntop} ref_labels={nright}")
    with (args.outdir / "manifest.tsv").open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(("original_png", "labelled_png", "query", "reference", "query_labels", "reference_labels", "query_bed", "reference_bed"))
        writer.writerows(rows)
    print(f"[OK] Wrote {len(rows)} relabelled dotplots to {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
