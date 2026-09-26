#!/usr/bin/env python3
"""Prepare native Circos inputs for the revised PZI/PPA/PSP assignments.

The script only consumes the revised FASTA/GFF3/EDTA package and JCVI
``blocks.simple`` files produced from the same revised protein/BED inputs.
It creates one self-contained Circos directory per species, including
chromosome ideograms, 1-Mb gene/TE/GC tracks and JCVI synteny ribbons.
"""

from __future__ import annotations

import argparse
import math
import os
from collections import defaultdict
from pathlib import Path
from typing import Iterable


SPECIES = {
    "pzi": ("A1", "A2", "B1", "B2"),
    "psp": ("A1", "A2", "B1", "B2"),
    "ppa": ("A", "B"),
}
HAP_COLORS = {"A1": "hapa1", "A2": "hapa2", "B1": "hapb1", "B2": "hapb2", "A": "hapa1", "B": "hapb1"}
HAP_RGB = {"hapa1": "55,95,145", "hapa2": "120,168,205", "hapb1": "186,85,40", "hapb2": "222,151,55"}
PAIR_RGB = (
    ("c01", "0,114,178"),
    ("c02", "230,159,0"),
    ("c03", "0,158,115"),
    ("c04", "213,94,0"),
    ("c05", "204,121,167"),
    ("c06", "86,180,233"),
    ("c07", "190,168,0"),
    ("c08", "80,80,80"),
)


def parse_attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";\r\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value.strip('"')
        elif " " in item:
            key, value = item.split(None, 1)
            result[key] = value.strip('"')
    return result


def fasta_sizes_and_windows(path: Path, window: int) -> tuple[dict[str, int], dict[str, list[tuple[int, int, float]]]]:
    """Stream FASTA once and return chromosome lengths plus GC fractions."""
    sizes: dict[str, int] = {}
    windows: dict[str, list[tuple[int, int, float]]] = {}
    name: str | None = None
    pos = 0
    # Per-window GC state is tracked below; the chromosome-level counters are
    # intentionally not retained because only the exact FASTA length matters.

    def finish() -> None:
        nonlocal name, pos, current_windows, current_start, current_gc, current_total
        if name is None:
            return
        if current_total:
            current_windows.append((current_start, pos, current_gc / current_total))
        sizes[name] = pos
        windows[name] = current_windows
        pos = 0
        current_windows = []
        current_start = 0
        current_gc = 0
        current_total = 0

    current_windows: list[tuple[int, int, float]] = []
    current_start = 0
    current_gc = 0
    current_total = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                finish()
                name = raw[1:].split()[0]
                current_windows = []
                current_start = 0
                current_gc = 0
                current_total = 0
                continue
            sequence = raw.strip().upper()
            offset = 0
            while offset < len(sequence):
                used = pos - current_start
                remaining = window - used
                chunk = sequence[offset : offset + remaining]
                current_gc += chunk.count("G") + chunk.count("C")
                current_total += len(chunk)
                pos += len(chunk)
                offset += len(chunk)
                if pos - current_start == window:
                    current_windows.append((current_start, pos, current_gc / current_total if current_total else 0.0))
                    current_start = pos
                    current_gc = 0
                    current_total = 0
    finish()
    if not sizes:
        raise ValueError(f"No FASTA sequences found: {path}")
    return sizes, windows


def gene_density(path: Path, sizes: dict[str, int], window: int) -> dict[str, list[int]]:
    result = {chrom: [0] * math.ceil(length / window) for chrom, length in sizes.items()}
    seen: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[2].lower() != "gene":
                continue
            chrom = fields[0]
            if chrom not in result:
                continue
            attrs = parse_attrs(fields[8])
            feature_id = attrs.get("ID") or f"{chrom}:{fields[3]}-{fields[4]}"
            if feature_id in seen:
                continue
            seen.add(feature_id)
            start = max(0, int(fields[3]) - 1)
            result[chrom][min(start // window, len(result[chrom]) - 1)] += 1
    return result


def te_class(fields: list[str], attrs: dict[str, str]) -> str:
    text = " ".join((fields[2], attrs.get("classification", ""), attrs.get("Name", ""))).lower()
    if "gypsy" in text:
        return "gypsy"
    if "copia" in text:
        return "copia"
    if "helitron" in text:
        return "helitron"
    if "tir" in text or any(token in text for token in ("dna/dta", "dna/dtc", "dna/dtm", "dna/dth")):
        return "tir"
    return "other"


def te_intervals(path: Path, sizes: dict[str, int]) -> dict[str, list[tuple[str, int, int]]]:
    result: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9 or fields[0] not in sizes:
                continue
            try:
                start = max(0, int(fields[3]) - 1)
                end = min(sizes[fields[0]], int(fields[4]))
            except ValueError:
                continue
            if end <= start:
                continue
            kind = te_class(fields, parse_attrs(fields[8]))
            result["all"].append((fields[0], start, end))
            result[kind].append((fields[0], start, end))
    return result


def merge_intervals(intervals: Iterable[tuple[str, int, int]]) -> list[tuple[str, int, int]]:
    grouped: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for chrom, start, end in intervals:
        grouped[chrom].append((start, end))
    merged: list[tuple[str, int, int]] = []
    for chrom, rows in grouped.items():
        rows.sort()
        current_start, current_end = rows[0]
        for start, end in rows[1:]:
            if start <= current_end:
                current_end = max(current_end, end)
            else:
                merged.append((chrom, current_start, current_end))
                current_start, current_end = start, end
        merged.append((chrom, current_start, current_end))
    return merged


def window_coverage(intervals: Iterable[tuple[str, int, int]], sizes: dict[str, int], window: int) -> dict[str, list[int]]:
    result = {chrom: [0] * math.ceil(length / window) for chrom, length in sizes.items()}
    for chrom, start, end in merge_intervals(intervals):
        first = start // window
        last = (end - 1) // window
        for index in range(first, last + 1):
            left = max(start, index * window)
            right = min(end, (index + 1) * window)
            if right > left:
                result[chrom][index] += right - left
    return result


def internal(hap: str, chrom: str) -> str:
    return f"{hap}_{chrom}"


def parse_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    result = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) >= 4:
                result[fields[3]] = (fields[0], int(fields[1]), int(fields[2]))
    if not result:
        raise ValueError(f"Empty BED: {path}")
    return result


def parse_blocks(path: Path, left_bed: dict, right_bed: dict):
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 6:
                continue
            a1, a2, b1, b2, score, orient = fields
            if any(item not in left_bed for item in (a1, a2)) or any(item not in right_bed for item in (b1, b2)):
                continue
            l1, l2 = left_bed[a1], left_bed[a2]
            r1, r2 = right_bed[b1], right_bed[b2]
            if l1[0] != l2[0] or r1[0] != r2[0]:
                continue
            yield (l1[0], min(l1[1], l2[1]), max(l1[2], l2[2]), r1[0], min(r1[1], r2[1]), max(r1[2], r2[2]), score, orient)


def pair_specs(haps: tuple[str, ...]) -> tuple[tuple[str, str, str], ...]:
    if haps == ("A1", "A2", "B1", "B2"):
        return (("A1", "A2", "background"), ("B1", "B2", "background"), ("A1", "B1", "highlight"), ("A2", "B2", "highlight"))
    return (("A", "B", "highlight"),)


def chromosome_order(haps: tuple[str, ...], sizes: dict[str, dict[str, int]]) -> list[tuple[str, str]]:
    rows: list[tuple[str, str]] = []
    for hap in haps:
        chroms = sorted(sizes[hap], key=lambda c: (c[-1], c))
        if hap in {"A2", "B2", "A", "B"}:
            chroms = list(reversed(chroms))
        rows.extend((hap, chrom) for chrom in chroms)
    return rows


def write_config(path: Path, species: str, haps: tuple[str, ...], order: list[tuple[str, str]], scales: dict[str, float], maxima: dict[str, float], gc_range: tuple[float, float]) -> None:
    colors = [f"{name} = {rgb}" for name, rgb in HAP_RGB.items()]
    colors.extend(f"{name} = {rgb}" for name, rgb in PAIR_RGB)
    # Keep the background-link colour as plain RGB.  Some Circos builds do
    # not accept the legacy ``a48`` alpha suffix in custom colour definitions.
    # Link z-order already keeps these anchors visually subordinate.
    colors += ["gene_track = 85,90,95", "te_track = 112,135,52", "gypsy_track = 196,91,82", "copia_track = 76,120,168", "tir_track = 126,88,160", "helitron_track = 40,150,145", "gc_track = 42,88,108", "within_gray = 165,165,165"]
    chrom_text = ";".join(internal(hap, chrom) for hap, chrom in order)
    scale_text = ",".join(f"{internal(hap, chrom)}={scales[hap]:.8f}" for hap, chrom in order)
    reverse_text = ",".join(internal(hap, chrom) for hap, chrom in order if hap in {"A2", "B2", "B", "A"})
    gap_left = internal(haps[0], "C1A")
    gap_right = internal(haps[-2] if len(haps) > 2 else haps[-1], "C1B")
    config = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
{os.linesep.join(colors)}
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {chrom_text}
chromosomes_scale = {scale_text}
chromosomes_reverse = {reverse_text}
<image>
background = white
dir = figures
file = {species}_repartition_circos.png
png = yes
svg = yes
radius = 1800p
angle_offset = -90
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.0025r
<pairwise {gap_left} {gap_right}>
spacing = 12.4r
</pairwise>
</spacing>
radius = 0.940r
thickness = 22p
fill = yes
stroke_color = white
stroke_thickness = 1p
show_label = yes
label_radius = 1.030r
label_size = 14p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 1r
color = grey
thickness = 1p
multiplier = 1e-6
format = %d
<tick>
spacing = 10u
size = 7p
show_label = yes
label_size = 9p
label_offset = 3p
</tick>
</ticks>
<plots>
<plot>
type = histogram
file = data/gene_density.txt
r0 = 0.932r
r1 = 0.972r
min = 0
max = {maxima['gene']:.6f}
fill_color = gene_track
color = gene_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/te_coverage.txt
r0 = 0.875r
r1 = 0.915r
min = 0
max = 1
fill_color = te_track
color = te_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/gypsy_coverage.txt
r0 = 0.818r
r1 = 0.858r
min = 0
max = {maxima['gypsy']:.6f}
fill_color = gypsy_track
color = gypsy_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/copia_coverage.txt
r0 = 0.761r
r1 = 0.801r
min = 0
max = {maxima['copia']:.6f}
fill_color = copia_track
color = copia_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/tir_coverage.txt
r0 = 0.704r
r1 = 0.744r
min = 0
max = {maxima['tir']:.6f}
fill_color = tir_track
color = tir_track
thickness = 0
</plot>
<plot>
type = histogram
file = data/helitron_coverage.txt
r0 = 0.647r
r1 = 0.687r
min = 0
max = {maxima['helitron']:.6f}
fill_color = helitron_track
color = helitron_track
thickness = 0
</plot>
<plot>
type = line
file = data/gc_fraction.txt
r0 = 0.590r
r1 = 0.630r
min = {gc_range[0]:.6f}
max = {gc_range[1]:.6f}
color = gc_track
thickness = 2p
</plot>
</plots>
<links>
<link>
file = data/links_background.txt
radius = 0.555r
bezier_radius = 0.065r
ribbon = yes
thickness = 1p
crest = 0.45
</link>
<link>
file = data/links_highlight.txt
radius = 0.545r
bezier_radius = 0.055r
ribbon = yes
thickness = 1p
crest = 0.45
</link>
</links>
"""
    path.write_text(config, encoding="ascii")


def prepare_species(root: Path, out_root: Path, jcvi_root: Path, species: str, window: int) -> dict[str, int]:
    haps = SPECIES[species]
    species_root = out_root / species
    data = species_root / "data"
    conf = species_root / "conf"
    figures = species_root / "figures"
    for directory in (data, conf, figures):
        directory.mkdir(parents=True, exist_ok=True)

    sizes: dict[str, dict[str, int]] = {}
    windows: dict[str, dict[str, list[tuple[int, int, float]]]] = {}
    genes: dict[str, dict[str, list[int]]] = {}
    coverages: dict[str, dict[str, dict[str, list[int]]]] = {}
    manifests = []
    for hap in haps:
        prefix = f"{species}_hap{hap}"
        genome = root / "haploid_repartition" / "inputs" / species / "genome" / f"{prefix}.fa"
        gff = root / "haploid_repartition" / "inputs" / species / "gff" / f"{prefix}.gff3"
        te = root / "haploid_repartition" / "inputs" / species / "edta" / hap / f"{prefix}.fa.mod.EDTA.TEanno.gff3"
        if not te.is_file():
            # Compatibility with packages generated before the canonical
            # inputs/{species}/edta/{hap} layout was introduced.
            te = root / "haploid_repartition" / "inputs" / species / "repeat" / f"{prefix}_edta" / f"{prefix}.fa.mod.EDTA.TEanno.gff3"
        for required in (genome, gff, te):
            if not required.is_file() or required.stat().st_size == 0:
                raise SystemExit(f"Missing revised input: {required}")
        sizes[hap], windows[hap] = fasta_sizes_and_windows(genome, window)
        genes[hap] = gene_density(gff, sizes[hap], window)
        raw_te = te_intervals(te, sizes[hap])
        coverages[hap] = {kind: window_coverage(rows, sizes[hap], window) for kind, rows in raw_te.items()}
        manifests.append((hap, genome, gff, te))

    order = chromosome_order(haps, sizes)
    with (data / "karyotype.txt").open("w", encoding="ascii") as handle:
        for hap, chrom in order:
            label = f"{chrom}_{hap}"
            handle.write(f"chr - {internal(hap, chrom)} {label} 0 {sizes[hap][chrom]} {HAP_COLORS[hap]}\n")
    with (data / "chromosome_map.tsv").open("w", encoding="ascii") as handle:
        handle.write("internal_chromosome\thaplotype\tshared_chromosome\tlength_bp\n")
        for hap, chrom in order:
            handle.write(f"{internal(hap, chrom)}\t{hap}\t{chrom}\t{sizes[hap][chrom]}\n")

    totals = {hap: sum(sizes[hap].values()) for hap in haps}
    mean_total = sum(totals.values()) / len(totals)
    scales = {hap: mean_total / total for hap, total in totals.items()}
    with (data / "layout_scaling.tsv").open("w", encoding="ascii") as handle:
        handle.write("haplotype\ttotal_bp\tangular_scale\n")
        for hap in haps:
            handle.write(f"{hap}\t{totals[hap]}\t{scales[hap]:.8f}\n")

    all_gene = []
    all_gc = []
    max_values = {"gene": 1.0, "gypsy": 0.01, "copia": 0.01, "tir": 0.01, "helitron": 0.01}
    with (data / "window_tracks.tsv").open("w", encoding="ascii") as raw, \
         (data / "gene_density.txt").open("w", encoding="ascii") as gene_out, \
         (data / "te_coverage.txt").open("w", encoding="ascii") as te_out, \
         (data / "gypsy_coverage.txt").open("w", encoding="ascii") as gypsy_out, \
         (data / "copia_coverage.txt").open("w", encoding="ascii") as copia_out, \
         (data / "tir_coverage.txt").open("w", encoding="ascii") as tir_out, \
         (data / "helitron_coverage.txt").open("w", encoding="ascii") as helitron_out, \
         (data / "gc_fraction.txt").open("w", encoding="ascii") as gc_out:
        raw.write("chrom\tstart\tend\tgene_count\tte_fraction\tgypsy_fraction\tcopia_fraction\ttir_fraction\thel itron_fraction\tgc_fraction\n".replace("hel itron", "helitron"))
        for hap, chrom in order:
            target = internal(hap, chrom)
            for index, (start, end, gc) in enumerate(windows[hap][chrom]):
                width = end - start
                gene_count = genes[hap][chrom][index]
                fractions = {}
                for kind in ("all", "gypsy", "copia", "tir", "helitron"):
                    fractions[kind] = coverages[hap].get(kind, {}).get(chrom, [0] * len(windows[hap][chrom]))[index] / width
                max_values["gene"] = max(max_values["gene"], gene_count)
                for kind in ("gypsy", "copia", "tir", "helitron"):
                    max_values[kind] = max(max_values[kind], fractions[kind])
                all_gene.append(gene_count)
                all_gc.append(gc)
                raw.write(f"{target}\t{start}\t{end}\t{gene_count}\t{fractions['all']:.6f}\t{fractions['gypsy']:.6f}\t{fractions['copia']:.6f}\t{fractions['tir']:.6f}\t{fractions['helitron']:.6f}\t{gc:.6f}\n")
                gene_out.write(f"{target} {start} {end} {gene_count}\n")
                te_out.write(f"{target} {start} {end} {fractions['all']:.6f}\n")
                gypsy_out.write(f"{target} {start} {end} {fractions['gypsy']:.6f}\n")
                copia_out.write(f"{target} {start} {end} {fractions['copia']:.6f}\n")
                tir_out.write(f"{target} {start} {end} {fractions['tir']:.6f}\n")
                helitron_out.write(f"{target} {start} {end} {fractions['helitron']:.6f}\n")
                gc_out.write(f"{target} {start} {end} {gc:.6f}\n")

    # JCVI link generation is deliberately strict: a missing block file is a
    # failed reproducibility check, not a reason to draw an incomplete figure.
    bed_cache = {}
    input_root = jcvi_root.parent / "jcvi_inputs"
    for hap in haps:
        bed = input_root / f"{species}_{hap}.bed"
        if not bed.is_file() or bed.stat().st_size == 0:
            raise SystemExit(f"Missing revised JCVI BED: {bed}")
        bed_cache[hap] = parse_bed(bed)
    summary = []
    with (data / "links_background.txt").open("w", encoding="ascii") as background, (data / "links_highlight.txt").open("w", encoding="ascii") as highlight:
        for left, right, layer in pair_specs(haps):
            # The JCVI job prefixes each haplotype with its species in the
            # directory name (for example psp_A1_psp_A2 and ppa_A_ppa_B).
            block_file = jcvi_root / f"{species}_{left}_{species}_{right}" / "blocks.simple"
            if not block_file.is_file() or block_file.stat().st_size == 0:
                raise SystemExit(f"Missing revised JCVI blocks: {block_file}")
            output = background if layer == "background" else highlight
            count = 0
            cross_chromosome = 0
            for lchr, ls, le, rchr, rs, re, score, orient in parse_blocks(block_file, bed_cache[left], bed_cache[right]):
                left_id = internal(left, lchr)
                right_id = internal(right, rchr)
                if lchr != rchr:
                    cross_chromosome += 1
                if layer == "background":
                    color = "within_gray"
                    z = 5
                else:
                    try:
                        color = dict(PAIR_RGB)[f"c{int(lchr[1:-1]):02d}"]
                    except (KeyError, ValueError):
                        color = "c01"
                    z = 100
                twist = ",twist=yes" if orient == "-" else ",twist=no"
                output.write(f"{left_id} {ls} {le} {right_id} {rs} {re} color={color},z={z}{twist}\n")
                count += 1
            summary.append((f"{left}-{right}", layer, count, cross_chromosome, str(block_file)))
    # Absolute paths may contain non-ASCII characters (for example a Chinese
    # Windows user directory mounted under WSL), so keep audit metadata UTF-8.
    with (data / "jcvi_link_summary.tsv").open("w", encoding="utf-8") as handle:
        handle.write("pair\tlayer\tblocks\tcross_chromosome_blocks\tfile\n")
        for row in summary:
            handle.write("\t".join(map(str, row)) + "\n")

    gc_low = min(all_gc) if all_gc else 0.0
    gc_high = max(all_gc) if all_gc else 1.0
    padding = max(0.005, (gc_high - gc_low) * 0.08)
    gc_range = (max(0.0, gc_low - padding), min(1.0, gc_high + padding))
    write_config(conf / "circos.conf", species, haps, order, scales, max_values, gc_range)
    with (data / "input_manifest.tsv").open("w", encoding="utf-8") as handle:
        handle.write("haplotype\tgenome\tgene_gff3\tte_gff3\n")
        for hap, genome, gff, te in manifests:
            handle.write(f"{hap}\t{genome}\t{gff}\t{te}\n")
    with (data / "track_layout.tsv").open("w", encoding="ascii") as handle:
        handle.write("track\tr0\tr1\tmin\tmax\n")
        handle.write(f"gene_density\t0.932r\t0.972r\t0\t{max_values['gene']:.6f}\n")
        handle.write("all_TE\t0.875r\t0.915r\t0\t1\n")
        for kind, r0, r1 in (("Gypsy", "0.818r", "0.858r"), ("Copia", "0.761r", "0.801r"), ("TIR_DNA", "0.704r", "0.744r"), ("Helitron", "0.647r", "0.687r")):
            key = kind.lower().replace("_dna", "").replace("tir", "tir")
            if kind == "TIR_DNA":
                key = "tir"
            handle.write(f"{kind}\t{r0}\t{r1}\t0\t{max_values[key]:.6f}\n")
        handle.write(f"GC\t0.590r\t0.630r\t{gc_range[0]:.6f}\t{gc_range[1]:.6f}\n")
    return {"windows": sum(len(v) for by_chrom in windows.values() for v in by_chrom.values()), "links": sum(row[2] for row in summary)}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--jcvi-root", type=Path, required=True)
    parser.add_argument("--window", type=int, default=1_000_000)
    args = parser.parse_args()
    manifest_rows = []
    for species in SPECIES:
        stats = prepare_species(args.root.resolve(), args.out_root.resolve(), args.jcvi_root.resolve(), species, args.window)
        manifest_rows.append((species, stats["windows"], stats["links"], f"{species}/figures/{species}_repartition_circos.png", f"{species}/figures/{species}_repartition_circos.svg"))
        print(f"[OK] {species}: {stats['windows']} windows, {stats['links']} JCVI links")
    manifest = args.out_root.resolve() / "figure_manifest.tsv"
    with manifest.open("w", encoding="ascii") as handle:
        handle.write("species\twindow_count\tjcvi_link_count\tpng\tsvg\n")
        for row in manifest_rows:
            handle.write("\t".join(map(str, row)) + "\n")
    (args.out_root.resolve() / "README_methods.txt").write_text(
        "Revised native Circos outputs generated from haploid_repartition/inputs.\n"
        "All TE and class-specific intervals were unioned within 1-Mb windows;\n"
        "JCVI ribbons were converted from blocks.simple generated from the\n"
        "revised BED/protein inputs. A1-A2 and B1-B2 are phase anchors;\n"
        "colored A-B ribbons are subgenome comparisons, not SV calls.\n",
        encoding="ascii",
    )


if __name__ == "__main__":
    main()
