#!/usr/bin/env python3
"""Prepare a three-sector, native Circos comparison of PZI, PSP and PPA.

PZI and PSP use their A1/B1 chromosomes; PPA uses its A/B chromosomes.  The
existing single-species tracks are filtered to those chromosomes and retained
at their original base-pair coordinates.  The two native JCVI pairwise block
files are used directly; PZI-PPA links are chained through shared PSP anchors
because no direct native PZI-PPA ``blocks.simple`` file is present.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


SPECIES = ("pzi", "psp", "ppa")
COMMON_ORDER = (
    "C1A", "C6B", "C7B", "C2A", "C5B", "C3A", "C8B", "C7A",
    "C4A", "C1B", "C5A", "C2B", "C6A", "C3B", "C8A", "C4B",
)
TRACKS = (
    "gene_density",
    "te_coverage",
    "gypsy_coverage",
    "copia_coverage",
    "tir_coverage",
    "helitron_coverage",
    "gc_fraction",
    "ncrna_density",
)
PAIR_RGB = {
    "C1": "0,114,178",
    "C2": "230,159,0",
    "C3": "0,158,115",
    "C4": "213,94,0",
    "C5": "204,121,167",
    "C6": "86,180,233",
    "C7": "190,168,0",
    "C8": "80,80,80",
}
LINK_RGB = {
    "pzi_psp": "0,114,178",
    "psp_ppa": "0,158,115",
    "ppa_pzi": "230,159,0",
}
SPECIES_RGB = {
    "pzi": "0,114,178",
    "psp": "0,158,115",
    "ppa": "230,159,0",
}
SOURCE_HAP = {"pzi": "A1", "psp": "A1", "ppa": "A"}
NATIVE_DEFAULT = Path("0903/native_multispecies_A1B1_pcam_pairwise")
CIRCOS_DEFAULT = Path("0903/haploid_repartition/circos")
OUT_DEFAULT = Path("0903/haploid_repartition/circos/pzi_psp_ppa")


def canonical(value: str) -> str | None:
    match = re.search(r"C0*(\d+)_?([AB])", value)
    if not match:
        return None
    number = int(match.group(1))
    side = match.group(2)
    if not 1 <= number <= 8:
        return None
    return f"C{number}{side}"


def namespaced(species: str, chrom: str) -> str:
    return f"{species}_{chrom}"


def parse_karyotype(path: Path) -> list[tuple[str, int, str]]:
    rows: list[tuple[str, int, str]] = []
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = raw.split()
        if len(fields) < 6 or fields[0] != "chr":
            continue
        try:
            length = int(fields[5])
        except ValueError:
            continue
        rows.append((fields[2], length, fields[6] if len(fields) > 6 else "pair1"))
    if not rows:
        raise ValueError(f"No chromosomes found in {path}")
    return rows


def selected_chromosomes(species: str, circos_root: Path) -> dict[str, tuple[str, int]]:
    rows = parse_karyotype(circos_root / species / "data" / "karyotype.txt")
    selected: dict[str, tuple[str, int]] = {}
    for source_id, length, _color in rows:
        chrom = canonical(source_id)
        if not chrom:
            continue
        if species in ("pzi", "psp"):
            keep = source_id.startswith(("A1_", "B1_"))
        else:
            keep = source_id.startswith(("A_", "B_"))
        if keep:
            selected[chrom] = (source_id, length)
    missing = [chrom for chrom in COMMON_ORDER if chrom not in selected]
    if missing:
        raise ValueError(f"{species}: missing selected chromosomes: {','.join(missing)}")
    return selected


def read_bed(path: Path) -> tuple[dict[str, tuple[str, int, int]], dict[str, int]]:
    genes: dict[str, tuple[str, int, int]] = {}
    lengths: dict[str, int] = defaultdict(int)
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            genes.setdefault(fields[3], (fields[0], start, end))
            lengths[fields[0]] = max(lengths[fields[0]], end)
    if not genes:
        raise ValueError(f"Empty BED: {path}")
    return genes, dict(lengths)


def parse_blocks(path: Path):
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if "*" in line:
                line = line.split("*", 1)[1].strip()
            fields = line.split()
            if len(fields) != 6:
                continue
            try:
                yield fields[0], fields[1], fields[2], fields[3], int(float(fields[4])), fields[5]
            except ValueError:
                continue


def block_records(path: Path, left_bed, right_bed):
    for a, b, c, d, score, orient in parse_blocks(path):
        if any(g not in left_bed for g in (a, b)) or any(g not in right_bed for g in (c, d)):
            continue
        left_a, left_b = left_bed[a], left_bed[b]
        right_a, right_b = right_bed[c], right_bed[d]
        if left_a[0] != left_b[0] or right_a[0] != right_b[0]:
            continue
        yield {
            "left_genes": (a, b),
            "right_genes": (c, d),
            "left_chrom": left_a[0],
            "left_start": min(left_a[1], left_b[1]),
            "left_end": max(left_a[2], left_b[2]),
            "right_chrom": right_a[0],
            "right_start": min(right_a[1], right_b[1]),
            "right_end": max(right_a[2], right_b[2]),
            "score": score,
            "orient": orient,
        }


def write_karyotype(data: Path, selected: dict[str, dict[str, tuple[str, int]]]) -> list[str]:
    order: list[str] = []
    with (data / "karyotype.txt").open("w", encoding="ascii", newline="\n") as handle:
        for species in SPECIES:
            for chrom in COMMON_ORDER:
                source_id, length = selected[species][chrom]
                internal = namespaced(species, chrom)
                handle.write(f"chr - {internal} {chrom} 0 {length} {SPECIES_RGB[species]}\n")
                order.append(internal)
    return order


def remap_tracks(data: Path, selected: dict[str, dict[str, tuple[str, int]]], circos_root: Path) -> dict[str, list[float]]:
    values: dict[str, list[float]] = defaultdict(list)
    for species in SPECIES:
        source_ids = {source_id: chrom for chrom, (source_id, _length) in selected[species].items()}
        for track in TRACKS:
            source = circos_root / species / "data" / f"{track}.txt"
            if not source.is_file():
                raise FileNotFoundError(source)
            rows: list[str] = []
            for raw in source.read_text(encoding="utf-8", errors="replace").splitlines():
                fields = raw.split()
                if len(fields) < 4 or fields[0] not in source_ids:
                    continue
                chrom = source_ids[fields[0]]
                try:
                    value = float(fields[3])
                except ValueError:
                    continue
                values[track].append(value)
                rows.append(f"{namespaced(species, chrom)} {fields[1]} {fields[2]} {fields[3]}")
            destination = data / f"{species}_{track}.txt"
            destination.write_text("\n".join(rows) + "\n", encoding="ascii")
            if not rows:
                raise ValueError(f"No selected rows in {source}")
    return values


def link_line(left_species: str, left_chrom: str, left_start: int, left_end: int,
              right_species: str, right_chrom: str, right_start: int, right_end: int,
              color: str, orient: str, z: int = 100) -> str:
    twist = "yes" if orient == "-" else "no"
    return (
        f"{namespaced(left_species, left_chrom)} {left_start} {left_end} "
        f"{namespaced(right_species, right_chrom)} {right_start} {right_end} "
        f"color={color}_a55,z={z},twist={twist}"
    )


def write_direct_links(data: Path, native: Path, beds, filename: str,
                       left_species: str, right_species: str, color: str) -> int:
    rows: list[str] = []
    left_bed, right_bed = beds[left_species][0], beds[right_species][0]
    path = native / "anchors" / filename
    for rec in block_records(path, left_bed, right_bed):
        if canonical(rec["left_chrom"]) not in COMMON_ORDER or canonical(rec["right_chrom"]) not in COMMON_ORDER:
            continue
        left_chrom = canonical(rec["left_chrom"])
        right_chrom = canonical(rec["right_chrom"])
        rows.append(link_line(left_species, left_chrom, rec["left_start"], rec["left_end"],
                              right_species, right_chrom, rec["right_start"], rec["right_end"], color, rec["orient"]))
    output = data / f"links_{left_species}_{right_species}.txt"
    output.write_text("\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_chained_pzi_ppa(data: Path, native: Path, beds) -> int:
    """Chain PZI-PSP and PSP-PPA blocks by shared PSP gene anchors."""
    pzi_bed, psp_bed, ppa_bed = beds["pzi"][0], beds["psp"][0], beds["ppa"][0]
    first = list(block_records(native / "anchors" / "pzi.psp.blocks.simple", pzi_bed, psp_bed))
    second = list(block_records(native / "anchors" / "psp.ppa.blocks.simple", psp_bed, ppa_bed))
    by_psp_gene: dict[str, list[dict]] = defaultdict(list)
    for rec in second:
        for gene in rec["left_genes"]:
            by_psp_gene[gene].append(rec)
    rows: list[str] = []
    seen: set[tuple] = set()
    for rec in first:
        candidates: list[tuple[int, int, dict]] = []
        for gene in rec["right_genes"]:
            for target in by_psp_gene.get(gene, ()):
                overlap = len(set(rec["right_genes"]) & set(target["left_genes"]))
                candidates.append((overlap, target["score"], target))
        if not candidates:
            continue
        _overlap, _score, target = max(candidates, key=lambda item: (item[0], item[1]))
        left_chrom = canonical(rec["left_chrom"])
        right_chrom = canonical(target["right_chrom"])
        if left_chrom not in COMMON_ORDER or right_chrom not in COMMON_ORDER:
            continue
        orient = "+" if rec["orient"] == target["orient"] else "-"
        key = (left_chrom, rec["left_start"], rec["left_end"], right_chrom, target["right_start"], target["right_end"])
        if key in seen:
            continue
        seen.add(key)
        rows.append(link_line("pzi", left_chrom, rec["left_start"], rec["left_end"],
                              "ppa", right_chrom, target["right_start"], target["right_end"], "link_ppa_pzi", orient))
    (data / "links_ppa_pzi.txt").write_text("\n".join(rows) + "\n", encoding="ascii")
    return len(rows)


def write_config(out_root: Path, order: list[str], selected: dict[str, dict[str, tuple[str, int]]], maxima: dict[str, float]) -> None:
    totals = {species: sum(length for _source, length in selected[species].values()) for species in SPECIES}
    sector_fraction = 110.0 / 360.0
    scale_lines: list[str] = []
    for species in SPECIES:
        for chrom in COMMON_ORDER:
            length = selected[species][chrom][1]
            scale = sector_fraction * length / totals[species]
            scale_lines.append(f"{namespaced(species, chrom)}={scale:.10f}r")
    scale_text = ",".join(scale_lines)
    pairwise = (
        f"<pairwise {namespaced('pzi', COMMON_ORDER[-1])} {namespaced('psp', COMMON_ORDER[0])}>\n"
        "spacing = 277.7778r\n</pairwise>\n"
        f"<pairwise {namespaced('psp', COMMON_ORDER[-1])} {namespaced('ppa', COMMON_ORDER[0])}>\n"
        "spacing = 277.7778r\n</pairwise>\n"
        f"<pairwise {namespaced('ppa', COMMON_ORDER[-1])} {namespaced('pzi', COMMON_ORDER[0])}>\n"
        "spacing = 277.7778r\n</pairwise>"
    )
    colors = [f"pair{i} = {PAIR_RGB[f'C{i}']}" for i in range(1, 9)]
    colors.extend(f"link_{name} = {rgb}" for name, rgb in LINK_RGB.items())
    colors.extend((
        "gene_track = 85,90,95",
        "te_track = 112,135,52",
        "gypsy_track = 196,91,82",
        "copia_track = 76,120,168",
        "tir_track = 126,88,160",
        "helitron_track = 40,150,145",
        "gc_track = 42,88,108",
        "ncrna_track = 111,158,158",
    ))
    plot_defs: list[str] = []
    layout = {
        "gene_density": ("heatmap", "0.932r", "0.972r", 0.0, maxima["gene_density"], "color = spectral-11-div\ncolor_mapping = 0"),
        "te_coverage": ("histogram", "0.875r", "0.915r", 0.0, 1.0, "fill_color = te_track\ncolor = te_track"),
        "gypsy_coverage": ("histogram", "0.818r", "0.858r", 0.0, maxima["gypsy_coverage"], "fill_color = gypsy_track\ncolor = gypsy_track"),
        "copia_coverage": ("histogram", "0.761r", "0.801r", 0.0, maxima["copia_coverage"], "fill_color = copia_track\ncolor = copia_track"),
        "tir_coverage": ("histogram", "0.704r", "0.744r", 0.0, maxima["tir_coverage"], "fill_color = tir_track\ncolor = tir_track"),
        "helitron_coverage": ("histogram", "0.647r", "0.687r", 0.0, maxima["helitron_coverage"], "fill_color = helitron_track\ncolor = helitron_track"),
        "gc_fraction": ("line", "0.590r", "0.630r", maxima["gc_min"], maxima["gc_max"], "color = gc_track\nthickness = 1.5p"),
        "ncrna_density": ("histogram", "0.535r", "0.575r", 0.0, maxima["ncrna_density"], "fill_color = ncrna_track\ncolor = ncrna_track"),
    }
    for species in SPECIES:
        for track in TRACKS:
            kind, r0, r1, minimum, maximum, options = layout[track]
            thickness = "" if kind == "line" else "\n thickness = 0"
            plot_defs.append(
                f"<plot>\n type = {kind}\n file = data/{species}_{track}.txt\n"
                f" r0 = {r0}\n r1 = {r1}\n min = {minimum:.6f}\n max = {maximum:.6f}\n"
                f" {options}{thickness}\n</plot>"
            )
    config = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
{chr(10).join(colors)}
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {';'.join(order)}
chromosomes_scale = {scale_text}
<image>
background = white
dir = figures
file = pzi_psp_ppa_circos.png
png = yes
svg = yes
pdf = yes
radius = 900p
angle_offset = -205
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.0001r
{pairwise}
</spacing>
radius = 0.900r
thickness = 7.06p
fill = yes
stroke_color = white
stroke_thickness = 1p
show_label = yes
label_radius = 1.045r
label_size = 15p
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
spacing = 50u
size = 5p
show_label = yes
label_size = 9p
label_offset = 2p
</tick>
</ticks>
<plots>
{chr(10).join(plot_defs)}
</plots>
<links>
<link>
file = data/links_pzi_psp.txt
radius = 0.520r
bezier_radius = 0.100r
ribbon = yes
thickness = 0.8p
crest = 0.45
</link>
<link>
file = data/links_psp_ppa.txt
radius = 0.520r
bezier_radius = 0.100r
ribbon = yes
thickness = 0.8p
crest = 0.45
</link>
<link>
file = data/links_ppa_pzi.txt
radius = 0.520r
bezier_radius = 0.100r
ribbon = yes
thickness = 0.8p
crest = 0.45
</link>
</links>
"""
    (out_root / "conf" / "circos.conf").write_text(config, encoding="ascii")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--native-dir", type=Path, default=NATIVE_DEFAULT)
    parser.add_argument("--circos-root", type=Path, default=CIRCOS_DEFAULT)
    parser.add_argument("--out-root", type=Path, default=OUT_DEFAULT)
    args = parser.parse_args()
    native = args.native_dir.resolve()
    circos_root = args.circos_root.resolve()
    out_root = args.out_root.resolve()
    data = out_root / "data"
    (out_root / "conf").mkdir(parents=True, exist_ok=True)
    data.mkdir(parents=True, exist_ok=True)
    (out_root / "figures").mkdir(parents=True, exist_ok=True)

    selected = {species: selected_chromosomes(species, circos_root) for species in SPECIES}
    order = write_karyotype(data, selected)
    values = remap_tracks(data, selected, circos_root)
    beds = {species: read_bed(native / "bed" / f"{species}.bed") for species in SPECIES}
    counts = {
        "pzi_psp": write_direct_links(data, native, beds, "pzi.psp.blocks.simple", "pzi", "psp", "link_pzi_psp"),
        "psp_ppa": write_direct_links(data, native, beds, "psp.ppa.blocks.simple", "psp", "ppa", "link_psp_ppa"),
        "ppa_pzi": write_chained_pzi_ppa(data, native, beds),
    }
    all_gc = values["gc_fraction"]
    maxima = {
        "gene_density": max(values["gene_density"] or [1.0]),
        "gypsy_coverage": max(values["gypsy_coverage"] or [0.01]),
        "copia_coverage": max(values["copia_coverage"] or [0.01]),
        "tir_coverage": max(values["tir_coverage"] or [0.01]),
        "helitron_coverage": max(values["helitron_coverage"] or [0.01]),
        "ncrna_density": max(values["ncrna_density"] or [1.0]),
        "gc_min": min(all_gc or [0.0]),
        "gc_max": max(all_gc or [1.0]),
    }
    padding = max(0.005, (maxima["gc_max"] - maxima["gc_min"]) * 0.08)
    maxima["gc_min"] = max(0.0, maxima["gc_min"] - padding)
    maxima["gc_max"] = min(1.0, maxima["gc_max"] + padding)
    write_config(out_root, order, selected, maxima)

    with (data / "chromosome_map.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("species", "internal_id", "source_id", "canonical", "length_bp"))
        for species in SPECIES:
            for chrom in COMMON_ORDER:
                source_id, length = selected[species][chrom]
                writer.writerow((species, namespaced(species, chrom), source_id, chrom, length))
    with (out_root / "figure_manifest.tsv").open("w", encoding="ascii", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("item", "value"))
        writer.writerow(("layout", "PZI upper-left, PSP upper-right, PPA lower; 110 degrees per species; 10 degree gaps"))
        writer.writerow(("chromosomes", "PZI/PSP A1+B1; PPA A+B; 16 chromosomes per species"))
        writer.writerow(("tracks", "gene_density, te_coverage, gypsy_coverage, copia_coverage, tir_coverage, helitron_coverage, gc_fraction, ncrna_density"))
        writer.writerow(("link_source_pzi_psp", "native JCVI pzi.psp.blocks.simple"))
        writer.writerow(("link_source_psp_ppa", "native JCVI psp.ppa.blocks.simple"))
        writer.writerow(("link_source_ppa_pzi", "chained native JCVI blocks through shared PSP anchors"))
        for name, count in counts.items():
            writer.writerow((f"links_{name}", count))
    (out_root / "README_methods.txt").write_text(
        "Three-sector native Circos comparison. PZI and PSP use A1/B1 and PPA uses A/B. "
        "Each sector occupies approximately 110 degrees and species gaps are 10 degrees. "
        "Tracks are the existing 1-Mb gene, TE-class, GC and ncRNA tracks filtered to the selected chromosomes. "
        "PZI-PSP and PSP-PPA links are native JCVI blocks.simple; PPA-PZI links are conservative transitive links "
        "from shared PSP anchors because no direct native PZI-PPA blocks.simple is available.\n",
        encoding="ascii",
    )
    print(f"[OK] prepared {out_root}")
    for name, count in counts.items():
        print(f"[OK] {name}: {count} links")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
