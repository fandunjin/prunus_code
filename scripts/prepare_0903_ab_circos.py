#!/usr/bin/env python3
"""Prepare three standalone CircOS A-B subgenome figures for the 0903 inputs.

Only existing anchors.simple files for A-B comparisons are consumed. A1-A2 and
B1-B2 links are intentionally not written. PZI and PSP use phase-matched A1-B1
and A2-B2 comparisons; PPA uses A-B.  Chromosomes follow the explicit
cross-subgenome order in ``CHROMOSOME_ORDER``.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


SPECIES = {
    "pzi": {"haps": ("A1", "A2", "B1", "B2"), "pairs": (("A1", "B1"), ("A2", "B2"))},
    "psp": {"haps": ("A1", "A2", "B1", "B2"), "pairs": (("A1", "B1"), ("A2", "B2"))},
    "ppa": {"haps": ("A", "B"), "pairs": (("A", "B"),)},
}

# Explicit chromosome order requested for the ribbon layout.  PZI/PSP expand
# each entry to its available phased copies (A1/A2 or B1/B2); PPA has one A/B
# copy per entry.
CHROMOSOME_ORDER = (
    "C1A", "C6B", "C7B", "C2A", "C5B", "C3A", "C8B", "C7A",
    "C4A", "C1B", "C5A", "C2B", "C6A", "C3B", "C8A", "C4B",
)

# Eight distinguishable colours used consistently for C1-C8 homologous pairs.
# The colour is assigned from the A-side chromosome, including rearranged blocks.
PAIR_COLORS = {
    1: "0,114,178",
    2: "230,159,0",
    3: "0,158,115",
    4: "213,94,0",
    5: "204,121,167",
    6: "86,180,233",
    7: "190,168,0",
    8: "80,80,80",
}

# CircOS uses ``p`` as a pixel unit.  Render at 600 dpi so a 45-mm square has
# a diameter of 45 / 25.4 * 600 = 1062.99 px.  The post-processor preserves
# the 45-mm physical size while retaining this high-resolution source.
CIRCOS_RADIUS_PX = 531.50
IDEOGRAM_THICKNESS_PX = 7.06
STROKE_THICKNESS_PX = 6.25
TICK_THICKNESS_PX = 6.25
TICK_SIZE_PX = 12.50
IDEOGRAM_LABEL_SIZE_P = 31.25
TICK_LABEL_SIZE_P = 25.00
TICK_LABEL_OFFSET_P = 12.50
GC_LINE_THICKNESS_PX = 12.50
LINK_THICKNESS_PX = 6.25
HAP_COLORS = {
    "A1": "55,95,145", "A2": "120,168,205",
    "B1": "186,85,40", "B2": "222,151,55",
    "A": "55,95,145", "B": "186,85,40",
}
TRACK_FILES = (
    "gene_density.txt",
    "te_coverage.txt",
    "gypsy_coverage.txt",
    "copia_coverage.txt",
    "tir_coverage.txt",
    "helitron_coverage.txt",
    "gc_fraction.txt",
)
TRACK_COLORS = {
    "gene_density": "85,90,95",
    "te": "112,135,52",
    "gypsy": "196,91,82",
    "copia": "76,120,168",
    "tir": "126,88,160",
    "helitron": "40,150,145",
    "gc": "42,88,108",
}


def fasta_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    name: str | None = None
    length = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    if length <= 0:
                        raise ValueError(f"Empty FASTA sequence: {name} in {path}")
                    if name in sizes:
                        raise ValueError(f"Duplicate FASTA sequence: {name} in {path}")
                    sizes[name] = length
                name = raw[1:].split()[0]
                length = 0
            else:
                length += len(raw.strip())
    if name is not None:
        if length <= 0:
            raise ValueError(f"Empty FASTA sequence: {name} in {path}")
        if name in sizes:
            raise ValueError(f"Duplicate FASTA sequence: {name} in {path}")
        sizes[name] = length
    if not sizes:
        raise ValueError(f"No FASTA sequences found: {path}")
    return sizes


def parse_bed(path: Path) -> dict[str, tuple[str, int, int]]:
    result: dict[str, tuple[str, int, int]] = {}
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_number, raw in enumerate(handle, 1):
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.rstrip("\r\n").split()
            if len(fields) < 4:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError as error:
                raise ValueError(f"Invalid BED coordinates at {path}:{line_number}") from error
            if end <= start:
                continue
            identifier = fields[3]
            if identifier in result:
                raise ValueError(f"Duplicate BED gene ID {identifier} in {path}")
            result[identifier] = (fields[0], max(0, start), end)
    if not result:
        raise ValueError(f"Empty BED: {path}")
    return result


def locate_simple_blocks(jcvi_root: Path, species: str, left: str, right: str) -> Path:
    pair_dir = jcvi_root / f"{species}_hap{left}_vs_{species}_hap{right}"
    if not pair_dir.is_dir():
        raise FileNotFoundError(f"Missing JCVI pair directory: {pair_dir}")
    candidates = sorted(pair_dir.glob("*.anchors.simple"))
    if not candidates:
        candidates = sorted(pair_dir.glob("blocks.simple"))
    if not candidates:
        raise FileNotFoundError(
            f"Missing JCVI A-B simple blocks for {species} {left}-{right}: {pair_dir}"
        )
    if len(candidates) > 1:
        raise ValueError(f"Multiple candidate simple block files in {pair_dir}: {candidates}")
    return candidates[0]


def parse_simple_blocks(
    path: Path,
    left_bed: dict[str, tuple[str, int, int]],
    right_bed: dict[str, tuple[str, int, int]],
) -> tuple[list[tuple[str, int, int, str, int, int, str]], int, int]:
    """Return valid links and counts of malformed/cross-chromosome rows."""
    links: list[tuple[str, int, int, str, int, int, str]] = []
    malformed = 0
    cross_chromosome = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if not raw.strip() or raw.startswith("#"):
                continue
            fields = raw.split()
            if len(fields) != 6:
                malformed += 1
                continue
            left_1, left_2, right_1, right_2, _score, orientation = fields
            if any(gene not in left_bed for gene in (left_1, left_2)) or any(
                gene not in right_bed for gene in (right_1, right_2)
            ):
                malformed += 1
                continue
            l1, l2 = left_bed[left_1], left_bed[left_2]
            r1, r2 = right_bed[right_1], right_bed[right_2]
            if l1[0] != l2[0] or r1[0] != r2[0]:
                cross_chromosome += 1
                continue
            orientation = "-" if orientation == "-" else "+"
            links.append((
                l1[0], min(l1[1], l2[1]), max(l1[2], l2[2]),
                r1[0], min(r1[1], r2[1]), max(r1[2], r2[2]), orientation,
            ))
    return links, malformed, cross_chromosome


def chrom_number(seqid: str) -> tuple[int, str]:
    match = re.search(r"C0*(\d+)[AB]$", seqid)
    return (int(match.group(1)), seqid) if match else (10**9, seqid)


def chromosome_group(seqid: str) -> str | None:
    match = re.search(r"C0*(\d+)([AB])$", seqid)
    if not match:
        return None
    return f"C{int(match.group(1))}{match.group(2)}"


def internal_id(species: str, hap: str, chrom: str) -> str:
    """Return the per-species CircOS ID used by the established track files."""
    # Each species has its own configuration directory, so the haplotype and
    # chromosome are sufficient to keep IDs unique.  This also preserves the
    # IDs used by the existing 1-Mb track files (for example A1_C1A).
    return f"{hap}_{chrom}"


def write_karyotype(
    path: Path, species: str, haps: tuple[str, ...],
    sizes: dict[str, dict[str, int]],
) -> list[tuple[str, str, str, str, int]]:
    by_group = {
        (hap, chromosome_group(chrom)): chrom
        for hap in haps
        for chrom in sizes[hap]
        if chromosome_group(chrom) is not None
    }
    entries: list[tuple[int, str, str]] = []
    used: set[tuple[str, str]] = set()
    for group in CHROMOSOME_ORDER:
        side = group[-1]
        number = int(group[1:-1])
        for hap in haps:
            if not hap.startswith(side):
                continue
            chrom = by_group.get((hap, group))
            if chrom is None:
                raise ValueError(f"Missing {group} chromosome for {species} {hap}")
            entries.append((number, hap, chrom))
            used.add((hap, group))
    expected = sum(len(sizes[hap]) for hap in haps)
    if len(entries) != expected:
        missing = [
            f"{hap}:{group}"
            for hap in haps
            for group in sorted({chromosome_group(chrom) for chrom in sizes[hap] if chromosome_group(chrom)})
            if (hap, group) not in used
        ]
        raise ValueError(
            f"Chromosome order does not cover all {species} records; missing {missing[:6]}"
        )
    rows: list[tuple[str, str, str, str, int]] = []
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for number, hap, chrom in entries:
            internal = internal_id(species, hap, chrom)
            chrom_match = re.search(r"C0*(\d+)([AB])$", chrom)
            if chrom_match:
                side = chrom_match.group(2)
                phase = "" if hap in {"A", "B"} else hap[-1]
                label = f"C{int(chrom_match.group(1))}{side}{phase}"
            else:
                label = f"{chrom}_{hap}"
            length = sizes[hap][chrom]
            # Ideograms sharing the same homologous C number use one colour,
            # independent of phase (A1/A2/B1/B2) or subgenome side.
            color = f"pair{number}" if number in PAIR_COLORS else f"hap{hap}"
            handle.write(f"chr - {internal} {label} 0 {length} {color}\n")
            rows.append((internal, species, hap, chrom, length))
    return rows


def write_config(
    path: Path, species: str, haps: tuple[str, ...],
    order: Iterable[tuple[str, str, str, str, int]],
    track_maxima: dict[str, float],
    gc_range: tuple[float, float],
) -> None:
    order_list = list(order)
    chromosome_text = ";".join(row[0] for row in order_list)
    color_lines = [f"hap{hap} = {HAP_COLORS[hap]}" for hap in haps]
    color_lines.extend(f"pair{i} = {rgb}" for i, rgb in PAIR_COLORS.items())
    color_lines.extend((
        "gene_track = 85,90,95",
        "te_track = 112,135,52",
        "gypsy_track = 196,91,82",
        "copia_track = 76,120,168",
        "tir_track = 126,88,160",
        "helitron_track = 40,150,145",
        "gc_track = 42,88,108",
    ))
    text = f"""<<include etc/colors_fonts_patterns.conf>>
<<include etc/housekeeping.conf>>
<colors>
{chr(10).join(color_lines)}
</colors>
karyotype = data/karyotype.txt
chromosomes_units = 1000000
chromosomes_display_default = no
chromosomes = {chromosome_text}
<image>
background = white
dir = figures
file = {species}_ab_circos.png
png = yes
svg = yes
radius = {CIRCOS_RADIUS_PX:.2f}p
angle_offset = -90
auto_alpha_colors = yes
auto_alpha_steps = 100
24bit = yes
</image>
<ideogram>
<spacing>
default = 0.004r
</spacing>
radius = 0.900r
thickness = {IDEOGRAM_THICKNESS_PX:.2f}p
fill = yes
stroke_color = white
stroke_thickness = {STROKE_THICKNESS_PX:.2f}p
show_label = yes
label_radius = 1.045r
label_size = {IDEOGRAM_LABEL_SIZE_P:.1f}p
label_parallel = yes
</ideogram>
show_ticks = yes
show_tick_labels = yes
<ticks>
radius = 1r
color = grey
thickness = {TICK_THICKNESS_PX:.2f}p
multiplier = 1e-6
format = %d
<tick>
spacing = 50u
size = {TICK_SIZE_PX:.1f}p
show_label = yes
label_size = {TICK_LABEL_SIZE_P:.1f}p
label_offset = {TICK_LABEL_OFFSET_P:.2f}p
</tick>
</ticks>
<plots>
<plot>
type = heatmap
file = data/gene_density.txt
r0 = 0.932r
r1 = 0.972r
min = 0
max = {track_maxima['gene']:.6f}
color = spectral-11-div
color_mapping = 0
stroke_thickness = 0
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
max = {track_maxima['gypsy']:.6f}
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
max = {track_maxima['copia']:.6f}
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
max = {track_maxima['tir']:.6f}
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
max = {track_maxima['helitron']:.6f}
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
thickness = {GC_LINE_THICKNESS_PX:.2f}p
</plot>
</plots>
<links>
<link>
file = data/links_ab.txt
radius = 0.520r
bezier_radius = 0.040r
ribbon = yes
thickness = {LINK_THICKNESS_PX:.2f}p
crest = 0.45
</link>
</links>
"""
    path.write_text(text, encoding="ascii")


def stage_tracks(
    track_root: Path,
    species: str,
    data: Path,
    order: list[tuple[str, str, str, str, int]],
) -> tuple[dict[str, float], tuple[float, float], list[tuple[str, str, int, int]]]:
    """Copy and validate the established 1-Mb tracks for one species.

    The track data were generated from the same repartitioned haplotype
    annotations and use ``<hap>_<chromosome>`` IDs.  Keeping them as ordinary
    CircOS text files makes the central tracks reproducible and editable.
    Coordinates are clipped only when the current karyotype is shorter than a
    source interval; no values are rescaled.
    """
    source_dir = track_root / species / "data"
    expected = {row[0]: row[4] for row in order}
    maxima = {"gene": 1.0, "gypsy": 0.01, "copia": 0.01, "tir": 0.01, "helitron": 0.01}
    gc_values: list[float] = []
    audit: list[tuple[str, str, int, int]] = []
    stem_to_key = {
        "gene_density.txt": "gene",
        "gypsy_coverage.txt": "gypsy",
        "copia_coverage.txt": "copia",
        "tir_coverage.txt": "tir",
        "helitron_coverage.txt": "helitron",
    }
    for filename in TRACK_FILES:
        source = source_dir / filename
        if not source.is_file() or source.stat().st_size == 0:
            raise FileNotFoundError(f"Missing central track for {species}: {source}")
        target = data / filename
        seen: set[str] = set()
        kept = 0
        clipped = 0
        with source.open("r", encoding="utf-8", errors="replace") as src, target.open(
            "w", encoding="ascii", newline="\n"
        ) as dst:
            for line_number, raw in enumerate(src, 1):
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.split()
                if len(fields) < 4:
                    raise ValueError(f"Malformed track row at {source}:{line_number}")
                chrom = fields[0]
                if chrom not in expected:
                    raise ValueError(
                        f"Track chromosome {chrom} in {source} is absent from the 0903 karyotype"
                    )
                try:
                    start, end = int(fields[1]), int(fields[2])
                    value = float(fields[3])
                except ValueError as error:
                    raise ValueError(f"Invalid track row at {source}:{line_number}") from error
                if start < 0 or end <= start:
                    continue
                limit = expected[chrom]
                if start >= limit:
                    continue
                if end > limit:
                    end = limit
                    clipped += 1
                if end <= start:
                    continue
                dst.write(f"{chrom} {start} {end} {value:.6f}\n")
                seen.add(chrom)
                kept += 1
                if filename in stem_to_key:
                    maxima[stem_to_key[filename]] = max(maxima[stem_to_key[filename]], value)
                elif filename == "gc_fraction.txt":
                    gc_values.append(value)
        missing = set(expected) - seen
        if missing:
            raise ValueError(
                f"Track {source} does not cover {len(missing)} karyotype chromosomes; "
                f"examples: {', '.join(sorted(missing)[:4])}"
            )
        audit.append((filename, str(source), kept, clipped))
    # The total-TE track is a fraction and therefore has a fixed unit scale.
    maxima["te"] = 1.0
    if gc_values:
        low, high = min(gc_values), max(gc_values)
        pad = max(0.005, (high - low) * 0.08)
        gc_range = (max(0.0, low - pad), min(1.0, high + pad))
    else:
        gc_range = (0.0, 1.0)
    with (data / "track_manifest.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("track", "source", "rows_written", "rows_clipped"))
        writer.writerows(audit)
    with (data / "track_layout.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("track", "r0", "r1", "min", "max"))
        writer.writerow(("gene_density", "0.932r", "0.972r", "0", f"{maxima['gene']:.6f}"))
        writer.writerow(("all_TE", "0.875r", "0.915r", "0", "1"))
        writer.writerow(("Gypsy", "0.818r", "0.858r", "0", f"{maxima['gypsy']:.6f}"))
        writer.writerow(("Copia", "0.761r", "0.801r", "0", f"{maxima['copia']:.6f}"))
        writer.writerow(("TIR_DNA", "0.704r", "0.744r", "0", f"{maxima['tir']:.6f}"))
        writer.writerow(("Helitron", "0.647r", "0.687r", "0", f"{maxima['helitron']:.6f}"))
        writer.writerow(("GC", "0.590r", "0.630r", f"{gc_range[0]:.6f}", f"{gc_range[1]:.6f}"))
    return maxima, gc_range, audit


def prepare_species(
    root: Path, out_root: Path, jcvi_root: Path,
    bed_root: Path, track_root: Path, species: str,
) -> dict[str, int]:
    spec = SPECIES[species]
    haps: tuple[str, ...] = spec["haps"]
    species_root = out_root / species
    data = species_root / "data"
    conf = species_root / "conf"
    figures = species_root / "figures"
    for directory in (data, conf, figures):
        directory.mkdir(parents=True, exist_ok=True)

    sizes: dict[str, dict[str, int]] = {}
    fasta_lengths: dict[str, dict[str, int]] = {}
    bed_maxima: dict[str, dict[str, int]] = {}
    coordinate_warnings: list[tuple[str, str, int, int, int, str]] = []
    beds: dict[str, dict[str, tuple[str, int, int]]] = {}
    input_rows: list[tuple[str, str, str]] = []
    for hap in haps:
        genome = root / "haploid_repartition" / "inputs" / species / "genome" / f"{species}_hap{hap}.fa"
        bed = bed_root / f"{species}_hap{hap}" / f"{species}_hap{hap}.bed"
        if not genome.is_file() or genome.stat().st_size == 0:
            raise FileNotFoundError(f"Missing genome FASTA: {genome}")
        if not bed.is_file() or bed.stat().st_size == 0:
            raise FileNotFoundError(f"Missing JCVI BED: {bed}")
        sizes[hap] = fasta_sizes(genome)
        fasta_lengths[hap] = dict(sizes[hap])
        beds[hap] = parse_bed(bed)
        unknown_bed_chroms = sorted({row[0] for row in beds[hap].values()} - set(sizes[hap]))
        if unknown_bed_chroms:
            raise ValueError(
                f"BED sequence IDs absent from {genome}: {','.join(unknown_bed_chroms[:5])}"
            )
        bed_maxima[hap] = {}
        for chrom, _start, end in beds[hap].values():
            bed_maxima[hap][chrom] = max(bed_maxima[hap].get(chrom, 0), end)
        for chrom, fasta_length in fasta_lengths[hap].items():
            bed_max = bed_maxima[hap].get(chrom, 0)
            used_length = max(fasta_length, bed_max)
            status = "pass" if bed_max <= fasta_length else "bed_exceeds_fasta"
            coordinate_warnings.append((hap, chrom, fasta_length, bed_max, used_length, status))
            sizes[hap][chrom] = used_length
        input_rows.append((hap, str(genome), str(bed)))

    order = write_karyotype(data / "karyotype.txt", species, haps, sizes)
    with (data / "chromosome_map.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "internal_id", "species", "haplotype", "source_chromosome",
            "fasta_length_bp", "bed_max_bp", "used_length_bp", "coordinate_status",
        ))
        for internal, row_species, hap, chrom, used_length in order:
            fasta_length = fasta_lengths[hap][chrom]
            bed_max = bed_maxima[hap].get(chrom, 0)
            status = "pass" if bed_max <= fasta_length else "bed_exceeds_fasta"
            writer.writerow((
                internal, row_species, hap, chrom, fasta_length, bed_max,
                used_length, status,
            ))
    with (data / "coordinate_audit.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "haplotype", "chromosome", "fasta_length_bp", "bed_max_bp",
            "used_length_bp", "status",
        ))
        writer.writerows(coordinate_warnings)
    with (data / "input_manifest.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("haplotype", "genome", "bed"))
        writer.writerows(input_rows)

    link_rows: list[tuple[str, int, int, str, int, int, str, str]] = []
    summary_rows: list[tuple[str, str, int, int, int, int]] = []
    for pair_index, (left, right) in enumerate(spec["pairs"], 1):
        block_file = locate_simple_blocks(jcvi_root, species, left, right)
        parsed, malformed, cross_chromosome = parse_simple_blocks(
            block_file, beds[left], beds[right]
        )
        for lchrom, lstart, lend, rchrom, rstart, rend, orientation in parsed:
            match = re.search(r"C0*(\d+)A$", lchrom)
            color_index = int(match.group(1)) if match else 1
            color_index = color_index if color_index in PAIR_COLORS else 1
            link_rows.append((
                internal_id(species, left, lchrom), lstart, lend,
                internal_id(species, right, rchrom), rstart, rend,
                f"pair{color_index}_a60", orientation,
            ))
        summary_rows.append((
            f"{left}-{right}", str(block_file), len(parsed),
            malformed, cross_chromosome, pair_index,
        ))

    with (data / "links_ab.txt").open("w", encoding="ascii", newline="\n") as handle:
        for lchrom, lstart, lend, rchrom, rstart, rend, color, orientation in link_rows:
            twist = "yes" if orientation == "-" else "no"
            handle.write(
                f"{lchrom} {lstart} {lend} {rchrom} {rstart} {rend} "
                f"color={color},z=100,twist={twist}\n"
            )
    with (data / "link_summary.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "pair", "block_file", "valid_links", "malformed_rows",
            "cross_chromosome_blocks", "pair_index",
        ))
        writer.writerows(summary_rows)
        writer.writerow(("TOTAL", "", len(link_rows), "", "", ""))
    track_maxima, gc_range, track_audit = stage_tracks(track_root, species, data, order)
    write_config(conf / "circos.conf", species, haps, order, track_maxima, gc_range)
    # Source paths may contain non-ASCII characters on the local workstation.
    with (data / "track_provenance.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("species", "track", "source_directory", "note"))
        for filename, source, _rows, _clipped in track_audit:
            writer.writerow((species, filename, str(Path(source).parent), "validated against current karyotype"))
    (data / "README_methods.txt").write_text(
        "A-B-only CircOS inputs for the 0903 repartition. PZI/PSP use A1-B1 "
        "and A2-B2 JCVI simple blocks; PPA uses A-B. No A1-A2 or B1-B2 "
        "links are included. Link colours are assigned by the A-side C1-C8 "
        "chromosome; a60 is the CircOS alpha suffix. The central tracks are "
        "1-Mb gene density, total TE, Gypsy, Copia, TIR/DNA, Helitron and GC; "
        "track_manifest.tsv records their source and coordinate checks. "
        "Chromosomes follow the explicit cross-subgenome order in the preparation "
        "script; each C1-C8 homologous group shares its pairN ideogram colour. "
        "Gene density is rendered as a 1-Mb heatmap using the native CircOS "
        "spectral colour scale; total TE, Gypsy, Copia, TIR/DNA and Helitron "
        "remain histogram tracks, and GC remains a line track. Coordinate "
        "ticks are placed every 50 Mb with labels in Mb; the A-B ribbons are drawn "
        "at 0.520r to keep them separate from the innermost track. "
        "Ideogram lengths use the maximum of FASTA length and JCVI BED extent; "
        "see coordinate_audit.tsv.\n",
        encoding="ascii",
    )
    return {"chromosomes": len(order), "links": len(link_rows)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--jcvi-root", type=Path, required=True)
    parser.add_argument("--bed-root", type=Path, required=True)
    parser.add_argument(
        "--track-root", type=Path,
        default=Path("downstream_analysis/circos/repartitioned"),
        help="Existing validated track root containing <species>/data/*.txt",
    )
    args = parser.parse_args()
    root = args.root.resolve()
    out_root = args.out_root.resolve()
    jcvi_root = args.jcvi_root.resolve()
    bed_root = args.bed_root.resolve()
    track_root = args.track_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    manifest_rows: list[tuple[str, int, int, str, str, str, str, str]] = []
    for species in SPECIES:
        stats = prepare_species(root, out_root, jcvi_root, bed_root, track_root, species)
        manifest_rows.append((
            species, stats["chromosomes"], stats["links"],
            f"{species}/figures/{species}_ab_circos.png",
            f"{species}/figures/{species}_ab_circos.svg",
            f"{species}/figures/{species}_ab_circos.pdf",
            "45",
            "50",
        ))
        print(f"[OK] {species}: {stats['chromosomes']} chromosomes, {stats['links']} A-B links")
    with (out_root / "figure_manifest.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow((
            "species", "chromosomes", "ab_links", "png", "svg", "pdf",
            "canvas_mm", "tick_spacing_mb",
        ))
        writer.writerows(manifest_rows)
    (out_root / "README_methods.txt").write_text(
        "0903 PZI/PSP/PPA CircOS figures. Existing JCVI anchors.simple files "
        "were converted to A-B ribbons only; within-haplotype links were omitted. "
        "Tick spacing is 50 Mb (0, 50, 100 Mb where chromosome length permits). "
        "The CircOS canvas is 45 mm x 45 mm (85.04-pixel radius at the native "
        "96-dpi unit). "
        "Each species includes seven central 1-Mb tracks: gene density, total TE, "
        "Gypsy, Copia, TIR/DNA, Helitron and GC. Track source and coverage checks "
        "are recorded in each species/data/track_manifest.tsv and "
        "track_provenance.tsv.\n",
        encoding="ascii",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
