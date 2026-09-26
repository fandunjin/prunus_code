#!/usr/bin/env python3
"""Add haplotype-resolved ncRNA density tracks to the 0903 CircOS inputs."""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


BIN_SIZE = 1_000_000
SPECIES_HAPS = {
    "pzi": ("A1", "A2", "B1", "B2"),
    "psp": ("A1", "A2", "B1", "B2"),
    "ppa": ("A", "B"),
}
NCRNA_COLOR = "166,54,95"


def parse_karyotype(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        fields = raw.split()
        if not fields or fields[0] != "chr":
            continue
        if len(fields) < 6:
            raise ValueError(f"Malformed karyotype row at {path}:{line_number}")
        lengths[fields[2]] = int(fields[5])
    if not lengths:
        raise ValueError(f"No chromosomes found in {path}")
    return lengths


def parse_attrs(text: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value.replace("%3B", ";").replace("%2C", ",")
    return attrs


def read_gff(path: Path, hap: str, lengths: dict[str, int]):
    records: list[tuple[str, int, int, str]] = []
    audit: list[tuple[str, int, str, str, int, int, int, int, str]] = []
    for line_number, raw in enumerate(path.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
        if not raw.strip() or raw.startswith("#"):
            continue
        fields = raw.rstrip("\r\n").split("\t")
        if len(fields) != 9:
            audit.append((str(path), line_number, "", "", 0, 0, 0, 0, "skip_invalid_columns"))
            continue
        source_seqid = fields[0]
        target = f"{hap}_{source_seqid}"
        if target not in lengths:
            audit.append((str(path), line_number, source_seqid, target, 0, 0, 0, 0, "skip_unknown_seqid"))
            continue
        try:
            old_start, old_end = int(fields[3]), int(fields[4])
        except ValueError:
            audit.append((str(path), line_number, source_seqid, target, 0, 0, 0, 0, "skip_invalid_coordinates"))
            continue
        start, end = sorted((old_start, old_end))
        limit = lengths[target]
        if end < 1 or start > limit:
            audit.append((str(path), line_number, source_seqid, target, old_start, old_end, 0, 0, "skip_outside_karyotype"))
            continue
        start = max(1, start)
        end = min(limit, end)
        attrs = parse_attrs(fields[8])
        nc_class = attrs.get("ncRNA_class", "other ncRNA")
        records.append((target, start, end, nc_class))
        status = "pass" if (start, end) == (old_start, old_end) else "clipped_to_karyotype"
        audit.append((str(path), line_number, source_seqid, target, old_start, old_end, start, end, status))
    return records, audit


def write_tracks(
    records: list[tuple[str, int, int, str]],
    lengths: dict[str, int],
    data_dir: Path,
    audit: list[tuple[str, int, str, str, int, int, int, int, str]],
    inputs: list[tuple[str, str]],
) -> int:
    counts: defaultdict[str, Counter[int]] = defaultdict(Counter)
    class_counts: Counter[str] = Counter()
    class_bp: Counter[str] = Counter()
    for target, start, end, nc_class in records:
        midpoint = (start + end) // 2
        bin_start = ((midpoint - 1) // BIN_SIZE) * BIN_SIZE
        counts[target][bin_start] += 1
        class_counts[nc_class] += 1
        class_bp[nc_class] += end - start + 1

    density = data_dir / "ncrna_density.txt"
    with density.open("w", encoding="ascii", newline="\n") as handle:
        for target, length in lengths.items():
            for start in range(0, length, BIN_SIZE):
                end = min(start + BIN_SIZE, length)
                handle.write(f"{target} {start} {end} {counts[target][start]}\n")

    with (data_dir / "ncrna_category_summary.tsv").open("w", encoding="ascii", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("class", "count", "total_bp"))
        for nc_class in sorted(class_counts):
            writer.writerow((nc_class, class_counts[nc_class], class_bp[nc_class]))

    with (data_dir / "ncrna_coordinate_audit.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gff3", "gff_line", "source_seqid", "circos_seqid", "original_start", "original_end", "used_start", "used_end", "status"))
        writer.writerows(audit)

    with (data_dir / "ncrna_track_manifest.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        for hap, path in inputs:
            writer.writerow((f"input_gff3_{hap}", path))
        writer.writerow(("karyotype", str(data_dir / "karyotype.txt")))
        writer.writerow(("bin_size_bp", BIN_SIZE))
        writer.writerow(("mapped_records", len(records)))
        writer.writerow(("skipped_records", sum(row[-1].startswith("skip_") for row in audit)))
        writer.writerow(("max_1Mb_bin_count", max(counts[target][start] for target, length in lengths.items() for start in range(0, length, BIN_SIZE))))
    return max(counts[target][start] for target, length in lengths.items() for start in range(0, length, BIN_SIZE))


def add_config(config: Path, output: Path, max_count: int, output_stem: str) -> None:
    text = config.read_text(encoding="utf-8")
    if "file = data/ncrna_density.txt" in text:
        raise ValueError(f"Config already contains an ncRNA plot: {config}")
    if "ncrna_track =" not in text:
        text = text.replace("gc_track = 42,88,108", f"gc_track = 42,88,108\nncrna_track = {NCRNA_COLOR}")
    plot = (
        "<plot>\n"
        "type = histogram\n"
        "file = data/ncrna_density.txt\n"
        "r0 = 0.535r\n"
        "r1 = 0.575r\n"
        "min = 0\n"
        f"max = {max(1, max_count)}\n"
        "fill_color = ncrna_track\n"
        "color = ncrna_track\n"
        "thickness = 0\n"
        "</plot>\n"
    )
    text = text.replace("</plots>", plot + "</plots>", 1)
    text = re.sub(r"file = [^\n]+", f"file = {output_stem}.png", text, count=1)
    output.write_text(text, encoding="ascii")


def prepare_species(root: Path, out_root: Path, ncrna_root: Path, species: str) -> None:
    species_root = out_root / species
    data_dir = species_root / "data"
    config = species_root / "conf" / "circos.conf"
    karyotype = data_dir / "karyotype.txt"
    lengths = parse_karyotype(karyotype)
    records: list[tuple[str, int, int, str]] = []
    audit: list[tuple[str, int, str, str, int, int, int, int, str]] = []
    inputs: list[tuple[str, str]] = []
    for hap in SPECIES_HAPS[species]:
        gff = ncrna_root / species / f"hap{hap}" / "final" / f"{species}_hap{hap}.ncRNA.final.gff3"
        if not gff.is_file() or gff.stat().st_size == 0:
            raise FileNotFoundError(f"Missing ncRNA GFF3: {gff}")
        current, current_audit = read_gff(gff, hap, lengths)
        records.extend(current)
        audit.extend(current_audit)
        inputs.append((hap, str(gff)))
    max_count = write_tracks(records, lengths, data_dir, audit, inputs)
    output_stem = f"{species}_ab_circos_ncrna"
    add_config(config, species_root / "conf" / "circos_ncrna.conf", max_count, output_stem)
    print(f"[OK] {species}: {len(records)} ncRNA records, max bin {max_count}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--ncrna-root", type=Path, required=True)
    parser.add_argument("--species", nargs="*", choices=tuple(SPECIES_HAPS), default=tuple(SPECIES_HAPS))
    args = parser.parse_args()
    for species in args.species:
        prepare_species(args.root.resolve(), args.out_root.resolve(), args.ncrna_root.resolve(), species)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
