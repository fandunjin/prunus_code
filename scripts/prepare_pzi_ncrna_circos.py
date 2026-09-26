#!/usr/bin/env python3
"""Prepare a PZI ncRNA density track for the existing A-B CircOS figure.

The ncRNA GFF3 is generated on the complete PZI genome (chr01-chr32), while
the CircOS karyotype uses A1_C1A ... B2_C8B identifiers.  The repartition map
is therefore used to transfer every locus to its unique haplotype chromosome
before 1-Mb binning.  The original CircOS files are never overwritten.
"""

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


BIN_SIZE = 1_000_000
CLASS_FILES = {
    "tRNA": "ncrna_trna_density.txt",
    "tRNA pseudogene": "ncrna_trna_pseudogene_density.txt",
    "rRNA": "ncrna_rrna_density.txt",
    "snRNA": "ncrna_snRNA_density.txt",
    "snoRNA": "ncrna_snoRNA_density.txt",
    "miRNA": "ncrna_miRNA_density.txt",
    "SRP RNA": "ncrna_srp_density.txt",
    "other ncRNA": "ncrna_other_density.txt",
}


def parse_karyotype(path):
    lengths = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        fields = line.split()
        if not fields or fields[0] != "chr" or len(fields) < 6:
            continue
        lengths[fields[2]] = int(fields[5])
    if not lengths:
        raise ValueError("No chromosome records found in karyotype: " + str(path))
    return lengths


def parse_map(path):
    mapping = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            hap = row["haplotype"]
            phase = row.get("phase", "1")
            source = row.get("copy1_source_chr32") if phase == "1" else row.get("copy2_source_chr32")
            if not source:
                source = row.get("original_chr32")
            if not source:
                raise ValueError("Map row has no source chromosome: " + repr(row))
            internal = hap.replace("hap", "", 1) + "_" + row["new_seqid"]
            if source in mapping:
                raise ValueError(f"Source chromosome maps more than once: {source}")
            mapping[source] = internal
    if len(mapping) != 32:
        raise ValueError(f"Expected 32 unique PZI source chromosomes, found {len(mapping)}")
    return mapping


def parse_attrs(text):
    result = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value.replace("%3B", ";").replace("%2C", ",")
    return result


def parse_gff(path, seq_map, lengths):
    records = []
    audit = []
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line_number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                audit.append((line_number, fields[0] if fields else "", "", "", "", "", "", "skip_invalid_columns"))
                continue
            source_seqid = fields[0]
            target = seq_map.get(source_seqid)
            if target is None:
                audit.append((line_number, source_seqid, "", "", "", "", "", "skip_unmapped_seqid"))
                continue
            try:
                start, end = sorted((int(fields[3]), int(fields[4])))
            except ValueError:
                audit.append((line_number, source_seqid, target, fields[3], fields[4], "", "", "skip_invalid_coordinates"))
                continue
            old_start, old_end = start, end
            limit = lengths[target]
            if end < 1 or start > limit:
                audit.append((line_number, source_seqid, target, old_start, old_end, "", "", "skip_outside_karyotype"))
                continue
            start = max(1, start)
            end = min(limit, end)
            attrs = parse_attrs(fields[8])
            nc_class = attrs.get("ncRNA_class", "other ncRNA")
            if nc_class not in CLASS_FILES:
                nc_class = "other ncRNA"
            records.append((target, start, end, nc_class, attrs.get("ID", f"line{line_number}")))
            status = "pass" if (start, end) == (old_start, old_end) else "clipped_to_karyotype"
            audit.append((line_number, source_seqid, target, old_start, old_end, start, end, status))
    if not records:
        raise ValueError("No mapped ncRNA records remain after coordinate validation")
    return records, audit


def write_tracks(records, lengths, out_dir):
    counts = {name: defaultdict(Counter) for name in ["total"] + list(CLASS_FILES)}
    class_total = Counter()
    class_bp = Counter()
    for target, start, end, nc_class, _identifier in records:
        midpoint = (start + end) // 2
        bin_start = ((midpoint - 1) // BIN_SIZE) * BIN_SIZE
        counts["total"][target][bin_start] += 1
        counts[nc_class][target][bin_start] += 1
        class_total[nc_class] += 1
        class_bp[nc_class] += end - start + 1

    def write_one(path, key):
        with path.open("w", encoding="ascii", newline="\n") as handle:
            for target in lengths:
                for start in range(0, lengths[target], BIN_SIZE):
                    end = min(start + BIN_SIZE, lengths[target])
                    handle.write(f"{target} {start} {end} {counts[key][target][start]}\n")

    total_path = out_dir / "ncrna_density.txt"
    write_one(total_path, "total")
    for nc_class, filename in CLASS_FILES.items():
        write_one(out_dir / filename, nc_class)

    with (out_dir / "ncrna_category_summary.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("class", "count", "total_bp"))
        for nc_class in CLASS_FILES:
            writer.writerow((nc_class, class_total[nc_class], class_bp[nc_class]))

    max_count = max(counts["total"][target][start] for target in lengths for start in range(0, lengths[target], BIN_SIZE))
    return total_path, max_count


def add_config(config, output_config, max_count, output_stem):
    text = config.read_text(encoding="utf-8")
    if "file = data/ncrna_density.txt" in text:
        raise ValueError("Config already contains an ncRNA track: " + str(config))
    text = text.replace("gc_track = 42,88,108", "gc_track = 42,88,108\nncrna_track = 166,54,95")
    plot = f'''<plot>\ntype = histogram\nfile = data/ncrna_density.txt\nr0 = 0.535r\nr1 = 0.575r\nmin = 0\nmax = {max(1, max_count)}\nfill_color = ncrna_track\ncolor = ncrna_track\nthickness = 0\n</plot>\n'''
    text = text.replace("</plots>", plot + "</plots>")
    text = re.sub(r"file = [^\n]+", f"file = {output_stem}.png", text, count=1)
    output_config.write_text(text, encoding="ascii")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff3", type=Path, required=True)
    parser.add_argument("--map", dest="map_path", type=Path, required=True)
    parser.add_argument("--karyotype", type=Path, required=True)
    parser.add_argument("--out-data", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-config", type=Path, required=True)
    parser.add_argument("--output-stem", default="pzi_ab_circos_ncrna")
    args = parser.parse_args()
    for path in (args.gff3, args.map_path, args.karyotype, args.config):
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(path)
    args.out_data.mkdir(parents=True, exist_ok=True)
    lengths = parse_karyotype(args.karyotype)
    seq_map = parse_map(args.map_path)
    records, audit = parse_gff(args.gff3, seq_map, lengths)
    total_path, max_count = write_tracks(records, lengths, args.out_data)
    with (args.out_data / "ncrna_coordinate_audit.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("gff_line", "source_seqid", "circos_seqid", "original_start", "original_end", "used_start", "used_end", "status"))
        writer.writerows(audit)
    add_config(args.config, args.output_config, max_count, args.output_stem)
    with (args.out_data / "ncrna_track_manifest.tsv").open("w", encoding="utf-8", newline="\n") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(("metric", "value"))
        writer.writerow(("input_gff3", str(args.gff3)))
        writer.writerow(("input_repartition_map", str(args.map_path)))
        writer.writerow(("karyotype", str(args.karyotype)))
        writer.writerow(("bin_size_bp", BIN_SIZE))
        writer.writerow(("source_gff_records", sum(1 for row in audit if row[0])))
        writer.writerow(("mapped_records", len(records)))
        writer.writerow(("skipped_records", sum(1 for row in audit if row[-1].startswith("skip_"))))
        for status, count in sorted(Counter(row[-1] for row in audit).items()):
            writer.writerow((f"audit_{status}", count))
        writer.writerow(("max_1Mb_bin_count", max_count))
        writer.writerow(("total_track", str(total_path)))
    print(f"mapped ncRNA records: {len(records)}")
    print(f"maximum 1-Mb ncRNA count: {max_count}")
    print(f"wrote: {total_path}")


if __name__ == "__main__":
    main()
