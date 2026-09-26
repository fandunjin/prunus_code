#!/usr/bin/env python3
"""Summarize standardized RepeatMasker annotations and library-family maps."""

import argparse
import collections
import csv
from pathlib import Path


def write_tsv(path, header, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def read_fasta_headers(path):
    headers = []
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            raw = line[1:].strip()
            token = raw.split()[0] if raw else ""
            family = token.split("#", 1)[0]
            class_family = token.split("#", 1)[1] if "#" in token else "unknown"
            headers.append((token, family, class_family, raw))
    return headers


def map_library(args):
    original_headers = read_fasta_headers(args.original_lib)
    paf_best = {}
    with open(args.paf) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            matches = int(fields[9])
            aln_len = int(fields[10])
            mapq = int(fields[11])
            if aln_len == 0:
                continue
            target_token = fields[5]
            target_family = target_token.split("#", 1)[0]
            target_class = target_token.split("#", 1)[1] if "#" in target_token else "unknown"
            record = (matches, aln_len, mapq, target_family, target_class)
            previous = paf_best.get(query)
            if previous is None or record[:3] > previous[:3]:
                paf_best[query] = record

    rows = []
    for original_token, original_family, original_class, raw_header in original_headers:
        hit = paf_best.get(original_token)
        if hit is None:
            rows.append([args.label, original_family, original_class, "", "", 0, 0, 0, 0.0, "unmapped", raw_header])
            continue
        matches, aln_len, mapq, pan_family, pan_class = hit
        rows.append([
            args.label, original_family, original_class, pan_family, pan_class,
            matches, aln_len, mapq, "{:.6f}".format(matches / aln_len), "mapped", raw_header,
        ])
    write_tsv(
        args.output,
        [
            "target", "original_family", "original_class_family", "pan_family", "pan_class_family",
            "matching_bp", "alignment_bp", "mapq", "identity", "mapping_status", "original_raw_header",
        ],
        rows,
    )


def fasta_length(path):
    total = 0
    with open(path) as handle:
        for line in handle:
            if not line.startswith(">"):
                total += len(line.strip())
    return total


def parse_rmout(path):
    class_intervals = collections.defaultdict(list)
    family_intervals = collections.defaultdict(list)
    total_intervals = collections.defaultdict(list)
    counts_class = collections.Counter()
    counts_family = collections.Counter()

    with open(path) as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 11 or not fields[0].isdigit():
                continue
            seq_id = fields[4]
            start = int(fields[5]) - 1
            end = int(fields[6])
            family = fields[9]
            class_family = fields[10]
            if end <= start:
                continue
            class_name = class_family.split("/", 1)[0]
            total_intervals[seq_id].append((start, end))
            class_intervals[(seq_id, class_name)].append((start, end))
            family_key = class_family + "|" + family
            family_intervals[(seq_id, family_key)].append((start, end))
            counts_class[class_name] += 1
            counts_family[(class_family, family)] += 1

    def union_length(intervals):
        if not intervals:
            return 0
        total = 0
        start, end = sorted(intervals)[0]
        for next_start, next_end in sorted(intervals)[1:]:
            if next_start > end:
                total += end - start
                start, end = next_start, next_end
            else:
                end = max(end, next_end)
        return total + end - start

    by_class = collections.defaultdict(int)
    by_family = collections.defaultdict(int)
    total_union = sum(union_length(intervals) for intervals in total_intervals.values())
    for (_, class_name), intervals in class_intervals.items():
        by_class[class_name] += union_length(intervals)
    for (_, family_key), intervals in family_intervals.items():
        by_family[family_key] += union_length(intervals)
    return by_class, by_family, counts_class, counts_family, total_union


def read_targets(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def summarize(args):
    output_dir = Path(args.out_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    class_rows = []
    family_rows = []
    total_rows = []
    for target in read_targets(args.targets_tsv):
        genome_bp = fasta_length(target["mask_input_fasta"])
        by_class, by_family, counts_class, counts_family, total_union = parse_rmout(target["repeatmasker_out"])
        total_rows.append([
            target["target"], target["group"], genome_bp, total_union,
            "{:.6f}".format(total_union / genome_bp if genome_bp else 0), target["repeatmasker_out"],
        ])
        for class_name in sorted(by_class):
            bp = by_class[class_name]
            class_rows.append([
                target["target"], target["group"], class_name, counts_class[class_name], bp,
                "{:.6f}".format(bp / genome_bp if genome_bp else 0),
            ])
        for combined in sorted(by_family):
            class_family, family = combined.split("|", 1)
            bp = by_family[combined]
            family_rows.append([
                target["target"], target["group"], class_family, family,
                counts_family[(class_family, family)], bp,
                "{:.6f}".format(bp / genome_bp if genome_bp else 0),
            ])
    write_tsv(
        output_dir / "prunus_pan_te_total_summary.tsv",
        ["target", "group", "genome_bp", "repeat_union_bp", "repeat_union_fraction", "repeatmasker_out"],
        total_rows,
    )
    write_tsv(
        output_dir / "prunus_pan_te_class_summary.tsv",
        ["target", "group", "repeat_class", "annotation_count", "union_bp", "genome_fraction"],
        class_rows,
    )
    write_tsv(
        output_dir / "prunus_pan_te_family_summary.tsv",
        ["target", "group", "class_family", "pan_family", "annotation_count", "union_bp", "genome_fraction"],
        family_rows,
    )


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    mapping = subparsers.add_parser("map-library")
    mapping.add_argument("--label", required=True)
    mapping.add_argument("--original-lib", required=True)
    mapping.add_argument("--paf", required=True)
    mapping.add_argument("--output", required=True)
    mapping.set_defaults(func=map_library)
    summary = subparsers.add_parser("summarize")
    summary.add_argument("--targets-tsv", required=True)
    summary.add_argument("--out-dir", required=True)
    summary.set_defaults(func=summarize)
    return parser.parse_args()


def main():
    args = parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
