#!/usr/bin/env python3
"""Split copy1 chromosomes at manually confirmed junctions into neutral components.

Only review-sheet rows with final_status=confirmed_rejoined_junction are used.
The output is a coordinate and external-fragment ledger, not an A/B, LF/MF,
or parental-origin assignment.
"""

from __future__ import print_function

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--review-sheet", required=True)
    parser.add_argument("--fragments", required=True)
    parser.add_argument("--query-fai", required=True)
    parser.add_argument("--copy-map", required=True)
    parser.add_argument("--outdir", required=True)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fai(path):
    lengths = {}
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(fields[1])
    return lengths


def read_copy_map(path):
    mapping = {}
    for row in read_tsv(path):
        normalized = row.get("normalized_chr") or row.get("copy1_chr")
        copy2 = row.get("copy2_chr") or normalized
        if normalized:
            mapping["copy1_" + normalized] = "copy2_" + copy2
    if not mapping:
        raise SystemExit("No copy1/copy2 chromosome pairs found in {}".format(path))
    return mapping


def natural_key(value):
    output = []
    number = ""
    for char in value:
        if char.isdigit():
            number += char
        else:
            if number:
                output.append((0, int(number)))
                number = ""
            output.append((1, char))
    if number:
        output.append((0, int(number)))
    return output


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def resolve_fai_chrom(query_bed_chrom, lengths):
    if query_bed_chrom in lengths:
        return query_bed_chrom
    if query_bed_chrom.startswith("copy1_"):
        candidate = query_bed_chrom[len("copy1_"):]
        if candidate in lengths:
            return candidate
    raise SystemExit("No copy1 FASTA sequence corresponds to {}".format(query_bed_chrom))


def resolve_bed_chrom(fasta_chrom, known_bed_chroms):
    for candidate in (fasta_chrom, "copy1_" + fasta_chrom):
        if candidate in known_bed_chroms:
            return candidate
    return "copy1_" + fasta_chrom


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    lengths = read_fai(args.query_fai)
    expected_copy2 = read_copy_map(args.copy_map)
    confirmed = [row for row in read_tsv(args.review_sheet)
                 if row.get("final_status", "") == "confirmed_rejoined_junction"]
    cuts = defaultdict(list)
    for row in confirmed:
        chrom = row["copy1_chr"]
        fasta_chrom = resolve_fai_chrom(chrom, lengths)
        coordinate = int(float(row["boundary_mid_bp"]))
        if coordinate <= 0 or coordinate >= lengths[fasta_chrom]:
            raise SystemExit("Confirmed junction outside {} bounds: {}".format(chrom, coordinate))
        cuts[chrom].append((coordinate, row["junction_id"]))
    for chrom in cuts:
        cuts[chrom].sort()
        if len({coordinate for coordinate, _identifier in cuts[chrom]}) != len(cuts[chrom]):
            raise SystemExit("Two confirmed junctions share one coordinate on {}".format(chrom))

    fragments = [row for row in read_tsv(args.fragments) if row.get("accepted") == "1"]
    by_chrom = defaultdict(list)
    for row in fragments:
        by_chrom[row["q_chr"]].append(row)

    known_bed_chroms = set(by_chrom).union(cuts)
    fields = [
        "component_id", "copy1_bed_chr", "copy1_fasta_chr", "expected_copy2_bed_chr", "start", "end", "length_bp", "left_confirmed_junction",
        "right_confirmed_junction", "accepted_external_fragment_count", "external_fragment_profile",
        "assignment_status", "claim_boundary",
    ]
    rows = []
    component_index = 0
    for fasta_chrom in sorted(lengths, key=natural_key):
        chrom = resolve_bed_chrom(fasta_chrom, known_bed_chroms)
        points = [(0, "chromosome_start")] + cuts.get(chrom, []) + [(lengths[fasta_chrom], "chromosome_end")]
        for (start, left_id), (end, right_id) in zip(points, points[1:]):
            if end <= start:
                continue
            component_index += 1
            overlap = []
            for fragment in by_chrom.get(chrom, []):
                frag_start = int(float(fragment["q_start"]))
                frag_end = int(float(fragment["q_end"]))
                if frag_end > start and frag_start < end:
                    overlap.append(fragment)
            profile = ";".join(sorted(
                "{}:{}:{}:{}".format(row["label"], row["r_chr"], row["orientation"], row["fragment_id"])
                for row in overlap)) or "NA"
            rows.append({
                "component_id": "component_{:03d}".format(component_index), "copy1_bed_chr": chrom,
                "copy1_fasta_chr": fasta_chrom, "expected_copy2_bed_chr": expected_copy2.get(chrom, "NA"),
                "start": start, "end": end, "length_bp": end - start,
                "left_confirmed_junction": left_id, "right_confirmed_junction": right_id,
                "accepted_external_fragment_count": len(overlap), "external_fragment_profile": profile,
                "assignment_status": "neutral_component_not_subgenome_assigned",
                "claim_boundary": "Component coordinates follow manually confirmed junctions; external profile does not identify a biological subgenome or parent.",
            })
    write_tsv(outdir / "neutral_component_ledger.tsv", fields, rows)
    write_tsv(outdir / "parameters.tsv", ["parameter", "value"], [
        {"parameter": "confirmed_status", "value": "confirmed_rejoined_junction"},
        {"parameter": "unconfirmed_junctions", "value": "do not split chromosomes"},
        {"parameter": "interpretation", "value": "neutral components only; no A/B, LF/MF or parental-origin assignment"},
    ])
    print("[OK] Confirmed junctions={}; neutral components={}".format(len(confirmed), len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
