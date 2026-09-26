#!/usr/bin/env python3
"""Prepare a complete PZI-versus-PSP chromosome heatmap matrix from PAF."""

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path


CHR_NAME = re.compile(r"^C(0[1-8])_([AB])([12])$")
PALETTE = [
    "#F6AD55",
    "#68D391",
    "#63B3ED",
    "#F687B3",
    "#4FD1C5",
    "#ED8936",
    "#9F7AEA",
    "#ECC94B",
]


def natural_key(text):
    parts = []
    buffer = ""
    digit = None
    for char in text:
        current = char.isdigit()
        if digit is not None and current != digit:
            parts.append(int(buffer) if digit else buffer)
            buffer = char
        else:
            buffer += char
        digit = current
    if buffer:
        parts.append(int(buffer) if digit else buffer)
    return parts


def read_fasta_ids(path):
    ids = []
    seen = set()
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seqid = line[1:].split()[0]
            if seqid in seen:
                raise SystemExit(f"ERROR: duplicate FASTA ID in {path}: {seqid}")
            seen.add(seqid)
            ids.append(seqid)
    return sorted(ids, key=natural_key)


def interval_union_bp(intervals):
    if not intervals:
        return 0
    ordered = sorted(intervals)
    start, end = ordered[0]
    total = 0
    for next_start, next_end in ordered[1:]:
        if next_start <= end:
            end = max(end, next_end)
        else:
            total += end - start
            start, end = next_start, next_end
    return total + end - start


def parse_paf(path, reference_ids, query_ids, min_mapq, min_aln_len):
    agg = defaultdict(
        lambda: {
            "aln_len": 0,
            "matches": 0,
            "target_intervals": [],
            "query_intervals": [],
            "plus_bp": 0,
            "minus_bp": 0,
        }
    )
    retained = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            query = fields[0]
            target = fields[5]
            if target not in reference_ids or query not in query_ids:
                continue
            aln_len = int(fields[10])
            mapq = int(fields[11])
            if mapq < min_mapq or aln_len < min_aln_len:
                continue
            key = (target, query)
            agg[key]["aln_len"] += aln_len
            agg[key]["matches"] += int(fields[9])
            agg[key]["target_intervals"].append((int(fields[7]), int(fields[8])))
            agg[key]["query_intervals"].append((int(fields[2]), int(fields[3])))
            if fields[4] == "+":
                agg[key]["plus_bp"] += aln_len
            else:
                agg[key]["minus_bp"] += aln_len
            retained += 1
    return agg, retained


def write_tsv(path, rows, fields):
    with open(path, "w", encoding="utf-8", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def chromosome_color(seqid):
    match = CHR_NAME.match(seqid)
    if not match:
        return "#A0AEC0", "unclassified"
    group = int(match.group(1))
    return PALETTE[group - 1], f"C{group:02d}_quartet"


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--paf", required=True)
    parser.add_argument("--reference-fasta", required=True)
    parser.add_argument("--query-fasta", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--axis-colors", required=True)
    parser.add_argument("--best-hits", required=True)
    parser.add_argument("--validation", required=True)
    parser.add_argument("--min-mapq", type=int, default=10)
    parser.add_argument("--min-aln-len", type=int, default=10000)
    return parser.parse_args()


def main():
    args = parse_args()
    reference_ids = read_fasta_ids(args.reference_fasta)
    query_ids = read_fasta_ids(args.query_fasta)
    if len(reference_ids) != 32 or len(query_ids) != 32:
        raise SystemExit(
            f"ERROR: expected 32 reference and 32 query chromosomes; observed "
            f"{len(reference_ids)} and {len(query_ids)}"
        )

    agg, retained = parse_paf(
        args.paf,
        set(reference_ids),
        set(query_ids),
        args.min_mapq,
        args.min_aln_len,
    )
    if retained == 0:
        raise SystemExit("ERROR: no PAF alignments passed the requested filters")

    summary_rows = []
    best_rows = []
    for target in reference_ids:
        for query in query_ids:
            stats = agg.get((target, query))
            if stats:
                aligned = stats["aln_len"]
                matches = stats["matches"]
                target_covered = interval_union_bp(stats["target_intervals"])
                query_covered = interval_union_bp(stats["query_intervals"])
                orientation = "+" if stats["plus_bp"] >= stats["minus_bp"] else "-"
            else:
                aligned = 0
                matches = 0
                target_covered = 0
                query_covered = 0
                orientation = ""
            summary_rows.append(
                {
                    "reference_chr": target,
                    "query_chr": query,
                    "aligned_block_bp": aligned,
                    "matching_bp": matches,
                    "aligned_covered_bp": target_covered,
                    "query_covered_bp": query_covered,
                    "dominant_orientation": orientation,
                }
            )

    for query in query_ids:
        candidates = [row for row in summary_rows if row["query_chr"] == query]
        ranked = sorted(
            candidates,
            key=lambda row: (-int(row["aligned_block_bp"]), natural_key(row["reference_chr"])),
        )
        best = ranked[0]
        total = sum(int(row["aligned_block_bp"]) for row in candidates)
        best_bp = int(best["aligned_block_bp"])
        best_fraction = best_bp / total if total else 0.0
        expected = query if query in reference_ids else ""
        best_rows.append(
            {
                "psp_chr": query,
                "best_pzi_chr": best["reference_chr"],
                "expected_same_name": expected,
                "best_aligned_block_bp": best_bp,
                "total_aligned_block_bp": total,
                "best_fraction": f"{best_fraction:.6f}",
                "dominant_orientation": best["dominant_orientation"],
                "status": "pass" if best["reference_chr"] == expected else "review",
            }
        )

    color_rows = []
    for seqid in reference_ids:
        color, label = chromosome_color(seqid)
        color_rows.append(
            {"axis": "target", "chr": seqid, "color": color, "label": label}
        )
    for seqid in query_ids:
        color, label = chromosome_color(seqid)
        color_rows.append(
            {"axis": "query", "chr": seqid, "color": color, "label": label}
        )

    same_name_count = sum(row["status"] == "pass" for row in best_rows)
    validation_rows = [
        {"check": "reference_chromosomes", "observed": len(reference_ids), "expected": 32, "status": "pass"},
        {"check": "query_chromosomes", "observed": len(query_ids), "expected": 32, "status": "pass"},
        {"check": "retained_paf_alignments", "observed": retained, "expected": ">0", "status": "pass"},
        {
            "check": "same_name_best_hits",
            "observed": same_name_count,
            "expected": 32,
            "status": "pass" if same_name_count == 32 else "review",
        },
        {
            "check": "reference_query_name_sets_equal",
            "observed": len(set(reference_ids) & set(query_ids)),
            "expected": 32,
            "status": "pass" if set(reference_ids) == set(query_ids) else "fail",
        },
    ]

    write_tsv(
        args.summary,
        summary_rows,
        [
            "reference_chr",
            "query_chr",
            "aligned_block_bp",
            "matching_bp",
            "aligned_covered_bp",
            "query_covered_bp",
            "dominant_orientation",
        ],
    )
    write_tsv(args.axis_colors, color_rows, ["axis", "chr", "color", "label"])
    write_tsv(
        args.best_hits,
        best_rows,
        [
            "psp_chr",
            "best_pzi_chr",
            "expected_same_name",
            "best_aligned_block_bp",
            "total_aligned_block_bp",
            "best_fraction",
            "dominant_orientation",
            "status",
        ],
    )
    write_tsv(args.validation, validation_rows, ["check", "observed", "expected", "status"])
    print(f"Retained PAF alignments: {retained}")
    print(f"Same-name best hits: {same_name_count}/32")


if __name__ == "__main__":
    main()

