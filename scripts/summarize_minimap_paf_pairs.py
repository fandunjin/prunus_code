#!/usr/bin/env python3
"""Summarize primary minimap2 PAF chains as non-overlapping coverage.

PAF query records are PZI copy1 chromosomes and targets are PPA sequences.
For every PPA-PZI chromosome pair, coverage is calculated as a union of PAF
intervals, avoiding the inflation caused by summing overlapping chains.
"""

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def union_length(intervals):
    if not intervals:
        return 0
    ordered = sorted(intervals)
    start, end = ordered[0]
    total = 0
    for next_start, next_end in ordered[1:]:
        if next_start > end:
            total += end - start
            start, end = next_start, next_end
        else:
            end = max(end, next_end)
    return total + end - start


def write_tsv(path, rows, columns):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--paf", type=Path, required=True)
    parser.add_argument("--ppa-map", type=Path, required=True)
    parser.add_argument("--out-filtered-paf", type=Path, required=True)
    parser.add_argument("--out-summary", type=Path, required=True)
    parser.add_argument("--out-best", type=Path, required=True)
    parser.add_argument("--min-mapq", type=int, required=True)
    parser.add_argument("--min-alignment-len", type=int, required=True)
    args = parser.parse_args()

    with args.ppa_map.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "ppa_source_seqid" not in reader.fieldnames:
            raise SystemExit(f"ERROR: invalid PPA chromosome map: {args.ppa_map}")
        ppa_chromosomes = {row["ppa_source_seqid"] for row in reader if row.get("ppa_source_seqid")}
    if len(ppa_chromosomes) != 16:
        raise SystemExit(f"ERROR: expected 16 PPA chromosomes, found {len(ppa_chromosomes)}")

    pairs = defaultdict(lambda: {
        "ppa_intervals": [], "pzi_intervals": [], "matching_bp_sum": 0,
        "block_bp_sum": 0, "records": 0, "ppa_len": 0, "pzi_len": 0,
    })
    total = retained = 0
    with args.paf.open("r", encoding="utf-8", errors="replace") as inp, args.out_filtered_paf.open("w", encoding="utf-8", newline="") as out:
        for line_no, raw in enumerate(inp, 1):
            if not raw.strip():
                continue
            total += 1
            fields = raw.rstrip("\n").split("\t")
            if len(fields) < 12:
                raise SystemExit(f"ERROR: malformed PAF at line {line_no}")
            pzi, pzi_len, pzi_start, pzi_end = fields[0], int(fields[1]), int(fields[2]), int(fields[3])
            ppa, ppa_len, ppa_start, ppa_end = fields[5], int(fields[6]), int(fields[7]), int(fields[8])
            matching_bp, block_bp, mapq = int(fields[9]), int(fields[10]), int(fields[11])
            if ppa not in ppa_chromosomes or mapq < args.min_mapq or block_bp < args.min_alignment_len or "tp:A:P" not in fields[12:]:
                continue
            out.write(raw if raw.endswith("\n") else raw + "\n")
            retained += 1
            item = pairs[(ppa, pzi)]
            item["ppa_intervals"].append((ppa_start, ppa_end))
            item["pzi_intervals"].append((pzi_start, pzi_end))
            item["matching_bp_sum"] += matching_bp
            item["block_bp_sum"] += block_bp
            item["records"] += 1
            item["ppa_len"] = ppa_len
            item["pzi_len"] = pzi_len
    if retained == 0:
        raise SystemExit("ERROR: no primary alignments passed filters")

    columns = ["reference_chr", "query_chr", "ppa_length_bp", "pzi_length_bp", "ppa_union_covered_bp", "pzi_union_covered_bp", "aligned_block_bp_sum", "matching_bp_sum", "alignment_records", "mean_identity_raw_blocks", "ppa_coverage_fraction", "pzi_coverage_fraction"]
    rows = []
    for (ppa, pzi), item in pairs.items():
        ppa_length = int(item["ppa_len"])
        pzi_length = int(item["pzi_len"])
        ppa_covered = union_length(item["ppa_intervals"])
        pzi_covered = union_length(item["pzi_intervals"])
        block_sum = int(item["block_bp_sum"])
        matching_sum = int(item["matching_bp_sum"])
        rows.append({
            "reference_chr": ppa, "query_chr": pzi,
            "ppa_length_bp": ppa_length, "pzi_length_bp": pzi_length,
            "ppa_union_covered_bp": ppa_covered, "pzi_union_covered_bp": pzi_covered,
            "aligned_block_bp_sum": block_sum, "matching_bp_sum": matching_sum,
            "alignment_records": int(item["records"]),
            "mean_identity_raw_blocks": f"{matching_sum / block_sum:.6f}",
            "ppa_coverage_fraction": f"{ppa_covered / ppa_length:.6f}",
            "pzi_coverage_fraction": f"{pzi_covered / pzi_length:.6f}",
        })
    rows.sort(key=lambda row: (str(row["reference_chr"]), -int(row["ppa_union_covered_bp"]), str(row["query_chr"])))
    write_tsv(args.out_summary, rows, columns)

    best_columns = ["ppa_primary_chr", "best_pzi_copy1_chr", "ppa_length_bp", "pzi_length_bp", "ppa_union_covered_bp", "pzi_union_covered_bp", "ppa_coverage_fraction", "pzi_coverage_fraction", "alignment_records", "mean_identity_raw_blocks", "second_best_ppa_union_covered_bp", "best_to_second_ratio"]
    by_ppa = defaultdict(list)
    for row in rows:
        by_ppa[str(row["reference_chr"])].append(row)
    best_rows = []
    for ppa in sorted(by_ppa):
        ranked = sorted(by_ppa[ppa], key=lambda row: (-int(row["ppa_union_covered_bp"]), str(row["query_chr"])))
        best = ranked[0]
        best_bp = int(best["ppa_union_covered_bp"])
        second_bp = int(ranked[1]["ppa_union_covered_bp"]) if len(ranked) > 1 else 0
        best_rows.append({
            "ppa_primary_chr": ppa, "best_pzi_copy1_chr": best["query_chr"],
            "ppa_length_bp": best["ppa_length_bp"], "pzi_length_bp": best["pzi_length_bp"],
            "ppa_union_covered_bp": best_bp, "pzi_union_covered_bp": best["pzi_union_covered_bp"],
            "ppa_coverage_fraction": best["ppa_coverage_fraction"], "pzi_coverage_fraction": best["pzi_coverage_fraction"],
            "alignment_records": best["alignment_records"], "mean_identity_raw_blocks": best["mean_identity_raw_blocks"],
            "second_best_ppa_union_covered_bp": second_bp,
            "best_to_second_ratio": "inf" if second_bp == 0 else f"{best_bp / second_bp:.6f}",
        })
    write_tsv(args.out_best, best_rows, best_columns)
    print(f"[OK] Retained {retained}/{total} primary PAF records")


if __name__ == "__main__":
    main()
