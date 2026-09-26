#!/usr/bin/env python3
"""Create a manually auditable junction-review sheet without making a call.

This joins the synteny validation table with optional HiFi support.  Hi-C
coordinate compatibility and local-map interpretation remain human-reviewed
fields because no available matrix has yet been verified as copy1-based.
"""

from __future__ import print_function

import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--validation", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--hifi", default="", help="Optional junction_hifi_support.tsv")
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    args = parse_args()
    validation = read_tsv(args.validation)
    hifi = {}
    if args.hifi:
        for row in read_tsv(args.hifi):
            hifi[row["junction_id"]] = row
    fields = [
        "junction_id", "copy1_chr", "boundary_mid_bp", "external_support_count", "supporting_references",
        "copy1_copy2_two_sided_support", "same_expected_copy2_target", "copy2_target_switch",
        "copy1_copy2_pattern", "haplotypes_with_two_sided_coverage", "synteny_status",
        "hifi_boundary_spanning_reads", "hifi_status", "hic_file_coordinate_system_confirmed",
        "hic_local_adjacency_support", "manual_reviewer", "manual_review_notes", "final_status",
    ]
    rows = []
    for row in validation:
        reads = hifi.get(row["junction_id"], {})
        rows.append({
            "junction_id": row["junction_id"], "copy1_chr": row["copy1_chr"],
            "boundary_mid_bp": row["boundary_mid_bp"],
            "external_support_count": row["external_support_count"],
            "supporting_references": row["supporting_references"],
            "copy1_copy2_two_sided_support": row["copy1_copy2_two_sided_support"],
            "same_expected_copy2_target": row.get("same_expected_copy2_target", "NA"),
            "copy2_target_switch": row.get("copy2_target_switch", "NA"),
            "copy1_copy2_pattern": row.get("copy1_copy2_pattern", "NA"),
            "haplotypes_with_two_sided_coverage": row["haplotypes_with_two_sided_coverage"],
            "synteny_status": row["status"],
            "hifi_boundary_spanning_reads": reads.get("boundary_spanning_reads", "NA"),
            "hifi_status": reads.get("status", "not_run"),
            "hic_file_coordinate_system_confirmed": "",
            "hic_local_adjacency_support": "",
            "manual_reviewer": "",
            "manual_review_notes": "",
            "final_status": "unreviewed",
        })
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print("[OK] Wrote {} manual junction-review rows: {}".format(len(rows), out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
