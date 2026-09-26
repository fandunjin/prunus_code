#!/usr/bin/env python3
"""Convert OrthoFinder gene counts to a filtered CAFE5 input table."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-counts", required=True, help="OrthoFinder Orthogroups.GeneCount.tsv")
    parser.add_argument("--species", required=True, help="Comma-separated species columns to keep.")
    parser.add_argument("--out", required=True)
    parser.add_argument("--max-family-size", type=int, default=100)
    parser.add_argument("--min-species-present", type=int, default=2)
    args = parser.parse_args()

    species = [x.strip() for x in args.species.split(",") if x.strip()]
    rows = []
    with Path(args.gene_counts).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            counts = []
            for sp in species:
                try:
                    counts.append(int(row.get(sp, "0") or 0))
                except ValueError:
                    counts.append(0)
            if sum(1 for x in counts if x > 0) < args.min_species_present:
                continue
            if max(counts) > args.max_family_size:
                continue
            rows.append([row.get("Orthogroup", row.get("HOG", ""))] + counts)

    with Path(args.out).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["Desc", "Family ID"] + species)
        for family_id, *counts in rows:
            writer.writerow(["(null)", family_id] + counts)
    print("families=%d" % len(rows))


if __name__ == "__main__":
    main()
