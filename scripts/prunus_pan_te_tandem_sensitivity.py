#!/usr/bin/env python3
"""Aggregate TRF-score sensitivity partitions without treating haplotypes as replicates."""

import argparse
import csv
import statistics
from pathlib import Path


GROUPS = ("PZI", "PSP", "PPA", "PCAM", "PPE")
FIELDS = ("genome_bp", "trf_tandem_bp", "pan_te_non_tandem_bp", "other_bp")


def read_tsv(path):
    with open(path, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, help="label:path to tandem_main_partition.tsv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    rows = []
    for spec in args.input:
        if ":" not in spec:
            raise SystemExit("ERROR: --input must be label:path")
        label, raw_path = spec.split(":", 1)
        path = Path(raw_path)
        if not path.is_file():
            raise SystemExit(f"ERROR: missing partition table: {path}")
        source = read_tsv(path)
        for group in GROUPS:
            selected = [row for row in source if row["group"] == group]
            if not selected:
                continue
            out = {"trf_score": label, "group": group, "genome_packages": len(selected)}
            for field in FIELDS:
                values = [int(row[field]) for row in selected]
                out[f"mean_{field}"] = round(statistics.mean(values))
                out[f"min_{field}"] = min(values)
                out[f"max_{field}"] = max(values)
            rows.append(out)
    fields = ["trf_score", "group", "genome_packages"] + [f"{stat}_{field}" for field in FIELDS for stat in ("mean", "min", "max")]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader(); writer.writerows(rows)
    print(f"Wrote {len(rows)} group-score rows to {args.output}")


if __name__ == "__main__":
    main()
