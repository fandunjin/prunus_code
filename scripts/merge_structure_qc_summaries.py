#!/usr/bin/env python3
"""Combine per-species gene_structure_summary.tsv files into one table."""

import argparse
import csv
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", action="append", required=True, metavar="SPECIES=PATH")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    metric_columns = ["species", "Gene number", "RNA number", "CDS number", "Exon number", "Intron number", "Mean exons per transcript", "Mean CDS length (bp)", "Mean exon length (bp)", "Mean intron length (bp)", "Long introns >=1 kb", "Long introns >=10 kb", "TE-overlapping introns"]
    for item in args.input:
        if "=" not in item:
            raise SystemExit(f"ERROR: --input must use SPECIES=PATH: {item}")
        species, path_text = item.split("=", 1)
        path = Path(path_text)
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: summary missing or empty for {species}: {path}")
        with path.open(encoding="utf-8", newline="") as handle:
            source_rows = list(csv.DictReader(handle, delimiter="\t"))
        if source_rows and "metric" in source_rows[0]:
            values = {row.get("metric", ""): row.get("value", "") for row in source_rows}
            rows.append({"species": species, **{key: values.get(key, "") for key in metric_columns[1:]}})
        else:
            row = source_rows[0] if source_rows else {}
            rows.append({"species": species, **{key: row.get(key, "") for key in metric_columns[1:]}})
    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=metric_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Combined structure summary: {args.output}")


if __name__ == "__main__":
    main()
