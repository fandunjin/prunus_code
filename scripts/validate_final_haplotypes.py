#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


HAPS = ["hapA1", "hapA2", "hapB1", "hapB2"]


def read_tsv(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: missing or empty table: %s" % path)
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def count_fasta(path):
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit("ERROR: missing or empty FASTA: %s" % path)
    count = 0
    bp = 0
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
            elif line.strip():
                bp += len(line.strip())
    return count, bp


def write_tsv(path, rows, fieldnames):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key, "") for key in fieldnames})


def parse_args():
    parser = argparse.ArgumentParser(description="Validate final x=8 PZI haplotype outputs.")
    parser.add_argument("--final-hap-dir", required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    root = Path(args.final_hap_dir)
    overview = read_tsv(root / "FINAL_HAPLOTYPE_OVERVIEW.tsv")
    rename = read_tsv(root / "pzi_A1A2B1B2.rename_map.tsv")
    rows = []

    rows.append({
        "metric": "overview_quartets",
        "value": str(len(overview)),
        "expected": "8",
        "status": "pass" if len(overview) == 8 else "fail",
    })
    rename_counts = Counter(row.get("hap_set", "") for row in rename)
    rows.append({
        "metric": "rename_map_chromosomes",
        "value": str(len(rename)),
        "expected": "32",
        "status": "pass" if len(rename) == 32 else "fail",
    })
    for hap in HAPS:
        rows.append({
            "metric": "%s_rename_count" % hap,
            "value": str(rename_counts.get(hap, 0)),
            "expected": "8",
            "status": "pass" if rename_counts.get(hap, 0) == 8 else "fail",
        })
        count, bp = count_fasta(root / ("pzi_%s.fa" % hap))
        rows.append({
            "metric": "%s_fasta_seq_count" % hap,
            "value": str(count),
            "expected": "8",
            "status": "pass" if count == 8 else "fail",
        })
        rows.append({
            "metric": "%s_fasta_bp" % hap,
            "value": str(bp),
            "expected": ">0",
            "status": "pass" if bp > 0 else "fail",
        })
    count, bp = count_fasta(root / "pzi_A1A2B1B2.renamed.chr32.fa")
    rows.append({
        "metric": "all_renamed_chr32_fasta_seq_count",
        "value": str(count),
        "expected": "32",
        "status": "pass" if count == 32 else "fail",
    })
    rows.append({
        "metric": "all_renamed_chr32_fasta_bp",
        "value": str(bp),
        "expected": ">0",
        "status": "pass" if bp > 0 else "fail",
    })

    write_tsv(root / "FINAL_VALIDATION.tsv", rows, ["metric", "value", "expected", "status"])
    failures = [row["metric"] for row in rows if row["status"] != "pass"]
    if failures:
        raise SystemExit("ERROR: final haplotype validation failed: %s" % ",".join(failures))
    print("Final haplotype validation passed:", root)


if __name__ == "__main__":
    main()
