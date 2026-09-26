#!/usr/bin/env python3
"""Validate cleaned haplotype ncRNA GFF3 files without changing them."""

import argparse
import csv
from collections import Counter
from pathlib import Path


ALLOWED_CLASSES = {
    "tRNA", "tRNA pseudogene", "rRNA", "snRNA", "snoRNA",
    "miRNA", "SRP RNA", "other ncRNA",
}


def attrs(text):
    result = {}
    for item in text.rstrip(";").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def validate(path: Path):
    seen_ids = set()
    bad_rows = bad_coords = bad_strand = duplicate_ids = bad_class = 0
    exon_rows = 0
    observed = Counter()
    seqids = set()
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                bad_rows += 1
                continue
            seqid, _source, feature, start, end, score, strand, _phase, attribute_text = fields
            seqids.add(seqid)
            try:
                start_i, end_i = int(start), int(end)
                if start_i < 1 or end_i < start_i:
                    bad_coords += 1
            except ValueError:
                bad_coords += 1
            if strand not in {"+", "-", "."}:
                bad_strand += 1
            if score not in {".", ""}:
                try:
                    float(score)
                except ValueError:
                    bad_rows += 1
            parsed = attrs(attribute_text)
            record_id = parsed.get("ID")
            if not record_id:
                bad_rows += 1
            elif record_id in seen_ids:
                duplicate_ids += 1
            else:
                seen_ids.add(record_id)
            feature_class = parsed.get("ncRNA_class", "")
            if feature == "exon":
                exon_rows += 1
            if feature_class not in ALLOWED_CLASSES:
                bad_class += 1
            observed[feature_class] += 1
    errors = bad_rows + bad_coords + bad_strand + duplicate_ids + bad_class + exon_rows
    return {
        "file": str(path),
        "records": sum(observed.values()),
        "seqids": len(seqids),
        "bad_rows": bad_rows,
        "bad_coords": bad_coords,
        "bad_strand": bad_strand,
        "duplicate_ids": duplicate_ids,
        "bad_class": bad_class,
        "exon_rows": exon_rows,
        "status": "PASS" if errors == 0 else "FAIL",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903/annotation/ncrna_haploid_clean"))
    args = parser.parse_args()
    files = sorted(args.root.glob("*/hap*/final/*.ncRNA.final.gff3"))
    if len(files) != 10:
        raise SystemExit(f"expected 10 haplotype GFF3 files, found {len(files)} under {args.root}")
    rows = [validate(path) for path in files]
    out = args.root / "ncrna.haploid.validation.tsv"
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    for row in rows:
        print(f"{row['status']} {row['file']} records={row['records']} seqids={row['seqids']}")
    print(f"validation: {out}")
    if any(row["status"] != "PASS" for row in rows):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
