#!/usr/bin/env python3
"""Translate reviewed PPA-to-PZI-copy1 matches into PZI-style C##_A/B labels."""

import argparse
import csv
import re
from pathlib import Path


PZI_LABEL_RE = re.compile(r"^C(0[1-8])_([AB])[12]$")


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best-hits", type=Path, required=True)
    parser.add_argument("--pzi-naming-map", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--review-status", choices=("approved", "pending"), default="approved")
    args = parser.parse_args()
    for label, path in (("best-hit table", args.best_hits), ("PZI copy-haplotype naming map", args.pzi_naming_map)):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing {label}: {path}")

    best_rows = read_tsv(args.best_hits)
    required_best = {"ppa_primary_chr", "best_pzi_copy1_chr", "ppa_union_covered_bp", "best_to_second_ratio"}
    if len(best_rows) != 16 or not best_rows or not required_best <= set(best_rows[0]):
        raise SystemExit("ERROR: best-hit table must contain the 16-row primary-chain coverage summary")

    copy1_to_label = {}
    for row in read_tsv(args.pzi_naming_map):
        copy1 = row.get("copy1", "").strip()
        if not copy1:
            continue
        labels = [row.get(column, "").strip() for column in ("hapA1", "hapA2", "hapB1", "hapB2") if row.get(column, "").strip()]
        if len(labels) != 1 or not PZI_LABEL_RE.fullmatch(labels[0]) or copy1 in copy1_to_label:
            raise SystemExit(f"ERROR: invalid or duplicate copy1 naming row: {row}")
        copy1_to_label[copy1] = labels[0]
    if len(copy1_to_label) != 16:
        raise SystemExit(f"ERROR: expected 16 PZI copy1 labels, found {len(copy1_to_label)}")

    output_rows = []
    seen_ppa = set()
    seen_labels = set()
    for row in best_rows:
        ppa = row["ppa_primary_chr"]
        pzi_copy1 = row["best_pzi_copy1_chr"]
        pzi_hap_label = copy1_to_label.get(pzi_copy1)
        if not pzi_hap_label:
            raise SystemExit(f"ERROR: cannot resolve {pzi_copy1} to a PZI copy1 haplotype label")
        match = PZI_LABEL_RE.fullmatch(pzi_hap_label)
        assert match is not None
        renamed = f"C{match.group(1)}_{match.group(2)}"
        if ppa in seen_ppa or renamed in seen_labels:
            raise SystemExit(f"ERROR: non-bijective PPA correspondence involving {ppa} or {renamed}")
        seen_ppa.add(ppa)
        seen_labels.add(renamed)
        output_rows.append({
            "ppa_source_seqid": ppa,
            "pzi_copy1_chr": pzi_copy1,
            "pzi_copy1_haplotype_label": pzi_hap_label,
            "subgenome": match.group(2),
            "renamed_seqid": renamed,
            "ppa_union_covered_bp": row["ppa_union_covered_bp"],
            "ppa_coverage_fraction": row.get("ppa_coverage_fraction", ""),
            "best_to_second_ratio": row["best_to_second_ratio"],
            "review_status": args.review_status,
        })
    expected = {f"C{i:02d}_{sg}" for i in range(1, 9) for sg in ("A", "B")}
    if seen_labels != expected:
        raise SystemExit(f"ERROR: PZI-derived labels do not form C01-C08 for A/B: {sorted(seen_labels)}")
    output_rows.sort(key=lambda row: row["renamed_seqid"])
    columns = ["ppa_source_seqid", "pzi_copy1_chr", "pzi_copy1_haplotype_label", "subgenome", "renamed_seqid", "ppa_union_covered_bp", "ppa_coverage_fraction", "best_to_second_ratio", "review_status"]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    print(f"[OK] Wrote PZI-style PPA assignment: {args.out}")


if __name__ == "__main__":
    main()
