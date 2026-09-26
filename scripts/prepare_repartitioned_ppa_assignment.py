#!/usr/bin/env python3
"""Create the PPA assignment table for the revised PZI copy1 A/B partition.

The PPA-to-PZI copy1 chromosome matches are unchanged by the revised
partition. Only the subgenome and chromosome labels are recalculated from
the current copy1 order.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


A_COPY1 = ("chr01", "chr02", "chr05", "chr06", "chr07", "chr10", "chr12", "chr14")
B_COPY1 = ("chr03", "chr04", "chr08", "chr09", "chr11", "chr13", "chr15", "chr16")


def norm_chr(value: str) -> str:
    value = value.strip()
    if value.lower().startswith("chr"):
        return f"chr{int(value[3:]):02d}"
    return f"chr{int(value):02d}"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"ERROR: missing header: {path}")
        return [{key: (value or "").strip() for key, value in row.items()} for row in reader]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--best-hits", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--review-status", choices=("approved", "pending"), default="approved")
    args = parser.parse_args()

    rows = read_tsv(args.best_hits)
    required = {
        "ppa_primary_chr",
        "best_pzi_copy1_chr",
        "ppa_union_covered_bp",
        "best_to_second_ratio",
    }
    if len(rows) != 16 or not required.issubset(rows[0]):
        raise SystemExit(
            "ERROR: PPA best-hit table must contain 16 rows and columns "
            f"{sorted(required)}: {args.best_hits}"
        )

    expected_copy1 = set(A_COPY1 + B_COPY1)
    seen_ppa: set[str] = set()
    seen_copy1: set[str] = set()
    output: list[dict[str, str]] = []
    for row in rows:
        ppa = row["ppa_primary_chr"]
        copy1 = norm_chr(row["best_pzi_copy1_chr"])
        if ppa in seen_ppa or copy1 in seen_copy1:
            raise SystemExit(f"ERROR: PPA assignment is not one-to-one: {ppa} / {copy1}")
        if copy1 not in expected_copy1:
            raise SystemExit(f"ERROR: unexpected PZI copy1 chromosome: {copy1}")
        seen_ppa.add(ppa)
        seen_copy1.add(copy1)
        if copy1 in A_COPY1:
            subgenome = "A"
            rank = A_COPY1.index(copy1) + 1
        else:
            subgenome = "B"
            rank = B_COPY1.index(copy1) + 1
        label = f"C{rank}{subgenome}"
        output.append(
            {
                "ppa_source_seqid": ppa,
                "pzi_copy1_chr": copy1,
                "subgenome": subgenome,
                "renamed_seqid": label,
                "pzi_copy1_haplotype_label": f"{label}_{subgenome}1",
                "ppa_union_covered_bp": row["ppa_union_covered_bp"],
                "ppa_coverage_fraction": row.get("ppa_coverage_fraction", ""),
                "best_to_second_ratio": row["best_to_second_ratio"],
                "review_status": args.review_status,
            }
        )

    if seen_copy1 != expected_copy1:
        raise SystemExit(
            f"ERROR: expected all copy1 chromosomes, missing={sorted(expected_copy1 - seen_copy1)}"
        )

    output.sort(key=lambda row: (row["subgenome"], int(row["renamed_seqid"][1:-1])))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "ppa_source_seqid",
        "pzi_copy1_chr",
        "subgenome",
        "renamed_seqid",
        "ppa_union_covered_bp",
        "best_to_second_ratio",
        "review_status",
        "ppa_coverage_fraction",
        "pzi_copy1_haplotype_label",
    ]
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output)
    print(f"[OK] Wrote revised PPA assignment: {args.out}")


if __name__ == "__main__":
    main()
