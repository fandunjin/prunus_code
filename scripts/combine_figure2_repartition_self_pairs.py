#!/usr/bin/env python3
"""Convert corrected i-ADHoRe self-synteny results to the common Ks pair schema."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "comparison", "comparison_type", "block_id", "query_label", "ref_label",
    "query_gene", "ref_gene", "query_seqid", "query_start0", "query_end",
    "query_strand", "ref_seqid", "ref_start0", "ref_end", "ref_strand",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise RuntimeError("Missing header: {}".format(path))
        return list(reader)


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-plan", type=Path, required=True)
    parser.add_argument("--self-root", type=Path, required=True)
    parser.add_argument("--out-pairs", type=Path, required=True)
    parser.add_argument("--out-status", type=Path, required=True)
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    status: list[dict[str, object]] = []
    seen: set[tuple[str, str, str]] = set()
    for plan in read_tsv(args.self_plan):
        label = plan.get("label", "")
        if not label:
            raise RuntimeError("Self plan has an empty label")
        path = args.self_root / label / "03_pairs" / "self_collinear_pairs.tsv"
        if not path.is_file() or path.stat().st_size == 0:
            status.append({"label": label, "status": "missing_self_pairs", "pairs": 0, "source": path})
            continue
        source = read_tsv(path)
        retained = 0
        for item in source:
            gene_a, gene_b = item.get("gene_1", ""), item.get("gene_2", "")
            seq_a, seq_b = item.get("seqid_1", ""), item.get("seqid_2", "")
            if not gene_a or not gene_b or not seq_a or not seq_b or gene_a == gene_b or seq_a == seq_b:
                raise RuntimeError("Invalid self-synteny pair in {}".format(path))
            pair_key = (label,) + tuple(sorted((gene_a, gene_b)))
            if pair_key in seen:
                continue
            seen.add(pair_key)
            rows.append({
                "comparison": "{}_self".format(label),
                "comparison_type": "self_synteny",
                "block_id": item.get("multiplicon", "") or item.get("family_id", ""),
                "query_label": label,
                "ref_label": label,
                "query_gene": gene_a,
                "ref_gene": gene_b,
                "query_seqid": seq_a,
                "query_start0": "",
                "query_end": "",
                "query_strand": "",
                "ref_seqid": seq_b,
                "ref_start0": "",
                "ref_end": "",
                "ref_strand": "",
            })
            retained += 1
        status.append({"label": label, "status": "ok" if retained else "no_pairs", "pairs": retained, "source": path})

    write_tsv(args.out_pairs, rows, FIELDS)
    write_tsv(args.out_status, status, ["label", "status", "pairs", "source"])
    if not rows or any(row["status"] != "ok" for row in status):
        return 2
    print("[OK] Combined {} real cross-chromosome self-synteny pairs across {} tracks.".format(len(rows), len(status)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
