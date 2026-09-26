#!/usr/bin/env python3
"""Combine the uniform YN Ks tables for the three-reference Figure 2a plot."""

from __future__ import annotations

import argparse
import csv
import math
import re
import statistics
from collections import Counter, defaultdict
from pathlib import Path


EXPECTED = {
    f"{query}_vs_{reference}"
    for query in (
        "pzi_hapA1", "pzi_hapB1", "psp_hapA1", "psp_hapB1", "ppa_hapA", "ppa_hapB"
    )
    for reference in ("pcam", "ppe", "gtr")
}
QUERY_INFO = {
    "pzi_hapA1": ("PZI", "A"),
    "pzi_hapB1": ("PZI", "B"),
    "psp_hapA1": ("PSP", "A"),
    "psp_hapB1": ("PSP", "B"),
    "ppa_hapA": ("PPA", "A"),
    "ppa_hapB": ("PPA", "B"),
}
REFERENCE_INFO = {"pcam": "PCAM", "ppe": "PPE", "gtr": "GTR"}
FIELDNAMES = [
    "comparison", "species", "subgenome", "reference", "source_result",
    "pair_id", "query_gene", "ref_gene", "block_id", "Ks", "method",
]


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtr", type=Path, required=True)
    parser.add_argument("--pzi-reference", type=Path, required=True)
    parser.add_argument("--new-reference", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def source_rows(path: Path, expected: set[str]) -> list[dict[str, object]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"[ERROR] Missing Ks table: {path}")
    rows: list[dict[str, object]] = []
    seen = Counter()
    for row in read_rows(path):
        comparison = row.get("comparison", "")
        if comparison not in expected:
            continue
        try:
            value = float(row.get("Ks", ""))
        except (TypeError, ValueError):
            continue
        if not math.isfinite(value) or not 0.0 < value < 5.0:
            continue
        query, reference = comparison.split("_vs_", 1)
        species, subgenome = QUERY_INFO[query]
        rows.append({
            "comparison": comparison,
            "species": species,
            "subgenome": subgenome,
            "reference": REFERENCE_INFO[reference],
            "source_result": str(path),
            "pair_id": row.get("pair_id", ""),
            "query_gene": row.get("query_gene", ""),
            "ref_gene": row.get("ref_gene", ""),
            "block_id": row.get("block_id", ""),
            "Ks": f"{value:.8f}",
            "method": row.get("method", ""),
        })
        seen[comparison] += 1
    missing = sorted(expected - set(seen))
    if missing:
        raise SystemExit(f"[ERROR] Missing comparisons in {path}: {', '.join(missing)}")
    return rows


def main() -> int:
    args = parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()) and not args.overwrite:
        raise SystemExit(f"[ERROR] Output directory is not empty: {args.outdir}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    gtr_expected = {f"{query}_vs_gtr" for query in QUERY_INFO}
    pzi_expected = {
        "pzi_hapA1_vs_pcam", "pzi_hapB1_vs_pcam",
        "pzi_hapA1_vs_ppe", "pzi_hapB1_vs_ppe",
    }
    new_expected = EXPECTED - gtr_expected - pzi_expected
    rows = (
        source_rows(args.gtr, gtr_expected)
        + source_rows(args.pzi_reference, pzi_expected)
        + source_rows(args.new_reference, new_expected)
    )
    found = {str(row["comparison"]) for row in rows}
    if found != EXPECTED:
        raise SystemExit(f"[ERROR] Combined comparisons do not match expected set: {sorted(EXPECTED - found)}")
    bad_methods = sorted({str(row["method"]) for row in rows if str(row["method"]) != "YN"})
    if bad_methods:
        raise SystemExit(f"[ERROR] Non-YN method values found: {', '.join(bad_methods)}")

    write_tsv(args.outdir / "subgenome_reference_ks.tsv", FIELDNAMES, rows)
    summary: list[dict[str, object]] = []
    grouped: dict[str, list[float]] = defaultdict(list)
    for row in rows:
        grouped[str(row["comparison"])].append(float(row["Ks"]))
    for comparison in sorted(grouped):
        values = grouped[comparison]
        query, reference = comparison.split("_vs_", 1)
        species, subgenome = QUERY_INFO[query]
        summary.append({
            "comparison": comparison,
            "species": species,
            "subgenome": subgenome,
            "reference": REFERENCE_INFO[reference],
            "n": len(values),
            "Ks_median": f"{statistics.median(values):.8f}",
            "Ks_mean": f"{statistics.mean(values):.8f}",
            "method": "YN",
            "filter": "0<Ks<5",
        })
    write_tsv(
        args.outdir / "subgenome_reference_summary.tsv",
        ["comparison", "species", "subgenome", "reference", "n", "Ks_median", "Ks_mean", "method", "filter"],
        summary,
    )
    write_tsv(
        args.outdir / "source_result_manifest.tsv",
        ["source_result", "n_rows", "comparisons"],
        [
            {"source_result": str(path), "n_rows": sum(1 for row in rows if row["source_result"] == str(path)),
             "comparisons": ",".join(sorted({str(row["comparison"]) for row in rows if row["source_result"] == str(path)}))}
            for path in (args.gtr, args.pzi_reference, args.new_reference)
        ],
    )
    (args.outdir / "method_contract.txt").write_text(
        "All comparisons use JCVI collinear pairs, MAFFT protein alignment, PAL2NAL codon back-translation, "
        "and KaKs_Calculator YN output with 0<Ks<5 retained for plotting.\n"
        "A and B are displayed as independent query subgenomes; PCAM, PPE and GTR are reference genomes.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Three-reference subgenome Ks source data completed.\n", encoding="ascii")
    print(f"[OK] Combined {len(rows)} Ks rows across {len(found)} comparisons: {args.outdir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
