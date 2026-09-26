#!/usr/bin/env python3
"""Combine genome-native self-syntenic Ks source tables for six taxa."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


UPPER = ("PZI", "PSP", "PPA")
LOWER = ("PPE", "PCAM", "GTR")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError("No header in {}".format(path))
        return list(reader)


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_ks(value: str, path: Path) -> float:
    try:
        ks = float(value)
    except (TypeError, ValueError):
        raise ValueError("Invalid Ks value in {}: {}".format(path, value))
    if not math.isfinite(ks) or not 0.0 < ks < 5.0:
        raise ValueError("Unfiltered Ks outside 0<Ks<5 in {}: {}".format(path, value))
    return ks


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-root", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    if args.outdir.exists() and any(args.outdir.iterdir()):
        raise SystemExit("[ERROR] Output directory is not empty: {}".format(args.outdir))
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, object]] = []
    for panel, labels in (("upper", UPPER), ("lower", LOWER)):
      for label in labels:
        path = args.self_root / label / "audit" / "{}_SELF.ks_0_5.tsv".format(label)
        complete = args.self_root / label / "COMPLETE"
        if not complete.is_file() or not path.is_file() or path.stat().st_size == 0:
            raise SystemExit("[ERROR] Missing completed native self-syntenic Ks for {}".format(label))
        for row in read_tsv(path):
            ks = parse_ks(row.get("dS", ""), path)
            rows.append({
                "panel": panel, "comparison": label,
                "comparison_class": "representative_A_B_composite_self_syntenic_paralogs" if label in UPPER else "genome_native_cross_chromosome_self_syntenic_paralogs",
                "family_id": row.get("family_id", ""),
                "gene_1": row.get("copy_1", row.get("gene1", "")),
                "gene_2": row.get("copy_2", row.get("gene2", "")),
                "Ks": "{:.8f}".format(ks), "source": "iADHoRe_self_anchor_CodeML",
            })

    fields = ["panel", "comparison", "comparison_class", "family_id", "gene_1", "gene_2", "Ks", "source"]
    write_tsv(args.outdir / "within_genome_ks.tsv", fields, rows)
    summary: list[dict[str, object]] = []
    for panel, labels in (("upper", UPPER), ("lower", LOWER)):
        for label in labels:
            values = [float(row["Ks"]) for row in rows if row["panel"] == panel and row["comparison"] == label]
            if not values:
                raise SystemExit("[ERROR] Empty source for {}".format(label))
            summary.append({"panel": panel, "comparison": label, "comparison_class": next(
                str(row["comparison_class"]) for row in rows if row["panel"] == panel and row["comparison"] == label),
                "n": len(values), "Ks_median": "{:.8f}".format(statistics.median(values)),
                "Ks_mean": "{:.8f}".format(statistics.mean(values)), "filter": "0<Ks<5"})
    write_tsv(args.outdir / "within_genome_ks_summary.tsv",
              ["panel", "comparison", "comparison_class", "n", "Ks_median", "Ks_mean", "filter"], summary)
    (args.outdir / "figure_contract.txt").write_text(
        "Core conclusion: Genome-native self-syntenic paralog divergence can be compared across six taxa after reducing the phased allopolyploids to one representative A+B component set.\n"
        "Archetype: quantitative two-row comparison.\n"
        "Backend: Python matplotlib.\n"
        "Upper evidence: i-ADHoRe self-syntenic anchors in representative PZI A1+B1, PSP A1+B1 and PPA A+B component sets.\n"
        "Lower evidence: genome-native, cross-chromosome i-ADHoRe self-syntenic anchors in PPE, PCAM and GTR.\n"
        "Review risk: a peak or its similarity alone does not prove a shared WGD.\n",
        encoding="ascii")
    (args.outdir / "method_contract.txt").write_text(
        "All plotted values are pairwise CodeML dS estimates filtered to 0<Ks<5.\n"
        "Upper: representative PZI A1+B1, PSP A1+B1 and PPA A+B genome self-syntenic i-ADHoRe anchors.\n"
        "Lower: PPE, PCAM and GTR cross-chromosome self-syntenic i-ADHoRe anchors.\n"
        "The lower pairs exclude same-list anchors, self-hits and duplicated anchor records.\n",
        encoding="ascii")
    (args.outdir / "COMPLETE").write_text("Within-genome Ks source data combined.\n", encoding="ascii")
    print("[OK] Combined within-genome Ks: {} rows".format(len(rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
