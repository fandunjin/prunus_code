#!/usr/bin/env python3
"""Retain i-ADHoRe A-B anchors whose genes occur in the same MCL family."""

from __future__ import print_function

import argparse
import csv
from collections import Counter
from pathlib import Path


def read_tsv(path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_mcl(path):
    gene_to_family = {}
    family_sizes = {}
    with path.open("r", encoding="utf-8") as handle:
        for index, line in enumerate(handle, 1):
            genes = line.rstrip("\n").split("\t")
            genes = [gene for gene in genes if gene]
            if not genes:
                continue
            family_id = "MCL_{:08d}".format(index)
            family_sizes[family_id] = len(genes)
            for gene in genes:
                if gene in gene_to_family:
                    raise ValueError("Gene occurs in multiple MCL families: {}".format(gene))
                gene_to_family[gene] = family_id
    if not family_sizes:
        raise ValueError("No MCL families found in {}".format(path))
    return gene_to_family, family_sizes


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, required=True)
    parser.add_argument("--mcl", type=Path, required=True)
    parser.add_argument("--dataset", choices=("PZI", "PSP", "PPA"), required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    return parser.parse_args()


def main():
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    pair_rows = read_tsv(args.pairs)
    gene_to_family, family_sizes = read_mcl(args.mcl)

    retained = []
    audit = []
    reasons = Counter()
    for row in pair_rows:
        gene_a = row["gene_a"]
        gene_b = row["gene_b"]
        family_a = gene_to_family.get(gene_a, "")
        family_b = gene_to_family.get(gene_b, "")
        if not family_a and not family_b:
            reason = "both_absent_from_mcl"
        elif not family_a:
            reason = "gene_a_absent_from_mcl"
        elif not family_b:
            reason = "gene_b_absent_from_mcl"
        elif family_a != family_b:
            reason = "different_mcl_families"
        else:
            reason = "retained_same_mcl_family"
        reasons[reason] += 1

        audit.append({
            "species": args.dataset,
            "pair_id": row["pair_id"],
            "gene_a": gene_a,
            "gene_b": gene_b,
            "mcl_family_a": family_a or "NA",
            "mcl_family_b": family_b or "NA",
            "mcl_family_size": family_sizes.get(family_a, "NA") if family_a == family_b else "NA",
            "status": reason,
            "representative_multiplicon": row["representative_multiplicon"],
            "duplicate_anchor_records": row["duplicate_anchor_records"],
        })
        if reason == "retained_same_mcl_family":
            retained.append({
                "family_id": row["pair_id"],
                "A": gene_a,
                "B": gene_b,
            })

    if not retained:
        raise SystemExit("[ERROR] No same-family i-ADHoRe A-B pairs retained for {}".format(args.dataset))

    audit_fields = [
        "species", "pair_id", "gene_a", "gene_b", "mcl_family_a", "mcl_family_b",
        "mcl_family_size", "status", "representative_multiplicon", "duplicate_anchor_records",
    ]
    write_tsv(args.outdir / (args.dataset + "_AB.mcl_pair_audit.tsv"), audit_fields, audit)
    write_tsv(
        args.outdir / (args.dataset + "_AB.families.tsv"),
        ["family_id", "A", "B"],
        retained,
    )
    summary = []
    for reason in (
        "retained_same_mcl_family", "different_mcl_families", "gene_a_absent_from_mcl",
        "gene_b_absent_from_mcl", "both_absent_from_mcl",
    ):
        summary.append({"dataset": args.dataset, "status": reason, "pairs": reasons[reason]})
    summary.extend([
        {"dataset": args.dataset, "status": "candidate_pairs_total", "pairs": len(pair_rows)},
        {"dataset": args.dataset, "status": "mcl_families_total", "pairs": len(family_sizes)},
        {"dataset": args.dataset, "status": "genes_assigned_to_mcl", "pairs": len(gene_to_family)},
    ])
    write_tsv(args.outdir / (args.dataset + "_AB.mcl_pair_summary.tsv"), ["dataset", "status", "pairs"], summary)
    print("[OK] {}: retained {}/{} i-ADHoRe A-B pairs in the same MCL family".format(
        args.dataset, len(retained), len(pair_rows)))


if __name__ == "__main__":
    main()
