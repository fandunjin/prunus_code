#!/usr/bin/env python3
"""Summarize OrthoFinder duplication events for a copyAB species tree."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--result", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    ap.add_argument("--support", type=float, default=0.5)
    args = ap.parse_args()

    dup = args.result / "Gene_Duplication_Events" / "Duplications.tsv"
    counts = (
        args.result
        / "Comparative_Genomics_Statistics"
        / "Duplications_per_Species_Tree_Node.tsv"
    )
    tree = args.result / "Species_Tree" / "SpeciesTree_rooted_node_labels.txt"
    for p in (dup, counts, tree):
        if not p.is_file():
            raise SystemExit(f"missing OrthoFinder output: {p}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    with dup.open(encoding="utf-8", errors="replace", newline="") as h:
        reader = csv.DictReader(h, delimiter="\t")
        fields = reader.fieldnames or []
        required = {"Orthogroup", "Species Tree Node", "Support", "Type"}
        if not required.issubset(fields):
            raise SystemExit(f"unexpected Duplications.tsv columns: {fields}")
        for row in reader:
            try:
                row["_support"] = float(row["Support"])
            except (TypeError, ValueError):
                row["_support"] = float("nan")
            rows.append(row)

    by_node = defaultdict(list)
    for row in rows:
        by_node[row["Species Tree Node"]].append(row)

    # Native OrthoFinder summary is retained alongside a recomputed family-level
    # table so the displayed values are directly auditable.
    native = {}
    with counts.open(encoding="utf-8", errors="replace", newline="") as h:
        for row in csv.DictReader(h, delimiter="\t"):
            native[row["Species Tree Node"]] = row

    out_rows = []
    for node in sorted(by_node):
        families = by_node[node]
        all_n = len(families)
        supported = sum(1 for r in families if r["_support"] >= args.support)
        nonterminal = sum(1 for r in families if r.get("Type") == "Non-Terminal")
        out_rows.append(
            {
                "species_tree_node": node,
                "duplication_events_all_rows": all_n,
                "duplication_events_support_ge_cutoff": supported,
                "nonterminal_events": nonterminal,
                "support_cutoff": args.support,
                "native_duplications_all": native.get(node, {}).get(
                    "Duplications (all)", ""
                ),
                "native_duplications_support_ge_50": native.get(node, {}).get(
                    "Duplications (50% support)", ""
                ),
            }
        )
    out = args.outdir / "copyAB_duplication_node_counts.tsv"
    with out.open("w", encoding="utf-8", newline="") as h:
        writer = (
            csv.DictWriter(
                h, fieldnames=list(out_rows[0]), delimiter="\t", lineterminator="\n"
            )
            if out_rows
            else None
        )
        if writer:
            writer.writeheader()
            writer.writerows(out_rows)
    (args.outdir / "copyAB_species_tree_rooted_node_labels.txt").write_text(
        tree.read_text(encoding="utf-8"), encoding="utf-8"
    )
    (args.outdir / "README.md").write_text(
        "# CopyAB duplication mapping\n\n"
        "Input is native OrthoFinder reconciliation for 14 terminals. PZI, PSP and PPA are each represented as an A1+B1 copyAB terminal.\n\n"
        "The table reports duplication rows mapped to each species-tree node and the number meeting the requested species-overlap support cutoff (not gene-tree bootstrap support). These are duplication-event counts, not direct WGD-event counts.\n",
        encoding="ascii",
    )
    print(
        f"[OK] Summarized {len(rows)} duplication rows across {len(out_rows)} nodes: {out}"
    )


if __name__ == "__main__":
    main()
