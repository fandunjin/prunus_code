#!/usr/bin/env python3
"""Summarize PUG mapped nodes at the pre-specified support cutoff."""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pug-root", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--cutoff", type=float, default=95.0)
    parser.add_argument("--mark-complete", action="store_true")
    args = parser.parse_args()

    rows = []
    counts = defaultdict(set)
    source_names = {
        "shared_candidate": "fig2b_shared_candidate_Paralog_Pairs_Nodes_Bootstraps.txt",
        "pzi_ab": "fig2b_pzi_ab_Paralog_Pairs_Nodes_Bootstraps.txt",
        "psp_ab": "fig2b_psp_ab_Paralog_Pairs_Nodes_Bootstraps.txt",
        "ppa_ab": "fig2b_ppa_ab_Paralog_Pairs_Nodes_Bootstraps.txt",
    }
    for event, filename in source_names.items():
        path = args.pug_root / event / filename
        if not path.is_file():
            raise SystemExit(f"Missing PUG result: {path}")
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            expected = {"Orthogroup", "Node", "Bootstrap", "Paralog1", "Paralog2"}
            if not expected.issubset(set(reader.fieldnames or [])):
                raise SystemExit(f"Unexpected PUG columns in {path}: {reader.fieldnames}")
            for record in reader:
                try:
                    support = float(record["Bootstrap"])
                except (TypeError, ValueError):
                    continue
                if support < args.cutoff:
                    continue
                og = Path(record["Orthogroup"]).stem
                node = record["Node"]
                key = (event, node)
                counts[key].add(og)
                rows.append({
                    "event": event,
                    "orthogroup": og,
                    "node": node,
                    "support": f"{support:.4f}",
                    "paralog1": record["Paralog1"],
                    "paralog2": record["Paralog2"],
                })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    mapped = args.out.with_name(args.out.stem + "_mapped_rows.tsv")
    with mapped.open("w", encoding="utf-8", newline="") as handle:
        fields = ["event", "orthogroup", "node", "support", "paralog1", "paralog2"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    with args.out.open("w", encoding="utf-8", newline="") as handle:
        fields = ["event", "node", "support_cutoff", "gene_tree_count"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for (event, node), orthogroups in sorted(counts.items()):
            writer.writerow({"event": event, "node": node, "support_cutoff": args.cutoff, "gene_tree_count": len(orthogroups)})
    if args.mark_complete:
        (args.pug_root / "software_and_parameters.txt").write_text(
            "support_cutoff=95\n"
            "support_note=FastTree_default local support in 0-1; PUG Bayesian mode used only for x100 scaling, not posterior interpretation\n"
            "outgroup=rch\n"
            "counting=default PUG unique LCA per gene tree; no --all_pairs\n",
            encoding="ascii",
        )
        (args.pug_root / "COMPLETE").touch()
    print(f"[OK] PUG rows >= {args.cutoff:g}: {len(rows)}; node summaries: {len(counts)}")


if __name__ == "__main__":
    main()
