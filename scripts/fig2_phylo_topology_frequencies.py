#!/usr/bin/env python3
"""Count decisive gene-tree quartet topologies and gene-bootstrap uncertainty."""

import argparse
import csv
import random
from collections import Counter
from pathlib import Path
from fig2_phylo_newick import parse_newick, quartet_topology


def percentile(values, probability):
    ordered = sorted(values)
    index = int(round((len(ordered) - 1) * probability))
    return ordered[index]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trees", required=True)
    parser.add_argument("--quartets", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260813)
    args = parser.parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=False)

    trees = [parse_newick(line) for line in Path(args.trees).read_text(encoding="ascii").splitlines() if line.strip()]
    with Path(args.quartets).open("r", encoding="utf-8", newline="") as handle:
        quartets = list(csv.DictReader(handle, delimiter="\t"))
    rng = random.Random(args.seed)
    detail_rows = []
    summary_rows = []

    for row in quartets:
        taxa = [row["taxon1"], row["taxon2"], row["taxon3"], row["taxon4"]]
        calls = []
        for index, tree in enumerate(trees, start=1):
            call = quartet_topology(tree, taxa)
            detail_rows.append({"contrast": row["contrast"], "gene_index": str(index), "topology": call})
            if call not in ["missing", "unresolved"]:
                calls.append(call)
        counts = Counter(calls)
        if not calls:
            raise SystemExit("ERROR: no decisive genes for quartet %s" % row["contrast"])
        topologies = [
            taxa[0] + "+" + taxa[1] + "|" + taxa[2] + "+" + taxa[3],
            taxa[0] + "+" + taxa[2] + "|" + taxa[1] + "+" + taxa[3],
            taxa[0] + "+" + taxa[3] + "|" + taxa[1] + "+" + taxa[2],
        ]
        boot = dict((topology, []) for topology in topologies)
        for _ in range(args.replicates):
            sampled = [calls[rng.randrange(len(calls))] for __ in range(len(calls))]
            sampled_counts = Counter(sampled)
            for topology in topologies:
                boot[topology].append(float(sampled_counts[topology]) / len(sampled))
        for topology in topologies:
            frequency = float(counts[topology]) / len(calls)
            summary_rows.append({
                "contrast": row["contrast"],
                "topology": topology,
                "decisive_genes": str(len(calls)),
                "gene_count": str(counts[topology]),
                "frequency": "%.6f" % frequency,
                "gene_bootstrap_ci_low": "%.6f" % percentile(boot[topology], 0.025),
                "gene_bootstrap_ci_high": "%.6f" % percentile(boot[topology], 0.975),
                "replicates": str(args.replicates),
                "seed": str(args.seed),
            })

    with (out_dir / "gene_topology_calls.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["contrast", "gene_index", "topology"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(detail_rows)
    fields = ["contrast", "topology", "decisive_genes", "gene_count", "frequency", "gene_bootstrap_ci_low", "gene_bootstrap_ci_high", "replicates", "seed"]
    with (out_dir / "topology_frequencies.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary_rows)
    (out_dir / "COMPLETE").write_text("Topology frequencies completed.\n", encoding="ascii")


if __name__ == "__main__":
    main()
