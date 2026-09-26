#!/usr/bin/env python3
"""Derive validated reduced Fig. 2 b/d backbones from the ASTRAL tree.

The script prunes the A and B subgenome histories independently.  It only
creates a species-labelled consensus tree when both pruned topologies agree,
so A and B leaves are never mechanically collapsed into one terminal.
"""

import argparse
import csv
from pathlib import Path

from fig2_phylo_newick import Node, parse_newick, serialize


SIX_A = {
    "pzi_A": "PZI",
    "psp_A": "PSP",
    "ppa_A": "PPA",
    "ppe": "PPE",
    "pcam": "PCAM",
    "gtr": "GTR",
}

SIX_B = {
    "pzi_B": "PZI",
    "psp_B": "PSP",
    "ppa_B": "PPA",
    "ppe": "PPE",
    "pcam": "PCAM",
    "gtr": "GTR",
}

FOUR_A = {key: value for key, value in SIX_A.items() if value not in {"PPE", "PCAM"}}
FOUR_B = {key: value for key, value in SIX_B.items() if value not in {"PPE", "PCAM"}}


def copy_pruned(node, retained):
    if not node.children:
        return Node(node.name, node.length, []) if node.name in retained else None
    children = []
    for child in node.children:
        pruned = copy_pruned(child, retained)
        if pruned is not None:
            children.append(pruned)
    if not children:
        return None
    if len(children) == 1:
        return children[0]
    return Node("", node.length, children)


def replace_labels(node, label_map):
    if not node.children:
        if node.name not in label_map:
            raise ValueError("Unexpected leaf in pruned tree: {}".format(node.name))
        return Node(label_map[node.name], node.length, [])
    return Node("", node.length, [replace_labels(child, label_map) for child in node.children])


def descendant_taxa(node):
    if not node.children:
        return frozenset([node.name])
    values = frozenset()
    for child in node.children:
        values = values.union(descendant_taxa(child))
    return values


def canonical_split(taxa, all_taxa):
    other = all_taxa.difference(taxa)
    if len(taxa) < 2 or len(other) < 2:
        return None
    left = tuple(sorted(taxa))
    right = tuple(sorted(other))
    return left if left < right else right


def internal_splits(node, all_taxa, result=None):
    if result is None:
        result = set()
    if node.children:
        split = canonical_split(descendant_taxa(node), all_taxa)
        if split is not None:
            result.add(split)
        for child in node.children:
            internal_splits(child, all_taxa, result)
    return result


def tree_text(tree):
    return serialize(tree) + ";\n"


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_pair(root, first_labels, second_labels, name, outdir):
    first = copy_pruned(root, set(first_labels))
    second = copy_pruned(root, set(second_labels))
    if first is None or second is None:
        raise ValueError("No retained terminals for {}".format(name))
    first_display = replace_labels(first, first_labels)
    second_display = replace_labels(second, second_labels)
    taxa = frozenset(first_labels.values())
    if descendant_taxa(first_display) != taxa or descendant_taxa(second_display) != taxa:
        raise ValueError("Retained taxon set failed for {}".format(name))
    first_splits = internal_splits(first_display, taxa)
    second_splits = internal_splits(second_display, taxa)
    if first_splits != second_splits:
        raise ValueError(
            "A and B reduced topologies disagree for {}: A={}, B={}".format(
                name, sorted(first_splits), sorted(second_splits)
            )
        )
    (outdir / (name + "_A_backbone.tre")).write_text(tree_text(first_display), encoding="ascii")
    (outdir / (name + "_B_backbone.tre")).write_text(tree_text(second_display), encoding="ascii")
    (outdir / (name + "_consensus_backbone.tre")).write_text(tree_text(first_display), encoding="ascii")
    return [{"backbone": name, "subgenome": "A", "splits": ";".join("|".join(split) for split in sorted(first_splits))},
            {"backbone": name, "subgenome": "B", "splits": ";".join("|".join(split) for split in sorted(second_splits))}]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--astral-tree", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    root = parse_newick(Path(args.astral_tree).read_text(encoding="ascii"))
    rows = []
    rows.extend(build_pair(root, SIX_A, SIX_B, "panel_b_six_taxa", outdir))
    rows.extend(build_pair(root, FOUR_A, FOUR_B, "panel_d_four_taxa", outdir))
    write_tsv(outdir / "backbone_topology_validation.tsv", ["backbone", "subgenome", "splits"], rows)
    (outdir / "README.txt").write_text(
        "A and B were pruned independently from the existing subgenome ASTRAL tree.\n"
        "Each species-labelled consensus is written only after exact agreement of the A and B unrooted split sets.\n"
        "Panel b uses the six-taxon consensus only as the scaffold for mapped WGD-family events.\n"
        "Panel d uses the four-taxon consensus only as the POInT display scaffold; GTR is an outgroup anchor, not a duplicated subgenome.\n",
        encoding="ascii",
    )
    print("[OK] Validated six- and four-taxon A/B-consensus backbones: {}".format(outdir))


if __name__ == "__main__":
    main()
