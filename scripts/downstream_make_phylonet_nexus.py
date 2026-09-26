#!/usr/bin/env python3
"""Build a PhyloNet InferNetwork_MP NEXUS file from single-copy gene trees.

The script is project glue only: IQ-TREE produces gene trees, and PhyloNet
performs network inference. This helper filters trees to the expected species
set, roots each gene tree with the configured outgroup where possible, falls
back to one outgroup taxon if needed, and writes a valid NEXUS command file.
"""

import argparse
import re
from pathlib import Path


class Node:
    def __init__(self, name="", length=None, support=None):
        self.name = name
        self.length = length
        self.support = support
        self.children = []
        self.parent = None

    def is_leaf(self):
        return not self.children


def strip_comments(text):
    return re.sub(r"\[[^\]]*\]", "", text)


def parse_label(text, idx):
    start = idx
    while idx < len(text) and text[idx] not in ",():;":
        idx += 1
    return text[start:idx].strip(), idx


def parse_length(text, idx):
    if idx < len(text) and text[idx] == ":":
        idx += 1
        start = idx
        while idx < len(text) and text[idx] not in ",();":
            idx += 1
        return text[start:idx].strip(), idx
    return None, idx


def parse_subtree(text, idx=0):
    if text[idx] == "(":
        idx += 1
        node = Node()
        while True:
            child, idx = parse_subtree(text, idx)
            child.parent = node
            node.children.append(child)
            if text[idx] == ",":
                idx += 1
                continue
            if text[idx] == ")":
                idx += 1
                break
            raise ValueError("unexpected character while parsing tree: %s" % text[idx])
        label, idx = parse_label(text, idx)
        if label:
            node.support = label
        node.length, idx = parse_length(text, idx)
        return node, idx

    label, idx = parse_label(text, idx)
    if not label:
        raise ValueError("empty leaf label")
    node = Node(name=label)
    node.length, idx = parse_length(text, idx)
    return node, idx


def parse_newick(text):
    text = strip_comments(text).strip()
    if "=" in text and text.lower().lstrip().startswith("tree"):
        text = text.split("=", 1)[1].strip()
    if text.endswith(";"):
        text = text[:-1]
    root, idx = parse_subtree(text, 0)
    if idx != len(text):
        trailing = text[idx:].strip()
        if trailing:
            raise ValueError("unexpected trailing text: %s" % trailing[:80])
    return root


def leaves(node):
    if node.is_leaf():
        return [node.name]
    out = []
    for child in node.children:
        out.extend(leaves(child))
    return out


def all_nodes(root):
    stack = [root]
    out = []
    while stack:
        node = stack.pop()
        out.append(node)
        stack.extend(node.children)
    return out


def leaf_set(node):
    return set(leaves(node))


def find_edge_for_split(root, target):
    all_taxa = leaf_set(root)
    for node in all_nodes(root):
        if node.parent is None:
            continue
        below = leaf_set(node)
        above = all_taxa - below
        if below == target or above == target:
            return node.parent, node, below, above
    return None


def find_leaf(root, name):
    for node in all_nodes(root):
        if node.is_leaf() and node.name == name:
            return node
    return None


def clone_away(current, previous=None):
    node = Node(name=current.name)
    neighbors = []
    if current.parent is not None:
        neighbors.append(current.parent)
    neighbors.extend(current.children)
    for neighbor in neighbors:
        if neighbor is previous:
            continue
        child = clone_away(neighbor, current)
        child.parent = node
        node.children.append(child)
    return node


def reroot_on_edge(a, b, b_side_is_target):
    if b_side_is_target:
        left = clone_away(b, a)
        right = clone_away(a, b)
    else:
        left = clone_away(a, b)
        right = clone_away(b, a)
    root = Node()
    left.parent = root
    right.parent = root
    root.children = [left, right]
    return root


def root_tree(root, outgroups, fallback_outgroup):
    split = find_edge_for_split(root, outgroups)
    if split:
        parent, child, below, _above = split
        return reroot_on_edge(parent, child, below == outgroups), "outgroup_clade"

    fallback = find_leaf(root, fallback_outgroup)
    if fallback is not None and fallback.parent is not None:
        return reroot_on_edge(fallback.parent, fallback, True), "fallback_outgroup"

    return root, "unrooted_kept"


def to_topology_newick(node, is_root=True):
    if node.is_leaf():
        text = node.name
    else:
        text = "(" + ",".join(to_topology_newick(child, False) for child in node.children) + ")"
    return text + (";" if is_root else "")


def read_gene_trees(path):
    trees = []
    for line in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        trees.append(line)
    return trees


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-trees", required=True)
    parser.add_argument("--out-nexus", required=True)
    parser.add_argument("--stats", required=True)
    parser.add_argument("--species", required=True, help="Comma-separated required species labels.")
    parser.add_argument("--outgroups", default="fve,rch", help="Comma-separated outgroup labels.")
    parser.add_argument("--fallback-outgroup", default="fve")
    parser.add_argument("--reticulations", type=int, default=1)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--max-trees", type=int, default=0, help="Use all trees when 0.")
    args = parser.parse_args()

    required = {x.strip() for x in args.species.split(",") if x.strip()}
    outgroups = {x.strip() for x in args.outgroups.split(",") if x.strip()}
    if not outgroups:
        raise SystemExit("ERROR: at least one outgroup is required")
    if args.fallback_outgroup not in outgroups:
        raise SystemExit("ERROR: fallback outgroup must be one of --outgroups")

    out_nexus = Path(args.out_nexus)
    stats_path = Path(args.stats)
    out_nexus.parent.mkdir(parents=True, exist_ok=True)
    stats_path.parent.mkdir(parents=True, exist_ok=True)

    records = []
    kept = []
    raw_trees = read_gene_trees(args.gene_trees)
    for idx, tree_text in enumerate(raw_trees, start=1):
        try:
            root = parse_newick(tree_text)
            taxa = set(leaves(root))
        except Exception as exc:
            records.append((idx, "skip_parse_error", str(exc), ""))
            continue
        missing = sorted(required - taxa)
        extra = sorted(taxa - required)
        if missing:
            records.append((idx, "skip_missing_required_taxa", ",".join(missing), ",".join(extra)))
            continue
        rooted, status = root_tree(root, outgroups, args.fallback_outgroup)
        tree_id = "gt%06d" % (len(kept) + 1)
        kept.append((tree_id, to_topology_newick(rooted)))
        records.append((idx, "used_" + status, "", ",".join(extra)))
        if args.max_trees and len(kept) >= args.max_trees:
            break

    if not kept:
        raise SystemExit("ERROR: no gene trees were retained for PhyloNet")

    tree_ids = ",".join(tree_id for tree_id, _tree in kept)
    with out_nexus.open("w", encoding="utf-8", newline="\n") as out:
        out.write("#NEXUS\n\n")
        out.write("BEGIN TREES;\n")
        for tree_id, tree in kept:
            out.write("Tree %s = [&R] %s\n" % (tree_id, tree))
        out.write("END;\n\n")
        out.write("BEGIN PHYLONET;\n")
        out.write("InferNetwork_MP (%s) %d -pl %d;\n" % (tree_ids, args.reticulations, args.threads))
        out.write("END;\n")

    with stats_path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("input_index\tstatus\tdetail\textra_taxa\n")
        for row in records:
            out.write("%s\t%s\t%s\t%s\n" % row)

    status_counts = {}
    for _idx, status, _detail, _extra in records:
        status_counts[status] = status_counts.get(status, 0) + 1
    print("input_trees=%d retained=%d out_nexus=%s stats=%s" % (
        len(raw_trees), len(kept), out_nexus, stats_path
    ))
    for status in sorted(status_counts):
        print("%s=%d" % (status, status_counts[status]))


if __name__ == "__main__":
    main()
