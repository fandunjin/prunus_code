#!/usr/bin/env python3
"""Outgroup affinity check for the five-species inflorescence run.

For every subfamily group in every gene tree, test where the single-copy
PCAM/PPE orthologs attach relative to the PZI/PSP/PPA A- and B-track clades.
A pre-WGD ortholog is expected to attach outside (above) the combined trio
clade; attachment inside the A or B clade would indicate subgenome-specific
affinity of the non-racemose lineage.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies5_20260912"

PHASE_RE = re.compile(r"^(pzi|psp|ppa|pcam|ppe)_hap(A1|A2|B1|B2|A|B|S)_")
TRACK = {"A1": "A", "A2": "A", "B1": "B", "B2": "B", "A": "A", "B": "B", "S": "S", "S1": "S"}
TRIO = {"pzi", "psp", "ppa"}


class Node:
    __slots__ = ("children", "length", "label")

    def __init__(self):
        self.children: list["Node"] = []
        self.length = 0.0
        self.label = ""


re_label = re.compile(r"[^:,()]+")
re_len = re.compile(r":\s*([0-9.eE+-]+)")


def parse_newick(text: str) -> Node:
    s = text.strip().rstrip(";")
    pos = 0

    def parse() -> Node:
        nonlocal pos
        node = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                node.children.append(parse())
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
            m = re_label.match(s, pos)
            if m:
                node.label = m.group(0)
                pos = m.end()
        else:
            m = re_label.match(s, pos)
            if m:
                node.label = m.group(0)
                pos = m.end()
        m = re_len.match(s, pos)
        if m:
            node.length = float(m.group(1))
            pos = m.end()
        return node

    return parse()


def leaves_of(node: Node) -> set[str]:
    if not node.children:
        return {node.label}
    out: set[str] = set()
    for c in node.children:
        out |= leaves_of(c)
    return out


def build_parent_map(root: Node) -> dict[int, Node]:
    parents: dict[int, Node] = {}

    def walk(n: Node) -> None:
        for c in n.children:
            parents[id(c)] = n
            walk(c)

    walk(root)
    return parents


def depths(root: Node) -> dict[int, float]:
    out: dict[int, float] = {}

    def walk(n: Node, d: float) -> None:
        out[id(n)] = d
        for c in n.children:
            walk(c, d + c.length)

    walk(root, 0.0)
    return out


def distance(a: str, b: str, leaf_index, parents, depth) -> float:
    m = mrca([a, b], leaf_index, parents)
    return depth[id(leaf_index[a])] + depth[id(leaf_index[b])] - 2 * depth[id(m)]


def mrca(tips: list[str], leaf_index: dict[str, Node], parents: dict[int, Node]) -> Node:
    chain: dict[int, int] = {}
    cur = leaf_index[tips[0]]
    step = 0
    while cur is not None:
        chain[id(cur)] = step
        cur = parents.get(id(cur))
        step += 1
    best = None
    best_cost = None
    for tip in tips[1:]:
        cur = leaf_index[tip]
        step = 0
        while cur is not None:
            if id(cur) in chain:
                cost = chain[id(cur)] + step
                if best_cost is None or cost < best_cost:
                    best, best_cost = cur, cost
                break
            cur = parents.get(id(cur))
            step += 1
    return best  # type: ignore[return-value]


def main() -> int:
    out = DEFAULT_OUT
    placements_path = out / "trees" / "tree_placements_crossspecies.tsv"
    rows = list(csv.DictReader(placements_path.open(encoding="utf-8"), delimiter="\t"))
    group_by_tip: dict[str, tuple[str, str, str]] = {}  # tip -> (family, group, species)
    for r in rows:
        if not r.get("subfamily_group"):
            continue
        group_by_tip[r["target_id"]] = (r["family"], r["subfamily_group"], r["species"])

    result_rows: list[dict[str, object]] = []
    for treefile in sorted((out / "trees").glob("*.treefile")):
        family = treefile.stem
        root = parse_newick(treefile.read_text(encoding="utf-8", errors="ignore"))
        leaf_index: dict[str, Node] = {}

        def collect(n: Node) -> None:
            if not n.children:
                leaf_index[n.label] = n
                return
            for c in n.children:
                collect(c)

        collect(root)
        parents = build_parent_map(root)
        depth = depths(root)
        groups: dict[str, dict[str, list[str]]] = {}
        for tip, (fam, group, species) in group_by_tip.items():
            if fam != family or tip not in leaf_index:
                continue
            track = "single" if species not in TRIO else TRACK.get(PHASE_RE.match(tip).group(2), "?")  # type: ignore[union-attr]
            groups.setdefault(group, {"A": [], "B": [], "single": []})[track].append(tip)
        for group, members in sorted(groups.items()):
            trio_tips = members["A"] + members["B"]
            singles = members["single"]
            if len(trio_tips) < 2 or not singles:
                continue
            ab_node = mrca(trio_tips, leaf_index, parents)
            ab_set = leaves_of(ab_node)
            a_node = mrca(members["A"], leaf_index, parents) if members["A"] else None
            b_node = mrca(members["B"], leaf_index, parents) if members["B"] else None
            a_set = leaves_of(a_node) if a_node is not None else set()
            b_set = leaves_of(b_node) if b_node is not None else set()
            for s in singles:
                s_species = group_by_tip[s][2]
                m = mrca([s] + trio_tips, leaf_index, parents)
                m_set = leaves_of(m)
                if m_set == ab_set:
                    cls = "at_AB_node"
                elif m_set >= ab_set:
                    cls = "outside_trio_clade"
                else:
                    if a_set and m_set <= a_set:
                        cls = "inside_A_clade"
                    elif b_set and m_set <= b_set:
                        cls = "inside_B_clade"
                    else:
                        cls = "inside_trio_internal"
                ranked = sorted(
                    trio_tips,
                    key=lambda t: distance(s, t, leaf_index, parents, depth),
                )[:3]
                nearest_tracks = "".join(
                    sorted({TRACK.get(PHASE_RE.match(t).group(2), "?") for t in ranked})  # type: ignore[union-attr]
                )
                nearest_species = ";".join(sorted({PHASE_RE.match(t).group(1) for t in ranked}))  # type: ignore[union-attr,arg-type]
                result_rows.append({
                    "family": family, "subfamily_group": group,
                    "single_tip": s, "single_species": s_species,
                    "n_trio_A": len(members["A"]), "n_trio_B": len(members["B"]),
                    "affinity_class": cls,
                    "nearest_tracks": nearest_tracks,
                    "nearest_species": nearest_species,
                    "nearest_distance": f"{distance(s, ranked[0], leaf_index, parents, depth):.5f}",
                })
    outdir = out / "stats"
    outdir.mkdir(parents=True, exist_ok=True)
    if result_rows:
        with (outdir / "outgroup_affinity_check.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(result_rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(result_rows)
    from collections import Counter
    summary = Counter((r["family"], r["subfamily_group"], r["affinity_class"]) for r in result_rows)
    print("family  group  class  n")
    for k, v in sorted(summary.items()):
        print("  ".join(str(x) for x in k), v)
    print()
    print("nearest-track summary (per single tip, top-3 nearest trio tips):")
    nt = Counter((r["family"], r["subfamily_group"], r["nearest_tracks"]) for r in result_rows)
    for k, v in sorted(nt.items()):
        print("  ".join(str(x) for x in k), v)
    return 0


if __name__ == "__main__":
    sys.exit(main())
