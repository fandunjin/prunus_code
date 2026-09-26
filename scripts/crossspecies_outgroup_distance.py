#!/usr/bin/env python3
"""A-vs-B outgroup distance test (five-species run).

For every homeolog anchor component in PZI/PSP/PPA, compare the patristic
distance of its A copies and B copies to the orthologous PCAM (corymbose) or
PPE (solitary-flowered) copies in the same subfamily clade.  Exact sign tests
per family and outgroup species; exploratory because the nearest outgroup tip
is a sequence-similarity proxy for orthology.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

from scipy.stats import binomtest

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "analysis" / "inflorescence_family_pipeline_crossspecies5_20260912"
PHASE_RE = re.compile(r"^(pzi|psp|ppa|pcam|ppe)_hap(A1|A2|B1|B2|A|B|S)_")
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

    def sub() -> Node:
        nonlocal pos
        n = Node()
        if pos < len(s) and s[pos] == "(":
            pos += 1
            while True:
                n.children.append(sub())
                if pos < len(s) and s[pos] == ",":
                    pos += 1
                    continue
                if pos < len(s) and s[pos] == ")":
                    pos += 1
                    break
            m = re_label.match(s, pos)
            if m:
                n.label = m.group(0)
                pos = m.end()
        else:
            m = re_label.match(s, pos)
            if m:
                n.label = m.group(0)
                pos = m.end()
        m = re_len.match(s, pos)
        if m:
            n.length = float(m.group(1))
            pos = m.end()
        return n

    return sub()


def parents_and_depth(root: Node):
    parents: dict[int, Node] = {}
    depth: dict[int, float] = {}

    def walk(n: Node, d: float) -> None:
        depth[id(n)] = d
        for c in n.children:
            parents[id(c)] = n
            walk(c, d + c.length)

    walk(root, 0.0)
    return parents, depth


def mrca(tips, leaf_index, parents):
    chain = {}
    cur = leaf_index[tips[0]]
    step = 0
    while cur is not None:
        chain[id(cur)] = step
        cur = parents.get(id(cur))
        step += 1
    best = best_cost = None
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
    return best


def distance(a, b, leaf_index, parents, depth):
    m = mrca([a, b], leaf_index, parents)
    return depth[id(leaf_index[a])] + depth[id(leaf_index[b])] - 2 * depth[id(m)]


def main() -> int:
    out = DEFAULT_OUT
    placements = list(csv.DictReader((out / "trees" / "tree_placements_crossspecies.tsv").open(encoding="utf-8"), delimiter="\t"))
    group_by_tip = {}
    tip_by_gene = {}
    for r in placements:
        if r.get("subfamily_group"):
            group_by_tip[r["target_id"]] = (r["family"], r["subfamily_group"], r["species"])
        if r.get("gene_key"):
            tip_by_gene[r["gene_key"]] = r["target_id"]

    components = list(csv.DictReader((out / "copy_fate" / "components_crossspecies.tsv").open(encoding="utf-8"), delimiter="\t"))
    rows = []
    for treefile in sorted((out / "trees").glob("*.treefile")):
        family = treefile.stem
        root = parse_newick(treefile.read_text(encoding="utf-8", errors="ignore"))
        leaf_index: dict[str, Node] = {}

        def collect(n):
            if not n.children:
                leaf_index[n.label] = n
                return
            for c in n.children:
                collect(c)

        collect(root)
        parents, depth = parents_and_depth(root)
        group_out: dict[str, list[str]] = defaultdict(list)
        for tip, (fam, grp, sp) in group_by_tip.items():
            if fam == family and sp in {"pcam", "ppe"} and tip in leaf_index:
                group_out[f"{sp}|{grp}"].append(tip)
        for comp in components:
            if family not in comp["families"].split(";"):
                continue
            members = [tip_by_gene[m] for m in comp["members"].split(";") if m in tip_by_gene]
            members = [m for m in members if m in leaf_index]
            a_members = [m for m in members if PHASE_RE.match(m) and PHASE_RE.match(m).group(2)[0] == "A"]
            b_members = [m for m in members if PHASE_RE.match(m) and PHASE_RE.match(m).group(2)[0] == "B"]
            if not a_members or not b_members:
                continue
            groups = {group_by_tip[m][1] for m in members if m in group_by_tip}
            for sp in ("pcam", "ppe"):
                outs = [t for grp in groups for t in group_out.get(f"{sp}|{grp}", [])]
                if not outs:
                    continue
                d_a = min(distance(a, o, leaf_index, parents, depth) for a in a_members for o in outs)
                d_b = min(distance(b, o, leaf_index, parents, depth) for b in b_members for o in outs)
                rows.append({
                    "component_id": comp["component_id"], "trio_species": comp["species"],
                    "family": family, "subfamily_group": ";".join(sorted(groups)),
                    "outgroup": sp, "d_A": f"{d_a:.5f}", "d_B": f"{d_b:.5f}",
                    "closer": "A" if d_a < d_b else ("B" if d_b < d_a else "tie"),
                })
    outdir = out / "stats"
    if rows:
        with (outdir / "outgroup_distance_test.tsv").open("w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
            w.writeheader()
            w.writerows(rows)

    print("=== per family x outgroup: closer to outgroup (A vs B) ===")
    summary = defaultdict(Counter)
    for r in rows:
        summary[(r["family"], r["outgroup"])][r["closer"]] += 1
    for key, cnt in sorted(summary.items()):
        n = cnt["A"] + cnt["B"]
        p = binomtest(cnt["A"], n, 0.5).pvalue if n else float("nan")
        print(f"{key[0]:10s} vs {key[1]:5s}  A-closer={cnt['A']:3d}  B-closer={cnt['B']:3d}  ties={cnt['tie']}  sign p={p:.4g}")
    print()
    print("=== per trio species x outgroup ===")
    summary2 = defaultdict(Counter)
    for r in rows:
        summary2[(r["trio_species"], r["outgroup"])][r["closer"]] += 1
    for key, cnt in sorted(summary2.items()):
        n = cnt["A"] + cnt["B"]
        p = binomtest(cnt["A"], n, 0.5).pvalue if n else float("nan")
        print(f"{key[0]:4s} vs {key[1]:5s}  A-closer={cnt['A']:3d}  B-closer={cnt['B']:3d}  ties={cnt['tie']}  sign p={p:.4g}")
    print()
    print("=== per subfamily group (pooled outgroups) ===")
    summary3 = defaultdict(Counter)
    for r in rows:
        summary3[r["subfamily_group"]][r["closer"]] += 1
    for key, cnt in sorted(summary3.items()):
        n = cnt["A"] + cnt["B"]
        p = binomtest(cnt["A"], n, 0.5).pvalue if n else float("nan")
        print(f"{key:25s}  A-closer={cnt['A']:3d}  B-closer={cnt['B']:3d}  ties={cnt['tie']}  sign p={p:.4g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
