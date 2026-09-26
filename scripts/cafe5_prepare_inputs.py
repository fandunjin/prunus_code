"""Build CAFE5 inputs for the racemose Prunus dataset (7 species).

1. Prune the MCMCtree dated species14 tree to pzi/psp/ppa/pcam/ppe/gtr/rch
   (collapse degree-2 nodes by summing branch lengths; verify ultrametric).
2. Build the gene-family count table from Orthogroups.txt (A+B tracks summed).
Input:  analysisphylogenetics/02_orthofinder/main/Results_Aug25/Orthogroups/Orthogroups.txt
       0903/analysisphylogenetics/09_mcmctree_species14_plant_outgroup_fossil55_35/03_mcmctree_chain1/FigTree.tre
Output: analysis/cafe5/main/{mcmctree_dated_ultrametric_tree.tre, gene_family_counts.tsv}
"""
from __future__ import annotations

import csv
import io
import re
import sys
from collections import Counter
from pathlib import Path

sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")

ROOT = Path(r"${PROJECT_ROOT}")
FIG = ROOT / "0903" / "analysisphylogenetics" / "09_mcmctree_species14_plant_outgroup_fossil55_35" / "03_mcmctree_chain1" / "FigTree.tre"
OG = ROOT / "analysisphylogenetics" / "02_orthofinder" / "main" / "Results_Aug25" / "Orthogroups" / "Orthogroups.txt"
OUT = ROOT / "analysis" / "cafe5" / "main"
OUT.mkdir(parents=True, exist_ok=True)

KEEP = {"pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"}
TRACK_TO_SPECIES = {
    "pzi_A": "pzi", "pzi_B": "pzi",
    "psp_A": "psp", "psp_B": "psp",
    "ppa_A": "ppa", "ppa_B": "ppa",
    "pcam": "pcam", "ppe": "ppe", "gtr": "gtr", "rch": "rch",
    "pmal": None, "pav": None, "pspe": None, "psal": None, "par": None, "pmu": None, "pdu": None,
}


# ---------------------------------------------------------------- tree parsing
class Node:
    __slots__ = ("name", "length", "children")

    def __init__(self, name="", length=None, children=None):
        self.name = name
        self.length = length
        self.children = children or []


def parse_newick(text: str):
    text = re.sub(r"\[[^\]]*\]", "", text)  # strip [&...] annotations
    pos = 0

    def parse():
        nonlocal pos
        node = Node()
        if text[pos] == "(":
            pos += 1
            while True:
                child = parse()
                node.children.append(child)
                if text[pos] == ",":
                    pos += 1
                elif text[pos] == ")":
                    pos += 1
                    break
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        node.name = text[start:pos].strip()
        if pos < len(text) and text[pos] == ":":
            pos += 1
            start = pos
            while pos < len(text) and text[pos] not in ",();":
                pos += 1
            node.length = float(text[start:pos])
        return node

    root = parse()
    return root


def to_newick(node: Node) -> str:
    if node.children:
        inner = ",".join(to_newick(c) for c in node.children)
        s = f"({inner}){node.name}"
    else:
        s = node.name
    if node.length is not None:
        s += f":{node.length:.6f}"
    return s


def prune(node: Node) -> Node | None:
    """Recursively prune; collapse degree-2 nodes summing branch lengths."""
    if not node.children:
        return node if node.name in KEEP else None
    kept = []
    for child in node.children:
        sub = prune(child)
        if sub is not None:
            kept.append(sub)
    if not kept:
        return None
    if len(kept) == 1 and node.name not in KEEP:
        # collapse: fold this node's length (if any) onto the single child
        child = kept[0]
        if node.length is not None:
            child.length = (child.length or 0.0) + node.length
        return child
    node.children = kept
    return node


def tip_depths(node: Node, acc=0.0):
    acc += node.length or 0.0
    if not node.children:
        yield node.name, acc
    else:
        for c in node.children:
            yield from tip_depths(c, acc)


raw = FIG.read_text(encoding="utf-8", errors="replace")
m = re.search(r"UTREE\s*1\s*=\s*([^;]+);", raw)
if not m:
    raise SystemExit("UTREE not found in FigTree.tre")
tree = parse_newick(m.group(1))
pruned = prune(tree)
nwk = to_newick(pruned) + ";"
(OUT / "mcmctree_dated_ultrametric_tree.tre").write_text(nwk + "\n", encoding="utf-8")

depths = list(tip_depths(pruned))
vals = [d for _n, d in depths]
print("pruned tree:", nwk)
print("tip depths (MCMCtree units):", {n: round(d, 6) for n, d in depths})
print("ultrametric spread:", round(max(vals) - min(vals), 8))
print("root age (Ma):", round(max(vals) * 100, 2))


# ---------------------------------------------------------------- gene counts
prefix = Counter()
rows = []
with OG.open(encoding="utf-8", errors="ignore") as fh:
    for line in fh:
        line = line.strip()
        if not line or ": " not in line:
            continue
        fam, members = line.split(": ", 1)
        counts = Counter()
        for mem in members.split():
            tr = mem.split("__", 1)[0]
            sp = TRACK_TO_SPECIES.get(tr, None)
            if sp is not None:
                counts[sp] += 1
        for tr in members.split():
            prefix[tr.split("__", 1)[0]] += 1
        vals7 = {sp: counts.get(sp, 0) for sp in ("pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch")}
        if sum(1 for v in vals7.values() if v > 0) < 2:
            continue
        if max(vals7.values()) > 100:
            continue
        rows.append((fam, vals7))

with (OUT / "gene_family_counts.tsv").open("w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n")
    w.writerow(["Desc", "Family ID", "pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch"])
    for fam, vals7 in rows:
        w.writerow(["(null)", fam] + [vals7[sp] for sp in ("pzi", "psp", "ppa", "pcam", "ppe", "gtr", "rch")])

print("families kept:", len(rows), "of 41780 (filter: >=2 taxa, max<=100)")
sizes = Counter()
for _f, v in rows:
    sizes[sum(v.values())] += 1
print("largest retained family size:", max(sizes))
