#!/usr/bin/env python3
"""Map corrected A/B homeolog LCAs from single-copy subgenome gene trees.

This is deliberately not a replacement for a full gene-family reconciliation.
Each accepted orthogroup has one representative at each labelled subgenome tip.
The script therefore maps the LCA of a *predefined* A/B homeolog pair from an
IQ-TREE gene tree onto the collapsed biological-species backbone.  It reports
the UFBoot value carried by that LCA and retains both all testable mappings and
the preregistered UFBoot >= 95 subset.
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Set


PAIR_DEFINITIONS = {
    "PZI": ("pzi_A", "pzi_B"),
    "PSP": ("psp_A", "psp_B"),
    "PPA": ("ppa_A", "ppa_B"),
}

COLLAPSED_TREE = (
    "(gtr,(rch,((ppa,(psp,pzi)),"
    "((pmal,(pav,(pspe,pcam))),((psal,(par,pmu)),(pdu,ppe))))));"
)


@dataclass
class Node:
    name: str = ""
    children: List["Node"] = field(default_factory=list)
    parent: Optional["Node"] = None
    leaves: Set[str] = field(default_factory=set)

    @property
    def is_leaf(self) -> bool:
        return not self.children


class NewickParser:
    def __init__(self, text: str):
        self.text = text.strip()
        self.i = 0

    def parse(self) -> Node:
        node = self._subtree()
        self._skip()
        if self.i < len(self.text) and self.text[self.i] == ";":
            self.i += 1
        self._annotate(node)
        return node

    def _skip(self) -> None:
        while self.i < len(self.text) and self.text[self.i].isspace():
            self.i += 1

    def _token(self, stops: str) -> str:
        start = self.i
        while self.i < len(self.text) and self.text[self.i] not in stops:
            self.i += 1
        return self.text[start:self.i].strip()

    def _branch(self) -> None:
        if self.i < len(self.text) and self.text[self.i] == ":":
            self.i += 1
            self._token(",();")

    def _subtree(self) -> Node:
        self._skip()
        if self.i >= len(self.text):
            raise ValueError("Unexpected end of Newick string")
        if self.text[self.i] == "(":
            self.i += 1
            children = [self._subtree()]
            self._skip()
            while self.i < len(self.text) and self.text[self.i] == ",":
                self.i += 1
                children.append(self._subtree())
                self._skip()
            if self.i >= len(self.text) or self.text[self.i] != ")":
                raise ValueError("Malformed Newick: missing ')'")
            self.i += 1
            label = self._token(":,();")
            self._branch()
            node = Node(name=label, children=children)
            for child in children:
                child.parent = node
            return node
        label = self._token(":,();")
        if not label:
            raise ValueError("Malformed Newick: unnamed leaf")
        self._branch()
        return Node(name=label)

    def _annotate(self, node: Node) -> Set[str]:
        if node.is_leaf:
            node.leaves = {node.name}
        else:
            node.leaves = set().union(*(self._annotate(child) for child in node.children))
        return node.leaves


def lca(root: Node, taxa: Iterable[str]) -> Optional[Node]:
    required = set(taxa)
    if not required.issubset(root.leaves):
        return None
    node = root
    while True:
        containing = [child for child in node.children if required.issubset(child.leaves)]
        if len(containing) != 1:
            return node
        node = containing[0]


def biological_taxon(label: str) -> str:
    return re.sub(r"_[AB]$", "", label)


def ufboot(label: str) -> Optional[float]:
    """IQ-TREE labels are commonly SH-aLRT/UFBoot, e.g. 83.1/91."""
    if not label:
        return None
    token = label.rsplit("/", 1)[-1]
    try:
        value = float(token)
    except ValueError:
        return None
    return value if 0.0 <= value <= 100.0 else None


def node_labels(root: Node) -> Dict[int, str]:
    """Name all collapsed-tree internal nodes reproducibly by descendant tips."""
    labels: Dict[int, str] = {}
    preferred = {
        frozenset({"ppa", "psp", "pzi"}): "Padus-Laurocerasus_crown",
        frozenset({"psp", "pzi"}): "PSP-PZI_crown",
        frozenset({"pmal", "pav", "pspe", "pcam"}): "Cerasus_crown",
        frozenset({"pav", "pspe", "pcam"}): "PAV-PSPE-PCAM_crown",
        frozenset({"pspe", "pcam"}): "PSPE-PCAM_crown",
        frozenset({"psal", "par", "pmu", "pdu", "ppe"}): "Amygdalus-Prunophora_crown",
        frozenset({"psal", "par", "pmu"}): "PSAL-PAR-PMU_crown",
        frozenset({"par", "pmu"}): "PAR-PMU_crown",
        frozenset({"pdu", "ppe"}): "PDU-PPE_crown",
    }

    def visit(node: Node) -> None:
        if node.is_leaf:
            labels[id(node)] = f"{node.name.upper()}_tip"
            return
        tipset = frozenset(node.leaves)
        if node.parent is None:
            labels[id(node)] = "Rosaceae_root"
        elif tipset == frozenset(root.leaves - {"gtr"}):
            labels[id(node)] = "Prunus-Rosa_crown"
        elif tipset == frozenset(root.leaves - {"gtr", "rch"}):
            labels[id(node)] = "Prunus_crown"
        else:
            labels[id(node)] = preferred.get(tipset, "_".join(sorted(tipset)) + "_crown")
        for child in node.children:
            visit(child)

    visit(root)
    return labels


def write_tsv(path: Path, rows: List[dict], fieldnames: List[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gene-trees", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--ufboot-cutoff", type=float, default=95.0)
    args = parser.parse_args()

    tree_paths = sorted(args.gene_trees.glob("*.treefile"))
    if not tree_paths:
        raise SystemExit(f"No '*.treefile' files under {args.gene_trees}")

    species_root = NewickParser(COLLAPSED_TREE).parse()
    species_node_names = node_labels(species_root)
    records: List[dict] = []
    failed: List[dict] = []
    for path in tree_paths:
        try:
            gene_root = NewickParser(path.read_text(encoding="utf-8")).parse()
        except Exception as exc:  # source record is retained for audit
            failed.append({"orthogroup": path.stem, "reason": str(exc)})
            continue
        for taxon, pair in PAIR_DEFINITIONS.items():
            duplication_lca = lca(gene_root, pair)
            if duplication_lca is None:
                records.append({
                    "orthogroup": path.stem, "taxon": taxon, "tip_a": pair[0], "tip_b": pair[1],
                    "status": "missing_pair", "gene_tree_lca_ufboot": "", "mapped_node": "",
                    "mapped_descendant_taxa": "", "high_support": "no",
                })
                continue
            collapsed_descendants = {biological_taxon(x) for x in duplication_lca.leaves}
            backbone_node = lca(species_root, collapsed_descendants)
            if backbone_node is None:
                raise RuntimeError(f"Cannot map {path.name} {taxon}: {collapsed_descendants}")
            support = ufboot(duplication_lca.name)
            records.append({
                "orthogroup": path.stem,
                "taxon": taxon,
                "tip_a": pair[0],
                "tip_b": pair[1],
                "status": "mapped" if support is not None else "mapped_no_ufboot",
                "gene_tree_lca_ufboot": "" if support is None else f"{support:.1f}",
                "mapped_node": species_node_names[id(backbone_node)],
                "mapped_descendant_taxa": ",".join(sorted(backbone_node.leaves)),
                "high_support": "yes" if support is not None and support >= args.ufboot_cutoff else "no",
            })

    fieldnames = [
        "orthogroup", "taxon", "tip_a", "tip_b", "status", "gene_tree_lca_ufboot",
        "mapped_node", "mapped_descendant_taxa", "high_support",
    ]
    write_tsv(args.outdir / "homeolog_pair_node_mappings.tsv", records, fieldnames)
    write_tsv(args.outdir / "unreadable_gene_trees.tsv", failed, ["orthogroup", "reason"])

    summary: List[dict] = []
    by_taxon_node: Dict[tuple, List[dict]] = defaultdict(list)
    for record in records:
        if record["status"].startswith("mapped"):
            by_taxon_node[(record["taxon"], record["mapped_node"])].append(record)
    for (taxon, node), values in sorted(by_taxon_node.items()):
        high = sum(row["high_support"] == "yes" for row in values)
        supported = [float(row["gene_tree_lca_ufboot"]) for row in values if row["gene_tree_lca_ufboot"]]
        summary.append({
            "taxon": taxon,
            "mapped_node": node,
            "all_testable_families": len(values),
            "ufboot_ge_cutoff_families": high,
            "high_support_fraction": f"{high / len(values):.6f}",
            "median_lca_ufboot": "" if not supported else f"{sorted(supported)[len(supported)//2]:.1f}",
            "ufboot_cutoff": f"{args.ufboot_cutoff:g}",
        })
    write_tsv(args.outdir / "homeolog_pair_node_summary.tsv", summary, list(summary[0]) if summary else ["taxon"])

    shared: List[dict] = []
    for orthogroup in sorted({row["orthogroup"] for row in records}):
        rows = [r for r in records if r["orthogroup"] == orthogroup and r["high_support"] == "yes"]
        node_to_taxa: Dict[str, Set[str]] = defaultdict(set)
        for row in rows:
            node_to_taxa[row["mapped_node"]].add(row["taxon"])
        for node, taxa in node_to_taxa.items():
            if len(taxa) >= 2:
                shared.append({"orthogroup": orthogroup, "mapped_node": node, "supporting_taxa": ",".join(sorted(taxa)), "n_taxa": len(taxa)})
    write_tsv(args.outdir / "shared_high_support_homeolog_mappings.tsv", shared,
              ["orthogroup", "mapped_node", "supporting_taxa", "n_taxa"])

    metadata = args.outdir / "README.md"
    metadata.write_text(
        "# Corrected-haplotype homeolog-pair node mapping\n\n"
        "Input: complete single-copy orthogroups reconstructed from the corrected A1/B1"
        " haplotype partition. For each orthogroup, the LCA of each predefined A/B"
        " pair (PZI, PSP, PPA) was mapped to the LCA of its collapsed biological"
        " descendant set on the fixed 14-tip species backbone. IQ-TREE UFBoot values"
        f" at that LCA were filtered at >= {args.ufboot_cutoff:g}.\n\n"
        "Interpretation boundary: this is a phylogenetic placement of retained A/B"
        " homeolog divergence. It is not a genome-wide count of duplication events"
        " and does not by itself distinguish allopolyploidy from other histories.\n",
        encoding="utf-8",
    )
    (args.outdir / "COMPLETE").touch()
    print(f"[OK] Mapped {len(records)} homeolog pairs from {len(tree_paths)} corrected gene trees: {args.outdir}")


if __name__ == "__main__":
    main()
