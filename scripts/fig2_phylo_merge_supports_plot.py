#!/usr/bin/env python3
"""Merge phylogenomic branch support by split and draw annotated ASTRAL trees."""

from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnnotationBbox, DrawingArea
from matplotlib.patches import Circle, Patch, Wedge


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "font.size": 7,
        "axes.linewidth": 0.8,
        "legend.frameon": False,
    }
)


Q_COLORS = ("#3366A8", "#D58A3A", "#9A9A9A")
TREE_COLOR = "#303030"
TIP_COLORS = {
    "A": "#A85432",
    "B": "#207A78",
    "outgroup": "#666666",
    "default": "#202020",
}

DISPLAY_ORDER = [
    "ppa_A",
    "psp_A",
    "pzi_A",
    "ppa_B",
    "psp_B",
    "pzi_B",
    "pmal",
    "pav",
    "pspe",
    "pcam",
    "psal",
    "par",
    "pmu",
    "pdu",
    "ppe",
    "gtr",
    "rch",
]


@dataclass
class Node:
    name: str = ""
    length: str = ""
    children: list["Node"] = field(default_factory=list)


@dataclass
class RootedNode:
    graph_id: int | None = None
    name: str = ""
    children: list["RootedNode"] = field(default_factory=list)
    parent: "RootedNode | None" = None
    x: float = 0.0
    y: float = 0.0
    descendant_taxa: frozenset[str] = field(default_factory=frozenset)
    support: dict[str, float] | None = None
    node_id: str = ""


def tokenize_newick(text: str) -> list[str]:
    return re.findall(r"\(|\)|,|:|;|'(?:[^']|'')*'|[^\s\(\),:;]+", text.strip())


def parse_newick(text: str) -> Node:
    tokens = tokenize_newick(text)
    position = [0]

    def parse_label_length() -> tuple[str, str]:
        name = ""
        length = ""
        if position[0] < len(tokens) and tokens[position[0]] not in ["(", ")", ",", ":", ";"]:
            name = tokens[position[0]]
            if name.startswith("'") and name.endswith("'"):
                name = name[1:-1].replace("''", "'")
            position[0] += 1
        if position[0] < len(tokens) and tokens[position[0]] == ":":
            position[0] += 1
            if position[0] >= len(tokens):
                raise ValueError("missing branch length")
            length = tokens[position[0]]
            position[0] += 1
        return name, length

    def parse_node() -> Node:
        if position[0] < len(tokens) and tokens[position[0]] == "(":
            position[0] += 1
            children = [parse_node()]
            while position[0] < len(tokens) and tokens[position[0]] == ",":
                position[0] += 1
                children.append(parse_node())
            if position[0] >= len(tokens) or tokens[position[0]] != ")":
                raise ValueError("unbalanced Newick parentheses")
            position[0] += 1
            name, length = parse_label_length()
            return Node(name=name, length=length, children=children)
        name, length = parse_label_length()
        if not name:
            raise ValueError("unnamed Newick leaf")
        return Node(name=name, length=length)

    root = parse_node()
    if position[0] < len(tokens) and tokens[position[0]] == ";":
        position[0] += 1
    if position[0] != len(tokens):
        raise ValueError("unexpected Newick tokens")
    return root


def read_tree(path: Path) -> Node:
    return parse_newick(path.read_text(encoding="utf-8"))


def iter_nodes(root: Node):
    yield root
    for child in root.children:
        yield from iter_nodes(child)


def terminal_names(root: Node) -> frozenset[str]:
    return frozenset(node.name for node in iter_nodes(root) if not node.children)


def descendant_map(root: Node) -> dict[int, frozenset[str]]:
    result: dict[int, frozenset[str]] = {}

    def visit(node: Node) -> frozenset[str]:
        taxa = frozenset([node.name]) if not node.children else frozenset().union(*(visit(c) for c in node.children))
        result[id(node)] = taxa
        return taxa

    visit(root)
    return result


def canonical_split(side: frozenset[str], all_taxa: frozenset[str]) -> tuple[str, ...] | None:
    other = all_taxa - side
    if len(side) < 2 or len(other) < 2:
        return None
    left = tuple(sorted(side))
    right = tuple(sorted(other))
    if len(left) < len(right):
        return left
    if len(right) < len(left):
        return right
    return min(left, right)


def split_text(key: tuple[str, ...], all_taxa: frozenset[str]) -> str:
    other = tuple(sorted(all_taxa - frozenset(key)))
    return ",".join(key) + "|" + ",".join(other)


def parse_astral_annotation(label: str) -> dict[str, float] | None:
    if not label.startswith("[") or not label.endswith("]"):
        return None
    values: dict[str, float] = {}
    for field_text in label[1:-1].split(";"):
        if "=" not in field_text:
            continue
        key, value = field_text.split("=", 1)
        try:
            values[key] = float(value)
        except ValueError:
            continue
    required = {"q1", "q2", "q3", "pp1"}
    return values if required.issubset(values) else None


def support_map(root: Node, kind: str) -> dict[tuple[str, ...], dict[str, float]]:
    all_taxa = terminal_names(root)
    descendants = descendant_map(root)
    result: dict[tuple[str, ...], dict[str, float]] = {}
    for node in iter_nodes(root):
        if not node.children:
            continue
        split = canonical_split(descendants[id(node)], all_taxa)
        if split is None:
            continue
        values: dict[str, float] | None = None
        if kind == "astral":
            values = parse_astral_annotation(node.name)
        elif kind == "bootstrap":
            try:
                values = {"bootstrap": float(node.name)}
            except ValueError:
                values = None
        elif kind == "gcf":
            try:
                values = {"gcf": float(node.name.split("/", 1)[0])}
            except ValueError:
                values = None
        elif kind == "lscf":
            try:
                values = {"lscf": float(node.name)}
            except ValueError:
                values = None
        else:
            raise ValueError("unsupported map kind: " + kind)
        if values is None:
            continue
        if split in result:
            raise ValueError("duplicate supported split in %s: %s" % (kind, split))
        result[split] = values
    return result


def graph_from_tree(root: Node):
    adjacency: dict[int, set[int]] = {}
    leaf_names: dict[int, str] = {}
    counter = [0]

    def visit(node: Node, parent: int | None = None):
        ident = counter[0]
        counter[0] += 1
        adjacency[ident] = set()
        if parent is not None:
            adjacency[ident].add(parent)
            adjacency[parent].add(ident)
        if not node.children:
            leaf_names[ident] = node.name
        for child in node.children:
            visit(child, ident)

    visit(root)
    return adjacency, leaf_names


def suppress_degree_two_nodes(adjacency: dict[int, set[int]], leaf_names: dict[int, str]) -> None:
    changed = True
    while changed:
        changed = False
        for ident in list(adjacency):
            if ident in leaf_names or len(adjacency[ident]) != 2:
                continue
            left, right = tuple(adjacency[ident])
            adjacency[left].remove(ident)
            adjacency[right].remove(ident)
            adjacency[left].add(right)
            adjacency[right].add(left)
            del adjacency[ident]
            changed = True
            break


def root_on_outgroup(backbone: Node, outgroup: str) -> RootedNode:
    adjacency, leaf_names = graph_from_tree(backbone)
    suppress_degree_two_nodes(adjacency, leaf_names)
    matches = [ident for ident, name in leaf_names.items() if name == outgroup]
    if len(matches) != 1:
        raise ValueError("outgroup must occur exactly once: " + outgroup)
    out_id = matches[0]
    if len(adjacency[out_id]) != 1:
        raise ValueError("outgroup is not a terminal")
    neighbor = next(iter(adjacency[out_id]))

    def build(current: int, parent_id: int) -> RootedNode:
        node = RootedNode(graph_id=current, name=leaf_names.get(current, ""))
        for child_id in sorted(adjacency[current] - {parent_id}):
            child = build(child_id, current)
            child.parent = node
            node.children.append(child)
        return node

    root = RootedNode()
    out_node = RootedNode(graph_id=out_id, name=outgroup, parent=root)
    ingroup_node = build(neighbor, out_id)
    ingroup_node.parent = root
    root.children = [out_node, ingroup_node]
    return root


def rooted_nodes(root: RootedNode):
    yield root
    for child in root.children:
        yield from rooted_nodes(child)


def attach_descendants(root: RootedNode) -> frozenset[str]:
    if not root.children:
        root.descendant_taxa = frozenset([root.name])
    else:
        root.descendant_taxa = frozenset().union(*(attach_descendants(child) for child in root.children))
    return root.descendant_taxa


def orient_for_display(root: RootedNode) -> None:
    rank = {taxon: index for index, taxon in enumerate(DISPLAY_ORDER)}

    def visit(node: RootedNode) -> frozenset[str]:
        if not node.children:
            node.descendant_taxa = frozenset([node.name])
            return node.descendant_taxa
        for child in node.children:
            visit(child)
        node.children.sort(
            key=lambda child: min(rank.get(taxon, len(rank)) for taxon in child.descendant_taxa)
        )
        node.descendant_taxa = frozenset().union(*(child.descendant_taxa for child in node.children))
        return node.descendant_taxa

    visit(root)


def assign_supports(
    rooted: RootedNode,
    maps: list[dict[tuple[str, ...], dict[str, float]]],
) -> list[RootedNode]:
    all_taxa = attach_descendants(rooted)
    supported: list[RootedNode] = []
    for node in rooted_nodes(rooted):
        if node is rooted or not node.children:
            continue
        split = canonical_split(node.descendant_taxa, all_taxa)
        if split is None:
            continue
        values: dict[str, float] = {}
        for support_values in maps:
            if split not in support_values:
                raise ValueError("support missing for split: " + split_text(split, all_taxa))
            values.update(support_values[split])
        required = {"q1", "q2", "q3", "pp1", "bootstrap", "gcf", "lscf"}
        missing = sorted(required - set(values))
        if missing:
            raise ValueError("support fields missing for split %s: %s" % (split, ",".join(missing)))
        node.support = values
        supported.append(node)
    expected = len(all_taxa) - 3
    if len(supported) != expected:
        raise ValueError("expected %d internal splits, found %d" % (expected, len(supported)))
    return supported


def assign_layout(root: RootedNode) -> list[RootedNode]:
    leaves = [node for node in rooted_nodes(root) if not node.children]
    for index, leaf in enumerate(leaves):
        leaf.y = float(len(leaves) - 1 - index)

    def set_y(node: RootedNode):
        if node.children:
            for child in node.children:
                set_y(child)
            node.y = (node.children[0].y + node.children[-1].y) / 2.0

    def set_x(node: RootedNode, depth: int = 0):
        node.x = float(depth)
        for child in node.children:
            set_x(child, depth + 1)

    set_y(root)
    set_x(root)
    max_depth = max(leaf.x for leaf in leaves)
    for leaf in leaves:
        leaf.x = max_depth + 0.4
    return leaves


def add_q_pie(ax, x: float, y: float, values: tuple[float, float, float]) -> None:
    size = 12.0
    radius = 5.4
    drawing = DrawingArea(size, size, 0, 0)
    total = sum(values)
    start = 90.0
    for value, color in zip(values, Q_COLORS):
        end = start + 360.0 * value / total
        drawing.add_artist(
            Wedge(
                (size / 2.0, size / 2.0),
                radius,
                start,
                end,
                facecolor=color,
                edgecolor="white",
                linewidth=0.25,
            )
        )
        start = end
    drawing.add_artist(
        Circle((size / 2.0, size / 2.0), radius, facecolor="none", edgecolor=TREE_COLOR, linewidth=0.45)
    )
    ax.add_artist(
        AnnotationBbox(
            drawing,
            (x, y),
            xybox=(0, 7.5),
            xycoords="data",
            boxcoords="offset points",
            frameon=False,
            pad=0,
            zorder=6,
        )
    )


def fmt_number(value: float, digits: int = 1) -> str:
    if digits == 0:
        return str(int(round(value)))
    return ("%.*f" % (digits, value))


def scientific_names(manifest: Path) -> dict[str, str]:
    labels: dict[str, str] = {}
    with manifest.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scientific_name = row["scientific_name"]
            if scientific_name.startswith("Prunus "):
                scientific_name = "P. " + scientific_name[len("Prunus ") :]
            labels[row["label"]] = scientific_name
    return labels


def tip_color(name: str) -> str:
    if name.endswith("_A"):
        return TIP_COLORS["A"]
    if name.endswith("_B"):
        return TIP_COLORS["B"]
    if name in {"rch", "gtr"}:
        return TIP_COLORS["outgroup"]
    return TIP_COLORS["default"]


def render_tree(
    root: RootedNode,
    labels: dict[str, str],
    dataset: str,
    out_prefix: Path,
    width: float,
    height: float,
) -> None:
    leaves = assign_layout(root)
    fig, ax = plt.subplots(figsize=(width, height))

    for node in rooted_nodes(root):
        if not node.children:
            continue
        child_y = [child.y for child in node.children]
        ax.plot([node.x, node.x], [min(child_y), max(child_y)], color=TREE_COLOR, lw=0.85, zorder=1)
        for child in node.children:
            ax.plot([node.x, child.x], [child.y, child.y], color=TREE_COLOR, lw=0.85, zorder=1)
            if child.support is None:
                continue
            midpoint = (node.x + child.x) / 2.0
            add_q_pie(
                ax,
                midpoint,
                child.y,
                (child.support["q1"], child.support["q2"], child.support["q3"]),
            )
            value_label = "%s/%s/%s/%s" % (
                fmt_number(child.support["bootstrap"], 0),
                fmt_number(child.support["pp1"], 2),
                fmt_number(child.support["gcf"], 1),
                fmt_number(child.support["lscf"], 1),
            )
            ax.annotate(
                value_label,
                xy=(midpoint, child.y),
                xytext=(0, -7.0),
                textcoords="offset points",
                ha="center",
                va="top",
                fontsize=5.1,
                color="#202020",
                bbox={"boxstyle": "square,pad=0.08", "facecolor": "white", "edgecolor": "none", "alpha": 0.88},
                zorder=7,
            )

    label_x = max(leaf.x for leaf in leaves) + 0.10
    for leaf in leaves:
        ax.text(
            label_x,
            leaf.y,
            labels.get(leaf.name, leaf.name),
            ha="left",
            va="center",
            fontsize=6.6,
            fontstyle="italic",
            color=tip_color(leaf.name),
        )

    legend_handles = [Patch(facecolor=color, edgecolor="none", label=label) for color, label in zip(Q_COLORS, ["q1", "q2", "q3"])]
    legend = ax.legend(
        handles=legend_handles,
        ncol=3,
        loc="lower left",
        bbox_to_anchor=(0.0, 1.002),
        borderaxespad=0,
        handlelength=1.0,
        handletextpad=0.35,
        columnspacing=0.9,
        fontsize=6.2,
    )
    ax.add_artist(legend)
    ax.text(
        1.0,
        1.012,
        "BS/PP/gCF/l-sCF",
        transform=ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=6.2,
        color="#303030",
    )

    ax.set_xlim(-0.15, label_x + 2.35)
    ax.set_ylim(-0.65, len(leaves) - 0.35)
    ax.set_axis_off()
    fig.subplots_adjust(left=0.02, right=0.98, bottom=0.02, top=0.955)

    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_prefix) + ".pdf", bbox_inches="tight")
    plt.close(fig)


def support_label(values: dict[str, float]) -> str:
    return "%s/%s/%s/%s" % (
        fmt_number(values["bootstrap"], 0),
        fmt_number(values["pp1"], 2),
        fmt_number(values["gcf"], 1),
        fmt_number(values["lscf"], 1),
    )


def annotated_newick(node: RootedNode, all_taxa: frozenset[str]) -> str:
    if not node.children:
        return node.name
    body = "(" + ",".join(annotated_newick(child, all_taxa) for child in node.children) + ")"
    if node.support is None:
        return body
    values = node.support
    comment = "[&q1=%.8f;q2=%.8f;q3=%.8f]" % (values["q1"], values["q2"], values["q3"])
    return body + support_label(values) + comment


def write_support_table(path: Path, root: RootedNode, supported: list[RootedNode], dataset: str) -> None:
    all_taxa = root.descendant_taxa
    fields = [
        "dataset",
        "node_id",
        "displayed_descendant_taxa",
        "bipartition",
        "q1",
        "q2",
        "q3",
        "astral_pp",
        "multilocus_bootstrap",
        "gcf",
        "likelihood_scf",
        "effective_genes",
        "quartet_count",
        "display_label",
    ]
    rows = []
    for index, node in enumerate(supported, start=1):
        node.node_id = "N%02d" % index
        split = canonical_split(node.descendant_taxa, all_taxa)
        assert split is not None and node.support is not None
        values = node.support
        rows.append(
            {
                "dataset": dataset,
                "node_id": node.node_id,
                "displayed_descendant_taxa": ",".join(sorted(node.descendant_taxa)),
                "bipartition": split_text(split, all_taxa),
                "q1": "%.8f" % values["q1"],
                "q2": "%.8f" % values["q2"],
                "q3": "%.8f" % values["q3"],
                "astral_pp": "%.8f" % values["pp1"],
                "multilocus_bootstrap": "%.4f" % values["bootstrap"],
                "gcf": "%.4f" % values["gcf"],
                "likelihood_scf": "%.4f" % values["lscf"],
                "effective_genes": "%.6f" % values.get("EN", float("nan")),
                "quartet_count": "%.6f" % values.get("QC", float("nan")),
                "display_label": support_label(values),
            }
        )
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def process_dataset(base: Path, dataset: str, out_dir: Path, labels: dict[str, str], outgroup: str, width: float, height: float):
    source = base / "06_species_trees" / dataset
    astral_path = source / "astral" / "species_tree.bs10.full_annotation.tre"
    bootstrap_path = source / "astral" / "gene_site_bootstrap.final_tree.tre"
    gcf_path = source / "concordance.cf.tree"
    lscf_path = source / "concordance_lscf.cf.tree"
    for path in [astral_path, bootstrap_path, gcf_path, lscf_path]:
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("missing input: " + str(path))

    astral_tree = read_tree(astral_path)
    trees = {
        "astral": astral_tree,
        "bootstrap": read_tree(bootstrap_path),
        "gcf": read_tree(gcf_path),
        "lscf": read_tree(lscf_path),
    }
    taxon_sets = {name: terminal_names(tree) for name, tree in trees.items()}
    if len(set(taxon_sets.values())) != 1:
        raise ValueError("taxon sets differ among support trees: %s" % taxon_sets)

    maps = [support_map(trees[name], name) for name in ["astral", "bootstrap", "gcf", "lscf"]]
    rooted = root_on_outgroup(astral_tree, outgroup)
    orient_for_display(rooted)
    supported = assign_supports(rooted, maps)
    dataset_dir = out_dir / dataset
    write_support_table(dataset_dir / "final_node_support.tsv", rooted, supported, dataset)
    (dataset_dir / "final_astral_rch_rooted.annotated.tre").write_text(
        annotated_newick(rooted, rooted.descendant_taxa) + ";\n", encoding="utf-8"
    )
    render_tree(rooted, labels, dataset, dataset_dir / "final_astral_support_tree", width, height)
    (dataset_dir / "COMPLETE").write_text("Final support tree completed.\n", encoding="ascii")
    return len(supported)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", required=True, help="phylogenetics output root")
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--datasets", nargs="+", default=["main", "sensitivity"])
    parser.add_argument("--outgroup", default="rch")
    parser.add_argument("--width", type=float, default=7.2)
    parser.add_argument("--height", type=float, default=7.0)
    args = parser.parse_args()

    base = Path(args.base).resolve()
    out_dir = Path(args.out_dir).resolve()
    if out_dir.exists():
        raise SystemExit("ERROR: output directory already exists: " + str(out_dir))
    labels = scientific_names(Path(args.manifest))
    out_dir.mkdir(parents=True)

    counts = {}
    for dataset in args.datasets:
        counts[dataset] = process_dataset(
            base=base,
            dataset=dataset,
            out_dir=out_dir,
            labels=labels,
            outgroup=args.outgroup,
            width=args.width,
            height=args.height,
        )
    with (out_dir / "merge_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["dataset", "taxa", "internal_splits", "outgroup", "node_label_order"])
        for dataset in args.datasets:
            writer.writerow([dataset, "17", counts[dataset], args.outgroup, "bootstrap/PP/gCF/likelihood-sCF"])
    (out_dir / "COMPLETE").write_text("Final support figures completed.\n", encoding="ascii")
    print("[OK] Final support trees completed: " + str(out_dir))


if __name__ == "__main__":
    main()
