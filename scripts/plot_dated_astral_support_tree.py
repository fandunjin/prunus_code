#!/usr/bin/env python3
"""Render the final ASTRAL topology on the matched MCMCtree time scale.

The script matches nodes strictly by their descendant tip sets.  This is
deliberate: the MCMCtree file is a subgenome-resolved analysis and its dates
must never be copied from a species-only tree onto both A and B lineages.
"""

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
from matplotlib.patches import Circle, Rectangle, Wedge


plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans", "Liberation Sans"]
plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams.update(
    {
        "font.size": 7,
        "axes.linewidth": 0.7,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
    }
)


TREE_COLOR = "#303030"
HPD_COLOR = "#5B7FCA"
MEAN_COLOR = "#0F4D92"
Q_COLORS = ("#3366A8", "#D58A3A", "#9A9A9A")
TIP_COLORS = {"A": "#A85432", "B": "#207A78", "outgroup": "#666666", "default": "#202020"}
DISPLAY_ORDER = [
    "ppa_A", "psp_A", "pzi_A", "ppa_B", "psp_B", "pzi_B", "pmal", "pav", "pspe", "pcam",
    "psal", "par", "pmu", "pdu", "ppe", "gtr", "rch",
]


@dataclass
class Node:
    name: str = ""
    length: float | None = None
    annotations: list[str] = field(default_factory=list)
    children: list["Node"] = field(default_factory=list)
    parent: "Node | None" = None
    descendants: frozenset[str] = field(default_factory=frozenset)
    age_ma: float | None = None
    hpd_low_ma: float | None = None
    hpd_high_ma: float | None = None
    support: dict[str, float] | None = None
    x: float = 0.0
    y: float = 0.0
    date_id: str = ""


def tokenize_newick(text: str) -> list[str]:
    """Tokenize labels and BEAST/FigTree comments without splitting comments."""
    return re.findall(r"\[&[^\]]*\]|\(|\)|,|:|;|'(?:[^']|'')*'|[^\s\(\),:;\[\]]+", text.strip())


def parse_newick(text: str) -> Node:
    tokens = tokenize_newick(text)
    pos = [0]

    def metadata(node: Node, leaf: bool = False) -> None:
        if leaf and pos[0] < len(tokens) and tokens[pos[0]] not in {"(", ")", ",", ":", ";"} and not tokens[pos[0]].startswith("[&"):
            node.name = tokens[pos[0]].strip("'")
            pos[0] += 1
        if not leaf and pos[0] < len(tokens) and tokens[pos[0]] not in {"(", ")", ",", ":", ";"} and not tokens[pos[0]].startswith("[&"):
            node.name = tokens[pos[0]].strip("'")
            pos[0] += 1
        while pos[0] < len(tokens) and tokens[pos[0]].startswith("[&"):
            node.annotations.append(tokens[pos[0]])
            pos[0] += 1
        if pos[0] < len(tokens) and tokens[pos[0]] == ":":
            pos[0] += 1
            if pos[0] >= len(tokens):
                raise ValueError("Missing branch length in Newick input")
            try:
                node.length = float(tokens[pos[0]])
            except ValueError as exc:
                raise ValueError("Non-numeric branch length: %s" % tokens[pos[0]]) from exc
            pos[0] += 1

    def parse_node() -> Node:
        if pos[0] >= len(tokens):
            raise ValueError("Unexpected end of Newick input")
        node = Node()
        if tokens[pos[0]] == "(":
            pos[0] += 1
            node.children.append(parse_node())
            while pos[0] < len(tokens) and tokens[pos[0]] == ",":
                pos[0] += 1
                node.children.append(parse_node())
            if pos[0] >= len(tokens) or tokens[pos[0]] != ")":
                raise ValueError("Unbalanced Newick parentheses")
            pos[0] += 1
            metadata(node, leaf=False)
            for child in node.children:
                child.parent = node
        else:
            metadata(node, leaf=True)
            if not node.name:
                raise ValueError("Unnamed Newick leaf")
        return node

    root = parse_node()
    if pos[0] < len(tokens) and tokens[pos[0]] == ";":
        pos[0] += 1
    if pos[0] != len(tokens):
        raise ValueError("Unexpected trailing Newick tokens: %s" % " ".join(tokens[pos[0]:]))
    return root


def read_astral(path: Path) -> Node:
    return parse_newick(path.read_text(encoding="utf-8"))


def read_mcmctree_figtree(path: Path) -> Node:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"UTREE\s+\d+\s*=\s*(.+?);", text, flags=re.DOTALL | re.IGNORECASE)
    if not match:
        # MCMCtree writes the posterior summary to *.mcmctree.out and the
        # FigTree/NEXUS tree to a companion file named FigTree.tre.  Accepting
        # the former as input is convenient, but the latter is the authoritative
        # source for HPD annotations and branch lengths.
        companion = path.with_name("FigTree.tre")
        if companion.is_file() and companion != path:
            companion_text = companion.read_text(encoding="utf-8")
            match = re.search(
                r"UTREE\s+\d+\s*=\s*(.+?);",
                companion_text,
                flags=re.DOTALL | re.IGNORECASE,
            )
            if match:
                path = companion
        if not match:
            raise ValueError(
                "Could not locate UTREE statement in %s or companion FigTree.tre" % path
            )
    return parse_newick(match.group(1) + ";")


def iter_nodes(root: Node):
    yield root
    for child in root.children:
        yield from iter_nodes(child)


def attach_descendants(node: Node) -> frozenset[str]:
    if not node.children:
        node.descendants = frozenset([node.name])
    else:
        node.descendants = frozenset().union(*(attach_descendants(child) for child in node.children))
    return node.descendants


def orient_for_display(node: Node) -> None:
    ranks = {tip: rank for rank, tip in enumerate(DISPLAY_ORDER)}
    for child in node.children:
        orient_for_display(child)
    node.children.sort(key=lambda child: min(ranks.get(tip, len(ranks)) for tip in child.descendants))


def parse_support(node: Node) -> dict[str, float] | None:
    if not node.name or not node.annotations:
        return None
    parts = node.name.split("/")
    if len(parts) != 4:
        return None
    try:
        bootstrap, pp, gcf, lscf = (float(value) for value in parts)
    except ValueError:
        return None
    annotation = ";".join(node.annotations)
    values = {}
    for key in ("q1", "q2", "q3"):
        match = re.search(key + r"=([0-9.eE+-]+)", annotation)
        if not match:
            return None
        values[key] = float(match.group(1))
    values.update({"bootstrap": bootstrap, "pp": pp, "gcf": gcf, "lscf": lscf})
    return values


def attach_supports(root: Node) -> None:
    for node in iter_nodes(root):
        if node.children:
            node.support = parse_support(node)


def get_hpd(node: Node) -> tuple[float, float] | None:
    for annotation in node.annotations:
        match = re.search(r"95%HPD=\{([0-9.eE+-]+),\s*([0-9.eE+-]+)\}", annotation)
        if match:
            return float(match.group(1)) * 100.0, float(match.group(2)) * 100.0
    return None


def attach_mcmctree_ages(root: Node) -> None:
    """Recover absolute node ages from an ultrametric MCMCtree branch tree."""
    root_to_tip = []

    def walk(node: Node, distance: float) -> None:
        if not node.children:
            root_to_tip.append(distance)
            return
        for child in node.children:
            if child.length is None:
                raise ValueError("Every non-root MCMCtree node needs a branch length")
            walk(child, distance + child.length)

    walk(root, 0.0)
    root_age_units = max(root_to_tip)
    if max(root_to_tip) - min(root_to_tip) > 1e-5:
        raise ValueError("MCMCtree FigTree tree is not ultrametric within 1e-5 time units")

    def set_age(node: Node, age_units: float) -> None:
        node.age_ma = age_units * 100.0
        hpd = get_hpd(node)
        if hpd is not None:
            node.hpd_low_ma, node.hpd_high_ma = hpd
        for child in node.children:
            assert child.length is not None
            set_age(child, age_units - child.length)

    set_age(root, root_age_units)


def map_dates(astral: Node, dated: Node) -> list[Node]:
    astral_tips = attach_descendants(astral)
    dated_tips = attach_descendants(dated)
    if astral_tips != dated_tips:
        raise ValueError("Tip sets differ between ASTRAL and MCMCtree: ASTRAL-only=%s; MCMCtree-only=%s" % (
            ",".join(sorted(astral_tips - dated_tips)), ",".join(sorted(dated_tips - astral_tips))))
    dated_nodes = {node.descendants: node for node in iter_nodes(dated) if node.children}
    unmapped = []
    dated_astral = []
    dated_leaf_ages = {node.name: node.age_ma for node in iter_nodes(dated) if not node.children}
    for node in iter_nodes(astral):
        if not node.children:
            # Posterior trees are ultrametric and present-day tips are 0 Ma.
            # Preserve a strict check so a malformed dated input cannot shift tips.
            if node.name not in dated_leaf_ages or dated_leaf_ages[node.name] is None:
                raise ValueError("MCMCtree age missing for tip: " + node.name)
            # FigTree branch lengths are rounded to six decimals, which gives
            # residual present-day ages around 0.0002 Ma after conversion.
            if abs(dated_leaf_ages[node.name]) > 1e-2:
                raise ValueError("MCMCtree tip is not present-day (0 Ma): " + node.name)
            node.age_ma = 0.0
            continue
        match = dated_nodes.get(node.descendants)
        if match is None:
            unmapped.append(",".join(sorted(node.descendants)))
            continue
        if match.age_ma is None or match.hpd_low_ma is None or match.hpd_high_ma is None:
            raise ValueError("MCMCtree date/95% HPD missing for split: " + ",".join(sorted(node.descendants)))
        node.age_ma = match.age_ma
        node.hpd_low_ma = match.hpd_low_ma
        node.hpd_high_ma = match.hpd_high_ma
        dated_astral.append(node)
    if unmapped:
        raise ValueError("No exact MCMCtree match for ASTRAL node(s): " + " | ".join(unmapped))
    return dated_astral


def assign_layout(root: Node) -> list[Node]:
    leaves = [node for node in iter_nodes(root) if not node.children]
    for index, leaf in enumerate(leaves):
        leaf.y = float(len(leaves) - 1 - index)

    def visit(node: Node) -> None:
        if node.children:
            for child in node.children:
                visit(child)
            node.y = (node.children[0].y + node.children[-1].y) / 2.0
        if node.age_ma is None:
            raise ValueError("Layout requested before every node received an age")
        node.x = -node.age_ma

    visit(root)
    return leaves


def add_q_pie(ax, x: float, y: float, values: tuple[float, float, float]) -> None:
    size = 9.6
    radius = 4.25
    drawing = DrawingArea(size, size, 0, 0)
    total = sum(values)
    start = 90.0
    for value, color in zip(values, Q_COLORS):
        end = start + 360.0 * value / total
        drawing.add_artist(Wedge((size / 2.0, size / 2.0), radius, start, end, facecolor=color, edgecolor="white", linewidth=0.22))
        start = end
    drawing.add_artist(Circle((size / 2.0, size / 2.0), radius, facecolor="none", edgecolor=TREE_COLOR, linewidth=0.4))
    ax.add_artist(AnnotationBbox(drawing, (x, y), xybox=(0, 5.3), xycoords="data", boxcoords="offset points", frameon=False, pad=0, zorder=7))


def fmt(value: float, digits: int = 1) -> str:
    return ("%.0f" if digits == 0 else "%.1f") % value


def support_label(values: dict[str, float]) -> str:
    return "%s/%s/%s/%s" % (fmt(values["bootstrap"], 0), "%.2f" % values["pp"], fmt(values["gcf"]), fmt(values["lscf"]))


def pretty_labels(manifest: Path) -> dict[str, str]:
    labels = {}
    with manifest.open(encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            scientific = row["scientific_name"]
            if scientific.startswith("Prunus "):
                scientific = "P. " + scientific[7:]
            labels[row["label"]] = scientific
    return labels


def tip_color(tip: str) -> str:
    if tip.endswith("_A"):
        return TIP_COLORS["A"]
    if tip.endswith("_B"):
        return TIP_COLORS["B"]
    if tip in {"rch", "gtr"}:
        return TIP_COLORS["outgroup"]
    return TIP_COLORS["default"]


def date_label(node: Node) -> str:
    assert node.age_ma is not None and node.hpd_low_ma is not None and node.hpd_high_ma is not None
    return "%s [%s-%s]" % (fmt(node.age_ma), fmt(node.hpd_low_ma), fmt(node.hpd_high_ma))


def write_source_table(path: Path, dated_nodes: list[Node]) -> None:
    rows = []
    ordered = sorted(dated_nodes, key=lambda node: (-len(node.descendants), tuple(sorted(node.descendants))))
    for index, node in enumerate(ordered, start=1):
        node.date_id = "D%02d" % index
        assert node.age_ma is not None and node.hpd_low_ma is not None and node.hpd_high_ma is not None
        record = {
            "node_id": node.date_id,
            "descendant_taxa": ",".join(sorted(node.descendants)),
            "mean_age_ma": "%.6f" % node.age_ma,
            "hpd_95_low_ma": "%.6f" % node.hpd_low_ma,
            "hpd_95_high_ma": "%.6f" % node.hpd_high_ma,
            "support_label_order": "multilocus_bootstrap/ASTRAL_PP/gCF/likelihood_sCF",
            "multilocus_bootstrap": "",
            "astral_pp": "",
            "gcf": "",
            "likelihood_scf": "",
            "q1": "",
            "q2": "",
            "q3": "",
        }
        if node.support:
            record.update({
                "multilocus_bootstrap": "%.4f" % node.support["bootstrap"], "astral_pp": "%.8f" % node.support["pp"],
                "gcf": "%.4f" % node.support["gcf"], "likelihood_scf": "%.4f" % node.support["lscf"],
                "q1": "%.8f" % node.support["q1"], "q2": "%.8f" % node.support["q2"], "q3": "%.8f" % node.support["q3"],
            })
        rows.append(record)
    fields = list(rows[0])
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def draw_geological_scale(ax, oldest_age: float, youngest_age: float) -> None:
    periods = [
        (130.0, 66.0, "Cretaceous", "#E9E6C4"),
        (66.0, 23.03, "Paleogene", "#DCE8F3"),
        (23.03, 2.58, "Neogene", "#E2EED8"),
        (2.58, 0.0, "Quaternary", "#F2E4D8"),
    ]
    for old, young, name, color in periods:
        left = min(old, oldest_age)
        right = max(young, youngest_age)
        if left <= right:
            continue
        ax.add_patch(Rectangle((-left, 0.0), left - right, 0.52, facecolor=color, edgecolor="white", linewidth=0.5))
        if left - right > 5:
            ax.text(-(left + right) / 2.0, 0.26, name, ha="center", va="center", fontsize=5.3, color="#4D4D4D")
    ax.set_ylim(0, 0.52)
    ax.set_yticks([])
    ax.set_xticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)


def render(root: Node, labels: dict[str, str], out_prefix: Path, title: str | None) -> None:
    leaves = assign_layout(root)
    dated_nodes = [node for node in iter_nodes(root) if node.children]
    max_age = 130.0
    min_age = 0.0
    fig = plt.figure(figsize=(7.22, 8.35))
    grid = fig.add_gridspec(nrows=2, ncols=1, height_ratios=[16, 0.8], hspace=0.035)
    ax = fig.add_subplot(grid[0, 0])
    ax_time = fig.add_subplot(grid[1, 0], sharex=ax)

    # Tree topology, positioned by posterior mean age rather than substitutions.
    for node in iter_nodes(root):
        if not node.children:
            continue
        ys = [child.y for child in node.children]
        ax.plot([node.x, node.x], [min(ys), max(ys)], color=TREE_COLOR, lw=0.72, solid_capstyle="round", zorder=1)
        for child in node.children:
            ax.plot([node.x, child.x], [child.y, child.y], color=TREE_COLOR, lw=0.72, solid_capstyle="round", zorder=1)
            if child.support is not None:
                midpoint = (node.x + child.x) / 2.0
                add_q_pie(ax, midpoint, child.y, (child.support["q1"], child.support["q2"], child.support["q3"]))
                ax.annotate(
                    support_label(child.support), xy=(midpoint, child.y), xytext=(0, -5.6), textcoords="offset points",
                    ha="center", va="top", fontsize=4.65, color="#202020",
                    bbox={"boxstyle": "square,pad=0.04", "facecolor": "white", "edgecolor": "none", "alpha": 0.91}, zorder=8,
                )

    # 95% HPD bars use their genuine posterior interval; the circle is the posterior mean.
    for node in dated_nodes:
        assert node.age_ma is not None and node.hpd_low_ma is not None and node.hpd_high_ma is not None
        x_low, x_high = -node.hpd_high_ma, -node.hpd_low_ma
        ax.plot([x_low, x_high], [node.y, node.y], color=HPD_COLOR, alpha=0.65, lw=2.0, solid_capstyle="round", zorder=4)
        ax.plot([-node.age_ma], [node.y], marker="o", ms=2.8, mfc=MEAN_COLOR, mec="white", mew=0.28, zorder=6)
        # Support labels sit below their branches, so dates consistently sit
        # above the node to keep the two independent evidence layers legible.
        offset = 4.8
        ax.annotate(date_label(node), xy=(-node.age_ma, node.y), xytext=(0, offset), textcoords="offset points",
                    ha="center", va="bottom" if offset > 0 else "top", fontsize=4.35, color=MEAN_COLOR,
                    bbox={"boxstyle": "square,pad=0.03", "facecolor": "white", "edgecolor": "none", "alpha": 0.76}, zorder=9)

    for leaf in leaves:
        ax.text(2.1, leaf.y, labels.get(leaf.name, leaf.name), ha="left", va="center", fontsize=6.25,
                fontstyle="italic", color=tip_color(leaf.name), clip_on=False)

    ax.text(-129.0, len(leaves) - 0.12, "Time-calibrated ASTRAL topology", fontsize=7.1, fontweight="bold", ha="left", va="bottom")
    if title:
        ax.text(-129.0, len(leaves) - 0.72, title, fontsize=5.5, color="#4D4D4D", ha="left", va="bottom")
    ax.text(-4.5, len(leaves) - 0.13, "Support: BS / PP / gCF / l-sCF", fontsize=5.35, ha="right", va="bottom")
    for x, text, color in [(-23.0, "q1", Q_COLORS[0]), (-15.0, "q2", Q_COLORS[1]), (-7.0, "q3", Q_COLORS[2])]:
        ax.plot([x], [len(leaves) - 0.37], marker="o", ms=4.0, color=color, clip_on=False)
        ax.text(x + 1.2, len(leaves) - 0.37, text, fontsize=5.25, va="center", ha="left")
    ax.text(-129.0, -1.25, "Blue bar: 95% HPD; blue point: posterior mean (Ma)", fontsize=5.45, color=MEAN_COLOR, ha="left", va="center")

    ax.set_xlim(-max_age, 36.5)
    ax.set_ylim(-1.72, len(leaves) + 0.25)
    ax.set_axis_off()

    ticks = [120, 100, 80, 60, 40, 20, 0]
    ax_time.set_xlim(-max_age, 36.5)
    ax_time.set_ylim(0, 0.52)
    ax_time.set_xticks([-tick for tick in ticks])
    ax_time.set_xticklabels([str(tick) for tick in ticks], fontsize=5.4)
    ax_time.tick_params(axis="x", length=2.2, width=0.5, pad=1)
    ax_time.set_xlabel("Time before present (Ma)", fontsize=6.1, labelpad=1.4)
    draw_geological_scale(ax_time, max_age, min_age)
    ax_time.set_xticks([-tick for tick in ticks])
    ax_time.set_xticklabels([str(tick) for tick in ticks], fontsize=5.4)
    ax_time.tick_params(axis="x", length=2.2, width=0.5, pad=1)
    ax_time.set_xlabel("Time before present (Ma)", fontsize=6.1, labelpad=1.4)

    fig.subplots_adjust(left=0.055, right=0.835, bottom=0.085, top=0.965)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(out_prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".tiff", dpi=600, bbox_inches="tight")
    fig.savefig(str(out_prefix) + ".png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--astral", required=True, type=Path)
    parser.add_argument("--mcmctree", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--title", default=None)
    args = parser.parse_args()
    for path in (args.astral, args.mcmctree, args.manifest):
        if not path.is_file():
            raise SystemExit("ERROR: Missing input: " + str(path))
    if args.outdir.exists():
        raise SystemExit("ERROR: Output directory already exists: " + str(args.outdir))

    astral = read_astral(args.astral)
    dated = read_mcmctree_figtree(args.mcmctree)
    attach_descendants(astral)
    attach_descendants(dated)
    orient_for_display(astral)
    attach_mcmctree_ages(dated)
    dated_nodes = map_dates(astral, dated)
    attach_supports(astral)
    labels = pretty_labels(args.manifest)

    args.outdir.mkdir(parents=True)
    write_source_table(args.outdir / "dated_node_support_and_hpd.tsv", dated_nodes)
    render(astral, labels, args.outdir / "final_astral_support_tree_time_calibrated", args.title)
    (args.outdir / "README.txt").write_text(
        "MCMCtree dates were mapped to the ASTRAL topology by exact descendant-tip sets.\n"
        "Time values were converted from MCMCtree units of 100 Myr to Ma.\n"
        "These dates represent bifurcating lineage divergences, not allopolyploidization ages.\n",
        encoding="ascii",
    )
    (args.outdir / "COMPLETE").write_text("Time-calibrated ASTRAL support tree completed.\n", encoding="ascii")
    print("[OK] Dated ASTRAL support tree written to " + str(args.outdir))


if __name__ == "__main__":
    main()
