#!/usr/bin/env python3
"""Plot the five-event Prunus ancestral-karyotype breakpoint model.

The drawing is a structural schematic, not a base-pair-scaled rearrangement
map. It follows the nine-unit model used for Fig. 3/4:
AK1=C2B, AK2=C5B, and AK3--AK9=GTR Chr03--Chr09. The event classes shown
are deliberately limited to one shared Prunus RCT, two A-associated RCTs plus
one A-associated EEJ, and one B-associated RCT plus one B-associated EEJ.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.path import Path as MplPath


COLORS = {
    "AK1": "#00AEEF",  # C2B
    "AK2": "#DF8442",  # C5B
    "AK3": "#5B9255",  # GTR03
    "AK4": "#C9595F",  # GTR04
    "AK5": "#9E7497",  # GTR05
    "AK6": "#6EAAA6",  # GTR06
    "AK7": "#D7B75A",  # GTR07
    "AK8": "#937366",  # GTR08
    "AK9": "#3A5C93",  # GTR09
}


EVENTS = [
    {
        "event_id": "P-RCT",
        "assignment": "Shared Prunus",
        "event": "RCT",
        "breakpoint": "AK2 x AK9",
        "products": "AKP8 / AKP9",
        "basis": "AK2 = C5B; AK9 = GTR09",
    },
    {
        "event_id": "A-RCT1",
        "assignment": "A-associated",
        "event": "RCT",
        "breakpoint": "AK5 x AKP8",
        "products": "A product 1 / reciprocal product",
        "basis": "AK5 = GTR05",
    },
    {
        "event_id": "A-RCT2",
        "assignment": "A-associated",
        "event": "RCT",
        "breakpoint": "AK1 x reciprocal product",
        "products": "C2A configuration / residual product",
        "basis": "AK1 = C2B",
    },
    {
        "event_id": "A-EEJ",
        "assignment": "A-associated",
        "event": "EEJ",
        "breakpoint": "AK8 x residual product",
        "products": "C7A configuration",
        "basis": "AK8 = GTR08",
    },
    {
        "event_id": "B-RCT",
        "assignment": "B-associated",
        "event": "RCT",
        "breakpoint": "AKP9 x AK5",
        "products": "B product 1 / reciprocal product",
        "basis": "AK5 = GTR05",
    },
    {
        "event_id": "B-EEJ",
        "assignment": "B-associated",
        "event": "EEJ",
        "breakpoint": "AKP3 x AKP8",
        "products": "B n = 8 configuration",
        "basis": "terminal joining",
    },
]


def add_text(ax, x: float, y: float, text: str, **kwargs: object) -> None:
    defaults = {"ha": "center", "va": "center", "fontsize": 7, "color": "#222222"}
    defaults.update(kwargs)
    ax.text(x, y, text, **defaults)


def draw_bar(ax, x: float, y: float, width: float, height: float, segments: list[tuple[float, str, str]], label: str = "") -> None:
    cursor = x
    total = sum(item[0] for item in segments)
    for fraction, color_key, text in segments:
        segment_width = width * fraction / total
        ax.add_patch(Rectangle((cursor, y), segment_width, height, facecolor=COLORS[color_key], edgecolor="white", linewidth=0.7))
        if segment_width >= 0.16:
            add_text(ax, cursor + segment_width / 2, y + height / 2, text, fontsize=6.5, color="white", fontweight="bold")
        cursor += segment_width
    ax.add_patch(Rectangle((x, y), width, height, fill=False, edgecolor="#333333", linewidth=0.65))
    if label:
        add_text(ax, x - 0.05, y + height / 2, label, ha="right", fontsize=6.5, color="#444444")


def unit_visual(label: str) -> list[tuple[float, str, str]]:
    """Return palette segments for an ancestral or derived unit label."""
    if label in COLORS:
        return [(1.0, label, label)]
    if label == "AKP8":
        return [(0.52, "AK2", "AK2"), (0.48, "AK9", "AK9")]
    if label == "AKP9":
        return [(0.52, "AK9", "AK9"), (0.48, "AK2", "AK2")]
    if label == "AKP3":
        return [(1.0, "AK3", "AKP3")]
    if label in {"reciprocal product", "residual product"}:
        return [(0.52, "AK5", "R1"), (0.48, "AK1", "R2")]
    return [(1.0, "AK5", label[:5])]


def draw_rct_card(ax, x: float, y: float, width: float, height: float, event: dict[str, str], second_pair: bool = False) -> None:
    """Draw two reciprocal exchanges with crossing breakpoint connectors."""
    ax.add_patch(Rectangle((x, y), width, height, facecolor="#FFFFFF", edgecolor="#C9CDD2", linewidth=0.8))
    add_text(ax, x + 0.08, y + height - 0.15, event["event_id"], ha="left", fontsize=7, fontweight="bold", color="#222222")
    add_text(ax, x + width - 0.08, y + height - 0.15, event["event"], ha="right", fontsize=7, fontweight="bold", color="#A23B4A")
    add_text(ax, x + width / 2, y + height - 0.34, event["breakpoint"], fontsize=6.2, color="#444444")

    top_y1 = y + height * 0.61
    top_y2 = y + height * 0.46
    bottom_y1 = y + height * 0.19
    bottom_y2 = y + height * 0.04
    bar_x = x + width * 0.17
    bar_w = width * 0.66
    bar_h = height * 0.075
    left, right = event["breakpoint"].split(" x ", 1)
    left_segments = unit_visual(left)
    right_segments = unit_visual(right)
    draw_bar(ax, bar_x, top_y1, bar_w, bar_h, left_segments, "input 1")
    draw_bar(ax, bar_x, top_y2, bar_w, bar_h, right_segments, "input 2")
    draw_bar(ax, bar_x, bottom_y1, bar_w, bar_h, left_segments + right_segments, "reciprocal 1")
    draw_bar(ax, bar_x, bottom_y2, bar_w, bar_h, right_segments + left_segments, "reciprocal 2")

    split_x = bar_x + bar_w * 0.52
    for y0, y1, direction in ((top_y1 + bar_h / 2, bottom_y2 + bar_h / 2, -1), (top_y2 + bar_h / 2, bottom_y1 + bar_h / 2, 1)):
        verts = [(split_x, y0), (split_x + width * 0.10 * direction, y0 - height * 0.07), (split_x - width * 0.10 * direction, y1 + height * 0.07), (split_x, y1)]
        path = MplPath(verts, [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
        ax.add_patch(FancyArrowPatch(path=path, arrowstyle="-", linewidth=1.3, color="#333333", alpha=0.85))
    ax.plot([split_x, split_x], [top_y2 - 0.01, top_y1 + bar_h + 0.01], color="#222222", linewidth=1.0, linestyle=(0, (2, 2)))
    add_text(ax, x + width / 2, y + 0.005, event["products"], fontsize=5.5, color="#666666")


def draw_eej_card(ax, x: float, y: float, width: float, height: float, event: dict[str, str]) -> None:
    ax.add_patch(Rectangle((x, y), width, height, facecolor="#FFFFFF", edgecolor="#C9CDD2", linewidth=0.8))
    add_text(ax, x + 0.08, y + height - 0.15, event["event_id"], ha="left", fontsize=7, fontweight="bold", color="#222222")
    add_text(ax, x + width - 0.08, y + height - 0.15, event["event"], ha="right", fontsize=7, fontweight="bold", color="#3D7A57")
    add_text(ax, x + width / 2, y + height - 0.34, event["breakpoint"], fontsize=6.2, color="#444444")
    bar_x = x + width * 0.18
    bar_w = width * 0.64
    bar_h = height * 0.10
    left, right = event["breakpoint"].split(" x ", 1)
    left_segments = unit_visual(left)
    right_segments = unit_visual(right)
    draw_bar(ax, bar_x, y + height * 0.56, bar_w * 0.48, bar_h, left_segments, "end 1")
    draw_bar(ax, bar_x + bar_w * 0.52, y + height * 0.56, bar_w * 0.48, bar_h, right_segments, "end 2")
    add_text(ax, x + width / 2, y + height * 0.46, "join", fontsize=6, color="#3D7A57", fontweight="bold")
    draw_bar(ax, bar_x, y + height * 0.19, bar_w, bar_h, left_segments + right_segments, "joined")
    left_fraction = sum(item[0] for item in left_segments) / (sum(item[0] for item in left_segments) + sum(item[0] for item in right_segments))
    ax.plot([bar_x + bar_w * left_fraction, bar_x + bar_w * left_fraction], [y + height * 0.18, y + height * 0.31], color="#3D7A57", linewidth=1.8)
    add_text(ax, x + width / 2, y + 0.005, event["products"], fontsize=5.5, color="#666666")


def draw_ancestor(ax: plt.Axes, x: float, y: float, width: float) -> None:
    add_text(ax, x, y + 0.34, "Ancestral nine-unit reference", ha="left", fontsize=8, fontweight="bold")
    add_text(ax, x, y + 0.18, "AK1 = C2B; AK2 = C5B; AK3–AK9 = GTR Chr03–Chr09", ha="left", fontsize=6.5, color="#555555")
    bar_x = x + 3.3
    unit_w = width / 9.0
    for index in range(1, 10):
        key = f"AK{index}"
        ax.add_patch(Rectangle((bar_x + unit_w * (index - 1), y + 0.09), unit_w, 0.19, facecolor=COLORS[key], edgecolor="white", linewidth=0.7))
        add_text(ax, bar_x + unit_w * (index - 0.5), y + 0.185, key, fontsize=6.2, color="white", fontweight="bold")
    ax.add_patch(Rectangle((bar_x, y + 0.09), unit_w * 9, 0.19, fill=False, edgecolor="#333333", linewidth=0.65))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures"))
    args = parser.parse_args()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "svg.fonttype": "none",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
    })

    fig = plt.figure(figsize=(13.2, 8.0), facecolor="white")
    grid = fig.add_gridspec(16, 12, left=0.035, right=0.985, top=0.96, bottom=0.055, hspace=0.65, wspace=0.65)
    ax = fig.add_subplot(grid[:, :])
    ax.set_xlim(0, 13.2)
    ax.set_ylim(0, 8.0)
    ax.axis("off")

    add_text(ax, 0.05, 7.78, "Fig. 4 | Ancestral Prunus breakpoint model", ha="left", fontsize=12, fontweight="bold")
    add_text(ax, 13.1, 7.78, "five focal rearrangements", ha="right", fontsize=7, color="#666666")
    draw_ancestor(ax, 0.12, 7.20, 9.3)

    # Shared stem event
    add_text(ax, 0.18, 6.63, "Shared Prunus stem", ha="left", fontsize=9, fontweight="bold", color="#7C3F4A")
    add_text(ax, 0.18, 6.43, "one RCT inherited by the sampled Prunus lineages", ha="left", fontsize=6.8, color="#555555")
    draw_rct_card(ax, 3.55, 5.63, 5.25, 1.22, EVENTS[0])
    add_text(ax, 9.15, 6.18, "AK2–AK9 exchange", ha="left", fontsize=7, fontweight="bold", color="#7C3F4A")
    add_text(ax, 9.15, 5.99, "AKP8 and AKP9 are the reciprocal products", ha="left", fontsize=6.5, color="#555555")
    ax.annotate("", xy=(6.25, 5.60), xytext=(6.25, 5.32), arrowprops={"arrowstyle": "-|>", "color": "#7C3F4A", "linewidth": 1.2})

    # A branch
    ax.add_patch(Rectangle((0.08, 2.35), 12.98, 2.72, facecolor="#F7F4FB", edgecolor="#D5CBE5", linewidth=0.8))
    add_text(ax, 0.22, 4.80, "A-associated pathway", ha="left", fontsize=9, fontweight="bold", color="#66518C")
    add_text(ax, 0.22, 4.59, "two RCTs followed by one EEJ", ha="left", fontsize=6.8, color="#66518C")
    draw_rct_card(ax, 2.45, 2.66, 3.15, 1.78, EVENTS[1])
    draw_rct_card(ax, 5.78, 2.66, 3.15, 1.78, EVENTS[2])
    draw_eej_card(ax, 9.11, 2.66, 3.15, 1.78, EVENTS[3])
    add_text(ax, 1.01, 3.36, "A branch", fontsize=6.8, color="#66518C", fontweight="bold")
    ax.annotate("", xy=(2.35, 3.55), xytext=(1.18, 3.55), arrowprops={"arrowstyle": "-|>", "color": "#66518C", "linewidth": 1.0})
    ax.annotate("", xy=(5.68, 3.55), xytext=(5.61, 3.55), arrowprops={"arrowstyle": "-|>", "color": "#66518C", "linewidth": 1.0})
    ax.annotate("", xy=(9.01, 3.55), xytext=(8.94, 3.55), arrowprops={"arrowstyle": "-|>", "color": "#66518C", "linewidth": 1.0})
    add_text(ax, 12.50, 3.55, "n = 8", fontsize=8, fontweight="bold", color="#66518C")

    # B branch
    ax.add_patch(Rectangle((0.08, 0.52), 12.98, 1.56, facecolor="#F3F8F4", edgecolor="#C7DDCE", linewidth=0.8))
    add_text(ax, 0.22, 1.82, "B-associated pathway", ha="left", fontsize=9, fontweight="bold", color="#3D7A57")
    add_text(ax, 0.22, 1.61, "one RCT followed by one EEJ", ha="left", fontsize=6.8, color="#3D7A57")
    draw_rct_card(ax, 3.75, 0.66, 3.15, 1.15, EVENTS[4])
    draw_eej_card(ax, 7.12, 0.66, 3.15, 1.15, EVENTS[5])
    ax.annotate("", xy=(3.65, 1.24), xytext=(1.18, 1.24), arrowprops={"arrowstyle": "-|>", "color": "#3D7A57", "linewidth": 1.0})
    ax.annotate("", xy=(7.02, 1.24), xytext=(6.98, 1.24), arrowprops={"arrowstyle": "-|>", "color": "#3D7A57", "linewidth": 1.0})
    add_text(ax, 10.65, 1.24, "n = 8", fontsize=8, fontweight="bold", color="#3D7A57")

    # Explicit event-count summary and scope note.
    ax.add_patch(Rectangle((0.12, 0.10), 12.9, 0.25, facecolor="#FFFFFF", edgecolor="#D9DDE2", linewidth=0.5))
    add_text(ax, 0.25, 0.225, "Event summary", ha="left", fontsize=6.7, fontweight="bold")
    add_text(ax, 1.55, 0.225, "Prunus shared: 1 RCT", ha="left", fontsize=6.7, color="#7C3F4A")
    add_text(ax, 4.05, 0.225, "A-associated: 2 RCT + 1 EEJ", ha="left", fontsize=6.7, color="#66518C")
    add_text(ax, 7.25, 0.225, "B-associated: 1 RCT + 1 EEJ", ha="left", fontsize=6.7, color="#3D7A57")
    add_text(ax, 10.45, 0.225, "schematic, not bp scale", ha="left", fontsize=6.7, color="#666666")

    fig.savefig(out / "figure4_karyotype_breakpoint_events.svg", bbox_inches="tight")
    fig.savefig(out / "figure4_karyotype_breakpoint_events.png", dpi=600, bbox_inches="tight")
    # Keep a machine-readable event manifest beside the figure.
    with (out / "figure4_karyotype_breakpoint_events.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(EVENTS[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(EVENTS)
    print(f"[OK] {out / 'figure4_karyotype_breakpoint_events.svg'}")
    print(f"[OK] {out / 'figure4_karyotype_breakpoint_events.png'}")


if __name__ == "__main__":
    main()
