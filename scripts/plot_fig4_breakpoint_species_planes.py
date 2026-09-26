#!/usr/bin/env python3
"""Plot model-projected ancestral rearrangement breakpoints by species.

Six species panels are ordered left-to-right as PZI, PSP, PPA, PCAM, PPE and
RCH. All panels share one horizontal reconstructed GTR-centred coordinate,
and all event rows are aligned, so an event is at the same position in every
species in which it is reconstructed. Markers are ancestral-unit footprints,
not single-nucleotide breakpoints.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle


plt.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.65,
        "legend.frameon": False,
    }
)

SPECIES = ["PZI", "PSP", "PPA", "PCAM", "PPE", "RCH"]
SPECIES_COLORS = {
    "PZI": "#4C78A8",
    "PSP": "#59A14F",
    "PPA": "#F28E2B",
    "PCAM": "#B07AA1",
    "PPE": "#9C755F",
    "RCH": "#8F9399",
}
EVENT_COLORS = {
    "RCT_shared": "#8E3B46",
    "RCT_A": "#6C4E8C",
    "EEJ_A": "#3D7A57",
    "RCT_B": "#4C78A8",
    "EEJ_B": "#2F7F78",
}

# Event footprints follow the reconstructed model used for the karyotype
# history. Derived products are represented by all ancestral units they contain.
EVENTS = [
    {
        "event_id": "P-RCT", "title": "Prunus-shared RCT", "class": "RCT_shared", "event": "RCT",
        "breakpoint": "AK2 x AK9", "units": ["C5B", "GTR09"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": True, "PPE": True, "RCH": False},
        "basis": "shared Prunus stem rearrangement; absent from RCH", "evidence_pair": "C5B|GTR09",
    },
    {
        "event_id": "A-RCT1", "title": "A-associated RCT 1", "class": "RCT_A", "event": "RCT",
        "breakpoint": "AK5 x AKP8", "units": ["C5B", "GTR05", "GTR09"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": False, "PPE": False, "RCH": False},
        "basis": "A-associated pathway; AKP8 is represented by AK2+AK9", "evidence_pair": "GTR05|C5B,GTR09;derived-product event",
    },
    {
        "event_id": "A-RCT2", "title": "A-associated RCT 2", "class": "RCT_A", "event": "RCT",
        "breakpoint": "AK1 x reciprocal product", "units": ["C2B", "C5B", "GTR05", "GTR09"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": False, "PPE": False, "RCH": False},
        "basis": "A-associated pathway; AK1 is C2B", "evidence_pair": "C2B|GTR05;derived-product event",
    },
    {
        "event_id": "A-EEJ", "title": "A-associated EEJ", "class": "EEJ_A", "event": "EEJ",
        "breakpoint": "AK8 x residual product", "units": ["C2B", "GTR05", "GTR08"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": False, "PPE": False, "RCH": False},
        "basis": "A-associated end-to-end joining; AK8 is GTR08", "evidence_pair": "GTR08|derived residual product",
    },
    {
        "event_id": "B-RCT", "title": "B-associated RCT", "class": "RCT_B", "event": "RCT",
        "breakpoint": "AKP9 x AK5", "units": ["C5B", "GTR05", "GTR09"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": True, "PPE": True, "RCH": False},
        "basis": "B-associated pathway shared with the PCAM/PPE structural lineage", "evidence_pair": "GTR05|C5B,GTR09;derived-product event",
    },
    {
        "event_id": "B-EEJ", "title": "B-associated EEJ", "class": "EEJ_B", "event": "EEJ",
        "breakpoint": "AKP3 x AKP8", "units": ["C5B", "GTR03", "GTR09"],
        "present": {"PZI": True, "PSP": True, "PPA": True, "PCAM": True, "PPE": True, "RCH": False},
        "basis": "B-associated terminal joining; AKP3=AK3 and AKP8=AK2+AK9", "evidence_pair": "GTR03|C5B,GTR09;derived-product event",
    },
]

UNIT_ORDER = ["C2B", "C5B", "GTR03", "GTR04", "GTR05", "GTR06", "GTR07", "GTR08", "GTR09"]
UNIT_LABELS = {
    "C2B": "AK1\nC2B", "C5B": "AK2\nC5B", "GTR03": "AK3\nGTR03", "GTR04": "AK4\nGTR04",
    "GTR05": "AK5\nGTR05", "GTR06": "AK6\nGTR06", "GTR07": "AK7\nGTR07",
    "GTR08": "AK8\nGTR08", "GTR09": "AK9\nGTR09",
}


def read_lens(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 2:
                lengths[fields[0]] = int(float(fields[1]))
    return lengths


def build_coordinates(lengths: dict[str, int]) -> tuple[dict[str, float], dict[str, float], float]:
    starts: dict[str, float] = {}
    centers: dict[str, float] = {}
    cursor = 0.0
    for unit in UNIT_ORDER:
        starts[unit] = cursor
        centers[unit] = cursor + float(lengths[unit]) / 2.0
        cursor += float(lengths[unit])
    return starts, centers, cursor


def add_reference_track(ax: plt.Axes, starts: dict[str, float], lengths: dict[str, int], total: float) -> None:
    for index, unit in enumerate(UNIT_ORDER):
        start = starts[unit]
        width = lengths[unit]
        face = "#F4F5F6" if index % 2 == 0 else "#FFFFFF"
        ax.add_patch(Rectangle((start, -0.39), width, 0.25, facecolor=face, edgecolor="#D3D7DB", linewidth=0.4, zorder=0))
        if index:
            ax.axvline(start, color="#C7CBD1", linewidth=0.42, zorder=0)
    ax.axvline(total, color="#C7CBD1", linewidth=0.42, zorder=0)


def plot_species_panel(ax: plt.Axes, species: str, centers: dict[str, float], starts: dict[str, float], lengths: dict[str, int], total: float) -> None:
    """Draw one species panel with the shared GTR x-axis and event rows."""
    ax.set_xlim(0, total)
    ax.set_ylim(-0.60, len(EVENTS) - 0.25)
    ax.set_yticks(range(len(EVENTS)))
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2.0, width=0.5, labelsize=4.6, pad=2)
    for row in range(len(EVENTS)):
        ax.axhline(row, color="#E5E7EA", linewidth=0.55, zorder=0)
    add_reference_track(ax, starts, lengths, total)

    # These offsets separate coincident event footprints without changing the
    # common ancestral coordinate used by every species panel.
    event_offset = {"P-RCT": -0.16, "A-RCT1": -0.095, "A-RCT2": -0.03, "A-EEJ": 0.035, "B-RCT": 0.10, "B-EEJ": 0.165}
    for row, event in enumerate(EVENTS):
        color = EVENT_COLORS[str(event["class"])]
        event_x = [centers[str(unit)] for unit in event["units"]]
        left, right = min(event_x), max(event_x)
        dx = event_offset[str(event["event_id"])] * total / 22.0
        if bool(event["present"][species]):
            ax.plot([left + dx, right + dx], [row, row], color=color, linewidth=1.18, alpha=0.92, zorder=2)
            marker = "o" if str(event["event"]) == "RCT" else "s"
            ax.scatter([x + dx for x in event_x], [row] * len(event_x), s=17, color=color, edgecolor="white", linewidth=0.4, zorder=3)
            ax.scatter([left + dx, right + dx], [row, row], s=27, marker=marker, facecolor=color, edgecolor="#FFFFFF", linewidth=0.4, zorder=4)
        else:
            ax.scatter([(left + right) / 2.0 + dx], [row], s=24, marker="x", color="#B9BDC1", linewidth=0.7, zorder=2)

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(species, fontsize=8.4, fontweight="bold", pad=5, color=SPECIES_COLORS[species])
    subtitle = "haplotype-resolved" if species in ("PZI", "PSP", "PPA") else ("B-like diploid" if species in ("PCAM", "PPE") else "outgroup")
    ax.text(0.5, 1.075, subtitle, transform=ax.transAxes, ha="center", va="bottom", fontsize=4.8, color="#666B70")


def write_manifest(path: Path, centers: dict[str, float], events: list[dict[str, object]]) -> None:
    fields = ["event_id", "title", "event", "breakpoint", "units", "unit_centres_bp", "species", "present", "basis", "evidence_pair"]
    rows: list[dict[str, object]] = []
    for event in events:
        for species in SPECIES:
            rows.append(
                {
                    "event_id": event["event_id"], "title": event["title"], "event": event["event"],
                    "breakpoint": event["breakpoint"], "units": ",".join(str(x) for x in event["units"]),
                    "unit_centres_bp": ",".join(f"{centers[str(x)]:.0f}" for x in event["units"]),
                    "species": species, "present": "yes" if event["present"][species] else "no",
                    "basis": event["basis"], "evidence_pair": event["evidence_pair"],
                }
            )
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903"))
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures"))
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    lens_path = root / "karyotype" / "fig4_breakpoint_candidates" / "reference" / "prunus9_reference.lens"
    lengths = read_lens(lens_path)
    missing = [unit for unit in UNIT_ORDER if unit not in lengths]
    if missing:
        raise SystemExit(f"Missing ancestral reference units: {missing}")
    starts, centers, total = build_coordinates(lengths)

    fig, axes_grid = plt.subplots(1, len(SPECIES), figsize=(14.2, 5.2), sharex=True, sharey=True, squeeze=False)
    axes = list(axes_grid[0])
    fig.subplots_adjust(left=0.14, right=0.985, top=0.79, bottom=0.22, wspace=0.055)
    fig.text(0.015, 0.965, "Fig. 4 | Breakpoint correspondence across Prunus and Rosaceae", ha="left", va="top", fontsize=12, fontweight="bold")
    fig.text(0.015, 0.936, "Species panels are ordered left to right; every panel uses the same reconstructed GTR-centred coordinate, so shared breakpoints occupy the same horizontal position.", ha="left", va="top", fontsize=6.8, color="#5E6266")
    fig.text(0.015, 0.918, "Markers show model-projected ancestral-unit event footprints. Crosses indicate that an event is not reconstructed in that lineage.", ha="left", va="top", fontsize=6.8, color="#5E6266")

    for ax, species in zip(axes, SPECIES):
        plot_species_panel(ax, species, centers, starts, lengths, total)
        ax.set_xticks([centers[unit] for unit in UNIT_ORDER])
        ax.set_xticklabels([UNIT_LABELS[unit] for unit in UNIT_ORDER], fontsize=3.9, linespacing=0.85)
        ax.set_xlabel("GTR coordinate", fontsize=5.8, labelpad=8)

    axes[0].set_yticks(range(len(EVENTS)))
    axes[0].set_yticklabels([f"{event['event_id']}\n{event['event']}" for event in EVENTS], fontsize=6.0, ha="right")
    axes[0].tick_params(axis="y", pad=8)
    fig.text(0.055, 0.51, "rearrangement event", rotation=90, ha="center", va="center", fontsize=6.8, color="#555A5E")

    legend_handles = [
        Line2D([0], [0], color=EVENT_COLORS["RCT_shared"], marker="o", linewidth=1.25, markersize=4.2, label="shared Prunus RCT"),
        Line2D([0], [0], color=EVENT_COLORS["RCT_A"], marker="o", linewidth=1.25, markersize=4.2, label="A-associated RCT"),
        Line2D([0], [0], color=EVENT_COLORS["EEJ_A"], marker="s", linewidth=1.25, markersize=4.2, label="A-associated EEJ"),
        Line2D([0], [0], color=EVENT_COLORS["RCT_B"], marker="o", linewidth=1.25, markersize=4.2, label="B-associated RCT"),
        Line2D([0], [0], color=EVENT_COLORS["EEJ_B"], marker="s", linewidth=1.25, markersize=4.2, label="B-associated EEJ"),
        Line2D([0], [0], color="#B9BDC1", marker="x", linewidth=0, markersize=4.2, label="not reconstructed in lineage"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", bbox_to_anchor=(0.58, 0.055), ncol=3, fontsize=6.0, handlelength=1.8, columnspacing=1.1)

    svg = out / "figure4_breakpoint_species_planes.svg"
    png = out / "figure4_breakpoint_species_planes.png"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    write_manifest(out / "figure4_breakpoint_species_planes.tsv", centers, EVENTS)
    plt.close(fig)
    print(f"[OK] {svg}")
    print(f"[OK] {png}")
    print(f"[OK] {out / 'figure4_breakpoint_species_planes.tsv'}")


if __name__ == "__main__":
    main()
