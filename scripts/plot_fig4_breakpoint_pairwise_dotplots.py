#!/usr/bin/env python3
"""Plot one GTR-to-target chromosome dot plot per event and species.

The six columns are PZI, PSP, PPA, PCAM, PPE and RCH. Each cell contains one
target chromosome only. The x-axis is the same GTR-centred nine-unit gene-order
coordinate in every cell; points are WGDI primary-block midpoints. Large
diamonds are candidate junctions whose source-unit pair falls within the
event's reconstructed unit set. The event assignments remain model-projected
and are not claimed to be single-nucleotide calls.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


mpl.rcParams.update(
    {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
        "font.size": 7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.55,
        "legend.frameon": False,
    }
)

SPECIES = ["PZI", "PSP", "PPA", "PCAM", "PPE", "RCH"]
SPECIES_TARGETS = {
    "PZI": ["pzi_hapA1", "pzi_hapA2", "pzi_hapB1", "pzi_hapB2"],
    "PSP": ["psp_hapA1", "psp_hapA2", "psp_hapB1", "psp_hapB2"],
    "PPA": ["ppa_hapA", "ppa_hapB"],
    "PCAM": ["pcam"],
    "PPE": ["ppe"],
    "RCH": ["rch"],
}

UNIT_ORDER = ["C2B", "C5B", "GTR03", "GTR04", "GTR05", "GTR06", "GTR07", "GTR08", "GTR09"]
UNIT_LABELS = {
    "C2B": "AK1/C2B", "C5B": "AK2/C5B", "GTR03": "AK3/GTR03", "GTR04": "AK4/GTR04",
    "GTR05": "AK5/GTR05", "GTR06": "AK6/GTR06", "GTR07": "AK7/GTR07",
    "GTR08": "AK8/GTR08", "GTR09": "AK9/GTR09",
}
UNIT_COLORS = {
    "C2B": "#00AEEF", "C5B": "#DF8442", "GTR03": "#5B9255", "GTR04": "#C9595F",
    "GTR05": "#9E7497", "GTR06": "#6EAAA6", "GTR07": "#D7B75A", "GTR08": "#937366", "GTR09": "#3A5C93",
}
EVENT_COLORS = {
    "shared": "#8E3B46", "A_RCT": "#6C4E8C", "A_EEJ": "#3D7A57", "B_RCT": "#4C78A8", "B_EEJ": "#2F7F78",
}

# Units are the source coordinates of the reconstructed event. For derived
# products, all constituent ancestral units are retained on the common x-axis.
EVENTS = [
    {
        "event_id": "P-RCT", "label": "Prunus-shared RCT", "kind": "shared", "event": "RCT",
        "units": ("C5B", "GTR09"), "presence": {s: s != "RCH" for s in SPECIES}, "allowed": "any",
    },
    {
        "event_id": "A-RCT1", "label": "A-associated RCT 1", "kind": "A_RCT", "event": "RCT",
        "units": ("C5B", "GTR05", "GTR09"), "presence": {s: s in ("PZI", "PSP", "PPA") for s in SPECIES}, "allowed": "A",
    },
    {
        "event_id": "A-RCT2", "label": "A-associated RCT 2", "kind": "A_RCT", "event": "RCT",
        "units": ("C2B", "C5B", "GTR05", "GTR09"), "presence": {s: s in ("PZI", "PSP", "PPA") for s in SPECIES}, "allowed": "A",
    },
    {
        "event_id": "A-EEJ", "label": "A-associated EEJ", "kind": "A_EEJ", "event": "EEJ",
        "units": ("C2B", "GTR05", "GTR08"), "presence": {s: s in ("PZI", "PSP", "PPA") for s in SPECIES}, "allowed": "A",
    },
    {
        "event_id": "B-RCT", "label": "B-associated RCT", "kind": "B_RCT", "event": "RCT",
        "units": ("C5B", "GTR05", "GTR09"), "presence": {s: s != "RCH" for s in SPECIES}, "allowed": "B_or_diploid",
    },
    {
        "event_id": "B-EEJ", "label": "B-associated EEJ", "kind": "B_EEJ", "event": "EEJ",
        "units": ("C5B", "GTR03", "GTR09"), "presence": {s: s != "RCH" for s in SPECIES}, "allowed": "B_or_diploid",
    },
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_lens(path: Path) -> dict[str, tuple[int, int]]:
    result: dict[str, tuple[int, int]] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) >= 3:
                try:
                    result[fields[0]] = (int(float(fields[1])), int(float(fields[2])))
                except ValueError:
                    continue
    return result


def read_reference_units(path: Path) -> dict[str, int]:
    return {row["reference_unit"]: int(float(row["source_gene_count"])) for row in read_tsv(path)}


def species_of_target(target: str) -> str:
    lower = target.lower()
    for species in ("pzi", "psp", "ppa", "pcam", "ppe", "rch"):
        if lower.startswith(species):
            return species.upper()
    return target.upper()


def allowed_target(target: str, species: str, allowed: str) -> bool:
    if target not in SPECIES_TARGETS[species]:
        return False
    if allowed == "any":
        return True
    if allowed == "A":
        return "_hapA" in target
    if allowed == "B_or_diploid":
        return ("_hapB" in target) or target in {"pcam", "ppe"}
    return True


def number(row: dict[str, str], field: str, default: float = 0.0) -> float:
    try:
        return float(row.get(field, ""))
    except (TypeError, ValueError):
        return default


def build_x(row: dict[str, str], unit_counts: dict[str, int], offsets: dict[str, float]) -> float | None:
    unit = row.get("source_unit", "")
    if unit not in offsets:
        return None
    count = max(1, unit_counts.get(unit, 1))
    start = number(row, "source_start_order", 1.0)
    end = number(row, "source_end_order", start)
    return offsets[unit] + (max(1.0, min(start, end)) - 1.0 + abs(end - start) / 2.0) / count


def unit_offsets(unit_counts: dict[str, int]) -> dict[str, float]:
    return {unit: float(index) for index, unit in enumerate(UNIT_ORDER)}


def choose_target(
    rows: list[dict[str, str]], species: str, event: dict[str, object]
) -> tuple[str | None, str | None, int, int]:
    """Select exactly one target chromosome using event-unit anchor support."""
    allowed = str(event["allowed"])
    units = set(str(item) for item in event["units"])
    counts: dict[tuple[str, str], tuple[int, int]] = defaultdict(lambda: (0, 0))
    for row in rows:
        target = row["target_label"]
        if not allowed_target(target, species, allowed) or row["source_unit"] not in units:
            continue
        key = (target, row["target_chromosome"])
        count, total_len = counts[key]
        counts[key] = (count + 1, total_len + int(number(row, "length_genes")))
    if not counts:
        return None, None, 0, 0
    target, chrom = max(counts, key=lambda key: (counts[key][0], counts[key][1], key[0], key[1]))
    count, total_len = counts[(target, chrom)]
    return target, chrom, count, total_len


def plot_pair(ax: plt.Axes, rows: list[dict[str, str]], target: str, chrom: str, lens_cache: dict[str, dict[str, tuple[int, int]]], unit_counts: dict[str, int], offsets: dict[str, float], event: dict[str, object], candidate_rows: list[dict[str, str]]) -> int:
    target_lens = lens_cache[target]
    target_count = max(1, target_lens.get(chrom, (0, 1))[1])
    selected = [row for row in rows if row["target_label"] == target and row["target_chromosome"] == chrom]
    colors: list[str] = []
    xs: list[float] = []
    ys: list[float] = []
    for row in selected:
        x = build_x(row, unit_counts, offsets)
        if x is None:
            continue
        start = number(row, "target_start_order", 1.0)
        end = number(row, "target_end_order", start)
        y = (min(start, end) - 1.0 + abs(end - start) / 2.0) / target_count
        xs.append(x)
        ys.append(y)
        colors.append(UNIT_COLORS.get(row["source_unit"], "#64707B"))
    if xs:
        ax.scatter(xs, ys, s=6.0, c=colors, alpha=0.62, linewidths=0, rasterized=False, zorder=2)

    event_units = set(str(item) for item in event["units"])
    pair_count = 0
    for row in candidate_rows:
        if row["target_label"] != target or row["target_chromosome"] != chrom:
            continue
        pair = [part for part in row.get("source_unit_pair", "").split("|") if part]
        if len(pair) != 2 or len(set(pair)) < 2 or not set(pair).issubset(event_units):
            continue
        pair_x = sum(offsets.get(unit, 0.0) + 0.5 for unit in pair) / len(pair)
        chrom_bp = number(row, "target_chromosome_length_bp", 1.0)
        y = max(0.0, min(1.0, number(row, "candidate_center_bp", 0.0) / max(1.0, chrom_bp)))
        ax.scatter([pair_x], [y], s=26, marker="D", facecolor=EVENT_COLORS[str(event["kind"])], edgecolor="white", linewidth=0.55, zorder=4)
        pair_count += 1

    ax.set_xlim(0, 9)
    ax.set_ylim(0, 1)
    ax.set_xticks([i + 0.5 for i in range(9)])
    ax.set_xticklabels([f"AK{i + 1}" for i in range(9)], fontsize=4.4)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_yticklabels(["0", "0.5", "1.0"], fontsize=4.1)
    ax.tick_params(axis="both", length=1.5, width=0.45, pad=1)
    for i in range(10):
        ax.axvline(i, color="#D1D6DB", linewidth=0.38, zorder=0)
    ax.grid(axis="y", color="#E5E8EB", linewidth=0.38, zorder=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_title(f"{target}:{chrom}", loc="left", fontsize=5.1, fontweight="bold", pad=2)
    return pair_count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("0903"))
    parser.add_argument("--output-dir", type=Path, default=Path("0903/figures"))
    parser.add_argument("--gtr-fasta", type=Path, default=Path("data/gtr/GCA_018257905.1_V1_genomic.fna"))
    args = parser.parse_args()
    root = args.root.resolve()
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    work = root / "karyotype" / "fig4_breakpoint_candidates"
    block_rows = read_tsv(work / "selected_primary_blocks.tsv")
    candidate_rows = read_tsv(work / "junction_candidates.tsv")
    unit_counts = read_reference_units(work / "reference" / "prunus9_reference_units.tsv")
    offsets = unit_offsets(unit_counts)
    lens_cache: dict[str, dict[str, tuple[int, int]]] = {}
    for target in sorted({row["target_label"] for row in block_rows}):
        lens_path = root / "karyotype" / "prepared" / target / f"{target}.lens"
        if lens_path.is_file():
            lens_cache[target] = read_lens(lens_path)

    fig = plt.figure(figsize=(16.5, 12.2))
    grid = fig.add_gridspec(len(EVENTS), len(SPECIES), left=0.07, right=0.995, top=0.86, bottom=0.14, hspace=0.72, wspace=0.12)
    fig.text(0.012, 0.975, "Fig. 4 | WGDI dot plots of shared and subgenome-associated breakpoints", ha="left", va="top", fontsize=12.5, fontweight="bold")
    fig.text(0.012, 0.949, "Each cell shows one GTR-to-target chromosome pair. The six species columns share the same GTR x-coordinate; large diamonds mark matching candidate junctions.", ha="left", va="top", fontsize=6.8, color="#5A6065")
    fig.text(0.012, 0.931, "PZI, PSP and PPA use the phased A/B targets; PCAM and PPE are included for B-associated events; RCH is blank for events not reconstructed in the outgroup.", ha="left", va="top", fontsize=6.8, color="#5A6065")

    manifest: list[dict[str, object]] = []
    for row_index, event in enumerate(EVENTS):
        for col_index, species in enumerate(SPECIES):
            ax = fig.add_subplot(grid[row_index, col_index])
            if col_index == 0:
                ax.set_ylabel("target gene order", fontsize=5.8, labelpad=5)
            event_present = bool(event["presence"][species])
            target, chrom, support_n, support_len = choose_target(block_rows, species, event)
            candidate_n = 0
            if event_present and target in lens_cache and chrom is not None:
                candidate_n = plot_pair(ax, block_rows, target, chrom, lens_cache, unit_counts, offsets, event, candidate_rows)
            else:
                ax.set_xlim(0, 9)
                ax.set_ylim(0, 1)
                ax.set_xticks([i + 0.5 for i in range(9)])
                ax.set_xticklabels([f"AK{i + 1}" for i in range(9)], fontsize=4.4)
                ax.set_yticks([])
                for i in range(10):
                    ax.axvline(i, color="#E1E4E7", linewidth=0.38, zorder=0)
                ax.text(0.5, 0.5, "not reconstructed", transform=ax.transAxes, ha="center", va="center", fontsize=5.0, color="#A0A5AA")
                for spine in ax.spines.values():
                    spine.set_visible(False)
                target, chrom = None, None
            if row_index == 0:
                ax.set_title(species, fontsize=8.0, fontweight="bold", pad=15, color={"PZI": "#4C78A8", "PSP": "#59A14F", "PPA": "#F28E2B", "PCAM": "#B07AA1", "PPE": "#9C755F", "RCH": "#8F9399"}[species])
            if col_index == 0:
                ax.text(-0.35, 0.5, str(event["event_id"]), transform=ax.transAxes, ha="right", va="center", fontsize=6.5, fontweight="bold", color=EVENT_COLORS[str(event["kind"])])
            manifest.append({
                "event_id": event["event_id"], "event": event["event"], "species": species,
                "present_in_model": "yes" if event_present else "no", "target_label": target or "",
                "target_chromosome": chrom or "", "anchor_block_count": support_n,
                "anchor_gene_count": support_len, "matching_candidate_count": candidate_n,
                "gtr_units": ",".join(str(x) for x in event["units"]),
            })

    # Apply the shared x-axis label only to the bottom row.
    for col_index in range(len(SPECIES)):
        bottom = fig.axes[(len(EVENTS) - 1) * len(SPECIES) + col_index]
        bottom.set_xlabel("GTR-centred ancestral gene order", fontsize=6.0, labelpad=7)

    unit_legend = [Patch(facecolor=UNIT_COLORS[u], edgecolor="none", label=UNIT_LABELS[u]) for u in UNIT_ORDER]
    event_legend = [
        Line2D([0], [0], marker="D", color="none", markerfacecolor=EVENT_COLORS[k], markeredgecolor="white", markersize=5.0, label=label)
        for k, label in (("shared", "Prunus-shared RCT"), ("A_RCT", "A-associated RCT"), ("A_EEJ", "A-associated EEJ"), ("B_RCT", "B-associated RCT"), ("B_EEJ", "B-associated EEJ"))
    ]
    fig.legend(handles=unit_legend, loc="lower center", bbox_to_anchor=(0.50, 0.080), ncol=9, fontsize=5.5, handlelength=1.0, columnspacing=0.7)
    fig.legend(handles=event_legend, loc="lower center", bbox_to_anchor=(0.50, 0.035), ncol=5, fontsize=5.7, handlelength=1.0, columnspacing=1.0)
    fig.text(0.012, 0.008, "Point positions use WGDI gene-order blocks and candidate junction intervals; they are not single-nucleotide breakpoint estimates.", fontsize=5.8, color="#686D72", ha="left")

    svg = out / "figure4_breakpoint_pairwise_dotplots.svg"
    png = out / "figure4_breakpoint_pairwise_dotplots.png"
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(png, dpi=600, bbox_inches="tight")
    plt.close(fig)

    with (out / "figure4_breakpoint_pairwise_dotplots.tsv").open("w", encoding="utf-8", newline="") as handle:
        fields = list(manifest[0]) if manifest else []
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"[OK] {svg}")
    print(f"[OK] {png}")
    print(f"[OK] {out / 'figure4_breakpoint_pairwise_dotplots.tsv'}")
    print(f"[INFO] GTR FASTA supplied: {args.gtr_fasta.resolve()} ({'present' if args.gtr_fasta.is_file() else 'not found'})")


if __name__ == "__main__":
    main()
