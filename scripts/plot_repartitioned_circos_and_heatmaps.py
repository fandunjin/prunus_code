#!/usr/bin/env python3
"""Reproducible Circos-style and shared-coordinate TE visualisations.

This workflow consumes the revised, MD5-verified EDTA extraction results under
``haploid_repartition/inputs``.  It uses fixed 1-Mb windows, reports the union
of overlapping EDTA intervals (nested insertions are not double-counted), and
keeps the PZI C1A-C8A/C1B-C8B coordinate vocabulary throughout.  Existing PAF
alignments may be supplied for optional PZI-to-PSP/PPA TE-density projection.
No TE discovery is performed here.
"""

from __future__ import annotations

import argparse
import csv
import math
import matplotlib
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Iterable

matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.patches import PathPatch, Wedge
from matplotlib.path import Path as MplPath


HAPS = {
    "pzi": ("hapA1", "hapA2", "hapB1", "hapB2"),
    "ppa": ("hapA", "hapB"),
    "psp": ("hapA1", "hapA2", "hapB1", "hapB2"),
}
SUBGENOME = {
    "hapA1": "A", "hapA2": "A", "hapB1": "B", "hapB2": "B",
    "hapA": "A", "hapB": "B",
}
HAP_COLORS = {
    "hapA1": "#386cb0", "hapA2": "#7fb8d8", "hapB1": "#e6550d", "hapB2": "#fdae6b",
    "hapA": "#386cb0", "hapB": "#e6550d",
}
CLASS_COLORS = {
    "all": "#4d4d4d", "Gypsy": "#c44e52", "Copia": "#4c78a8",
    "TIR/DNA": "#59a14f", "Helitron": "#b279a2", "Other TE": "#8c8c8c",
}
CLASS_KEYS = ("Gypsy", "Copia", "TIR/DNA", "Helitron")
COMPOSITION_KEYS = (*CLASS_KEYS, "Other TE")
WINDOW = 1_000_000


def configure() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.linewidth": 0.7,
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "savefig.facecolor": "white",
    })


def save_pub(fig: plt.Figure, stem: Path, dpi: int = 600) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(f"{stem}.svg", bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.tiff", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def fasta_lengths(path: Path) -> tuple[list[str], dict[str, int]]:
    order: list[str] = []
    lengths: dict[str, int] = {}
    name: str | None = None
    length = 0
    with path.open(encoding="ascii", errors="replace") as handle:
        for raw in handle:
            if raw.startswith(">"):
                if name is not None:
                    lengths[name] = length
                name = raw[1:].split(None, 1)[0]
                order.append(name)
                length = 0
            else:
                length += len(raw.strip())
    if name is not None:
        lengths[name] = length
    if not order or any(lengths[x] <= 0 for x in order):
        raise RuntimeError(f"Invalid FASTA: {path}")
    return order, lengths


def parse_attrs(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for item in text.rstrip(";\n").split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value
    return result


def classify(attrs_text: str) -> tuple[str, str, str]:
    attrs = parse_attrs(attrs_text)
    value = attrs.get("classification") or attrs.get("Classification") or attrs.get("Class") or ""
    if not value and "#" in attrs.get("Name", ""):
        value = attrs["Name"].rsplit("#", 1)[1]
    value = value.replace(" ", "")
    parts = [x for x in re.split(r"[/|]", value) if x]
    parent = parts[0] if parts else "unknown"
    child = parts[1] if len(parts) > 1 else "unknown"
    p = parent.lower()
    if p == "ltr":
        key = "Gypsy" if child.lower() == "gypsy" else "Copia" if child.lower() == "copia" else "Other LTR"
    elif p in {"dna", "tir", "mite"}:
        key = "TIR/DNA"
    elif p in {"rc", "helitron"}:
        key = "Helitron"
    else:
        key = "Other TE"
    return key, parent, child


def union_length(intervals: Iterable[tuple[int, int]], clip_start: int | None = None,
                 clip_end: int | None = None) -> int:
    clipped: list[tuple[int, int]] = []
    for start, end in intervals:
        if clip_start is not None:
            start = max(start, clip_start)
        if clip_end is not None:
            end = min(end, clip_end)
        if start <= end:
            clipped.append((start, end))
    if not clipped:
        return 0
    clipped.sort()
    total = 0
    start, end = clipped[0]
    for left, right in clipped[1:]:
        if left <= end + 1:
            end = max(end, right)
        else:
            total += end - start + 1
            start, end = left, right
    return total + end - start + 1


def parse_gff(path: Path) -> tuple[dict[str, list[tuple[int, int, str, str]]], dict[str, dict[str, list[tuple[int, int]]]]]:
    """Return all features and class-specific interval lists by chromosome."""
    records: dict[str, list[tuple[int, int, str, str]]] = defaultdict(list)
    by_class: dict[str, dict[str, list[tuple[int, int]]]] = {
        "all": defaultdict(list), **{key: defaultdict(list) for key in COMPOSITION_KEYS}
    }
    with path.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            if raw.startswith("#") or not raw.strip():
                continue
            fields = raw.rstrip("\r\n").split("\t")
            if len(fields) != 9:
                continue
            try:
                start, end = int(fields[3]), int(fields[4])
            except ValueError:
                continue
            if start < 1 or end < start:
                continue
            key, parent, child = classify(fields[8])
            records[fields[0]].append((start, end, key, parse_attrs(fields[8]).get("Name", "")))
            by_class["all"][fields[0]].append((start, end))
            if key in COMPOSITION_KEYS:
                by_class[key][fields[0]].append((start, end))
    return records, by_class


def disjoint_category_lengths(records: dict[str, list[tuple[int, int, str, str]]]) -> dict[str, int]:
    """Assign each covered base to one TE category, avoiding nested double counts.

    EDTA can report nested features.  For a composition bar, each base is
    assigned to the first active category in the fixed priority order; this
    makes the category stack exactly equal to the union of all TE intervals.
    """
    priority = {key: i for i, key in enumerate(COMPOSITION_KEYS)}
    totals = {key: 0 for key in COMPOSITION_KEYS}
    for features in records.values():
        events: dict[int, list[tuple[str, str]]] = defaultdict(list)
        for start, end, key, _ in features:
            if key not in priority:
                key = "Other TE"
            events[start].append(("add", key))
            events[end + 1].append(("remove", key))
        active: dict[str, int] = defaultdict(int)
        positions = sorted(events)
        if not positions:
            continue
        previous = positions[0]
        for position in positions:
            if position > previous and active:
                selected = min((key for key, count in active.items() if count > 0), key=priority.get)
                totals[selected] += position - previous
            for action, key in events[position]:
                active[key] += 1 if action == "add" else -1
                if active[key] <= 0:
                    active.pop(key, None)
            previous = position
    return totals


def make_metrics(species: str, hap: str, fasta: Path, gff: Path) -> list[dict[str, object]]:
    order, lengths = fasta_lengths(fasta)
    _, by_class = parse_gff(gff)
    rows: list[dict[str, object]] = []
    for chrom in order:
        length = lengths[chrom]
        for start0 in range(0, length, WINDOW):
            end0 = min(length, start0 + WINDOW)
            # GFF3 is one-based inclusive; convert the half-open window to
            # inclusive coordinates for exact union-length calculations.
            wstart, wend = start0 + 1, end0
            all_bp = union_length(by_class["all"].get(chrom, []), wstart, wend)
            row: dict[str, object] = {
                "species": species, "haplotype": hap, "chromosome": chrom,
                "start_bp": wstart, "end_bp": wend, "window_bp": end0 - start0,
                "all_te_bp": all_bp, "all_te_fraction": all_bp / (end0 - start0),
            }
            for key in CLASS_KEYS:
                bp = union_length(by_class[key].get(chrom, []), wstart, wend)
                row[f"{key}_bp"] = bp
                row[f"{key}_fraction"] = bp / (end0 - start0)
            rows.append(row)
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"No rows to write: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)


def input_paths(input_root: Path, species: str, hap: str) -> tuple[Path, Path]:
    base = input_root / species / "edta" / hap
    if not (base / f"{species}_{hap}.fa.mod.EDTA.TEanno.gff3").is_file():
        base = input_root / species / "repeat" / f"{species}_{hap}_edta"
    return base / f"{species}_{hap}.fa", base / f"{species}_{hap}.fa.mod.EDTA.TEanno.gff3"


def load_all_metrics(input_root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for species, haps in HAPS.items():
        for hap in haps:
            fasta, gff = input_paths(input_root, species, hap)
            if not fasta.is_file() or not gff.is_file():
                raise RuntimeError(f"Missing revised EDTA input for {species}/{hap}: {fasta} or {gff}")
            print(f"[READ] {species}/{hap}: {gff}")
            rows.extend(make_metrics(species, hap, fasta, gff))
    return rows


def group_rows(rows: list[dict[str, object]], species: str, hap: str) -> list[dict[str, object]]:
    return [x for x in rows if x["species"] == species and x["haplotype"] == hap]


def draw_circos(species: str, rows: list[dict[str, object]], out_root: Path) -> None:
    haps = HAPS[species]
    lengths: dict[tuple[str, str], int] = {}
    chroms: dict[tuple[str, str], list[str]] = {}
    for hap in haps:
        sub = group_rows(rows, species, hap)
        chroms[(species, hap)] = list(dict.fromkeys(str(x["chromosome"]) for x in sub))
        for x in sub:
            chroms_key = (species, hap)
            lengths[(hap, str(x["chromosome"]))] = lengths.get((hap, str(x["chromosome"])), 0) + int(x["window_bp"])
    total = sum(lengths.values())
    gap = math.radians(1.4)
    usable = 2 * math.pi - gap * sum(len(chroms[(species, h)]) for h in haps)
    spans: dict[tuple[str, str], tuple[float, float]] = {}
    angle = math.radians(90)
    for hap in haps:
        for chrom in chroms[(species, hap)]:
            delta = usable * lengths[(hap, chrom)] / total
            spans[(hap, chrom)] = (angle, angle + delta)
            angle += delta + gap

    fig, ax = plt.subplots(figsize=(8.0, 8.0), subplot_kw={"projection": "polar"})
    ax.set_theta_direction(-1); ax.set_theta_offset(0); ax.set_axis_off(); ax.set_ylim(0, 1.25)
    # Ideogram and chromosome labels.
    for hap in haps:
        for chrom in chroms[(species, hap)]:
            a0, a1 = spans[(hap, chrom)]
            color = HAP_COLORS[hap]
            ax.bar((a0 + a1) / 2, 0.045, width=a1 - a0, bottom=1.025,
                   color=color, edgecolor="white", linewidth=0.35, align="center")
            if a1 - a0 > math.radians(4):
                mid_deg = math.degrees((a0 + a1) / 2)
                rotation = mid_deg - 90
                if rotation < -90:
                    rotation += 180
                if rotation > 90:
                    rotation -= 180
                ax.text((a0 + a1) / 2, 1.10, chrom, fontsize=6, rotation=rotation,
                        rotation_mode="anchor", ha="center", va="center")
        # Hap label at the centre of its chromosome span.
        a0 = spans[(hap, chroms[(species, hap)][0])][0]
        a1 = spans[(hap, chroms[(species, hap)][-1])][1]
        ax.text((a0 + a1) / 2, 1.19, hap, color=HAP_COLORS[hap], fontsize=8, fontweight="bold",
                ha="center", va="center")

    # All-TE and class tracks.  Every track is 1-Mb and therefore directly
    # comparable across haplotypes; class tracks use fixed global maxima.
    track_specs = [("all_te_fraction", "All TE", 0.945, 0.055, "viridis", 1.0)]
    track_specs += [(f"{key}_fraction", key, 0.0, 0.0, None, 0.35) for key in CLASS_KEYS]
    for index, hap in enumerate(haps):
        base = 0.38 + (len(haps) - index - 1) * 0.13
        sub = group_rows(rows, species, hap)
        for key, label, _, _, cmap_name, vmax in track_specs[:1]:
            values = {(str(x["chromosome"]), int(x["start_bp"])): float(x[key]) for x in sub}
            for chrom in chroms[(species, hap)]:
                a0, a1 = spans[(hap, chrom)]
                chr_rows = [x for x in sub if x["chromosome"] == chrom]
                for x in chr_rows:
                    frac0 = (int(x["start_bp"]) - 1) / lengths[(hap, chrom)]
                    frac1 = int(x["end_bp"]) / lengths[(hap, chrom)]
                    aa0, aa1 = a0 + (a1-a0)*frac0, a0 + (a1-a0)*frac1
                    val = min(1.0, max(0.0, float(x[key])))
                    ax.bar((aa0+aa1)/2, 0.046 * val, width=aa1-aa0, bottom=base,
                           color=plt.get_cmap(cmap_name)(val), edgecolor="none", align="center")
        for j, key in enumerate(CLASS_KEYS):
            r = base - 0.015 - j * 0.014
            submax = 0.35
            for chrom in chroms[(species, hap)]:
                a0, a1 = spans[(hap, chrom)]
                for x in [x for x in sub if x["chromosome"] == chrom]:
                    frac0 = (int(x["start_bp"]) - 1) / lengths[(hap, chrom)]
                    frac1 = int(x["end_bp"]) / lengths[(hap, chrom)]
                    aa0, aa1 = a0 + (a1-a0)*frac0, a0 + (a1-a0)*frac1
                    value = min(submax, max(0.0, float(x[f"{key}_fraction"]))) / submax
                    ax.bar((aa0+aa1)/2, 0.009 * value, width=aa1-aa0, bottom=r,
                           color=CLASS_COLORS[key], edgecolor="none", align="center", alpha=0.92)

    # Link homologous chromosomes only within phase pairs. These are visual
    # anchors, not inferred rearrangement calls.
    if len(haps) >= 4:
        for left, right in (("hapA1", "hapA2"), ("hapB1", "hapB2")):
            for chrom in set(chroms[(species, left)]) & set(chroms[(species, right)]):
                a0 = sum(spans[(left, chrom)]) / 2; a1 = sum(spans[(right, chrom)]) / 2
                verts = [(a0, 0.34), (a0, 0.26), (a1, 0.26), (a1, 0.34)]
                path = MplPath(verts, [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4])
                ax.add_patch(PathPatch(path, transform=ax.transData, facecolor="none", edgecolor="#8c8c8c", lw=0.35, alpha=0.28))
    # Compact legend.
    handles = [plt.Line2D([0], [0], color=HAP_COLORS[h], lw=4, label=h) for h in haps]
    handles += [plt.Line2D([0], [0], color=CLASS_COLORS[k], lw=4, label=k) for k in CLASS_KEYS]
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.05), ncol=4, frameon=False, fontsize=7)
    ax.set_title(f"{species.upper()} revised haplotype genome overview", fontsize=11, pad=20, fontweight="bold")
    save_pub(fig, out_root / f"{species}_circos")


def draw_species_composition(rows: list[dict[str, object]], out_root: Path, input_root: Path) -> None:
    species_order = ("pzi", "psp", "ppa")
    labels = {"pzi": "PZI", "psp": "PSP", "ppa": "PPA"}
    means: dict[str, dict[str, float]] = {}
    for species in species_order:
        pooled_bp: dict[str, int] = defaultdict(int)
        pooled_denom = 0
        for hap in HAPS[species]:
            fasta, gff = input_paths(input_root, species, hap)
            _, lengths = fasta_lengths(fasta)
            records, _ = parse_gff(gff)
            category_bp = disjoint_category_lengths(records)
            denom = sum(lengths.values())
            for key in COMPOSITION_KEYS:
                pooled_bp[key] += category_bp[key]
            pooled_bp["Non-TE"] += max(0, denom - sum(category_bp.values()))
            pooled_denom += denom
        means[species] = {key: pooled_bp[key] / pooled_denom for key in (*COMPOSITION_KEYS, "Non-TE")}
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    x = np.arange(len(species_order)); bottom = np.zeros(len(x))
    order = ["Gypsy", "Copia", "TIR/DNA", "Helitron", "Other TE", "Non-TE"]
    colors = {**CLASS_COLORS, "Non-TE": "#e6e6e6"}
    for key in order:
        vals = np.array([means[s][key] for s in species_order])
        ax.bar(x, vals, bottom=bottom, color=colors[key], width=0.62, label=key)
        bottom += vals
    ax.set_xticks(x, [labels[s] for s in species_order]); ax.set_ylim(0, 1.02)
    ax.set_ylabel("Genome fraction"); ax.set_title("TE composition from revised independent EDTA annotations", fontsize=9)
    ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0))
    ax.legend(frameon=False, ncol=3, fontsize=6, loc="upper center", bbox_to_anchor=(0.5, -0.12))
    save_pub(fig, out_root / "species_te_composition")

    comp_rows = []
    for species in species_order:
        for key in order:
            comp_rows.append({"species": species, "category": key, "pooled_genome_fraction": f"{means[species][key]:.8f}"})
    write_tsv(out_root / "species_te_composition.tsv", comp_rows)


def shared_axis(rows: list[dict[str, object]], pzi_haps: tuple[str, ...]) -> tuple[list[str], dict[str, int], np.ndarray, list[str]]:
    chromosomes = [f"C{i}A" for i in range(1, 9)] + [f"C{i}B" for i in range(1, 9)]
    max_windows: dict[str, int] = {}
    for chrom in chromosomes:
        max_windows[chrom] = max([sum(1 for x in rows if x["species"] == "pzi" and x["haplotype"] == h and x["chromosome"] == chrom) for h in pzi_haps] + [0])
    offsets = {}; cursor = 0
    for chrom in chromosomes:
        offsets[chrom] = cursor; cursor += max_windows[chrom]
    matrix = np.full((len(pzi_haps), cursor), np.nan)
    for i, hap in enumerate(pzi_haps):
        for chrom in chromosomes:
            vals = [x for x in rows if x["species"] == "pzi" and x["haplotype"] == hap and x["chromosome"] == chrom]
            for j, x in enumerate(vals): matrix[i, offsets[chrom] + j] = float(x["all_te_fraction"])
    return chromosomes, offsets, matrix, [h for h in pzi_haps]


def draw_shared_heatmap(rows: list[dict[str, object]], out_root: Path) -> None:
    pzi_haps = HAPS["pzi"]
    chromosomes, offsets, all_matrix, labels = shared_axis(rows, pzi_haps)
    class_matrices = {"All TE": all_matrix}
    for key in CLASS_KEYS:
        matrix = np.full_like(all_matrix, np.nan)
        for i, hap in enumerate(pzi_haps):
            for chrom in chromosomes:
                vals = [x for x in rows if x["species"] == "pzi" and x["haplotype"] == hap and x["chromosome"] == chrom]
                for j, x in enumerate(vals): matrix[i, offsets[chrom] + j] = float(x[f"{key}_fraction"])
        class_matrices[key] = matrix
    fig, axes = plt.subplots(5, 1, figsize=(12.5, 6.3), sharex=True, gridspec_kw={"hspace": 0.10})
    for ax, (label, matrix) in zip(axes, class_matrices.items()):
        im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="YlOrRd", vmin=0, vmax=1 if label == "All TE" else 0.35)
        ax.set_yticks(range(len(labels)), ["PZI-" + x.replace("hap", "") for x in labels])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=24, fontsize=7)
        ax.tick_params(axis="y", length=0)
        for chrom in chromosomes[1:]: ax.axvline(offsets[chrom] - 0.5, color="white", lw=0.5)
        fig.colorbar(im, ax=ax, fraction=0.012, pad=0.006, aspect=30)
    centers = [offsets[c] + max(0, int(np.ceil(max(1, max([sum(1 for x in rows if x["species"] == "pzi" and x["haplotype"] == h and x["chromosome"] == c) for h in pzi_haps])) / 2)) - 0.5) for c in chromosomes]
    axes[-1].set_xticks(centers, chromosomes, rotation=45, ha="right", fontsize=7)
    axes[-1].set_xlabel("Shared PZI chromosome coordinate; each cell is a fixed 1-Mb window")
    fig.suptitle("PZI chromosome-resolved TE density across revised haplotypes", y=0.995, fontsize=11, fontweight="bold")
    save_pub(fig, out_root / "pzi_shared_coordinate_te_heatmap")

    heat_rows = []
    for label, matrix in class_matrices.items():
        for i, hap in enumerate(pzi_haps):
            for chrom in chromosomes:
                vals = [x for x in rows if x["species"] == "pzi" and x["haplotype"] == hap and x["chromosome"] == chrom]
                for j, x in enumerate(vals):
                    heat_rows.append({"track": label, "haplotype": hap, "chromosome": chrom, "window_index": j, "start_bp": x["start_bp"], "end_bp": x["end_bp"], "fraction": x["all_te_fraction"] if label == "All TE" else x[f"{label}_fraction"]})
    write_tsv(out_root / "pzi_shared_coordinate_te_heatmap.tsv", heat_rows)


def draw_pzi_profile(rows: list[dict[str, object]], out_root: Path) -> None:
    # A compact chromosome-level profile is useful beside the dense heatmap.
    pzi_rows = [x for x in rows if x["species"] == "pzi"]
    groups = [("A", ("hapA1", "hapA2")), ("B", ("hapB1", "hapB2"))]
    fig, axes = plt.subplots(2, 1, figsize=(10.2, 4.6), sharex=False, gridspec_kw={"hspace": 0.36})
    for ax, (sg, haps) in zip(axes, groups):
        for key in ("all_te_fraction", "Gypsy_fraction", "Copia_fraction", "TIR/DNA_fraction", "Helitron_fraction"):
            means = []
            labels = []
            for i in range(1, 9):
                chrom = f"C{i}{sg}"
                vals = [float(x[key]) for x in pzi_rows if x["chromosome"] == chrom and x["haplotype"] in haps]
                means.append(float(np.mean(vals)) if vals else np.nan); labels.append(f"C{i}{sg}")
            if key == "all_te_fraction":
                ax.plot(labels, means, color="#252525", lw=1.8, marker="o", ms=3, label="All TE")
            else:
                ax.plot(labels, means, lw=1.0, marker="o", ms=2.4, label=key.replace("_fraction", ""))
        ax.set_ylabel(f"PZI-{sg}\nmean fraction", rotation=0, ha="right", va="center", labelpad=42)
        ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(1.0)); ax.grid(axis="y", color="#dddddd", lw=0.4)
        ax.legend(frameon=False, ncol=5, fontsize=6, loc="upper left")
    axes[-1].set_xlabel("Shared chromosome rank")
    fig.suptitle("Chromosome-level TE density profile", fontsize=10, fontweight="bold")
    save_pub(fig, out_root / "pzi_chromosome_te_profile")


def pzi_coordinate_map(input_root: Path, hap: str) -> dict[str, str]:
    """Map legacy expansion-QC sequence IDs to revised PZI chromosomes."""
    path = input_root / "pzi" / "edta" / hap / "chromosome_map.tsv"
    if not path.is_file():
        path = input_root / "pzi" / "repeat" / f"pzi_{hap}_edta" / "chromosome_map.tsv"
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8", errors="replace") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            target = row.get("target_seqid", "")
            source = row.get("map_source_id", "")
            if target and source:
                mapping[source] = target
    return mapping


def load_expansion_overlay(input_root: Path, expansion_root: Path, rows: list[dict[str, object]],
                           out_root: Path) -> list[dict[str, object]]:
    """Project non-aligned candidates and ONT-supported candidates to 1-Mb PZI windows."""
    overlay: list[dict[str, object]] = []
    metric_lookup = {(str(x["haplotype"]), str(x["chromosome"]), int(x["start_bp"])): x for x in rows if x["species"] == "pzi"}
    for hap in HAPS["pzi"]:
        bed = expansion_root / hap / "comparative_alignment" / f"pzi_{hap}.non_aligned_candidates.bed"
        supported = expansion_root / hap / "candidate_support" / f"pzi_{hap}.candidate_ont_supported.tsv"
        if not bed.is_file():
            continue
        supported_ids: set[str] = set()
        if supported.is_file():
            with supported.open(encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for item in reader:
                    if (item.get("passes_ont_support", "").lower() == "yes" or
                            item.get("passes_ont_support", "").lower() in {"true", "1"}):
                        supported_ids.add(item.get("candidate_id", ""))
        mapping = pzi_coordinate_map(input_root, hap)
        intervals: dict[tuple[str, bool], list[tuple[int, int]]] = defaultdict(list)
        with bed.open(encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if not raw.strip() or raw.startswith("#"):
                    continue
                fields = raw.rstrip("\r\n").split("\t")
                if len(fields) < 3:
                    continue
                legacy, start_text, end_text = fields[:3]
                try:
                    start0, end0 = int(start_text), int(end_text)
                except ValueError:
                    continue
                if end0 <= start0:
                    continue
                source = legacy.rsplit("_", 1)[0] if legacy.rsplit("_", 1)[-1] in {"A1", "A2", "B1", "B2"} else legacy
                chrom = mapping.get(source)
                if chrom is None:
                    continue
                candidate_id = fields[3] if len(fields) >= 4 else f"{legacy}_{start0}_{end0}"
                intervals[(chrom, False)].append((start0 + 1, end0))
                if candidate_id in supported_ids:
                    intervals[(chrom, True)].append((start0 + 1, end0))
        for chrom in [f"C{i}{'A' if hap in {'hapA1', 'hapA2'} else 'B'}" for i in range(1, 9)]:
            metric_rows = [x for x in rows if x["species"] == "pzi" and x["haplotype"] == hap and x["chromosome"] == chrom]
            for x in metric_rows:
                start, end = int(x["start_bp"]), int(x["end_bp"])
                candidate_bp = union_length(intervals.get((chrom, False), []), start, end)
                ont_bp = union_length(intervals.get((chrom, True), []), start, end)
                overlay.append({
                    "haplotype": hap, "chromosome": chrom, "start_bp": start, "end_bp": end,
                    "window_bp": int(x["window_bp"]), "unmapped_candidate_bp": candidate_bp,
                    "ont_supported_bp": ont_bp,
                    "unmapped_candidate_fraction": candidate_bp / int(x["window_bp"]),
                    "ont_supported_fraction": ont_bp / int(x["window_bp"]),
                })
    if overlay:
        write_tsv(out_root / "pzi_candidate_ont_overlay_1Mb.tsv", overlay)
    return overlay


def draw_overlay_heatmap(overlay: list[dict[str, object]], rows: list[dict[str, object]], out_root: Path) -> None:
    """Draw candidate and ONT support tracks on the shared PZI coordinate axis."""
    if not overlay:
        return
    pzi_haps = HAPS["pzi"]
    chromosomes, offsets, _, labels = shared_axis(rows, pzi_haps)
    tracks = {
        "Non-aligned candidates": "unmapped_candidate_fraction",
        "ONT-supported candidates": "ont_supported_fraction",
    }
    matrices: dict[str, np.ndarray] = {}
    for label, key in tracks.items():
        matrix = np.full((len(pzi_haps), sum(max(1, max([sum(1 for z in rows if z["species"] == "pzi" and z["haplotype"] == h and z["chromosome"] == c) for h in pzi_haps] + [0])) for c in chromosomes)), np.nan)
        for i, hap in enumerate(pzi_haps):
            for chrom in chromosomes:
                vals = [x for x in overlay if x["haplotype"] == hap and x["chromosome"] == chrom]
                for j, x in enumerate(vals):
                    matrix[i, offsets[chrom] + j] = float(x[key])
        matrices[label] = matrix
    fig, axes = plt.subplots(2, 1, figsize=(12.5, 2.9), sharex=True, gridspec_kw={"hspace": 0.14})
    for ax, (label, matrix) in zip(axes, matrices.items()):
        im = ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap="YlGnBu", vmin=0, vmax=1)
        ax.set_yticks(range(len(labels)), ["PZI-" + x.replace("hap", "") for x in labels])
        ax.set_ylabel(label, rotation=0, ha="right", va="center", labelpad=36, fontsize=7)
        ax.tick_params(axis="y", length=0)
        for chrom in chromosomes[1:]:
            ax.axvline(offsets[chrom] - 0.5, color="white", lw=0.5)
        fig.colorbar(im, ax=ax, fraction=0.012, pad=0.006, aspect=30)
    centers = []
    for chrom in chromosomes:
        count = max([sum(1 for x in rows if x["species"] == "pzi" and x["haplotype"] == h and x["chromosome"] == chrom) for h in pzi_haps] + [1])
        centers.append(offsets[chrom] + (count - 1) / 2)
    axes[-1].set_xticks(centers, chromosomes, rotation=45, ha="right", fontsize=7)
    axes[-1].set_xlabel("Shared PZI chromosome coordinate; each cell is a fixed 1-Mb window")
    fig.suptitle("PZI expansion candidates and ONT support", y=0.995, fontsize=10, fontweight="bold")
    save_pub(fig, out_root / "pzi_shared_coordinate_expansion_qc_heatmap")


def parse_paf_summary(paf: Path, target_fasta: Path, query_fasta: Path, out: Path) -> None:
    """Write an auditable chromosome alignment matrix from an existing PAF."""
    t_order, t_len = fasta_lengths(target_fasta); q_order, q_len = fasta_lengths(query_fasta)
    aligned: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    with paf.open(encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            f = raw.rstrip().split("\t")
            if len(f) < 12: continue
            q, t = f[0], f[5]
            if q not in q_len or t not in t_len: continue
            try: qs, qe, ts, te, mapq = int(f[2]), int(f[3]), int(f[7]), int(f[8]), int(f[11])
            except ValueError: continue
            if qe - qs < 10000 or mapq < 20: continue
            aligned[(t, q)].append((ts, te))
    rows = []
    for t in t_order:
        for q in q_order:
            bp = union_length(aligned.get((t, q), []), 0, t_len[t])
            rows.append({"target_chr": t, "query_chr": q, "aligned_union_bp": bp, "target_length_bp": t_len[t], "query_length_bp": q_len[q], "target_coverage": f"{bp/t_len[t]:.8f}", "query_coverage": f"{bp/q_len[q]:.8f}"})
    write_tsv(out, rows)


def write_readme(out_root: Path, pafs: dict[str, Path]) -> None:
    text = """# Revised PZI/PPA/PSP Circos and TE coordinate analyses

## Method contract

The plots use the revised, MD5-verified EDTA extraction results. EDTA was not
rerun. GFF3 intervals are treated as one-based inclusive coordinates and are
unioned within each 1-Mb window, so nested or overlapping annotations are not
double-counted. `All TE` includes every EDTA feature; class tracks report
LTR/Gypsy, LTR/Copia, TIR/DNA and Helitron. Genome fractions use the complete
FASTA chromosome length as denominator.

The species composition stack additionally reports `Other TE` (including
unclassified LTR and other EDTA classes) and `Non-TE`. For this stack only,
nested intervals are assigned to one category using the fixed order
Gypsy > Copia > TIR/DNA > Helitron > Other TE, so all six segments sum to
exactly 100% of each genome. This prevents the unreported TE classes from
being incorrectly counted as non-TE sequence.

The PZI heatmap uses the shared chromosome vocabulary C1A-C8A and C1B-C8B.
Each cell is a fixed 1-Mb window; missing phase/chromosome combinations are
shown as blank rather than imputed. The Circos-style links between A1-A2 and
B1-B2 are homologous chromosome anchors for orientation only, not structural
variant calls.

Existing PAF alignment matrices are summarized with aligned length >=10 kb
and MAPQ >=20. They are audit tables, not TE-density evidence. ONT-supported
candidate intervals are included in a separate QC overlay after coordinate
projection. SyRI boundaries were not included because the available SyRI file
is a legacy local test rather than a complete revised-coordinate set; they can
be added once all revised haplotype boundaries are available.

Independent EDTA family names are retained for exploratory family-level work.
Cross-species family claims should use a pan-TE mask after the pan-TE GFF is
repartitioned to the same revised chromosome IDs.

PZI expansion-QC candidate BED files were projected through the chromosome
maps when available. The candidate track represents regions not aligned to
the comparator genomes; the ONT track is restricted to candidates marked as
passing both-boundary spanning-read and depth criteria in the supplied QC
tables. These tracks are independent validation evidence and are not added to
the TE union length.
"""
    if pafs:
        text += "\nPAF inputs summarized:\n" + "\n".join(f"- {key}" for key in sorted(pafs)) + "\n"
    (out_root / "README_methods.txt").write_text(text, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("haploid_repartition/inputs"))
    parser.add_argument("--out-root", type=Path, default=Path("haploid_repartition/circos"))
    parser.add_argument("--expansion-root", type=Path, default=Path("analysis/pzi_expansion_qc"))
    parser.add_argument("--pzi-psp-paf", type=Path, default=Path("haploid_repartition/figure/work/pzi_vs_psp.asm10.paf"))
    parser.add_argument("--pzi-ppa-paf", type=Path, default=Path("haploid_repartition/figure/work/pzi_vs_ppa.asm10.paf"))
    parser.add_argument("--psp-ppa-paf", type=Path, default=Path("haploid_repartition/figure/work/psp_vs_ppa.asm10.paf"))
    args = parser.parse_args()
    configure(); args.out_root.mkdir(parents=True, exist_ok=True)
    rows = load_all_metrics(args.input_root)
    write_tsv(args.out_root / "window_metrics_1Mb.tsv", rows)
    for species in HAPS:
        draw_circos(species, rows, args.out_root)
    draw_species_composition(rows, args.out_root, args.input_root)
    draw_shared_heatmap(rows, args.out_root)
    draw_pzi_profile(rows, args.out_root)
    overlay = load_expansion_overlay(args.input_root, args.expansion_root, rows, args.out_root)
    draw_overlay_heatmap(overlay, rows, args.out_root)

    pafs: dict[str, Path] = {}
    if args.pzi_psp_paf.is_file():
        pafs["pzi_vs_psp"] = args.pzi_psp_paf
        parse_paf_summary(args.pzi_psp_paf, args.input_root / "pzi" / "genome" / "pzi_hapA1.fa", args.input_root / "psp" / "genome" / "psp_hapA1.fa", args.out_root / "pzi_vs_psp_alignment_matrix.tsv")
    if args.pzi_ppa_paf.is_file():
        pafs["pzi_vs_ppa"] = args.pzi_ppa_paf
        parse_paf_summary(args.pzi_ppa_paf, args.input_root / "pzi" / "genome" / "pzi_hapA1.fa", args.input_root / "ppa" / "genome" / "ppa_hapA.fa", args.out_root / "pzi_vs_ppa_alignment_matrix.tsv")
    if args.psp_ppa_paf.is_file():
        pafs["psp_vs_ppa"] = args.psp_ppa_paf
        parse_paf_summary(args.psp_ppa_paf, args.input_root / "psp" / "genome" / "psp_hapA1.fa", args.input_root / "ppa" / "genome" / "ppa_hapA.fa", args.out_root / "psp_vs_ppa_alignment_matrix.tsv")
    write_readme(args.out_root, pafs)
    manifest = [{"file": p.name, "path": str(p), "status": "generated"} for p in sorted(args.out_root.iterdir()) if p.is_file()]
    write_tsv(args.out_root / "figure_manifest.tsv", manifest)
    print(f"[DONE] Outputs written to {args.out_root}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (RuntimeError, OSError, ValueError) as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise SystemExit(2)
