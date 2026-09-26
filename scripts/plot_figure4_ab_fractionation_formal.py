#!/usr/bin/env python3
"""Build a formal, source-data-traceable Figure 4 for Prunus A/B fractionation.

The script intentionally separates observed/anchored states from biological
gene-loss calls.  It uses the strict POInT posterior pillar table, the formal
FraBias 100-gene windows, the POInT tree-mapped expectations, and the existing
pan-TE composition table.  All plotted output is generated with matplotlib.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.image import imread
import numpy as np


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
    "svg.fonttype": "none",
    "pdf.fonttype": 42,
    "font.size": 7,
    "axes.linewidth": 0.75,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "legend.frameon": False,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
})

SPECIES = ("PZI", "PSP", "PPA")
COPIES = ("A", "B")
TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B")
FRA_BIAS_TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")
STATE_ORDER = ("AB_anchor", "A_only_anchor", "B_only_anchor", "neither_anchor", "unresolved")
STATE_LABELS = {
    "AB_anchor": "A+B anchored",
    "A_only_anchor": "A-only anchored",
    "B_only_anchor": "B-only anchored",
    "neither_anchor": "neither anchored",
    "unresolved": "unresolved",
}
STATE_COLORS = {
    "AB_anchor": "#4C78A8",
    "A_only_anchor": "#59A14F",
    "B_only_anchor": "#E15759",
    "neither_anchor": "#9D9D9D",
    "unresolved": "#D8C7A1",
}
SPECIES_COLORS = {"PZI": "#2F6B9A", "PSP": "#B24C63", "PPA": "#4C8C6B"}
FRA_TRACK_COLORS = {
    "PZI_A": "#D64A4A", "PZI_B": "#D64A4A",
    "PSP_A": "#3478B9", "PSP_B": "#3478B9",
    "PPA_A": "#3C965D", "PPA_B": "#3C965D",
    "PCAM": "#8064A2", "PPE": "#666666",
}
FRA_TRACK_STYLES = {
    "PZI_A": "-", "PZI_B": "--", "PSP_A": "-", "PSP_B": "--",
    "PPA_A": "-", "PPA_B": "--", "PCAM": "-.", "PPE": ":",
}
REFERENCE_TARGETS = ("GTR__Chr09", "GTR__Chr05", "GTR__Chr08")
REFERENCE_TARGET_LABELS = {"GTR__Chr09": "Chr09", "GTR__Chr05": "Chr05", "GTR__Chr08": "Chr08"}
TE_COLORS = {
    "LTR/Gypsy": "#4C78A8",
    "LTR/Copia": "#F28E2B",
    "Other LTR": "#B279A2",
    "TIR/DNA": "#59A14F",
    "Helitron": "#E15759",
    "LINE": "#79706E",
    "SINE": "#BAB0AC",
    "Other TE": "#9C755F",
}


def read_tsv(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    if fields is None:
        fields = list(rows[0].keys()) if rows else []
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def natural_key(value):
    nums = re.findall(r"\d+", str(value))
    return tuple(int(x) for x in nums) if nums else (10 ** 9, str(value))


def wilson(k, n, z=1.959963984540054):
    if n <= 0:
        return 0.0, 0.0, 0.0
    p = float(k) / float(n)
    den = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / den
    half = z * math.sqrt((p * (1.0 - p) / n) + z * z / (4.0 * n * n)) / den
    return p, max(0.0, centre - half), min(1.0, centre + half)


def bootstrap_mean(values, groups, seed=13, n_boot=2000):
    """Block bootstrap a mean using chromosome/target groups."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    unique = np.unique(groups)
    if len(unique) < 2:
        value = float(np.nanmean(values)) if len(values) else float("nan")
        return value, float("nan"), float("nan")
    grouped = [values[groups == group] for group in unique]
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(n_boot):
        picked = rng.integers(0, len(grouped), size=len(grouped))
        estimates.append(float(np.nanmean(np.concatenate([grouped[i] for i in picked]))))
    return float(np.nanmean(values)), float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5))


def bootstrap_paired_delta(samples, seed=13, n_boot=2000):
    samples = np.asarray(samples, dtype=float)
    samples = samples[np.isfinite(samples)]
    if len(samples) < 2:
        value = float(np.nanmean(samples)) if len(samples) else float("nan")
        return value, float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    estimates = [float(np.mean(rng.choice(samples, size=len(samples), replace=True))) for _ in range(n_boot)]
    return float(np.mean(samples)), float(np.percentile(estimates, 2.5)), float(np.percentile(estimates, 97.5))


def presence(row, species, copy):
    return bool(row.get(f"{species}_{copy}_gene", "").strip()) and row.get(f"{species}_{copy}_gene", "").strip().upper() != "NONE"


def anchor_state(row, species, posterior_cutoff):
    posterior = float(row.get("max_posterior", "nan"))
    if not np.isfinite(posterior) or posterior < posterior_cutoff:
        return "unresolved"
    a = presence(row, species, "A")
    b = presence(row, species, "B")
    if a and b:
        return "AB_anchor"
    if a:
        return "A_only_anchor"
    if b:
        return "B_only_anchor"
    return "neither_anchor"


def lineage_event(row, copy, posterior_cutoff):
    posterior = [float(row.get("max_posterior", "nan"))]
    if not np.isfinite(posterior[0]) or posterior[0] < posterior_cutoff:
        return "unresolved"
    flags = tuple(presence(row, species, copy) for species in SPECIES)
    return {
        (True, True, True): "retained_in_all_three",
        (False, True, True): "PZI-specific anchor loss",
        (True, False, True): "PSP-specific anchor loss",
        (True, True, False): "PPA-specific anchor loss",
        (False, False, True): "PZI+PSP shared anchor loss",
        (False, True, False): "PZI+PPA shared anchor loss",
        (True, False, False): "PSP+PPA shared anchor loss",
        (False, False, False): "older/shared anchor loss",
    }[flags]


def build_panel_a(posterior_rows, outdir, cutoff):
    summary_rows = []
    detail_rows = []
    for row in posterior_rows:
        out = dict(row)
        for species in SPECIES:
            state = anchor_state(row, species, cutoff)
            out[f"{species}_state"] = state
        for copy in COPIES:
            out[f"{copy}_lineage_event"] = lineage_event(row, copy, cutoff)
        detail_rows.append(out)
    for species in SPECIES:
        n = len(detail_rows)
        for state in STATE_ORDER:
            k = sum(row[f"{species}_state"] == state for row in detail_rows)
            fraction, low, high = wilson(k, n)
            summary_rows.append({
                "species": species,
                "category": state,
                "label": STATE_LABELS[state],
                "count": k,
                "n_pillars": n,
                "fraction": f"{fraction:.8f}",
                "ci_low": f"{low:.8f}",
                "ci_high": f"{high:.8f}",
                "posterior_cutoff": f"{cutoff:.3f}",
                "interpretation": "strict-POInT pillar anchor state; not independently validated gene loss",
            })
    event_rows = []
    for copy in COPIES:
        for event in sorted({row[f"{copy}_lineage_event"] for row in detail_rows}):
            k = sum(row[f"{copy}_lineage_event"] == event for row in detail_rows)
            fraction, low, high = wilson(k, len(detail_rows))
            event_rows.append({
                "copy": copy,
                "event": event,
                "count": k,
                "n_pillars": len(detail_rows),
                "fraction": f"{fraction:.8f}",
                "ci_low": f"{low:.8f}",
                "ci_high": f"{high:.8f}",
                "interpretation": "parsimony-style cross-species anchor-state event; not a validated gene-loss call",
            })
    write_tsv(outdir / "panel4a_anchor_state_by_pillar.tsv", detail_rows)
    write_tsv(outdir / "panel4a_anchor_state_summary.tsv", summary_rows)
    write_tsv(outdir / "panel4a_lineage_event_summary.tsv", event_rows)
    return detail_rows, summary_rows, event_rows


def query_code(seqid):
    match = re.search(r"C(\d+)", str(seqid))
    return match.group(1) if match else str(seqid)


def build_panel_b(window_rows, outdir):
    # Preserve the reference-style track trajectories: one row per formal
    # track, target chromosome and overlapping 100-gene window.  A/B summaries
    # are derived separately from paired windows and are not used to redraw the
    # trajectories.
    target_order = [t for t in REFERENCE_TARGETS if any(r["target_seqid"] == t for r in window_rows)]
    target_order += [t for t in sorted({r["target_seqid"] for r in window_rows}, key=natural_key) if t not in target_order]
    direct_rows = []
    for row in window_rows:
        if row["track_id"] not in FRA_BIAS_TRACKS:
            continue
        direct_rows.append({
            "track_id": row["track_id"],
            "target_seqid": row["target_seqid"],
            "query_seqid": row["query_seqid"],
            "window_index": row["window_index"],
            "window_start_gene": row["window_start_gene"],
            "window_end_gene": row["window_end_gene"],
            "window_center_gene": row["window_center_gene"],
            "syntenic_genes_in_window": row["syntenic_genes_in_window"],
            "window_size": row["window_size"],
            "retention_percent": row["retention_percent"],
            "retention_proportion": row["retention_proportion"],
            "delta_A_minus_B": "",
            "method": "formal quota-filtered syntenic anchor retention; not validated gene retention/loss",
        })
    pair_map = {(r["target_seqid"], r["window_index"], r["track_id"]): float(r["retention_proportion"]) for r in direct_rows}
    for row in direct_rows:
        track = row["track_id"]
        if track.endswith("_A") or track.endswith("_B"):
            species, copy = track.split("_")
            other = "B" if copy == "A" else "A"
            counterpart = pair_map.get((row["target_seqid"], row["window_index"], f"{species}_{other}"))
            if counterpart is not None:
                current = float(row["retention_proportion"])
                row["delta_A_minus_B"] = f"{(current - counterpart if copy == 'A' else counterpart - current):.8f}"
    direct_rows.sort(key=lambda r: (target_order.index(r["target_seqid"]), float(r["window_center_gene"]), FRA_BIAS_TRACKS.index(r["track_id"])))

    species_summary = []
    for species in SPECIES:
        a_rows = [r for r in direct_rows if r["track_id"] == f"{species}_A"]
        b_rows = [r for r in direct_rows if r["track_id"] == f"{species}_B"]
        a_key = {(r["target_seqid"], r["window_index"]): r for r in a_rows}
        b_key = {(r["target_seqid"], r["window_index"]): r for r in b_rows}
        common = sorted(set(a_key) & set(b_key), key=lambda x: (target_order.index(x[0]), int(x[1])))
        vals_a = np.asarray([float(a_key[k]["retention_proportion"]) for k in common])
        vals_b = np.asarray([float(b_key[k]["retention_proportion"]) for k in common])
        groups = np.asarray([k[0] for k in common])
        mean_a, low_a, high_a = bootstrap_mean(vals_a, groups, seed=11 + SPECIES.index(species))
        mean_b, low_b, high_b = bootstrap_mean(vals_b, groups, seed=21 + SPECIES.index(species))
        mean_d, low_d, high_d = bootstrap_mean(vals_a - vals_b, groups, seed=31 + SPECIES.index(species))
        species_summary.append({
            "species": species,
            "n_windows": len(common),
            "retention_A": f"{mean_a:.8f}", "retention_A_ci_low": f"{low_a:.8f}", "retention_A_ci_high": f"{high_a:.8f}",
            "retention_B": f"{mean_b:.8f}", "retention_B_ci_low": f"{low_b:.8f}", "retention_B_ci_high": f"{high_b:.8f}",
            "delta_A_minus_B": f"{mean_d:.8f}", "delta_ci_low": f"{low_d:.8f}", "delta_ci_high": f"{high_d:.8f}",
            "interpretation": "descriptive 100-gene-window summary; block bootstrap by GTR target chromosome",
        })
    write_tsv(outdir / "panel4b_window_summary.tsv", direct_rows)
    write_tsv(outdir / "panel4b_species_summary.tsv", species_summary)
    return direct_rows, species_summary, target_order


def build_panel_c(point_rows, point_tree_path, outdir):
    rows = []
    for row in point_rows:
        if row.get("copy_track") not in ("C1", "C2"):
            continue
        n = float(row["initial_pillars"])
        loss = float(row["expected_loss"])
        retained = float(row["expected_retained"])
        rows.append({
            "copy_track": row["copy_track"],
            "copy": row.get("copy") or row.get("state_retained", ""),
            "branch": row["branch"],
            "branch_length": row["branch_length"],
            "initial_pillars": row["initial_pillars"],
            "expected_loss": row["expected_loss"],
            "expected_retained": row["expected_retained"],
            "retention_ratio": f"{retained / n:.8f}",
            "interpretation": "POInT model-expected loss/retention; not observed gene-loss counts",
        })
    write_tsv(outdir / "panel4c_point_tree_data.tsv", rows)
    params = {"FixRate": "NA", "Copy2Bias": "NA", "SwitchProb": "NA", "Final_ln_likelihood": "NA"}
    if point_tree_path and point_tree_path.is_file():
        text = point_tree_path.read_text(encoding="utf-8", errors="replace")
        for key in params:
            match = re.search(re.escape(key) + r"\s*[:=]\s*([-+0-9.eE]+)", text)
            if match:
                params[key] = match.group(1)
        match = re.search(r"Final ln likelihood:\s*([-+0-9.eE]+)", text)
        if match:
            params["Final_ln_likelihood"] = match.group(1)
    param_rows = [{"parameter": key, "value": value} for key, value in params.items()]
    write_tsv(outdir / "panel4c_point_parameters.tsv", param_rows)
    return rows, params


def build_panel_d(window_rows, outdir):
    # Convert the long window table into one matrix with a stable chromosome
    # and window order. Empty cells are uninformative, not losses.
    targets = sorted({r["target_seqid"] for r in window_rows}, key=natural_key)
    target_rank = {x: i for i, x in enumerate(targets)}
    keys = sorted({(r["target_seqid"], int(r["window_index"])) for r in window_rows}, key=lambda x: (target_rank[x[0]], x[1]))
    aggregated = defaultdict(list)
    for row in window_rows:
        track = row["track_id"]
        if track in TRACKS:
            aggregated[(track, row["target_seqid"], int(row["window_index"]))].append(float(row["retention_proportion"]))
    rows = []
    for track in TRACKS:
        for target, win in keys:
            values = aggregated.get((track, target, win), [])
            rows.append({
                "track_id": track,
                "target_seqid": target,
                "window_index": win,
                "retention_proportion": "" if not values else f"{np.mean(values):.8f}",
                "n_query_chromosomes": len(values),
                "interpretation": "empty cell means no informative quota-filtered anchor, not validated loss",
            })
    write_tsv(outdir / "panel4d_retention_heatmap.tsv", rows)
    return rows, targets, keys


def build_panel_e(te_rows, outdir):
    # Average the two haplotypes within each A/B subgenome. This is deliberately
    # a genome-level composition context panel because coordinate-level TE
    # intervals for all six formal tracks are not available in this workspace.
    selected = [r for r in te_rows if r.get("group") in SPECIES and r.get("subgenome") in COPIES and r.get("category") != "Non-TE"]
    values = defaultdict(list)
    for row in selected:
        values[(row["group"], row["subgenome"], row["category"])].append(float(row["genome_fraction"]))
    categories = [x for x in TE_COLORS if any((s, c, x) in values for s in SPECIES for c in COPIES)]
    rows = []
    for species in SPECIES:
        for copy in COPIES:
            for category in categories:
                vals = values.get((species, copy, category), [])
                rows.append({
                    "species": species, "subgenome": copy, "category": category,
                    "mean_genome_fraction": "" if not vals else f"{np.mean(vals):.8f}",
                    "n_haplotypes": len(vals),
                    "interpretation": "genome-level TE composition context; no window-level TE–fractionation association tested",
                })
    write_tsv(outdir / "panel4e_te_composition.tsv", rows)
    return rows, categories


def draw_panel_a(ax, summary_rows, event_rows):
    # A 100% stacked bar makes the dominant unresolved and A+B states visible,
    # while the small right-hand text block reports the A-only/B-only contrast.
    plot_states = ("AB_anchor", "A_only_anchor", "B_only_anchor", "neither_anchor", "unresolved")
    y = np.arange(len(SPECIES))
    left = np.zeros(len(SPECIES))
    for state in plot_states:
        values = np.asarray([float(next(r for r in summary_rows if r["species"] == s and r["category"] == state)["fraction"]) for s in SPECIES])
        ax.barh(y, values, left=left, height=0.52, color=STATE_COLORS[state], edgecolor="white", linewidth=0.45, label=STATE_LABELS[state])
        for yi, start, value, species in zip(y, left, values, SPECIES):
            row = next(r for r in summary_rows if r["species"] == species and r["category"] == state)
            if value >= 0.045:
                text_color = "white" if state in ("AB_anchor", "unresolved") else "#1F1F1F"
                label = f"{100 * value:.1f}%"
                if value >= 0.12:
                    label += f"\n{int(row['count']):,}"
                ax.text(start + value / 2.0, yi, label, ha="center", va="center", fontsize=5.0, color=text_color, fontweight="bold", linespacing=0.9)
        left += values
    ax.set_xlim(0, 1.0)
    ax.set_yticks(y, SPECIES)
    ax.invert_yaxis()
    ax.set_xticks(np.linspace(0, 1, 6), ["0", "20", "40", "60", "80", "100"])
    ax.set_xlabel("strict pillar states (%)")
    ax.set_title("A/B anchor-state composition", loc="left", fontsize=8, pad=4)
    ax.text(-0.08, 1.08, "a", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    ax.grid(axis="x", color="#E5E5E5", linewidth=0.5)
    ax.legend(fontsize=5.0, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.16), handlelength=1.2, columnspacing=0.8)
    ax.text(1.02, 0.72, "n = 11,309 strict pillars\nposterior cutoff = 0.95\nUnresolved is an anchor-state\nclassification, not validated gene loss.",
            transform=ax.transAxes, fontsize=5.0, color="#555555", va="top", ha="left")
    contrast = []
    for species in SPECIES:
        a = next(r for r in summary_rows if r["species"] == species and r["category"] == "A_only_anchor")
        b = next(r for r in summary_rows if r["species"] == species and r["category"] == "B_only_anchor")
        contrast.append(f"{species}: {int(a['count']):,} ({100*float(a['fraction']):.1f}%) | {int(b['count']):,} ({100*float(b['fraction']):.1f}%)")
    ax.text(1.02, 0.02, "A-only vs B-only (count; %)\n" + "\n".join(contrast), transform=ax.transAxes, fontsize=4.9, color="#444444", va="bottom", ha="left")


def _label_positions(values, min_gap=6.0):
    order = sorted(range(len(values)), key=lambda i: values[i], reverse=True)
    positions = {}
    previous = 98.0
    for idx in order:
        positions[idx] = min(float(values[idx]), previous - min_gap)
        previous = positions[idx]
    if positions and min(positions.values()) < 4.0:
        for rank, idx in enumerate(order):
            positions[idx] = 96.0 - rank * (88.0 / max(1, len(order) - 1))
    return positions


def _track_label(track, query_seqid):
    suffix = str(query_seqid).split("__", 1)[-1]
    if track in ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B"):
        return f"{track.replace('_', '-') } {suffix}"
    return f"{track} {suffix}"


def draw_panel_b(axes, rows, targets):
    # Match the supplied FraBias reference style: independent chromosome
    # tracks, percentage scale, and direct query-chromosome labels.
    for axis_index, (ax, target) in enumerate(zip(axes, targets)):
        target_rows = [r for r in rows if r["target_seqid"] == target]
        max_x = max([float(r["window_center_gene"]) for r in target_rows], default=1.0)
        label_x = max_x * 1.035
        ax.set_xlim(0, max_x * 1.23)
        ax.set_ylim(0, 100)
        ax.set_yticks([0, 50, 100])
        ax.grid(axis="y", color="#E5E5E5", linewidth=0.5)
        ax.text(0.99, 0.95, REFERENCE_TARGET_LABELS.get(target, target.replace("GTR__", "")), transform=ax.transAxes,
                ha="right", va="top", fontsize=7.0, color="#555555")
        endpoints = []
        for track in FRA_BIAS_TRACKS:
            subset = sorted([r for r in target_rows if r["track_id"] == track], key=lambda r: float(r["window_center_gene"]))
            if not subset:
                continue
            x = np.asarray([float(r["window_center_gene"]) for r in subset])
            y = np.asarray([100.0 * float(r["retention_proportion"]) for r in subset])
            color = FRA_TRACK_COLORS[track]
            ax.plot(x, y, color=color, linestyle=FRA_TRACK_STYLES[track], lw=0.85, alpha=0.95, solid_capstyle="round")
            endpoints.append((track, x[-1], y[-1], _track_label(track, subset[-1]["query_seqid"])))
        labels = _label_positions([item[2] for item in endpoints])
        for i, (track, x_end, y_end, label) in enumerate(endpoints):
            y_label = labels[i]
            ax.plot([x_end, label_x], [y_end, y_label], color=FRA_TRACK_COLORS[track], lw=0.55, clip_on=False)
            ax.text(label_x * 1.01, y_label, label, color=FRA_TRACK_COLORS[track], fontsize=5.1, va="center", ha="left", clip_on=False,
                    bbox=dict(facecolor="white", edgecolor="none", pad=0.25))
        if axis_index < len(axes) - 1:
            ax.set_xticklabels([])
        else:
            ax.set_xlabel("GTR target chromosome gene order", labelpad=4)
        if axis_index == 0:
            ax.set_ylabel("Gene retention (%)")
            ax.set_title("Quota-filtered syntenic anchor retention", loc="left", fontsize=8, pad=5)
            ax.text(-0.12, 1.18, "b", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
        else:
            ax.set_ylabel("")


def draw_panel_c(ax, point_image_path):
    ax.axis("off")
    image = imread(point_image_path)
    ax.imshow(image, aspect="auto")
    # The reused source panel carries its own "d" label. Cover that source
    # identifier and add the Figure 4 panel label in the outer layout.
    ax.add_patch(plt.Rectangle((0.0, 0.92), 0.10, 0.08, transform=ax.transAxes, facecolor="white", edgecolor="none", zorder=3))
    ax.text(-0.04, 1.02, "c", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")


def draw_panel_d(ax, heat_rows, targets):
    tracks = list(TRACKS)
    keys = sorted({(r["target_seqid"], int(r["window_index"])) for r in heat_rows}, key=lambda x: (natural_key(x[0]), x[1]))
    target_rank = {t: i for i, t in enumerate(targets)}
    keys = sorted(keys, key=lambda x: (target_rank[x[0]], x[1]))
    matrix = np.full((len(tracks), len(keys)), np.nan)
    for i, track in enumerate(tracks):
        for j, key in enumerate(keys):
            row = next((r for r in heat_rows if r["track_id"] == track and r["target_seqid"] == key[0] and int(r["window_index"]) == key[1]), None)
            if row and row["retention_proportion"]:
                matrix[i, j] = float(row["retention_proportion"])
    cmap = mpl.colors.LinearSegmentedColormap.from_list("ret", ["#DCEAF5", "#BFD7EA", "#2166AC"])
    cmap.set_bad("#FFFFFF")
    im = ax.imshow(matrix, aspect="auto", vmin=0, vmax=1, cmap=cmap, interpolation="nearest")
    ax.set_yticks(np.arange(len(tracks)), tracks, fontsize=5.2)
    ax.set_xticks([])
    ax.set_title("Chromosome-level retention map", loc="left", fontsize=8, pad=4)
    ax.text(-0.08, 1.08, "d", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    starts = []
    for j, key in enumerate(keys):
        if j == 0 or key[0] != keys[j - 1][0]:
            starts.append((j, key[0].replace("GTR__", "")))
            ax.axvline(j - 0.5, color="white", lw=1.0)
    for start, label in starts:
        end = next((j for j in range(start + 1, len(keys)) if keys[j][0] != keys[start][0]), len(keys))
        ax.text((start + end - 1) / 2, 0.98, label, transform=ax.get_xaxis_transform(), ha="center", va="top", fontsize=5.2)
    cbar = ax.figure.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("retention", fontsize=5.5)
    cbar.ax.tick_params(labelsize=5)
    ax.text(0.0, -0.18, "White cells are uninformative windows, not validated losses.\nAK breakpoint overlay requires a separately exported high-confidence breakpoint table.",
            transform=ax.transAxes, fontsize=5.1, color="#555555")


def draw_panel_e(ax, te_rows, categories):
    x = np.arange(len(SPECIES) * 2)
    labels = [f"{s}-{c}" for s in SPECIES for c in COPIES]
    bottom = np.zeros(len(x))
    for category in categories:
        vals = []
        for species in SPECIES:
            for copy in COPIES:
                row = next((r for r in te_rows if r["species"] == species and r["subgenome"] == copy and r["category"] == category), None)
                vals.append(float(row["mean_genome_fraction"]) if row and row["mean_genome_fraction"] else 0.0)
        vals = np.asarray(vals)
        ax.bar(x, vals, bottom=bottom, width=0.72, color=TE_COLORS[category], label=category, linewidth=0)
        bottom += vals
    ax.set_xticks(x, labels, rotation=35, ha="right", fontsize=5.5)
    ax.set_ylim(0, 0.75)
    ax.set_ylabel("fraction of genome")
    ax.set_title("TE composition context", loc="left", fontsize=8, pad=4)
    ax.text(-0.08, 1.08, "e", transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")
    ax.legend(fontsize=4.8, ncol=4, loc="upper left", bbox_to_anchor=(0, 1.01), handlelength=1.0, columnspacing=0.7)
    ax.text(0.0, -0.34, "Genome-level composition only; no window-level TE–fractionation association was tested.\nUse associated/consistent with language, not causal TE claims.",
            transform=ax.transAxes, fontsize=5.1, color="#555555")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    root_default = Path(__file__).resolve().parents[1]
    parser.add_argument("--root", type=Path, default=root_default)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--posterior-cutoff", type=float, default=0.95)
    args = parser.parse_args()
    root = args.root
    outdir = args.outdir
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit(f"ERROR: output directory is not empty: {outdir}")
    outdir.mkdir(parents=True, exist_ok=True)

    posterior_path = root / "downstream_analysis/figure2_polyploid/strict_wgd_mapping/05_point_official_reference_exact_r13/03_model_runs/reference_WGD_bias_fix/normalized/point_posterior_summary.tsv"
    window_path = root / "downstream_analysis/figure2_polyploid/strict_wgd_mapping/13_fractbias_repartition_reference_exact/05_fractbias/fractbias_windows.tsv"
    point_path = root / "downstream_analysis/figure2_polyploid/strict_wgd_mapping/18_point_repartition_tree_mapped/panel_d_tree_mapped_source_data.tsv"
    point_tree_path = root / "downstream_analysis/figure2_polyploid/strict_wgd_mapping/15_point_repartition_official_reference_exact/03_model_runs/reference_WGD_bias_fix/point_species_tree.nex.out"
    point_image_path = root / "downstream_analysis/figure2_polyploid/strict_wgd_mapping/18_point_repartition_tree_mapped/figure2_panel_d_point_tree_mapped.png"
    te_path = root / "analysis/prunus_pan_te/postprocess_v2/te_composition.tsv"
    for path in (posterior_path, window_path, point_path, point_tree_path, point_image_path, te_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"ERROR: missing required input: {path}")

    posterior_rows = read_tsv(posterior_path)
    window_rows = read_tsv(window_path)
    point_rows = read_tsv(point_path)
    te_rows = read_tsv(te_path)
    a_detail, a_summary, a_events = build_panel_a(posterior_rows, outdir, args.posterior_cutoff)
    b_rows, b_summary, targets = build_panel_b(window_rows, outdir)
    c_rows, c_params = build_panel_c(point_rows, point_tree_path, outdir)
    d_rows, d_targets, _ = build_panel_d(window_rows, outdir)
    e_rows, e_categories = build_panel_e(te_rows, outdir)

    fig = plt.figure(figsize=(13.6, 10.5), constrained_layout=False)
    gs = fig.add_gridspec(3, 2, width_ratios=[1.0, 1.38], height_ratios=[1.65, 1.20, 1.15], wspace=0.30, hspace=0.58)
    ax_a = fig.add_subplot(gs[0, 0])
    sub_b = gs[0, 1].subgridspec(3, 1, hspace=0.10)
    ax_b = [fig.add_subplot(sub_b[i, 0]) for i in range(3)]
    ax_c = fig.add_subplot(gs[1:, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 1])
    draw_panel_a(ax_a, a_summary, a_events)
    draw_panel_b(ax_b, b_rows, targets)
    draw_panel_c(ax_c, point_image_path)
    draw_panel_d(ax_d, d_rows, d_targets)
    draw_panel_e(ax_e, e_rows, e_categories)
    fig.suptitle("A/B subgenome retention and fractionation in Prunus", fontsize=11, y=0.995)
    fig.text(0.01, 0.005, "Evidence labels: anchor states and quota-filtered syntenic retention are not independently validated gene loss; POInT values are model expectations.", fontsize=5.5, color="#4D4D4D")
    prefix = outdir / "figure4_ab_fractionation_formal"
    fig.savefig(str(prefix) + ".svg", bbox_inches="tight")
    fig.savefig(str(prefix) + ".pdf", bbox_inches="tight")
    fig.savefig(str(prefix) + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(str(prefix) + ".tiff", dpi=600, bbox_inches="tight")
    plt.close(fig)

    (outdir / "analysis_boundary.txt").write_text(
        "Panel 4a: strict POInT pillar anchor states with max_posterior >= 0.95; unresolved otherwise.\n"
        "Panel 4b: formal FraBias quota-filtered syntenic anchor retention in overlapping 100-gene windows; plotted as independent Chr09/Chr05/Chr08 trajectories with direct query-track labels.\n"
        "Panel 4c: reused POInT tree-mapped result image; values are model-expected losses and retained pillars, not observed gene-loss counts.\n"
        "Panel 4d: chromosome-level retention heatmap from the same FraBias windows; blank means uninformative.\n"
        "Panel 4e: genome-level TE composition context. Coordinate-level TE association with fractionation was not tested because complete interval tables for all six formal tracks were unavailable.\n"
        "No panel identifies an extant parental species.\n",
        encoding="ascii",
    )
    (outdir / "software_and_parameters.txt").write_text(
        "backend=Python matplotlib only\n"
        "posterior_cutoff={:.3f}\n".format(args.posterior_cutoff)
        + "panel4a=point_posterior_summary.tsv\n"
        + "panel4b=13_fractbias_repartition_reference_exact/05_fractbias/fractbias_windows.tsv; window_size=100; reference-style chromosome trajectories\n"
        + "panel4c=18_point_repartition_tree_mapped/figure2_panel_d_point_tree_mapped.png; source-data=panel_d_tree_mapped_source_data.tsv\n"
        + "panel4e=te_composition.tsv\n"
        + "confidence_intervals=Wilson for category proportions; target-chromosome block bootstrap for window summaries\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Formal Figure 4 analysis and visualization completed.\n", encoding="ascii")
    print(f"[OK] Wrote formal Figure 4 and source tables to {outdir}")


if __name__ == "__main__":
    main()
