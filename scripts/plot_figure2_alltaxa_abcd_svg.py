#!/usr/bin/env python3
"""Build the all-taxon Figure 2 a-d composite as editable SVG."""

from __future__ import annotations

import argparse
import csv
import html
import math
from collections import Counter, defaultdict
from pathlib import Path


WIDTH = 2400
HEIGHT = 1800

COLORS = {
    "PZI": "#B84A46",
    "PSP": "#2C64A7",
    "PPA": "#777777",
    "A": "#D85A55",
    "B": "#3F8F9D",
    "PP": "#202020",
    "gCF": "#D08A3C",
    "sCF": "#5B7FCA",
    "both_anchored_in_both": "#5B7FCA",
    "same_A_only_candidate": "#D9544D",
    "same_B_only_candidate": "#42949E",
    "opposite_single_anchor": "#B89BD9",
    "same_ambiguous_state": "#7B7B7B",
    "discordant_or_uninformative": "#D1D1D1",
}

CROSS_SPECIES_COLORS = {
    "PZI": "#202020",
    "PSP": "#3B6FB6",
    "PCAM": "#3F8C75",
    "PPA": "#D4554B",
    "PPE": "#D08A3C",
}

SELECTED_PZI_CHROMOSOMES = [
    {"copy1": "chr01", "pzi_label": "pzi_hapA1", "hap_seqid": "C06_A1", "psp_label": "psp_hapA1", "ppa_label": "ppa_hapA"},
    {"copy1": "chr02", "pzi_label": "pzi_hapA1", "hap_seqid": "C02_A1", "psp_label": "psp_hapA1", "ppa_label": "ppa_hapA"},
    {"copy1": "chr03", "pzi_label": "pzi_hapA1", "hap_seqid": "C03_A1", "psp_label": "psp_hapA1", "ppa_label": "ppa_hapA"},
    {"copy1": "chr12", "pzi_label": "pzi_hapB1", "hap_seqid": "C03_B1", "psp_label": "psp_hapB1", "ppa_label": "ppa_hapB"},
    {"copy1": "chr13", "pzi_label": "pzi_hapB1", "hap_seqid": "C06_B1", "psp_label": "psp_hapB1", "ppa_label": "ppa_hapB"},
    {"copy1": "chr15", "pzi_label": "pzi_hapB1", "hap_seqid": "C01_B1", "psp_label": "psp_hapB1", "ppa_label": "ppa_hapB"},
]

JOINT_ORDER = [
    "both_anchored_in_both",
    "same_A_only_candidate",
    "same_B_only_candidate",
    "opposite_single_anchor",
    "same_ambiguous_state",
    "discordant_or_uninformative",
]

JOINT_LABELS = {
    "both_anchored_in_both": "Both in both",
    "same_A_only_candidate": "Shared A-only",
    "same_B_only_candidate": "Shared B-only",
    "opposite_single_anchor": "Opposite",
    "same_ambiguous_state": "Same multi",
    "discordant_or_uninformative": "Discordant / uninformative",
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ks", required=True)
    parser.add_argument("--topologies", required=True)
    parser.add_argument("--support", required=True)
    parser.add_argument("--reference-loci", required=True)
    parser.add_argument("--fractionation-states", required=True)
    parser.add_argument("--pairs", required=True, help="All-taxon GTR syntenic-pair table")
    parser.add_argument("--jcvi-inputs", required=True, help="JCVI BED/seqids input directory")
    parser.add_argument("--joint-states", required=True)
    parser.add_argument("--candidates", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--ks-max", type=float, default=3.0)
    parser.add_argument("--ks-bandwidth", type=float, default=0.04)
    parser.add_argument("--window-genes", type=int, default=101)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def esc(value):
    return html.escape(str(value), quote=True)


def svg_start():
    return [
        '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">' % (WIDTH, HEIGHT, WIDTH, HEIGHT),
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>'
        'text{font-family:Arial,Helvetica,sans-serif;letter-spacing:0;fill:#202020}'
        '.panel{font-size:38px;font-weight:700}.title{font-size:29px;font-weight:700}'
        '.subtitle{font-size:18px;fill:#555}.axis{font-size:18px}.tick{font-size:16px;fill:#444}'
        '.small{font-size:15px;fill:#444}.tiny{font-size:13px;fill:#555}'
        '.tablehead{font-size:15px;font-weight:700}.tablebody{font-size:14px}'
        '</style>',
    ]


def text(lines, x, y, value, cls="small", anchor="start", weight=None, fill=None):
    attrs = ['class="%s"' % cls, 'x="%.2f"' % x, 'y="%.2f"' % y]
    if anchor != "start":
        attrs.append('text-anchor="%s"' % anchor)
    if weight:
        attrs.append('font-weight="%s"' % weight)
    if fill:
        attrs.append('style="fill:%s"' % fill)
    lines.append("<text %s>%s</text>" % (" ".join(attrs), esc(value)))


def line(lines, x1, y1, x2, y2, stroke="#333333", width=1.5, dash=None, opacity=None):
    attrs = [
        'x1="%.2f"' % x1, 'y1="%.2f"' % y1, 'x2="%.2f"' % x2, 'y2="%.2f"' % y2,
        'stroke="%s"' % stroke, 'stroke-width="%.2f"' % width,
    ]
    if dash:
        attrs.append('stroke-dasharray="%s"' % dash)
    if opacity is not None:
        attrs.append('opacity="%.3f"' % opacity)
    lines.append("<line %s/>" % " ".join(attrs))


def rect(lines, x, y, w, h, fill="none", stroke=None, width=1.0, opacity=None):
    attrs = ['x="%.2f"' % x, 'y="%.2f"' % y, 'width="%.2f"' % w, 'height="%.2f"' % h, 'fill="%s"' % fill]
    if stroke:
        attrs.extend(['stroke="%s"' % stroke, 'stroke-width="%.2f"' % width])
    if opacity is not None:
        attrs.append('opacity="%.3f"' % opacity)
    lines.append("<rect %s/>" % " ".join(attrs))


def circle(lines, x, y, radius, fill, stroke="#333333", width=1.0):
    lines.append('<circle cx="%.2f" cy="%.2f" r="%.2f" fill="%s" stroke="%s" stroke-width="%.2f"/>' % (x, y, radius, fill, stroke, width))


def polyline(lines, points, stroke, width=2.0, fill="none", opacity=None):
    if not points:
        return
    attrs = ['points="%s"' % " ".join("%.2f,%.2f" % point for point in points), 'fill="%s"' % fill, 'stroke="%s"' % stroke, 'stroke-width="%.2f"' % width, 'stroke-linejoin="round"', 'stroke-linecap="round"']
    if opacity is not None:
        attrs.append('opacity="%.3f"' % opacity)
    lines.append("<polyline %s/>" % " ".join(attrs))


def polygon(lines, points, fill, opacity=1.0, stroke="none", width=0.0):
    lines.append('<polygon points="%s" fill="%s" fill-opacity="%.3f" stroke="%s" stroke-width="%.2f"/>' % (
        " ".join("%.2f,%.2f" % point for point in points), fill, opacity, stroke, width,
    ))


def panel_heading(lines, label, title_value, x, y):
    text(lines, x, y, label, "panel")
    text(lines, x + 70, y, title_value, "title")


def chromosome_key(seqid):
    digits = "".join(character for character in seqid if character.isdigit())
    return int(digits or 999)


def smooth_density(values, xmax, bins, bandwidth):
    counts = [0.0] * bins
    bin_width = xmax / bins
    for value in values:
        if 0 < value < xmax:
            index = min(bins - 1, int(value / bin_width))
            counts[index] += 1.0
    sigma = max(0.5, bandwidth / bin_width)
    radius = max(2, int(math.ceil(4.0 * sigma)))
    kernel = [math.exp(-0.5 * (offset / sigma) ** 2) for offset in range(-radius, radius + 1)]
    kernel_total = sum(kernel)
    kernel = [value / kernel_total for value in kernel]
    density = []
    for index in range(bins):
        value = 0.0
        for offset, weight in zip(range(-radius, radius + 1), kernel):
            source = index + offset
            if 0 <= source < bins:
                value += counts[source] * weight
        density.append(value)
    maximum = max(density) or 1.0
    return [((index + 0.5) * bin_width, value / maximum) for index, value in enumerate(density)]


def draw_panel_a(lines, rows, ks_max, bandwidth):
    x0, y0, w, h = 95, 105, 920, 535
    panel_heading(lines, "a", "Syntenic A-B divergence", 35, 70)
    comparisons = [
        ("pzi_hapA1_vs_pzi_hapB1", "PZI A-B", "PZI"),
        ("psp_hapA1_vs_psp_hapB1", "PSP A-B", "PSP"),
        ("ppa_hapA_vs_ppa_hapB", "PPA A-B", "PPA"),
    ]
    grouped = defaultdict(list)
    source = []
    for row in rows:
        if row.get("comparison") not in {item[0] for item in comparisons} or row.get("kaks_status") != "ok":
            continue
        try:
            value = float(row["Ks"])
        except (KeyError, ValueError):
            continue
        if 0 < value < ks_max:
            grouped[row["comparison"]].append(value)
            source.append(row)
    plot_left, plot_right = x0 + 60, x0 + w - 15
    axis_y = y0 + h - 45
    baselines = [y0 + 135, y0 + 275, y0 + 415]
    for index, (comparison, label, color_key) in enumerate(comparisons):
        values = grouped[comparison]
        if not values:
            raise ValueError("No Ks values for " + comparison)
        density = smooth_density(values, ks_max, 430, bandwidth)
        base = baselines[index]
        points = []
        maximum_point = max(density, key=lambda item: item[1])
        for x_value, y_value in density:
            x = plot_left + (plot_right - plot_left) * x_value / ks_max
            y = base - 82 * y_value
            points.append((x, y))
        polygon(lines, [(plot_left, base)] + points + [(plot_right, base)], COLORS[color_key], opacity=0.25)
        polyline(lines, points, COLORS[color_key], width=3.0)
        line(lines, plot_left, base, plot_right, base, stroke="#A0A0A0", width=1.0)
        mode_x = plot_left + (plot_right - plot_left) * maximum_point[0] / ks_max
        line(lines, mode_x, base - 90, mode_x, base + 4, stroke=COLORS[color_key], width=1.5, dash="6,5")
        text(lines, plot_right - 8, base - 12, "%s  n=%d; mode=%.2f" % (label, len(values), maximum_point[0]), "small", anchor="end", fill=COLORS[color_key])
    line(lines, plot_left, axis_y, plot_right, axis_y, width=2.0)
    for tick in range(7):
        value = tick * 0.5
        x = plot_left + (plot_right - plot_left) * value / ks_max
        line(lines, x, axis_y, x, axis_y + 9, width=1.4)
        text(lines, x, axis_y + 31, "%.1f" % value, "tick", anchor="middle")
    text(lines, (plot_left + plot_right) / 2, axis_y + 68, "Synonymous substitutions per synonymous site (Ks)", "axis", anchor="middle")
    text(lines, plot_left, y0 + 15, "Common bandwidth = %.2f; dashed lines mark density modes" % bandwidth, "subtitle")
    return source


def select_support(rows):
    selected = {}
    for row in rows:
        taxa = set(filter(None, row.get("displayed_descendant_taxa", "").split(",")))
        if taxa == {"psp_A", "pzi_A"}:
            selected["A"] = row
        elif taxa == {"psp_B", "pzi_B"}:
            selected["B"] = row
    if set(selected) != {"A", "B"}:
        raise ValueError("A and B support nodes were not both found")
    return selected


def draw_panel_b(lines, topology_rows, support_rows):
    panel_heading(lines, "b", "Concordant subgenome histories", 1055, 70)
    support = select_support(support_rows)
    sx0, sy0, sw, sh = 1095, 145, 475, 390
    text(lines, sx0, sy0 - 30, "ASTRAL and concordance support", "subtitle", weight="700")
    plot_left, plot_right = sx0 + 85, sx0 + sw - 20
    for tick in [0, 25, 50, 75, 100]:
        x = plot_left + (plot_right - plot_left) * tick / 100.0
        line(lines, x, sy0 + 20, x, sy0 + 250, stroke="#E2E2E2", width=1.0)
        text(lines, x, sy0 + 278, str(tick), "tick", anchor="middle")
    metrics = [("PP", "astral_pp"), ("gCF", "gcf"), ("sCF", "likelihood_scf")]
    for row_index, subgenome in enumerate(["A", "B"]):
        base_y = sy0 + 85 + row_index * 115
        text(lines, sx0 + 20, base_y + 7, subgenome + " node", "axis", weight="700")
        for metric_index, (metric, field) in enumerate(metrics):
            value = float(support[subgenome][field])
            if field == "astral_pp":
                value *= 100.0
            x = plot_left + (plot_right - plot_left) * value / 100.0
            y = base_y + (metric_index - 1) * 25
            circle(lines, x, y, 8, COLORS[metric], stroke="white", width=1.2)
            text(lines, x + 13, y + 6, "%.1f" % value, "tiny")
    legend_x = sx0 + 95
    for index, (metric, _) in enumerate(metrics):
        circle(lines, legend_x + index * 105, sy0 + 335, 7, COLORS[metric], stroke="white")
        text(lines, legend_x + 14 + index * 105, sy0 + 341, metric, "small")
    text(lines, (plot_left + plot_right) / 2, sy0 + 310, "Support (%)", "axis", anchor="middle")

    selected = [row for row in topology_rows if row.get("contrast") == "pzi_psp_AB"]
    if len(selected) != 3:
        raise ValueError("Expected three pzi_psp_AB topology rows")
    tx0, ty0, tw, th = 1640, 145, 680, 390
    text(lines, tx0, ty0 - 30, "PZI/PSP four-taxon topology frequencies", "subtitle", weight="700")
    desired = [
        "pzi_A+psp_A|pzi_B+psp_B",
        "pzi_A+pzi_B|psp_A+psp_B",
        "pzi_A+psp_B|psp_A+pzi_B",
    ]
    lookup = {row["topology"]: row for row in selected}
    bar_bottom = ty0 + 275
    max_height = 240
    bar_width = 145
    gap = 70
    labels = ["A-A | B-B", "Within species", "Crossed"]
    bar_colors = [COLORS["PSP"], "#B9B9B9", "#B89BD9"]
    for index, topology in enumerate(desired):
        row = lookup[topology]
        frequency = 100.0 * float(row["frequency"])
        x = tx0 + 55 + index * (bar_width + gap)
        height = max_height * frequency / 100.0
        rect(lines, x, bar_bottom - height, bar_width, height, fill=bar_colors[index], stroke="#4A4A4A", width=1.3)
        text(lines, x + bar_width / 2, bar_bottom - height - 15, "%.2f%%" % frequency, "axis", anchor="middle", weight="700")
        text(lines, x + bar_width / 2, bar_bottom + 28, labels[index], "small", anchor="middle")
        text(lines, x + bar_width / 2, bar_bottom + 50, "%s/%s trees" % (row["gene_count"], row["decisive_genes"]), "tiny", anchor="middle")
    line(lines, tx0 + 30, bar_bottom, tx0 + tw - 20, bar_bottom, width=1.8)
    text(lines, tx0 + tw / 2, ty0 + 372, "Gene-tree bootstrap 95% CIs are retained in Source Data", "tiny", anchor="middle")
    return selected, support


def moving_average(values, window):
    radius = window // 2
    prefix = [0.0]
    for value in values:
        prefix.append(prefix[-1] + value)
    result = []
    for index in range(len(values)):
        left = max(0, index - radius)
        right = min(len(values), index + radius + 1)
        result.append((prefix[right] - prefix[left]) / float(right - left))
    return result


def build_retention_profiles(reference_rows, state_rows, window):
    main_rows = [row for row in state_rows if row.get("dataset") == "main" and row.get("species") in {"PZI", "PSP", "PPA"}]
    state_lookup = {(row["species"], row["reference_gene"]): row for row in main_rows}
    chromosomes = sorted({row["reference_seqid"] for row in reference_rows}, key=chromosome_key)
    loci_by_chr = defaultdict(list)
    for row in reference_rows:
        loci_by_chr[row["reference_seqid"]].append(row)
    for seqid in chromosomes:
        loci_by_chr[seqid].sort(key=lambda row: (int(row["reference_start0"]), row["reference_gene"]))
    gap = 55
    offsets = {}
    cursor = 0
    for seqid in chromosomes:
        offsets[seqid] = cursor
        cursor += len(loci_by_chr[seqid]) + gap
    total = max(1, cursor - gap)
    profiles = defaultdict(lambda: defaultdict(dict))
    source_rows = []
    for species in ["PZI", "PSP", "PPA"]:
        for seqid in chromosomes:
            loci = loci_by_chr[seqid]
            a_values = []
            b_values = []
            for locus in loci:
                row = state_lookup[(species, locus["reference_gene"])]
                a_values.append(1 if int(row["a_anchor_count"]) == 1 else 0)
                b_values.append(1 if int(row["b_anchor_count"]) == 1 else 0)
            a_smoothed = moving_average(a_values, window)
            b_smoothed = moving_average(b_values, window)
            profiles[species][seqid]["A"] = a_smoothed
            profiles[species][seqid]["B"] = b_smoothed
            for index, locus in enumerate(loci):
                source_rows.append({
                    "species": species,
                    "reference_gene": locus["reference_gene"],
                    "reference_seqid": seqid,
                    "reference_start0": locus["reference_start0"],
                    "reference_gene_order": offsets[seqid] + index + 1,
                    "window_genes": window,
                    "A_unique_anchor_fraction": "%.8f" % a_smoothed[index],
                    "B_unique_anchor_fraction": "%.8f" % b_smoothed[index],
                })
    return chromosomes, loci_by_chr, offsets, total, profiles, source_rows


def draw_panel_c(lines, reference_rows, state_rows, window):
    panel_heading(lines, "c", "GTR-anchored A/B retention landscapes", 35, 735)
    text(lines, 105, 780, "13,849 GTR loci retained as unique 1:1 anchors in both PCAM and PPE; %d-gene rolling windows" % window, "subtitle")
    chromosomes, loci_by_chr, offsets, total, profiles, source = build_retention_profiles(reference_rows, state_rows, window)
    x_left, x_right = 125, 1370
    top = 835
    track_height = 190
    track_gap = 35
    y_min, y_max = 0.45, 1.0
    for chr_index, seqid in enumerate(chromosomes):
        start_x = x_left + (x_right - x_left) * offsets[seqid] / total
        end_x = x_left + (x_right - x_left) * (offsets[seqid] + len(loci_by_chr[seqid])) / total
        if chr_index % 2 == 0:
            rect(lines, start_x, top - 15, end_x - start_x, 3 * track_height + 2 * track_gap + 30, fill="#F6F7F8")
        midpoint = (start_x + end_x) / 2
        text(lines, midpoint, top + 3 * track_height + 2 * track_gap + 55, seqid.replace("Chr", ""), "tick", anchor="middle")
        if chr_index:
            line(lines, start_x, top - 15, start_x, top + 3 * track_height + 2 * track_gap + 15, stroke="#D8D8D8", width=1.0)
    for species_index, species in enumerate(["PZI", "PSP", "PPA"]):
        y_top = top + species_index * (track_height + track_gap)
        y_bottom = y_top + track_height
        for tick in [0.5, 0.75, 1.0]:
            y = y_bottom - track_height * (tick - y_min) / (y_max - y_min)
            line(lines, x_left, y, x_right, y, stroke="#DDDDDD", width=1.0)
            text(lines, x_left - 12, y + 5, str(int(tick * 100)), "tick", anchor="end")
        text(lines, x_left - 62, y_top + track_height / 2 + 7, species, "axis", anchor="middle", weight="700")
        line(lines, x_left, y_bottom, x_right, y_bottom, width=1.5)
        for seqid in chromosomes:
            length = len(loci_by_chr[seqid])
            sample_step = max(1, length // 180)
            for subgenome in ["A", "B"]:
                points = []
                values = profiles[species][seqid][subgenome]
                indexes = list(range(0, length, sample_step))
                if indexes[-1] != length - 1:
                    indexes.append(length - 1)
                for index in indexes:
                    x = x_left + (x_right - x_left) * (offsets[seqid] + index) / total
                    clipped = min(y_max, max(y_min, values[index]))
                    y = y_bottom - track_height * (clipped - y_min) / (y_max - y_min)
                    points.append((x, y))
                polyline(lines, points, COLORS[subgenome], width=2.6)
    text(lines, (x_left + x_right) / 2, top + 3 * track_height + 2 * track_gap + 92, "GTR chromosome gene order", "axis", anchor="middle")
    text(lines, 58, top + 305, "Unique-anchor retention (%)", "small", anchor="middle")
    line(lines, 1110, 798, 1150, 798, stroke=COLORS["A"], width=3.0)
    text(lines, 1160, 804, "A", "small")
    line(lines, 1215, 798, 1255, 798, stroke=COLORS["B"], width=3.0)
    text(lines, 1265, 804, "B", "small")
    return source


def read_bed_gene_order(path, seqid):
    genes = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line_value in handle:
            fields = line_value.rstrip("\n").split("\t")
            if len(fields) < 4 or fields[0] != seqid:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            genes.append((start, end, fields[3]))
    genes.sort(key=lambda item: (item[0], item[1], item[2]))
    return genes


def read_primary_seqids(inputs, label):
    path = Path(inputs) / label / (label + ".seqids")
    seqids = [value for value in path.read_text(encoding="utf-8").strip().split(",") if value]
    if label == "ppe":
        return [value for value in seqids if value.startswith("NC_0340")]
    return seqids[:8]


def dominant_corresponding_seqid(ref_genes, label, pair_index, allowed_seqids):
    counter = Counter()
    allowed = set(allowed_seqids)
    for ref_gene in ref_genes:
        for seqid in pair_index.get((label, ref_gene), set()):
            if seqid in allowed:
                counter[seqid] += 1
    if not counter:
        return "", 0, 0
    seqid, count = counter.most_common(1)[0]
    return seqid, count, sum(counter.values())


def draw_panel_c_selected(lines, pair_rows, inputs, window):
    panel_heading(lines, "c", "Cross-species retention along selected PZI chromosomes", 35, 735)
    text(lines, 105, 780, "PZI copy1 targets are relabelled with haplotype chromosome IDs; all GTR chromosomes are eligible bridge anchors", "subtitle")

    pair_index = defaultdict(set)
    pzi_rows_by_label_seqid = defaultdict(list)
    for row in pair_rows:
        label = row["query_label"]
        pair_index[(label, row["ref_gene"])].add(row["query_seqid"])
        if label in {"pzi_hapA1", "pzi_hapB1"}:
            pzi_rows_by_label_seqid[(label, row["query_seqid"])].append(row)

    plot_left, plot_right = 115, 1390
    top = 845
    column_gap = 55
    row_gap = 58
    subplot_w = (plot_right - plot_left - column_gap) / 2.0
    subplot_h = 185
    legend_x = 730
    for index, species in enumerate(["PZI", "PSP", "PCAM", "PPA", "PPE"]):
        x = legend_x + index * 135
        line(lines, x, 802, x + 34, 802, stroke=CROSS_SPECIES_COLORS[species], width=3.2)
        text(lines, x + 43, 808, species, "small")

    source_rows = []
    mapping_rows = []
    for panel_index, target in enumerate(SELECTED_PZI_CHROMOSOMES):
        column = panel_index % 2
        row_index = panel_index // 2
        x0 = plot_left + column * (subplot_w + column_gap)
        y0 = top + row_index * (subplot_h + row_gap)
        bed_path = Path(inputs) / target["pzi_label"] / (target["pzi_label"] + ".bed")
        pzi_genes = read_bed_gene_order(bed_path, target["hap_seqid"])
        if not pzi_genes:
            raise ValueError("No PZI genes for " + target["hap_seqid"])

        pzi_pair_rows = pzi_rows_by_label_seqid[(target["pzi_label"], target["hap_seqid"])]
        gene_to_refs = defaultdict(set)
        gtr_chr_counts = Counter()
        for row in pzi_pair_rows:
            gene_to_refs[row["query_gene"]].add(row["ref_gene"])
            gtr_chr_counts[row["ref_seqid"]] += 1
        bridge_ref_genes = set()
        for values in gene_to_refs.values():
            bridge_ref_genes.update(values)

        target_labels = {
            "PSP": target["psp_label"],
            "PCAM": "pcam",
            "PPA": target["ppa_label"],
            "PPE": "ppe",
        }
        corresponding = {}
        for species, label in target_labels.items():
            seqid, dominant_count, all_count = dominant_corresponding_seqid(
                bridge_ref_genes, label, pair_index, read_primary_seqids(inputs, label)
            )
            corresponding[species] = (label, seqid)
            mapping_rows.append({
                "copy1_chromosome": target["copy1"],
                "pzi_haplotype_chromosome": target["hap_seqid"],
                "comparison_species": species,
                "comparison_label": label,
                "corresponding_chromosome": seqid,
                "dominant_shared_GTR_anchors": dominant_count,
                "all_primary_chromosome_shared_GTR_anchors": all_count,
                "GTR_chromosomes_used": ",".join("%s:%d" % item for item in sorted(gtr_chr_counts.items(), key=lambda item: chromosome_key(item[0]))),
            })

        raw = {species: [] for species in ["PZI", "PSP", "PCAM", "PPA", "PPE"]}
        per_gene_refs = []
        for _, _, gene in pzi_genes:
            refs = gene_to_refs.get(gene, set())
            per_gene_refs.append(refs)
            raw["PZI"].append(1 if refs else 0)
            for species in ["PSP", "PCAM", "PPA", "PPE"]:
                label, seqid = corresponding[species]
                present = any(seqid in pair_index.get((label, ref_gene), set()) for ref_gene in refs) if seqid else False
                raw[species].append(1 if present else 0)
        profiles = {species: moving_average(values, window) for species, values in raw.items()}

        for gene_index, (_, _, gene) in enumerate(pzi_genes):
            row = {
                "copy1_chromosome": target["copy1"],
                "pzi_haplotype_chromosome": target["hap_seqid"],
                "pzi_gene_order": gene_index + 1,
                "pzi_gene": gene,
                "GTR_bridge_genes": ",".join(sorted(per_gene_refs[gene_index])),
                "window_genes": window,
            }
            for species in ["PZI", "PSP", "PCAM", "PPA", "PPE"]:
                row[species + "_retention_fraction"] = "%.8f" % profiles[species][gene_index]
            source_rows.append(row)

        rect(lines, x0, y0, subplot_w, subplot_h, fill="#FBFBFB", stroke="#777777", width=0.8)
        for tick in [0, 50, 100]:
            y = y0 + subplot_h - subplot_h * tick / 100.0
            line(lines, x0, y, x0 + subplot_w, y, stroke="#E2E2E2", width=0.8)
            if column == 0:
                text(lines, x0 - 10, y + 5, tick, "tiny", anchor="end")
        text(lines, x0 + 8, y0 + 22, "%s (%s)" % (target["hap_seqid"], target["copy1"]), "small", weight="700")
        correspondence_text = "PSP %s | PCAM %s | PPA %s | PPE %s" % tuple(
            corresponding[species][1] or "NA" for species in ["PSP", "PCAM", "PPA", "PPE"]
        )
        text(lines, x0 + 8, y0 + 42, correspondence_text, "tiny")
        sample_step = max(1, len(pzi_genes) // 280)
        indexes = list(range(0, len(pzi_genes), sample_step))
        if indexes[-1] != len(pzi_genes) - 1:
            indexes.append(len(pzi_genes) - 1)
        for species in ["PZI", "PSP", "PCAM", "PPA", "PPE"]:
            points = []
            for gene_index in indexes:
                x = x0 + subplot_w * gene_index / max(1, len(pzi_genes) - 1)
                y = y0 + subplot_h - subplot_h * profiles[species][gene_index]
                points.append((x, y))
            polyline(lines, points, CROSS_SPECIES_COLORS[species], width=2.0, opacity=0.92)
        for tick_index, tick in enumerate([1, max(1, len(pzi_genes) // 2), len(pzi_genes)]):
            x = x0 + subplot_w * (tick - 1) / max(1, len(pzi_genes) - 1)
            text(lines, x, y0 + subplot_h + 20, tick, "tiny", anchor="middle")
    text(lines, 55, top + 350, "Retention across PZI gene-order windows (%)", "small", anchor="middle")
    text(lines, (plot_left + plot_right) / 2, 1645, "PZI target chromosome gene order", "axis", anchor="middle")
    return source_rows, mapping_rows


def candidate_chromosome_summary(reference_rows, candidate_rows):
    main = [row for row in candidate_rows if row.get("dataset") == "main"]
    chromosomes = sorted({row["reference_seqid"] for row in reference_rows}, key=chromosome_key)
    loci_counts = Counter(row["reference_seqid"] for row in reference_rows)
    counts = defaultdict(Counter)
    for row in main:
        seqid = row["reference_seqid"]
        state = row["joint_anchor_state"]
        if state == "same_A_only_candidate":
            counts[seqid]["shared_A_only"] += 1
            if row.get("ppa_anchor_state") == "A_only_anchored":
                counts[seqid]["ppa_same_A"] += 1
        elif state == "same_B_only_candidate":
            counts[seqid]["shared_B_only"] += 1
            if row.get("ppa_anchor_state") == "B_only_anchored":
                counts[seqid]["ppa_same_B"] += 1
    rows = []
    for seqid in chromosomes:
        a_count = counts[seqid]["shared_A_only"]
        b_count = counts[seqid]["shared_B_only"]
        rows.append({
            "reference_seqid": seqid,
            "reference_loci": loci_counts[seqid],
            "shared_A_only_candidates": a_count,
            "PPA_same_A_only": counts[seqid]["ppa_same_A"],
            "shared_B_only_candidates": b_count,
            "PPA_same_B_only": counts[seqid]["ppa_same_B"],
        })
    rows.append({
        "reference_seqid": "Total",
        "reference_loci": sum(loci_counts.values()),
        "shared_A_only_candidates": sum(row["shared_A_only_candidates"] for row in rows),
        "PPA_same_A_only": sum(row["PPA_same_A_only"] for row in rows),
        "shared_B_only_candidates": sum(row["shared_B_only_candidates"] for row in rows),
        "PPA_same_B_only": sum(row["PPA_same_B_only"] for row in rows),
    })
    return rows


def draw_panel_d(lines, reference_rows, joint_rows, candidate_rows):
    panel_heading(lines, "d", "Joint PZI/PSP anchor states", 1450, 735)
    selected = [row for row in joint_rows if row.get("dataset") == "main"]
    if not selected:
        raise ValueError("No main joint states")
    counts = Counter(row["joint_anchor_state"] for row in selected)
    total = len(selected)
    bar_x, bar_y, bar_w, bar_h = 1515, 825, 790, 66
    cursor = bar_x
    for state in JOINT_ORDER:
        fraction = counts.get(state, 0) / float(total)
        segment_width = bar_w * fraction
        rect(lines, cursor, bar_y, segment_width, bar_h, fill=COLORS[state], stroke="white", width=1.0)
        if segment_width > 75:
            text(lines, cursor + segment_width / 2, bar_y + 42, "%.1f%%" % (100 * fraction), "small", anchor="middle", weight="700", fill="#FFFFFF" if state != "discordant_or_uninformative" else "#333333")
        cursor += segment_width
    text(lines, bar_x, bar_y - 18, "Main haplotypes: n=%d conservative GTR loci" % total, "subtitle")
    for index, state in enumerate(JOINT_ORDER):
        column = index % 2
        row_index = index // 2
        lx = bar_x + column * 390
        ly = 930 + row_index * 35
        rect(lines, lx, ly - 14, 20, 20, fill=COLORS[state], stroke="#666666", width=0.6)
        text(lines, lx + 30, ly + 2, "%s: %d" % (JOINT_LABELS[state], counts.get(state, 0)), "small")

    summary_rows = candidate_chromosome_summary(reference_rows, candidate_rows)
    table_x, table_y = 1500, 1065
    widths = [105, 115, 145, 115, 145]
    headers = ["GTR chr", "A-only", "PPA same A", "B-only", "PPA same B"]
    table_width = sum(widths)
    row_height = 43
    rect(lines, table_x, table_y, table_width, row_height, fill="#ECEFF2")
    x = table_x
    for header, width_value in zip(headers, widths):
        text(lines, x + width_value / 2, table_y + 27, header, "tablehead", anchor="middle")
        x += width_value
    for row_index, row in enumerate(summary_rows):
        y = table_y + row_height * (row_index + 1)
        fill = "#F1F1F1" if row["reference_seqid"] == "Total" else ("#FAFAFA" if row_index % 2 == 0 else "#FFFFFF")
        rect(lines, table_x, y, table_width, row_height, fill=fill)
        values = [
            row["reference_seqid"].replace("Chr", ""),
            row["shared_A_only_candidates"],
            row["PPA_same_A_only"],
            row["shared_B_only_candidates"],
            row["PPA_same_B_only"],
        ]
        x = table_x
        for value, width_value in zip(values, widths):
            text(lines, x + width_value / 2, y + 27, value, "tablebody", anchor="middle", weight="700" if row["reference_seqid"] == "Total" else None)
            x += width_value
    rect(lines, table_x, table_y, table_width, row_height * (len(summary_rows) + 1), fill="none", stroke="#777777", width=1.0)
    text(lines, table_x, 1610, "Candidate shared fractionation states; absence of an anchor is not a verified gene loss.", "tiny")
    text(lines, table_x, 1640, "PPA columns provide an additional same-state comparison, not parental evidence.", "tiny")
    return selected, summary_rows


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: " + str(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    lines = svg_start()
    ks_source = draw_panel_a(lines, read_tsv(args.ks), args.ks_max, args.ks_bandwidth)
    topology_source, support_source = draw_panel_b(lines, read_tsv(args.topologies), read_tsv(args.support))
    reference_rows = read_tsv(args.reference_loci)
    retention_source, correspondence_source = draw_panel_c_selected(lines, read_tsv(args.pairs), args.jcvi_inputs, args.window_genes)
    joint_source, candidate_summary = draw_panel_d(lines, reference_rows, read_tsv(args.joint_states), read_tsv(args.candidates))
    lines.append("</svg>")
    (outdir / "figure2_alltaxa_abcd.svg").write_text("\n".join(lines) + "\n", encoding="utf-8")

    write_tsv(outdir / "panel_a_ks_source_data.tsv", list(ks_source[0].keys()), ks_source)
    write_tsv(outdir / "panel_b_topology_source_data.tsv", list(topology_source[0].keys()), topology_source)
    support_rows = []
    for subgenome in ["A", "B"]:
        row = dict(support_source[subgenome])
        row["subgenome_node"] = subgenome
        support_rows.append(row)
    write_tsv(outdir / "panel_b_support_source_data.tsv", list(support_rows[0].keys()), support_rows)
    write_tsv(outdir / "panel_c_retention_profiles.tsv", list(retention_source[0].keys()), retention_source)
    write_tsv(outdir / "panel_c_corresponding_chromosomes.tsv", list(correspondence_source[0].keys()), correspondence_source)
    write_tsv(outdir / "panel_d_joint_source_data.tsv", list(joint_source[0].keys()), joint_source)
    write_tsv(outdir / "panel_d_candidate_by_chromosome.tsv", list(candidate_summary[0].keys()), candidate_summary)
    manifest = [
        {"panel": "a", "question": "A-B syntenic Ks distributions", "boundary": "Ks gives relative divergence, not proof of one shared polyploidization."},
        {"panel": "b", "question": "PZI/PSP A-A and B-B phylogenetic concordance", "boundary": "Does not identify parental lineages."},
        {"panel": "c", "question": "Five-species retention along six selected PZI chromosomes using all GTR chromosomes as bridge anchors", "boundary": "Corresponding target chromosomes are inferred from dominant shared GTR syntenic anchors; absence is not verified loss."},
        {"panel": "d", "question": "PZI/PSP joint states and candidate shared fractionation", "boundary": "Candidates require local synteny, sequence search, gap and read-depth validation."},
    ]
    write_tsv(outdir / "figure_manifest.tsv", ["panel", "question", "boundary"], manifest)
    (outdir / "README.txt").write_text(
        "Editable master: figure2_alltaxa_abcd.svg\n"
        "Panel c uses all GTR chromosomes as bridge anchors for six selected PZI copy1 chromosomes.\n"
        "The plot labels use haplotype chromosome IDs rather than copy1 numbering.\n"
        "Panel d categories are candidate anchor states, not validated gene losses.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("Figure 2 all-taxon a-d SVG completed.\n", encoding="ascii")
    print("[OK] Wrote Figure 2 all-taxon SVG bundle: " + str(outdir))


if __name__ == "__main__":
    main()
