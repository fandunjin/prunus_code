#!/usr/bin/env python3
"""Cluster-bootstrap filtered Ks medians and histogram modes by comparison."""

from __future__ import annotations

import argparse
import csv
import math
import random
from collections import defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--out-summary", required=True)
    parser.add_argument("--out-contrasts", required=True)
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20260814)
    parser.add_argument("--min-ks", type=float, default=0.0)
    parser.add_argument("--max-ks", type=float, default=3.0)
    parser.add_argument("--bin-width", type=float, default=0.02)
    parser.add_argument(
        "--cluster-field",
        default="block_id",
        help="Resampling unit. Use an empty value only for gene-pair bootstrap.",
    )
    return parser.parse_args()


def quantile(values, probability):
    values = sorted(values)
    if not values:
        return float("nan")
    position = (len(values) - 1) * probability
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return values[lower]
    fraction = position - lower
    return values[lower] * (1.0 - fraction) + values[upper] * fraction


def median(values):
    return quantile(values, 0.5)


def histogram_mode(values, min_ks, max_ks, bin_width):
    n_bins = int(math.ceil((max_ks - min_ks) / bin_width))
    counts = [0] * n_bins
    for value in values:
        index = int((value - min_ks) / bin_width)
        if 0 <= index < n_bins:
            counts[index] += 1
    maximum = max(counts) if counts else 0
    if maximum == 0:
        return float("nan")
    winners = [i for i, count in enumerate(counts) if count == maximum]
    index = winners[len(winners) // 2]
    return min_ks + (index + 0.5) * bin_width


def read_values(path, min_ks, max_ks, cluster_field):
    grouped = defaultdict(lambda: defaultdict(list))
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        for row_index, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=1):
            try:
                value = float(row.get("Ks", ""))
            except ValueError:
                continue
            if math.isfinite(value) and min_ks < value < max_ks:
                cluster = row.get(cluster_field, "") if cluster_field else ""
                if not cluster:
                    cluster = "pair_{:09d}".format(row_index)
                grouped[row["comparison"]][cluster].append(value)
    return grouped


def bootstrap_cluster_sample(clusters, rng):
    cluster_values = list(clusters.values())
    sampled = []
    for _ in cluster_values:
        sampled.extend(cluster_values[rng.randrange(len(cluster_values))])
    return sampled


def write_tsv(path, rows, fields):
    out = Path(path)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    rng = random.Random(args.seed)
    grouped = read_values(args.input, args.min_ks, args.max_ks, args.cluster_field)
    boot = {}
    summary_rows = []
    for comparison, clusters in sorted(grouped.items()):
        values = [value for cluster_values in clusters.values() for value in cluster_values]
        medians = []
        modes = []
        for _ in range(args.replicates):
            sample = bootstrap_cluster_sample(clusters, rng)
            medians.append(median(sample))
            modes.append(histogram_mode(sample, args.min_ks, args.max_ks, args.bin_width))
        boot[comparison] = medians
        summary_rows.append({
            "comparison": comparison,
            "n_valid_ks": len(values),
            "n_bootstrap_clusters": len(clusters),
            "ks_median": "{:.8f}".format(median(values)),
            "ks_median_ci025": "{:.8f}".format(quantile(medians, 0.025)),
            "ks_median_ci975": "{:.8f}".format(quantile(medians, 0.975)),
            "ks_histogram_mode": "{:.8f}".format(histogram_mode(values, args.min_ks, args.max_ks, args.bin_width)),
            "ks_mode_ci025": "{:.8f}".format(quantile(modes, 0.025)),
            "ks_mode_ci975": "{:.8f}".format(quantile(modes, 0.975)),
            "bootstrap_replicates": args.replicates,
            "bootstrap_unit": args.cluster_field or "gene_pair",
            "ks_filter": "{}<Ks<{}".format(args.min_ks, args.max_ks),
        })

    write_tsv(
        args.out_summary,
        summary_rows,
        [
            "comparison", "n_valid_ks", "n_bootstrap_clusters", "ks_median", "ks_median_ci025",
            "ks_median_ci975", "ks_histogram_mode", "ks_mode_ci025",
            "ks_mode_ci975", "bootstrap_replicates", "bootstrap_unit", "ks_filter",
        ],
    )

    focal = [
        "pzi_hapA1_vs_pzi_hapB1",
        "psp_hapA1_vs_psp_hapB1",
        "ppa_hapA_vs_ppa_hapB",
    ]
    contrast_rows = []
    for index, left in enumerate(focal):
        for right in focal[index + 1 :]:
            if left not in boot or right not in boot:
                continue
            differences = [a - b for a, b in zip(boot[left], boot[right])]
            contrast_rows.append({
                "comparison_left": left,
                "comparison_right": right,
                "median_difference_left_minus_right": "{:.8f}".format(median(differences)),
                "ci025": "{:.8f}".format(quantile(differences, 0.025)),
                "ci975": "{:.8f}".format(quantile(differences, 0.975)),
                "bootstrap_replicates": args.replicates,
            })
    write_tsv(
        args.out_contrasts,
        contrast_rows,
        [
            "comparison_left", "comparison_right", "median_difference_left_minus_right",
            "ci025", "ci975", "bootstrap_replicates",
        ],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
