#!/usr/bin/env python3
"""Create 100-GTR-gene retention windows with all PZI A/B chromosomes.

The calculation reuses the formal quota-filtered anchors.  It keeps every
PZI-A and PZI-B chromosome in every GTR panel, including chromosomes with no
anchor (reported as zero retention), while retaining one maximum-anchor
chromosome for PSP, PPA, PCAM and PPE.
"""

from __future__ import print_function

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")
ALL_PZI = ("PZI_A", "PZI_B")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def chromosome_key(seqid):
    match = re.search(r"C(\d+)", seqid)
    if match:
        return (int(match.group(1)), seqid)
    match = re.search(r"(\d+)$", seqid)
    return (int(match.group(1)) if match else 10 ** 9, seqid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", required=True)
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window-size", type=int, default=100)
    args = parser.parse_args()
    if args.window_size != 100:
        raise SystemExit("ERROR: the reference comparison uses a 100-gene window")
    if any(token in str(args.anchors).lower() for token in ("jcvi", "dollo", "candidate")):
        raise SystemExit("ERROR: only formal quota-filtered anchors are valid input")

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    anchors = read_tsv(args.anchors)
    catalog = Path(args.catalog)
    gtr_order = read_tsv(catalog / "gtr_gene_order.tsv")
    query_map = read_tsv(catalog / "query_gene_track_map.tsv")
    if not anchors or not gtr_order or not query_map:
        raise SystemExit("ERROR: empty formal anchor, GTR-order or query-track input")

    genes_by_target = defaultdict(list)
    for row in gtr_order:
        genes_by_target[row["gtr_seqid"]].append(row)

    anchor_presence = defaultdict(set)
    anchor_counts = Counter()
    chromosomes_by_track = defaultdict(set)
    for row in anchors:
        track = row["track_id"]
        if track not in TRACKS:
            continue
        key = (row["gtr_seqid"], row["gtr_gene"], row["query_seqid"])
        anchor_presence[key].add(row["query_gene"])
        anchor_counts[(track, row["query_seqid"], row["gtr_seqid"])] += 1
        chromosomes_by_track[track].add(row["query_seqid"])

    # The catalog is authoritative for the complete PZI chromosome set.
    catalog_chromosomes = defaultdict(set)
    for row in query_map:
        if row["track_id"] in TRACKS:
            catalog_chromosomes[row["track_id"]].add(row["query_seqid"])

    selected_rows = []
    window_rows = []
    gene_rows = []
    for target_seqid in sorted(genes_by_target, key=chromosome_key):
        target_genes = genes_by_target[target_seqid]
        if len(target_genes) < args.window_size:
            raise SystemExit("ERROR: {} has fewer than {} genes".format(target_seqid, args.window_size))
        selected = {}
        for track in TRACKS:
            if track in ALL_PZI:
                query_ids = sorted(catalog_chromosomes[track], key=chromosome_key)
            else:
                ranked = sorted(
                    ((anchor_counts[(track, query_seqid, target_seqid)], query_seqid)
                     for query_seqid in chromosomes_by_track[track]),
                    key=lambda item: (-item[0], item[1]),
                )
                query_ids = [ranked[0][1]] if ranked and ranked[0][0] else []
            selected[track] = query_ids
            for query_seqid in query_ids:
                selected_rows.append({
                    "target_seqid": target_seqid,
                    "track_id": track,
                    "query_seqid": query_seqid,
                    "n_quota_filtered_anchors": anchor_counts[(track, query_seqid, target_seqid)],
                    "selection_rule": "all_eight_PZI_chromosomes" if track in ALL_PZI else "maximum_quota_filtered_anchors",
                })

        for track in TRACKS:
            for query_seqid in selected[track]:
                presence = []
                for gene_order, gene in enumerate(target_genes, 1):
                    observed = sorted(anchor_presence[(target_seqid, gene["gtr_gene"], query_seqid)])
                    present = 1 if observed else 0
                    presence.append(present)
                    gene_rows.append({
                        "target_seqid": target_seqid,
                        "target_gene_order": gene_order,
                        "gtr_gene": gene["gtr_gene"],
                        "track_id": track,
                        "query_seqid": query_seqid,
                        "anchored": present,
                        "status": "no_quota_filtered_anchor" if not anchor_counts[(track, query_seqid, target_seqid)] else "observed_or_partial",
                        "query_genes": ",".join(observed),
                    })
                running = sum(presence[:args.window_size])
                for start in range(0, len(presence) - args.window_size + 1):
                    if start:
                        running += presence[start + args.window_size - 1] - presence[start - 1]
                    retention = int(round(100.0 * running / args.window_size))
                    window_rows.append({
                        "target_seqid": target_seqid,
                        "window_index": start + 1,
                        "window_start_gene": start + 1,
                        "window_end_gene": start + args.window_size,
                        "window_center_gene": start + (args.window_size + 1) / 2.0,
                        "track_id": track,
                        "query_seqid": query_seqid,
                        "syntenic_genes_in_window": running,
                        "window_size": args.window_size,
                        "retention_percent": retention,
                        "retention_proportion": "{:.6f}".format(retention / 100.0),
                        "status": "no_quota_filtered_anchor" if not anchor_counts[(track, query_seqid, target_seqid)] else "observed",
                        "method": "FraBias_local_IO_port_all_gtr_all_PZI",
                    })

    write_tsv(outdir / "all_gtr_all_pzi_selected_tracks.tsv", list(selected_rows[0]), selected_rows)
    write_tsv(outdir / "all_gtr_all_pzi_fractbias_genes.tsv", list(gene_rows[0]), gene_rows)
    write_tsv(outdir / "all_gtr_all_pzi_fractbias_windows.tsv", list(window_rows[0]), window_rows)
    (outdir / "method_status.txt").write_text(
        "Formal quota-filtered anchors were projected to every GTR chromosome.\n"
        "PZI_A and PZI_B retain all eight catalog chromosomes in every GTR panel, including zero-anchor chromosomes.\n"
        "PSP/PPA/PCAM/PPE retain the maximum-anchor query chromosome per GTR target.\n"
        "FraBias denominator is every ordered GTR gene in an overlapping 100-gene window.\n"
        "Anchor absence is not a validated biological gene-loss call.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("All-GTR windows with all PZI A/B chromosomes completed.\n", encoding="ascii")
    print("[OK] Wrote {} windows for {} GTR chromosomes (all PZI A/B chromosomes retained): {}".format(
        len(window_rows), len(genes_by_target), outdir))


if __name__ == "__main__":
    main()
