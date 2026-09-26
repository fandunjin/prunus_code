#!/usr/bin/env python3
"""Calculate eight-track FraBias windows for every chromosome-scale GTR target.

This is an audit companion to the three-row formal panel c.  For every GTR
chromosome, each query track contributes the one physical chromosome with the
most quota-filtered anchors to that target.  The calculation therefore compares
the best chromosome-scale correspondence in each taxon; it does not call gene
loss or force chromosome names to match across assemblies.
"""

from __future__ import print_function

import argparse
import csv
import re
from collections import Counter, defaultdict
from pathlib import Path


TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def chromosome_key(seqid):
    match = re.search(r"(\d+)$", seqid)
    return (int(match.group(1)) if match else 10 ** 9, seqid)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--anchors", required=True, help="Formal quota_filtered_anchors.tsv")
    parser.add_argument("--catalog", required=True, help="Formal frozen Fig. 2c catalog")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window-size", type=int, default=100)
    args = parser.parse_args()
    if args.window_size != 100:
        raise SystemExit("ERROR: the reference Fig. 3c comparison uses a 100-gene window")
    if any(token in str(args.anchors).lower() for token in ("jcvi", "dollo", "candidate", "alltaxa_fractionation_states")):
        raise SystemExit("ERROR: only formal quota-filtered anchors are valid FraBias input")

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    catalog = Path(args.catalog)
    anchors = read_tsv(args.anchors)
    gtr_order = read_tsv(catalog / "gtr_gene_order.tsv")
    if not anchors or not gtr_order:
        raise SystemExit("ERROR: empty formal anchor or GTR-order input")

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

    selected_rows = []
    window_rows = []
    gene_rows = []
    for target_seqid in sorted(genes_by_target, key=chromosome_key):
        target_genes = genes_by_target[target_seqid]
        if len(target_genes) < args.window_size:
            raise SystemExit("ERROR: {} has fewer than {} genes".format(target_seqid, args.window_size))
        selected = {}
        for track in TRACKS:
            ranked = sorted(
                ((anchor_counts[(track, query_seqid, target_seqid)], query_seqid)
                 for query_seqid in chromosomes_by_track[track]),
                key=lambda item: (-item[0], item[1]),
            )
            n_anchors, query_seqid = ranked[0] if ranked else (0, "NA")
            if n_anchors == 0:
                query_seqid = "NA"
            selected[track] = query_seqid
            selected_rows.append({
                "target_seqid": target_seqid,
                "track_id": track,
                "query_seqid": query_seqid,
                "n_quota_filtered_anchors": n_anchors,
                "selection_rule": (
                    "maximum_quota_filtered_anchors_per_track"
                    if n_anchors else "no_quota_filtered_anchor"
                ),
            })

        for track in TRACKS:
            query_seqid = selected[track]
            no_anchor = query_seqid == "NA"
            presence = []
            for gene_order, gene in enumerate(target_genes, 1):
                observed_genes = [] if no_anchor else sorted(
                    anchor_presence[(target_seqid, gene["gtr_gene"], query_seqid)]
                )
                present = 1 if observed_genes else 0
                presence.append(present)
                gene_rows.append({
                    "target_seqid": target_seqid,
                    "target_gene_order": gene_order,
                    "gtr_gene": gene["gtr_gene"],
                    "track_id": track,
                    "query_seqid": query_seqid,
                    "anchored": present,
                    "status": "no_quota_filtered_anchor" if no_anchor else "observed",
                    "query_genes": ",".join(observed_genes),
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
                    "retention_percent": "NA" if no_anchor else retention,
                    "retention_proportion": "NA" if no_anchor else "{:.6f}".format(retention / 100.0),
                    "status": "no_quota_filtered_anchor" if no_anchor else "observed",
                    "method": "FraBias_local_IO_port_all_gtr_best_track",
                })

    write_tsv(outdir / "all_gtr_selected_tracks.tsv", list(selected_rows[0].keys()), selected_rows)
    write_tsv(outdir / "all_gtr_fractbias_genes.tsv", list(gene_rows[0].keys()), gene_rows)
    write_tsv(outdir / "all_gtr_fractbias_windows.tsv", list(window_rows[0].keys()), window_rows)
    (outdir / "method_status.txt").write_text(
        "Formal FraBias anchors were projected to every GTR chromosome.\n"
        "Each of eight tracks uses its chromosome with the maximum quota-filtered anchor count to that GTR target.\n"
        "FraBias denominator is every ordered GTR gene in an overlapping 100-gene window.\n"
        "Anchor absence is not a validated biological gene-loss call.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("All-GTR FraBias windows completed.\n", encoding="ascii")
    print("[OK] Wrote {} windows for {} GTR chromosomes: {}".format(
        len(window_rows), len(genes_by_target), outdir))


if __name__ == "__main__":
    main()
