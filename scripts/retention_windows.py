#!/usr/bin/env python3
"""Run the FraBias 100-gene sliding-window calculation on local inputs.

This is an I/O port of CoGe's ``fractionation_bias.py``.  The inferential
algorithm is unchanged: quota-filtered syntenic anchors are projected onto all
ordered genes of each target chromosome, presence/absence vectors are built
for every query chromosome, and retention is calculated in every overlapping
window as ``round(100 * sum(presence) / window_size)``.
"""

from __future__ import print_function

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


TRACKS = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, delimiter="\t", fieldnames=fields, lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def prefixed(track, chromosome):
    return track + "__" + chromosome


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--anchors", required=True, help="Formal quota_filtered_anchors.tsv"
    )
    parser.add_argument("--catalog", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window-size", type=int, default=100)
    args = parser.parse_args()
    if args.window_size != 100:
        raise SystemExit(
            "ERROR: the reference Figure 3 method requires a 100-gene window"
        )
    if "jcvi" in str(args.anchors).lower() or "alltaxa_fractionation_states" in str(
        args.anchors
    ):
        raise SystemExit(
            "ERROR: JCVI/candidate anchors are not valid formal FraBias input"
        )

    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)
    catalog = Path(args.catalog)
    anchors = read_tsv(args.anchors)
    gtr_order = read_tsv(catalog / "gtr_gene_order.tsv")
    pairs = read_tsv(catalog / "pzi_display_chromosome_pairs.tsv")
    if not anchors or not gtr_order or not pairs:
        raise SystemExit("ERROR: empty anchor, GTR-order, or display-pair input")

    gene_rows_by_chr = defaultdict(list)
    for row in gtr_order:
        gene_rows_by_chr[row["gtr_seqid"]].append(row)

    anchor_presence = defaultdict(set)
    anchor_counts = Counter()
    query_chromosomes_by_track = defaultdict(set)
    for row in anchors:
        track = row["track_id"]
        if track not in TRACKS:
            continue
        key = (row["gtr_seqid"], row["gtr_gene"], row["query_seqid"])
        anchor_presence[key].add(row["query_gene"])
        anchor_counts[(track, row["query_seqid"], row["gtr_seqid"])] += 1
        query_chromosomes_by_track[track].add(row["query_seqid"])

    selected_rows = []
    used_targets = set()
    for pair in pairs:
        pair_label = pair["pair_label"]
        fixed = {
            "PZI_A": prefixed("PZI_A", pair["pzi_A_chromosome"]),
            "PZI_B": prefixed("PZI_B", pair["pzi_B_chromosome"]),
        }
        candidates = []
        for target_chr in gene_rows_by_chr:
            score = sum(
                anchor_counts.get((track, query_chr, target_chr), 0)
                for track, query_chr in fixed.items()
            )
            candidates.append((score, target_chr))
        candidates.sort(key=lambda item: (-item[0], item[1]))
        available = [
            item for item in candidates if item[1] not in used_targets and item[0] > 0
        ]
        if not available:
            raise SystemExit(
                "ERROR: no unique GTR target chromosome found for {}".format(pair_label)
            )
        target_score, target_chr = available[0]
        used_targets.add(target_chr)
        for track in TRACKS:
            if track in fixed:
                query_chr = fixed[track]
            else:
                ranked = sorted(
                    (
                        (anchor_counts.get((track, qchr, target_chr), 0), qchr)
                        for qchr in query_chromosomes_by_track[track]
                    ),
                    key=lambda item: (-item[0], item[1]),
                )
                if not ranked or ranked[0][0] == 0:
                    raise SystemExit(
                        "ERROR: no quota-filtered {} chromosome anchors target {}".format(
                            track, target_chr
                        )
                    )
                query_chr = ranked[0][1]
            n_anchors = anchor_counts.get((track, query_chr, target_chr), 0)
            if n_anchors == 0:
                raise SystemExit(
                    "ERROR: selected chromosome {} has no {} anchors to {}".format(
                        query_chr, track, target_chr
                    )
                )
            selected_rows.append(
                {
                    "pair_label": pair_label,
                    "target_seqid": target_chr,
                    "track_id": track,
                    "query_seqid": query_chr,
                    "n_quota_filtered_anchors": n_anchors,
                    "pzi_pair_anchor_score": target_score,
                }
            )

    window_rows = []
    gene_presence_rows = []
    for selected in selected_rows:
        pair_label = selected["pair_label"]
        target_chr = selected["target_seqid"]
        track = selected["track_id"]
        query_chr = selected["query_seqid"]
        target_genes = gene_rows_by_chr[target_chr]
        presence = []
        for local_index, gene_row in enumerate(target_genes, 1):
            query_genes = sorted(
                anchor_presence.get(
                    (target_chr, gene_row["gtr_gene"], query_chr), set()
                )
            )
            observed = 1 if query_genes else 0
            presence.append(observed)
            gene_presence_rows.append(
                {
                    "pair_label": pair_label,
                    "target_seqid": target_chr,
                    "target_gene_order": local_index,
                    "gtr_gene": gene_row["gtr_gene"],
                    "track_id": track,
                    "query_seqid": query_chr,
                    "anchored": observed,
                    "query_genes": ",".join(query_genes),
                }
            )
        if len(presence) < args.window_size:
            raise SystemExit(
                "ERROR: {} has fewer than 100 target genes".format(target_chr)
            )
        running = sum(presence[: args.window_size])
        for start in range(0, len(presence) - args.window_size + 1):
            if start > 0:
                running += presence[start + args.window_size - 1] - presence[start - 1]
            percent = int(round(100.0 * running / float(args.window_size)))
            window_rows.append(
                {
                    "pair_label": pair_label,
                    "target_seqid": target_chr,
                    "window_index": start + 1,
                    "window_start_gene": start + 1,
                    "window_end_gene": start + args.window_size,
                    "window_center_gene": start + (args.window_size + 1) / 2.0,
                    "track_id": track,
                    "query_seqid": query_chr,
                    "syntenic_genes_in_window": running,
                    "window_size": args.window_size,
                    "retention_percent": percent,
                    "retention_proportion": "{:.6f}".format(percent / 100.0),
                    "method": "FraBias_local_IO_port",
                }
            )

    write_tsv(
        outdir / "selected_chromosome_tracks.tsv",
        list(selected_rows[0].keys()),
        selected_rows,
    )
    write_tsv(
        outdir / "fractbias_genes.tsv",
        list(gene_presence_rows[0].keys()),
        gene_presence_rows,
    )
    write_tsv(
        outdir / "fractbias_windows.tsv", list(window_rows[0].keys()), window_rows
    )
    (outdir / "method_status.txt").write_text(
        "Algorithm: CoGe FraBias presence vectors and every-overlapping-window retention.\n"
        "Window: 100 ordered genes; denominator includes all target genes.\n"
        "Only the CoGe API/file reader was replaced by a local TSV reader.\n"
        "Anchor absence is not by itself a validated biological gene-loss call.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text(
        "FraBias 100-gene calculation completed.\n", encoding="ascii"
    )
    print(
        "[OK] FraBias windows={} selections={}: {}".format(
            len(window_rows), len(selected_rows), outdir
        )
    )


if __name__ == "__main__":
    main()
