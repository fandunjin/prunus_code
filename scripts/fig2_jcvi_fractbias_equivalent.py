#!/usr/bin/env python3
"""Calculate JCVI-adapted 100-gene retention curves for Fig. 2c.

For each selected PZI A/B chromosome pair, GTR loci anchored by either selected
PZI chromosome define the displayed homologous target regions.  The eight query
tracks are scored for screened JCVI anchor presence in a centered 100-gene GTR
window.  The result is a reference-style fractionation curve, but it is not an
execution of the FractBias software and must be described as JCVI-adapted.
"""

from __future__ import print_function

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import numpy as np


TRACK_LABELS = {
    "pzi_hapA1": "PZI_A", "pzi_hapB1": "PZI_B",
    "psp_hapA1": "PSP_A", "psp_hapB1": "PSP_B",
    "ppa_hapA": "PPA_A", "ppa_hapB": "PPA_B",
    "pcam": "PCAM", "ppe": "PPE",
}
PAIRS = {
    "C06": {"PZI_A": "C06_A1", "PZI_B": "C06_B1"},
    "C02": {"PZI_A": "C02_A1", "PZI_B": "C02_B1"},
    "C03": {"PZI_A": "C03_A1", "PZI_B": "C03_B1"},
}
TRACK_ORDER = ("PZI_A", "PZI_B", "PSP_A", "PSP_B", "PPA_A", "PPA_B", "PCAM", "PPE")


def read_bed(path):
    chromosomes = []
    by_chromosome = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            seqid, gene = fields[0], fields[3]
            if seqid not in by_chromosome:
                chromosomes.append(seqid)
            by_chromosome[seqid].append(gene)
    return chromosomes, by_chromosome


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def rolling_mean(values, window):
    values = np.asarray(values, dtype=float)
    left = window // 2
    right = window - left
    cumulative = np.concatenate(([0.0], np.cumsum(values)))
    result = np.zeros(len(values), dtype=float)
    for index in range(len(values)):
        start = max(0, index - left)
        end = min(len(values), index + right)
        result[index] = (cumulative[end] - cumulative[start]) / float(end - start)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", required=True, help="gtr_alltaxa_syntenic_pairs.tsv")
    parser.add_argument("--gtr-bed", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--window", type=int, default=100)
    parser.add_argument("--step", type=int, default=5)
    args = parser.parse_args()
    if args.window != 100:
        raise SystemExit("ERROR: the registered Fig. 2c window is fixed at 100 genes")
    if args.step < 1:
        raise SystemExit("ERROR: --step must be positive")
    outdir = Path(args.outdir)
    if outdir.exists() and any(outdir.iterdir()):
        raise SystemExit("ERROR: output directory is not empty: {}".format(outdir))
    outdir.mkdir(parents=True, exist_ok=True)

    anchors = defaultdict(lambda: defaultdict(set))
    selected_pzi = defaultdict(lambda: defaultdict(set))
    with Path(args.pairs).open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            source = row.get("query_label", "")
            if source not in TRACK_LABELS or row.get("ref_label") != "gtr":
                continue
            track = TRACK_LABELS[source]
            gtr_gene = row.get("ref_gene", "")
            query_gene = row.get("query_gene", "")
            query_seqid = row.get("query_seqid", "")
            if not gtr_gene or not query_gene:
                continue
            anchors[track][gtr_gene].add(query_gene)
            if track in ("PZI_A", "PZI_B"):
                for pair, selected in PAIRS.items():
                    if query_seqid == selected[track]:
                        selected_pzi[pair][track].add(gtr_gene)
    for pair in PAIRS:
        if not selected_pzi[pair]["PZI_A"] or not selected_pzi[pair]["PZI_B"]:
            raise SystemExit("ERROR: selected PZI chromosomes have no GTR anchors for {}".format(pair))

    chromosomes, genes_by_chromosome = read_bed(args.gtr_bed)
    global_index = 0
    output_rows = []
    summary_rows = []
    for pair in ("C06", "C02", "C03"):
        target_region = selected_pzi[pair]["PZI_A"] | selected_pzi[pair]["PZI_B"]
        per_track_values = defaultdict(list)
        pair_global_index = 0
        for seqid in chromosomes:
            genes = genes_by_chromosome[seqid]
            for track in TRACK_ORDER:
                if track in ("PZI_A", "PZI_B"):
                    retained_genes = selected_pzi[pair][track]
                else:
                    retained_genes = set(anchors[track])
                binary = [1.0 if gene in target_region and gene in retained_genes else 0.0 for gene in genes]
                retention = rolling_mean(binary, args.window)
                for local_index in range(0, len(genes), args.step):
                    value = float(retention[local_index])
                    output_rows.append({
                        "pair_label": pair,
                        "track_id": track,
                        "gtr_global_index": pair_global_index + local_index,
                        "gtr_seqid": seqid,
                        "gtr_gene": genes[local_index],
                        "retention_proportion": "{:.8f}".format(value),
                        "window_genes": args.window,
                        "sampling_step_genes": args.step,
                        "method": "JCVI_screened_anchors_centered_rolling_window",
                    })
                    per_track_values[track].append(value)
            pair_global_index += len(genes)
        for track in TRACK_ORDER:
            summary_rows.append({
                "pair_label": pair,
                "track_id": track,
                "selected_gtr_anchor_loci": len(target_region),
                "mean_genomewide_window_retention": "{:.8f}".format(float(np.mean(per_track_values[track]))),
                "maximum_window_retention": "{:.8f}".format(float(np.max(per_track_values[track]))),
            })
        global_index = max(global_index, pair_global_index)
    write_tsv(outdir / "jcvi_fractbias_windows.tsv", list(output_rows[0].keys()), output_rows)
    write_tsv(outdir / "jcvi_fractbias_summary.tsv", list(summary_rows[0].keys()), summary_rows)
    (outdir / "method_status.txt").write_text(
        "JCVI-adapted fractionation curves from screened GTR anchors, 100-gene centered rolling windows, sampled every {} genes.\n".format(args.step)
        + "This is not a FractBias software run. Unanchored loci are candidate anchor attrition, not verified gene loss.\n",
        encoding="ascii",
    )
    (outdir / "COMPLETE").write_text("JCVI-adapted 100-gene curves completed.\n", encoding="ascii")
    print("[OK] Wrote {} curve rows across {} GTR genes: {}".format(len(output_rows), global_index, outdir))


if __name__ == "__main__":
    main()
