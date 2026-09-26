#!/usr/bin/env python3
"""Convert screened JCVI anchors into conservative syntenic fragments.

This implements the structural part of Huang et al. (2026, Nat. Commun.) in
a tool-agnostic way: anchors are first turned into non-overlapping, continuous
syntenic fragments; only then are adjacent fragments with different external
chromosome targets reported as *candidate* rejoined junctions.  A candidate is
not called a PZI fusion/fission by this script.
"""

from __future__ import print_function

import argparse
import csv
import math
from collections import Counter, defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-bed", required=True)
    parser.add_argument("--comparison", action="append", nargs=3, required=True,
                        metavar=("LABEL", "ANCHORS_NEW", "REFERENCE_BED"),
                        help="Repeat for each external reference")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-anchors", type=int, default=10,
                        help="Minimum unique anchors per retained fragment")
    parser.add_argument("--max-gap-genes", type=int, default=50,
                        help="Project discovery threshold for merging adjacent JCVI blocks")
    parser.add_argument("--max-gap-bp", type=int, default=300000,
                        help="Project discovery threshold for merging adjacent JCVI blocks")
    parser.add_argument("--junction-tolerance-bp", type=int, default=300000,
                        help="Coordinate tolerance when clustering support across references")
    parser.add_argument("--min-external-support", type=int, default=2)
    return parser.parse_args()


def natural_key(value):
    result = []
    number = ""
    for char in value:
        if char.isdigit():
            number += char
        else:
            if number:
                result.append((0, int(number)))
                number = ""
            result.append((1, char))
    if number:
        result.append((0, int(number)))
    return result


def read_bed(path):
    records = {}
    per_chrom = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                continue
            gene = fields[3]
            records[gene] = {"chrom": fields[0], "start": start, "end": end}
            per_chrom[fields[0]].append(gene)
    rank = {}
    for chrom, genes in per_chrom.items():
        genes.sort(key=lambda gene: (records[gene]["start"], records[gene]["end"], gene))
        for index, gene in enumerate(genes):
            rank[gene] = index
    return records, rank


def parse_blocks(path, qbed, qrank, rbed, rrank, label):
    blocks = []
    pairs = []
    block_index = 0

    def flush():
        nonlocal pairs, block_index
        if not pairs:
            return
        block_index += 1
        qchrom_counts = Counter(item[2]["chrom"] for item in pairs)
        rchrom_counts = Counter(item[3]["chrom"] for item in pairs)
        qchrom, qmode = qchrom_counts.most_common(1)[0]
        rchrom, rmode = rchrom_counts.most_common(1)[0]
        selected = [item for item in pairs if item[2]["chrom"] == qchrom and item[3]["chrom"] == rchrom]
        qgenes = sorted(set(item[0] for item in selected), key=lambda gene: qrank[gene])
        rgenes = sorted(set(item[1] for item in selected), key=lambda gene: rrank[gene])
        if not qgenes or not rgenes:
            pairs = []
            return
        slope = orientation([(qrank[item[0]], rrank[item[1]]) for item in selected])
        blocks.append({
            "label": label,
            "source_block": "{}_raw{:05d}".format(label, block_index),
            "q_chr": qchrom,
            "r_chr": rchrom,
            "orientation": slope,
            "anchor_pairs": len(selected),
            "unique_q_genes": len(qgenes),
            "unique_r_genes": len(rgenes),
            "q_start": min(qbed[gene]["start"] for gene in qgenes),
            "q_end": max(qbed[gene]["end"] for gene in qgenes),
            "r_start": min(rbed[gene]["start"] for gene in rgenes),
            "r_end": max(rbed[gene]["end"] for gene in rgenes),
            "q_rank_min": min(qrank[gene] for gene in qgenes),
            "q_rank_max": max(qrank[gene] for gene in qgenes),
            "r_rank_min": min(rrank[gene] for gene in rgenes),
            "r_rank_max": max(rrank[gene] for gene in rgenes),
            "q_chrom_mode_fraction": float(qmode) / len(pairs),
            "r_chrom_mode_fraction": float(rmode) / len(pairs),
        })
        pairs = []

    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text:
                continue
            if text.startswith("#"):
                flush()
                continue
            fields = text.split()
            if len(fields) < 2:
                continue
            qgene, rgene = fields[0], fields[1]
            if qgene in qbed and rgene in rbed:
                pairs.append((qgene, rgene, qbed[qgene], rbed[rgene]))
        flush()
    return blocks


def orientation(points):
    """Return + or - from rank covariance; a tie is unresolved."""
    by_query = defaultdict(list)
    for qrank, rrank in points:
        by_query[qrank].append(rrank)
    values = []
    for qrank, refs in by_query.items():
        refs.sort()
        values.append((qrank, refs[len(refs) // 2]))
    values.sort()
    if len(values) < 3:
        return "?"
    mean_q = sum(item[0] for item in values) / float(len(values))
    mean_r = sum(item[1] for item in values) / float(len(values))
    covariance = sum((item[0] - mean_q) * (item[1] - mean_r) for item in values)
    if math.isclose(covariance, 0.0):
        return "?"
    return "+" if covariance > 0 else "-"


def close_by_rank_or_bp(left, right, max_genes, max_bp):
    q_gap_genes = max(0, right["q_rank_min"] - left["q_rank_max"] - 1)
    q_gap_bp = max(0, right["q_start"] - left["q_end"])
    r_gap_genes = min(abs(right["r_rank_min"] - left["r_rank_max"]),
                      abs(left["r_rank_min"] - right["r_rank_max"]))
    r_gap_bp = min(abs(right["r_start"] - left["r_end"]),
                   abs(left["r_start"] - right["r_end"]))
    query_close = q_gap_genes < max_genes or q_gap_bp < max_bp
    ref_close = r_gap_genes < max_genes or r_gap_bp < max_bp
    return query_close and ref_close, q_gap_genes, q_gap_bp, r_gap_genes, r_gap_bp


def merge_blocks(raw_blocks, max_genes, max_bp, min_anchors):
    fragments = []
    grouped = defaultdict(list)
    for row in raw_blocks:
        grouped[(row["label"], row["q_chr"])].append(row)
    fragment_index = 0
    for key in sorted(grouped, key=lambda value: (value[0], natural_key(value[1]))):
        rows = sorted(grouped[key], key=lambda row: (row["q_start"], row["q_end"], row["source_block"]))
        active = None
        for row in rows:
            if active is None:
                active = dict(row)
                active["source_blocks"] = [row["source_block"]]
                continue
            merge_ok = False
            if active["r_chr"] == row["r_chr"] and active["orientation"] == row["orientation"]:
                merge_ok, _qg, _qb, _rg, _rb = close_by_rank_or_bp(active, row, max_genes, max_bp)
            if merge_ok:
                active["source_blocks"].append(row["source_block"])
                active["anchor_pairs"] += row["anchor_pairs"]
                active["unique_q_genes"] += row["unique_q_genes"]
                active["unique_r_genes"] += row["unique_r_genes"]
                active["q_start"] = min(active["q_start"], row["q_start"])
                active["q_end"] = max(active["q_end"], row["q_end"])
                active["r_start"] = min(active["r_start"], row["r_start"])
                active["r_end"] = max(active["r_end"], row["r_end"])
                active["q_rank_min"] = min(active["q_rank_min"], row["q_rank_min"])
                active["q_rank_max"] = max(active["q_rank_max"], row["q_rank_max"])
                active["r_rank_min"] = min(active["r_rank_min"], row["r_rank_min"])
                active["r_rank_max"] = max(active["r_rank_max"], row["r_rank_max"])
                active["q_chrom_mode_fraction"] = min(active["q_chrom_mode_fraction"], row["q_chrom_mode_fraction"])
                active["r_chrom_mode_fraction"] = min(active["r_chrom_mode_fraction"], row["r_chrom_mode_fraction"])
            else:
                fragment_index += 1
                active["fragment_id"] = "frag{:05d}".format(fragment_index)
                active["source_blocks"] = ",".join(active["source_blocks"])
                active["accepted"] = int(active["unique_q_genes"] >= min_anchors)
                fragments.append(active)
                active = dict(row)
                active["source_blocks"] = [row["source_block"]]
        if active is not None:
            fragment_index += 1
            active["fragment_id"] = "frag{:05d}".format(fragment_index)
            active["source_blocks"] = ",".join(active["source_blocks"])
            active["accepted"] = int(active["unique_q_genes"] >= min_anchors)
            fragments.append(active)
    return fragments


def make_junctions(fragments, max_genes, max_bp):
    junctions = []
    grouped = defaultdict(list)
    for row in fragments:
        if row["accepted"]:
            grouped[(row["label"], row["q_chr"])].append(row)
    for key in sorted(grouped, key=lambda value: (value[0], natural_key(value[1]))):
        rows = sorted(grouped[key], key=lambda row: (row["q_start"], row["q_end"], row["fragment_id"]))
        for left, right in zip(rows, rows[1:]):
            if left["r_chr"] == right["r_chr"] or right["q_rank_min"] <= left["q_rank_max"]:
                continue
            close, q_gap_genes, q_gap_bp, r_gap_genes, r_gap_bp = close_by_rank_or_bp(left, right, max_genes, max_bp)
            if not close:
                continue
            junctions.append({
                "label": key[0], "q_chr": key[1],
                "left_fragment": left["fragment_id"], "right_fragment": right["fragment_id"],
                "left_ref_chr": left["r_chr"], "right_ref_chr": right["r_chr"],
                "left_orientation": left["orientation"], "right_orientation": right["orientation"],
                "left_anchor_genes": left["unique_q_genes"], "right_anchor_genes": right["unique_q_genes"],
                "left_q_end": left["q_end"], "right_q_start": right["q_start"],
                "boundary_mid_bp": int((left["q_end"] + right["q_start"]) / 2.0),
                "query_gap_genes": q_gap_genes, "query_gap_bp": q_gap_bp,
                "reference_gap_genes": r_gap_genes, "reference_gap_bp": r_gap_bp,
            })
    return junctions


def cluster_junctions(junctions, tolerance, min_support):
    out = []
    grouped = defaultdict(list)
    for row in junctions:
        grouped[row["q_chr"]].append(row)
    index = 0
    for qchr in sorted(grouped, key=natural_key):
        rows = sorted(grouped[qchr], key=lambda row: row["boundary_mid_bp"])
        cluster = []
        center = None
        for row in rows:
            if not cluster or abs(row["boundary_mid_bp"] - center) <= tolerance:
                cluster.append(row)
                center = int(sum(item["boundary_mid_bp"] for item in cluster) / float(len(cluster)))
            else:
                index += 1
                out.append(finalize_cluster("junction{:04d}".format(index), qchr, cluster, center, min_support))
                cluster = [row]
                center = row["boundary_mid_bp"]
        if cluster:
            index += 1
            out.append(finalize_cluster("junction{:04d}".format(index), qchr, cluster, center, min_support))
    return out


def finalize_cluster(identifier, qchr, rows, center, min_support):
    labels = sorted(set(row["label"] for row in rows))
    details = ";".join("{}:{}+{}".format(row["label"], row["left_ref_chr"], row["right_ref_chr"]) for row in rows)
    return {
        "junction_id": identifier,
        "q_chr": qchr,
        "boundary_mid_bp": center,
        "external_support_count": len(labels),
        "supporting_references": ",".join(labels),
        "support_details": details,
        "status": "multi_reference_candidate" if len(labels) >= min_support else "single_reference_candidate",
        "interpretation": "requires copy1-copy2, phased-haplotype and read-level validation; not a fusion/fission call",
    }


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    qbed, qrank = read_bed(args.query_bed)
    raw = []
    for label, anchor_path, ref_bed_path in args.comparison:
        rbed, rrank = read_bed(ref_bed_path)
        blocks = parse_blocks(anchor_path, qbed, qrank, rbed, rrank, label)
        raw.extend(blocks)
    fragments = merge_blocks(raw, args.max_gap_genes, args.max_gap_bp, args.min_anchors)
    direct_junctions = make_junctions(fragments, args.max_gap_genes, args.max_gap_bp)
    clustered = cluster_junctions(direct_junctions, args.junction_tolerance_bp, args.min_external_support)

    raw_fields = ["label", "source_block", "q_chr", "r_chr", "orientation", "anchor_pairs", "unique_q_genes", "unique_r_genes",
                  "q_start", "q_end", "r_start", "r_end", "q_rank_min", "q_rank_max", "r_rank_min", "r_rank_max",
                  "q_chrom_mode_fraction", "r_chrom_mode_fraction"]
    frag_fields = ["fragment_id", "label", "q_chr", "r_chr", "orientation", "anchor_pairs", "unique_q_genes", "unique_r_genes",
                   "q_start", "q_end", "r_start", "r_end", "q_rank_min", "q_rank_max", "r_rank_min", "r_rank_max",
                   "q_chrom_mode_fraction", "r_chrom_mode_fraction", "source_blocks", "accepted"]
    junction_fields = ["label", "q_chr", "left_fragment", "right_fragment", "left_ref_chr", "right_ref_chr", "left_orientation",
                       "right_orientation", "left_anchor_genes", "right_anchor_genes", "left_q_end", "right_q_start",
                       "boundary_mid_bp", "query_gap_genes", "query_gap_bp", "reference_gap_genes", "reference_gap_bp"]
    cluster_fields = ["junction_id", "q_chr", "boundary_mid_bp", "external_support_count", "supporting_references",
                      "support_details", "status", "interpretation"]
    write_tsv(outdir / "raw_anchor_blocks.tsv", raw_fields, raw)
    write_tsv(outdir / "syntenic_fragments.tsv", frag_fields, fragments)
    write_tsv(outdir / "reference_junction_observations.tsv", junction_fields, direct_junctions)
    write_tsv(outdir / "candidate_junctions.tsv", cluster_fields, clustered)
    summary_rows = []
    for qchr in sorted(set(row["q_chr"] for row in fragments), key=natural_key):
        subset = [row for row in fragments if row["q_chr"] == qchr]
        candidates = [row for row in clustered if row["q_chr"] == qchr]
        summary_rows.append({"copy1_chr": qchr, "fragments_total": len(subset),
                             "fragments_accepted": sum(row["accepted"] for row in subset),
                             "candidate_junctions": len(candidates),
                             "multi_reference_candidates": sum(row["external_support_count"] >= args.min_external_support for row in candidates)})
    write_tsv(outdir / "fragment_summary_by_copy1_chr.tsv",
              ["copy1_chr", "fragments_total", "fragments_accepted", "candidate_junctions", "multi_reference_candidates"], summary_rows)
    write_tsv(outdir / "parameters.tsv", ["parameter", "value"], [
        {"parameter": "min_anchors", "value": args.min_anchors},
        {"parameter": "max_gap_genes", "value": args.max_gap_genes},
        {"parameter": "max_gap_bp", "value": args.max_gap_bp},
        {"parameter": "junction_tolerance_bp", "value": args.junction_tolerance_bp},
        {"parameter": "min_external_support", "value": args.min_external_support},
        {"parameter": "method_reference", "value": "Huang et al. 2026 Nat Commun 17:8244; JCVI-based discovery is a project-specific adaptation, while 20-gene/>10-anchor validation is applied separately"},
    ])
    print("[OK] Syntenic fragments={} accepted={} candidate_junctions={}".format(
        len(fragments), sum(row["accepted"] for row in fragments), len(clustered)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
