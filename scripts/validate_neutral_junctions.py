#!/usr/bin/env python3
"""Validate external syntenic-junction candidates with neutral PZI anchors.

The script applies the two-sided, 20-gene flank rule used by Huang et al.
(2026) to copy1-copy2 anchors, then records corresponding support in the four
anonymous PZI haplotype assemblies.  It intentionally reports a structural
*candidate*, never a final fusion/fission or subgenome call.
"""

from __future__ import print_function

import argparse
import bisect
import csv
from collections import Counter, defaultdict
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-junctions", required=True)
    parser.add_argument("--query-bed", required=True)
    parser.add_argument("--copy2-bed", required=True)
    parser.add_argument("--copy-map", required=True)
    parser.add_argument("--copy1-copy2-anchors", required=True)
    parser.add_argument("--query-fai", default="",
                        help="Optional copy1 FASTA index. When supplied, writes 100-kb candidate windows for read/Hi-C review.")
    parser.add_argument("--hap-evidence", action="append", nargs=3, required=True,
                        metavar=("LABEL", "ANCHORS_NEW", "HAP_BED"))
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--flank-genes", type=int, default=20)
    parser.add_argument("--min-side-anchors", type=int, default=10)
    parser.add_argument("--min-haplotype-support", type=int, default=2)
    return parser.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path, fields, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def read_bed(path):
    genes = {}
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
            genes[fields[3]] = {"chrom": fields[0], "start": start, "end": end}
            per_chrom[fields[0]].append(fields[3])
    starts = {}
    rank = {}
    for chrom, values in per_chrom.items():
        values.sort(key=lambda gene: (genes[gene]["start"], genes[gene]["end"], gene))
        starts[chrom] = [genes[gene]["start"] for gene in values]
        for index, gene in enumerate(values):
            rank[gene] = index
    return genes, per_chrom, starts, rank


def read_anchor_lookup(path, qgenes, rgenes):
    lookup = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split()
            if len(fields) < 2:
                continue
            query, ref = fields[0], fields[1]
            if query in qgenes and ref in rgenes:
                lookup[query].append(ref)
    return lookup


def suffix_chrom(value):
    return value.split("_", 1)[1] if "_" in value else value


def expected_copy2_map(path):
    out = {}
    for row in read_tsv(path):
        normalized = (row.get("normalized_chr") or "").strip()
        copy1 = (row.get("copy1_chr") or normalized).strip()
        copy2 = (row.get("copy2_chr") or normalized).strip()
        if copy1 and copy2:
            out["copy1_" + copy1] = "copy2_" + copy2
    if not out:
        raise SystemExit("No copy1/copy2 chromosome pairs found in {}".format(path))
    return out


def flank_gene_sets(query_genes, chrom_genes, chrom_starts, chrom, coordinate, flank_size):
    values = chrom_genes.get(chrom, [])
    starts = chrom_starts.get(chrom, [])
    if not values:
        return [], []
    boundary = bisect.bisect_left(starts, coordinate)
    upstream = values[max(0, boundary - flank_size):boundary]
    downstream = values[boundary:min(len(values), boundary + flank_size)]
    return upstream, downstream


def side_summary(gene_list, lookup, ref_genes, expected_chrom=None):
    target_sets = defaultdict(set)
    anchored = set()
    for gene in gene_list:
        for ref in lookup.get(gene, []):
            anchored.add(gene)
            target_sets[ref_genes[ref]["chrom"]].add(gene)
    ordered = sorted(target_sets.items(), key=lambda item: (-len(item[1]), item[0]))
    top_target = ordered[0][0] if ordered else "NA"
    top_count = len(ordered[0][1]) if ordered else 0
    expected_count = len(target_sets.get(expected_chrom, set())) if expected_chrom else "NA"
    return {
        "window_genes": len(gene_list),
        "anchored_query_genes": len(anchored),
        "top_target": top_target,
        "top_target_anchor_genes": top_count,
        "expected_target_anchor_genes": expected_count,
        "target_count": len(target_sets),
    }


def read_fai(path):
    """Return sequence lengths keyed by FASTA record name."""
    lengths = {}
    if not path:
        return lengths
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                try:
                    lengths[fields[0]] = int(fields[1])
                except ValueError:
                    pass
    return lengths


def resolve_fai_chrom(query_bed_chrom, lengths):
    """Map a neutral BED label (copy1_chr01) to its source FASTA name (chr01)."""
    if query_bed_chrom in lengths:
        return query_bed_chrom
    if query_bed_chrom.startswith("copy1_"):
        candidate = query_bed_chrom[len("copy1_"):]
        if candidate in lengths:
            return candidate
    raise SystemExit("Query FAI has no sequence corresponding to {}".format(query_bed_chrom))


def main():
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    candidates = read_tsv(args.candidate_junctions)
    qgenes, qchroms, qstarts, _qrank = read_bed(args.query_bed)
    c2genes, _c2chroms, _c2starts, _c2rank = read_bed(args.copy2_bed)
    query_lengths = read_fai(args.query_fai)
    expected = expected_copy2_map(args.copy_map)
    copy_lookup = read_anchor_lookup(args.copy1_copy2_anchors, qgenes, c2genes)
    hap_data = []
    for label, anchor_path, bed_path in args.hap_evidence:
        hgenes, _hchroms, _hstarts, _hrank = read_bed(bed_path)
        hap_data.append((label, hgenes, read_anchor_lookup(anchor_path, qgenes, hgenes)))

    copy_rows = []
    hap_rows = []
    final_rows = []
    flank_rows = []
    window_rows = []
    for candidate in candidates:
        chrom = candidate["q_chr"]
        coordinate = int(float(candidate["boundary_mid_bp"]))
        upstream, downstream = flank_gene_sets(qgenes, qchroms, qstarts, chrom, coordinate, args.flank_genes)
        for side, genes in (("upstream", upstream), ("downstream", downstream)):
            for gene in genes:
                flank_rows.append({
                    "junction_id": candidate["junction_id"], "copy1_chr": chrom,
                    "boundary_mid_bp": coordinate, "side": side, "gene_id": gene,
                    "gene_start": qgenes[gene]["start"], "gene_end": qgenes[gene]["end"],
                })
        if query_lengths:
            fasta_chrom = resolve_fai_chrom(chrom, query_lengths)
            chrom_length = query_lengths[fasta_chrom]
            window_rows.append({
                "copy1_bed_chr": chrom,
                "copy1_fasta_chr": fasta_chrom,
                "start": max(0, coordinate - 100000),
                "end": min(chrom_length, coordinate + 100000),
                "junction_id": candidate["junction_id"],
                "boundary_mid_bp": coordinate,
                "window_definition": "100 kb on each side of candidate boundary; review only, not evidence by itself",
            })
        expected_target = expected.get(chrom, "NA")
        up = side_summary(upstream, copy_lookup, c2genes, expected_target)
        down = side_summary(downstream, copy_lookup, c2genes, expected_target)
        # Preserve both two-sided patterns.  A candidate component boundary
        # can have a strong copy2 target switch, so requiring one unchanged
        # expected target would bias discovery against the event being tested.
        internal_supported = int(
            up["top_target_anchor_genes"] >= args.min_side_anchors and
            down["top_target_anchor_genes"] >= args.min_side_anchors
        )
        same_expected_target = int(
            internal_supported and up["top_target"] == expected_target and down["top_target"] == expected_target
        )
        target_switch = int(
            internal_supported and up["top_target"] != "NA" and down["top_target"] != "NA" and
            up["top_target"] != down["top_target"]
        )
        if same_expected_target:
            internal_pattern = "same_expected_primary_homolog"
        elif target_switch:
            internal_pattern = "copy2_target_switch"
        elif internal_supported:
            internal_pattern = "two_sided_nonprimary_or_ambiguous"
        else:
            internal_pattern = "insufficient_two_sided_copy2_support"
        copy_row = {
            "junction_id": candidate["junction_id"], "copy1_chr": chrom, "boundary_mid_bp": coordinate,
            "expected_copy2_chr": expected_target,
            "upstream_window_genes": up["window_genes"], "upstream_anchored_query_genes": up["anchored_query_genes"],
            "upstream_top_copy2_chr": up["top_target"], "upstream_top_anchor_genes": up["top_target_anchor_genes"],
            "upstream_expected_anchor_genes": up["expected_target_anchor_genes"],
            "downstream_window_genes": down["window_genes"], "downstream_anchored_query_genes": down["anchored_query_genes"],
            "downstream_top_copy2_chr": down["top_target"], "downstream_top_anchor_genes": down["top_target_anchor_genes"],
            "downstream_expected_anchor_genes": down["expected_target_anchor_genes"],
            "copy1_copy2_two_sided_support": internal_supported,
            "same_expected_copy2_target": same_expected_target,
            "copy2_target_switch": target_switch,
            "copy1_copy2_pattern": internal_pattern,
        }
        copy_rows.append(copy_row)
        supported_haps = 0
        for label, hgenes, lookup in hap_data:
            hup = side_summary(upstream, lookup, hgenes)
            hdown = side_summary(downstream, lookup, hgenes)
            coverage = int(hup["top_target_anchor_genes"] >= args.min_side_anchors and
                           hdown["top_target_anchor_genes"] >= args.min_side_anchors)
            supported_haps += coverage
            hap_rows.append({
                "junction_id": candidate["junction_id"], "haplotype": label, "copy1_chr": chrom,
                "boundary_mid_bp": coordinate,
                "upstream_top_target": hup["top_target"], "upstream_top_anchor_genes": hup["top_target_anchor_genes"],
                "downstream_top_target": hdown["top_target"], "downstream_top_anchor_genes": hdown["top_target_anchor_genes"],
                "top_target_switch": int(hup["top_target"] != "NA" and hdown["top_target"] != "NA" and hup["top_target"] != hdown["top_target"]),
                "two_sided_anchor_coverage": coverage,
            })
        if candidate["external_support_count"] and int(candidate["external_support_count"]) >= 2 and internal_supported and supported_haps >= args.min_haplotype_support:
            status = "structural_candidate_needs_read_validation"
        elif candidate["external_support_count"] and int(candidate["external_support_count"]) >= 2 and internal_supported:
            status = "external_and_copy_support_needs_haplotype_validation"
        elif candidate["external_support_count"] and int(candidate["external_support_count"]) >= 2:
            status = "external_multiblock_only"
        else:
            status = "single_reference_only"
        final_rows.append({
            "junction_id": candidate["junction_id"], "copy1_chr": chrom, "boundary_mid_bp": coordinate,
            "external_support_count": candidate["external_support_count"], "supporting_references": candidate["supporting_references"],
            "copy1_copy2_two_sided_support": internal_supported,
            "same_expected_copy2_target": same_expected_target,
            "copy2_target_switch": target_switch,
            "copy1_copy2_pattern": internal_pattern,
            "haplotypes_with_two_sided_coverage": supported_haps,
            "status": status,
            "claim_boundary": "No fusion/fission or component ancestry is called without read-level validation and subsequent block-ancestry analysis.",
        })

    copy_fields = ["junction_id", "copy1_chr", "boundary_mid_bp", "expected_copy2_chr", "upstream_window_genes",
                   "upstream_anchored_query_genes", "upstream_top_copy2_chr", "upstream_top_anchor_genes",
                   "upstream_expected_anchor_genes", "downstream_window_genes", "downstream_anchored_query_genes",
                   "downstream_top_copy2_chr", "downstream_top_anchor_genes", "downstream_expected_anchor_genes",
                    "copy1_copy2_two_sided_support", "same_expected_copy2_target", "copy2_target_switch",
                    "copy1_copy2_pattern"]
    hap_fields = ["junction_id", "haplotype", "copy1_chr", "boundary_mid_bp", "upstream_top_target",
                  "upstream_top_anchor_genes", "downstream_top_target", "downstream_top_anchor_genes",
                  "top_target_switch", "two_sided_anchor_coverage"]
    final_fields = ["junction_id", "copy1_chr", "boundary_mid_bp", "external_support_count", "supporting_references",
                    "copy1_copy2_two_sided_support", "same_expected_copy2_target", "copy2_target_switch",
                    "copy1_copy2_pattern", "haplotypes_with_two_sided_coverage", "status", "claim_boundary"]
    write_tsv(outdir / "copy1_copy2_junction_support.tsv", copy_fields, copy_rows)
    write_tsv(outdir / "haplotype_junction_support.tsv", hap_fields, hap_rows)
    write_tsv(outdir / "junction_validation_summary.tsv", final_fields, final_rows)
    write_tsv(outdir / "junction_flank_genes.tsv",
              ["junction_id", "copy1_chr", "boundary_mid_bp", "side", "gene_id", "gene_start", "gene_end"], flank_rows)
    if query_lengths:
        write_tsv(outdir / "junction_100kb_review_windows.bed",
                  ["copy1_bed_chr", "copy1_fasta_chr", "start", "end", "junction_id", "boundary_mid_bp", "window_definition"], window_rows)
    write_tsv(outdir / "parameters.tsv", ["parameter", "value"], [
        {"parameter": "flank_genes", "value": args.flank_genes},
        {"parameter": "min_side_anchors", "value": args.min_side_anchors},
        {"parameter": "min_haplotype_support", "value": args.min_haplotype_support},
        {"parameter": "method_reference", "value": "Huang et al. 2026 Nat Commun 17:8244; 20 genes per flank and >10 conserved anchors per side"},
        {"parameter": "read_review_window", "value": "100 kb on each side of a validated candidate boundary, following Huang et al. 2026" if query_lengths else "not written; no --query-fai supplied"},
    ])
    print("[OK] Junction candidates={} with internal two-sided support={}".format(
        len(final_rows), sum(int(row["copy1_copy2_two_sided_support"]) for row in final_rows)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
