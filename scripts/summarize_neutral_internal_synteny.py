#!/usr/bin/env python3
"""Summarize copy1-copy2 and neutral haplotype JCVI anchors.

This script reports reciprocal chromosome specificity, multi-target mapping and
order continuity.  For copy1-copy2 it additionally tests the supplied
assembly-level homolog map.  No A/B labels are used in the output tables.
"""

from __future__ import annotations

import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pair-plan", required=True)
    p.add_argument("--pair-root", required=True)
    p.add_argument("--input-root", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--copy-map", required=True,
                   help="Copy1/copy2 map (normalized map or corrected_homolog_pairs.tsv)")
    p.add_argument("--same-chr-alignment", required=True,
                   help="Whole-genome copy1/copy2 normalized chromosome alignment summary")
    p.add_argument("--gtr-profiles", required=True)
    p.add_argument("--min-cell-anchors", type=int, default=8)
    p.add_argument("--min-cell-fraction", type=float, default=0.05)
    return p.parse_args()


def read_tsv(path):
    with Path(path).open("r", encoding="utf-8", errors="replace", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, fields, rows):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)


def natural_key(value):
    chunks = []; buf = ""
    for ch in value:
        if ch.isdigit(): buf += ch
        else:
            if buf: chunks.append(int(buf)); buf = ""
            chunks.append(ch)
    if buf: chunks.append(int(buf))
    return chunks


def read_bed(path):
    gene = {}; order = defaultdict(list)
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"): continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 4: continue
            try: start, end = int(f[1]), int(f[2])
            except ValueError: continue
            gene[f[3]] = {"chr": f[0], "start": start, "end": end}
            order[f[0]].append(f[3])
    rank = {}
    for chrom, genes in order.items():
        for i, value in enumerate(genes): rank[value] = i
    return gene, rank, order


def parse_anchors(path):
    pairs = []; seen = set()
    with Path(path).open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            text = line.strip()
            if not text or text.startswith("#"): continue
            f = text.split()
            if len(f) < 2: continue
            pair = (f[0], f[1])
            if pair not in seen:
                seen.add(pair); pairs.append(pair)
    return pairs


def _average_ranks(values):
    ranked = sorted(enumerate(values), key=lambda item: item[1])
    out = [0.0] * len(values); i = 0
    while i < len(ranked):
        j = i + 1
        while j < len(ranked) and ranked[j][1] == ranked[i][1]: j += 1
        rank = (i + 1 + j) / 2.0
        for k in range(i, j): out[ranked[k][0]] = rank
        i = j
    return out


def _pearson(values_a, values_b):
    if len(values_a) < 2: return 0.0
    mean_a = sum(values_a) / len(values_a); mean_b = sum(values_b) / len(values_b)
    numerator = sum((a - mean_a) * (b - mean_b) for a, b in zip(values_a, values_b))
    denom_a = math.sqrt(sum((a - mean_a) ** 2 for a in values_a))
    denom_b = math.sqrt(sum((b - mean_b) ** 2 for b in values_b))
    return numerator / (denom_a * denom_b) if denom_a and denom_b else 0.0


def monotonicity(records):
    """Order continuity for one chromosome pair, allowing inversion."""
    by_q = defaultdict(list)
    for qrank, rrank in records: by_q[qrank].append(rrank)
    points = [(qrank, sorted(values)[len(values) // 2]) for qrank, values in by_q.items()]
    points.sort()
    if len(points) < 3: return 0.0
    qrank = [q for q, _ in points]; rrank = [r for _, r in points]
    rho = abs(_pearson(_average_ranks(qrank), _average_ranks(rrank)))
    steps = [b - a for a, b in zip(rrank, rrank[1:]) if b != a]
    if not steps: return 0.0
    positive = sum(step > 0 for step in steps)
    direction = max(positive, len(steps) - positive) / len(steps)
    return (rho + direction) / 2.0


def weighted_monotonicity(records_by_target):
    total = sum(len(records) for records in records_by_target.values())
    if not total: return 0.0
    return sum(monotonicity(records) * len(records) for records in records_by_target.values()) / total


def summarize_direction(matrix, totals, min_anchors, min_fraction):
    rows = {}
    for chrom, values in matrix.items():
        total = totals.get(chrom, 0)
        threshold = max(min_anchors, min_fraction * total) if total else float("inf")
        supported = [target for target, n in values.items() if n >= threshold]
        supported.sort(key=natural_key)
        retained = sum(values[x] for x in supported)
        ordered = sorted(values.items(), key=lambda x: x[1], reverse=True)
        top = ordered[0] if ordered else ("NA", 0)
        rows[chrom] = {"supported": supported, "retained": retained,
                       "top": top[0], "top_fraction": top[1] / total if total else 0.0,
                       "multi": len(supported) > 1}
    return rows


def suffix_chrom(value):
    # Neutral names are label_chrNN. Return chrNN for assembly-map joins.
    return value.split("_", 1)[1] if "_" in value else value


def read_expected_copy2_map(path):
    """Read either the packaged normalized map or the corrected pair table."""
    rows = read_tsv(Path(path))
    mapping = {}
    for row in rows:
        normalized = row.get("normalized_chr", "").strip()
        if normalized:
            # The neutral inputs use copy1_chrNN/copy2_chrNN IDs, so the
            # normalized legacy map is an identity mapping between copies.
            mapping[normalized] = normalized
            continue
        copy1 = row.get("copy1_chr", "").strip()
        copy2 = row.get("copy2_chr", "").strip()
        if copy1 and copy2:
            mapping[copy1] = copy2
    if not mapping:
        raise SystemExit("No copy1/copy2 chromosome pairs found in {}".format(path))
    return mapping


def main():
    a = parse_args(); out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    pair_root = Path(a.pair_root); input_root = Path(a.input_root)
    expected_copy2 = read_expected_copy2_map(a.copy_map)
    same_chr_alignment = {row["chr"]: row for row in read_tsv(Path(a.same_chr_alignment))}
    gtr_profiles = {row["copy1_chr"]: row for row in read_tsv(Path(a.gtr_profiles))}
    plans = read_tsv(Path(a.pair_plan))
    pair_rows = []; matrix_rows = []; pair_cache = {}

    for plan in plans:
        comp, query, ref = plan["comparison"], plan["query_label"], plan["ref_label"]
        work = pair_root / comp
        anchor_path = work / f"{query}.{ref}.anchors.new"
        qbed_path = input_root / query / f"{query}.bed"
        rbed_path = input_root / ref / f"{ref}.bed"
        if not anchor_path.is_file() or not qbed_path.is_file() or not rbed_path.is_file():
            pair_rows.append({"comparison": comp, "query": query, "ref": ref, "status": "missing_input"})
            continue
        qgene, qrank, qorder = read_bed(qbed_path); rgene, rrank, rorder = read_bed(rbed_path)
        anchors = []
        qmatrix = defaultdict(lambda: defaultdict(int)); rmatrix = defaultdict(lambda: defaultdict(int))
        qrecords = defaultdict(lambda: defaultdict(list)); rrecords = defaultdict(lambda: defaultdict(list))
        q_unique = defaultdict(set); r_unique = defaultdict(set)
        for qid, rid in parse_anchors(anchor_path):
            if qid not in qgene or rid not in rgene: continue
            qchr, rchr = qgene[qid]["chr"], rgene[rid]["chr"]
            qmatrix[qchr][rchr] += 1; rmatrix[rchr][qchr] += 1
            qrecords[qchr][rchr].append((qrank[qid], rrank[rid]))
            rrecords[rchr][qchr].append((rrank[rid], qrank[qid]))
            q_unique[qchr].add(qid); r_unique[rchr].add(rid)
            anchors.append((qid, rid, qchr, rchr))
        if not anchors:
            pair_rows.append({"comparison": comp, "query": query, "ref": ref, "status": "no_usable_anchors"})
            continue
        qtotals = {x: sum(v.values()) for x, v in qmatrix.items()}; rtotals = {x: sum(v.values()) for x, v in rmatrix.items()}
        qsum = summarize_direction(qmatrix, qtotals, a.min_cell_anchors, a.min_cell_fraction)
        rsum = summarize_direction(rmatrix, rtotals, a.min_cell_anchors, a.min_cell_fraction)
        q_order = {x: weighted_monotonicity(qrecords[x]) for x in qmatrix}
        r_order = {x: weighted_monotonicity(rrecords[x]) for x in rmatrix}
        qret = sum(qsum[x]["retained"] for x in qsum); rret = sum(rsum[x]["retained"] for x in rsum)
        row = {"comparison": comp, "query": query, "ref": ref, "status": "ok", "anchor_pairs": len(anchors),
               "query_chromosomes": len(qorder), "reference_chromosomes": len(rorder),
               "query_gene_coverage": sum(len(v) for v in q_unique.values()) / max(1, len(qgene)),
               "reference_gene_coverage": sum(len(v) for v in r_unique.values()) / max(1, len(rgene)),
               "query_top_specificity": sum(qsum[x]["top_fraction"] for x in qsum) / max(1, len(qsum)),
               "reference_top_specificity": sum(rsum[x]["top_fraction"] for x in rsum) / max(1, len(rsum)),
               "query_multi_target_chromosomes": sum(qsum[x]["multi"] for x in qsum),
               "reference_multi_target_chromosomes": sum(rsum[x]["multi"] for x in rsum),
               "query_order_continuity": sum(q_order.values()) / max(1, len(q_order)),
               "reference_order_continuity": sum(r_order.values()) / max(1, len(r_order))}
        if comp == "copy1_vs_copy2":
            expected = 0; correct = 0
            for qid, _rid, qchr, rchr in anchors:
                qbase, rbase = suffix_chrom(qchr), suffix_chrom(rchr)
                expected += int(qbase in expected_copy2)
                correct += int(expected_copy2.get(qbase) == rbase)
            row["expected_homolog_anchor_fraction"] = correct / expected if expected else "NA"
        else:
            row["expected_homolog_anchor_fraction"] = "NA"
        pair_rows.append(row)
        # Matrix rows are written in long form for all comparisons.
        for qchr, values in qmatrix.items():
            for rchr, count in values.items():
                matrix_rows.append({"comparison": comp, "query_chr": qchr, "reference_chr": rchr,
                                    "anchor_pairs": count, "query_total": qtotals[qchr],
                                    "reference_total": rtotals[rchr], "query_fraction": count / qtotals[qchr],
                                    "reference_fraction": count / rtotals[rchr]})
        pair_cache[comp] = {"qsum": qsum, "qmatrix": qmatrix, "qorder": q_order, "anchors": anchors,
                            "qgene": qgene, "rgene": rgene}

    pair_fields = ["comparison", "query", "ref", "status", "anchor_pairs", "query_chromosomes", "reference_chromosomes",
                   "query_gene_coverage", "reference_gene_coverage", "query_top_specificity", "reference_top_specificity",
                   "query_multi_target_chromosomes", "reference_multi_target_chromosomes", "query_order_continuity",
                   "reference_order_continuity", "expected_homolog_anchor_fraction"]
    write_tsv(out / "pair_summary.tsv", pair_fields, pair_rows)
    write_tsv(out / "chromosome_matrix.tsv",
              ["comparison", "query_chr", "reference_chr", "anchor_pairs", "query_total", "reference_total", "query_fraction", "reference_fraction"], matrix_rows)

    # Copy1-copy2 resolution against the previous neutral GTR profile.
    resolution = []
    internal = pair_cache.get("copy1_vs_copy2")
    if internal:
        qsum = internal["qsum"]
        for qchr in sorted(qsum, key=natural_key):
            base = suffix_chrom(qchr); profile = gtr_profiles.get(base, {})
            gtr_labels = [x for x in profile.get("supported_gtr_blocks", "").split(",") if x and x != "NA"]
            targets = qsum[qchr]["supported"]
            expected_target = expected_copy2.get(base, "NA")
            expected_fraction = "NA"
            expected_is_top = "NA"
            if base in expected_copy2:
                values = internal["qmatrix"].get(qchr, {})
                total = sum(values.values())
                expected_fraction = values.get(f"copy2_{expected_target}", 0) / total if total else 0.0
                expected_is_top = int(qsum[qchr]["top"] == f"copy2_{expected_target}")
            aln = same_chr_alignment.get(base, {})
            same_bp = int(aln.get("same_chr_plus_aligned_bp", 0)) + int(aln.get("same_chr_minus_aligned_bp", 0))
            other_bp = int(aln.get("other_chr_aligned_bp", 0))
            same_fraction = same_bp / (same_bp + other_bp) if same_bp + other_bp else 0.0
            breakpoint_candidate = int((expected_is_top != 1 or expected_fraction < 0.5) and same_fraction < 0.5)
            if expected_is_top == 1 and expected_fraction >= 0.5 and same_fraction >= 0.5:
                if len(gtr_labels) > 1:
                    interpretation = "GTR_multiblock_stable_primary_internal_homolog"
                else:
                    interpretation = "stable_primary_internal_homolog"
            elif breakpoint_candidate:
                interpretation = "candidate_copy_specific_breakpoint"
            else:
                interpretation = "ambiguous_internal_structure"
            resolution.append({"copy1_chr": base, "gtr_supported_blocks": ",".join(gtr_labels) or "NA",
                               "copy2_supported_targets": ",".join(targets) or "NA", "copy2_target_count": len(targets),
                               "assembly_expected_copy2_chr": expected_target, "expected_target_is_top": expected_is_top,
                               "expected_target_anchor_fraction": expected_fraction,
                               "same_chr_genome_alignment_fraction": same_fraction,
                               "same_chr_aligned_bp": same_bp, "other_chr_aligned_bp": other_bp,
                               "query_top_specificity": qsum[qchr]["top_fraction"],
                               "copy_specific_breakpoint_candidate": breakpoint_candidate,
                               "interpretation": interpretation})
    write_tsv(out / "copy1_copy2_gtr_resolution.tsv",
              ["copy1_chr", "gtr_supported_blocks", "copy2_supported_targets", "copy2_target_count",
               "assembly_expected_copy2_chr", "expected_target_is_top", "expected_target_anchor_fraction",
               "same_chr_genome_alignment_fraction", "same_chr_aligned_bp", "other_chr_aligned_bp",
               "query_top_specificity", "copy_specific_breakpoint_candidate", "interpretation"], resolution)

    independent = [r for r in pair_rows if r.get("query", "").startswith("hap") and r.get("ref", "").startswith("hap")]
    write_tsv(out / "independent_haplotype_pairs.tsv", pair_fields, independent)
    cross = [r for r in pair_rows if (r.get("query", "") == "copy1" and r.get("ref", "").startswith("hap")) or
             (r.get("query", "") == "copy2" and r.get("ref", "").startswith("hap"))]
    write_tsv(out / "copy_to_haplotype_pairs.tsv", pair_fields, cross)
    write_tsv(out / "summary_parameters.tsv", ["parameter", "value"], [
        {"parameter": "min_cell_anchors", "value": a.min_cell_anchors},
        {"parameter": "min_cell_fraction", "value": a.min_cell_fraction},
        {"parameter": "copy_normalization_map", "value": str(a.copy_map)},
        {"parameter": "same_chr_genome_alignment", "value": str(a.same_chr_alignment)},
        {"parameter": "order_continuity_definition", "value": "anchor-weighted within-chromosome-pair score: mean(abs(Spearman rho), dominant step direction)"},
        {"parameter": "breakpoint_rule", "value": "candidate only if expected copy2 target is not dominant or <0.50 anchors AND same-chromosome genome alignment is <0.50"},
        {"parameter": "interpretation", "value": "internal multi-target mapping is not a standalone breakpoint call"},
    ])
    (out / "README.txt").write_text(
        "Neutral internal synteny summary. copy1/copy2/hap01-hap04 are analysis labels only.\n"
        "expected_homolog_anchor_fraction uses the assembly-level copy1-copy2 homolog map.\n"
        "A GTR-multiblock chromosome with a stable primary copy2 homolog is not a PZI\n"
        "fusion/fission call. A copy-specific breakpoint candidate requires both a disrupted\n"
        "primary internal target and weak same-chromosome genome alignment, then requires\n"
        "coordinate-level inspection.\n", encoding="ascii")
    print(f"[OK] Neutral internal synteny summarized: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
